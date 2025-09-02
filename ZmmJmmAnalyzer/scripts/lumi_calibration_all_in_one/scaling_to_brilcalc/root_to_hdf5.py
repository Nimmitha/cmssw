#!/usr/bin/env python3
"""
Stage 1: ROOT to HDF5 Converter (Memory-Optimized)
Processes large ROOT files with minimal memory by:
1. Reading in chunks
2. Applying cuts
3. Writing data immediately to HDF5 without accumulation
"""

import uproot as up
import numpy as np
import h5py
import os
# import gc
import json
import tempfile
import shutil
from tqdm import tqdm
import argparse

class IncrementalHDF5Writer:
    """
    Handles incremental writing to HDF5 file without keeping data in memory.
    Uses resizable datasets that grow as needed.
    """
    
    def __init__(self, output_path, bin_width=50, compression='lzf'):
        self.output_path = output_path
        self.bin_width = bin_width
        self.compression = compression
        self.compression_level = 4 if compression == 'gzip' else None
        self.hf = None
        self.bin_datasets = {}
        self.bin_counts = {}
        self.chunk_size = 50000  # Initial chunk size for HDF5 datasets
        self.event_dtype = np.dtype([
            ("mass", np.float32),
            ("npv", np.int32),
        ])
        
    def __enter__(self):
        self.hf = h5py.File(self.output_path, 'w')
        # Create groups
        self.metadata_group = self.hf.create_group('metadata')
        self.bins_group = self.hf.create_group('bins')
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.hf:
            # Resize datasets to actual size
            for bin_key in self.bin_datasets:
                self._resize_bin_datasets(bin_key, final=True)
            self.hf.close()
    
    def _get_or_create_bin(self, bin_key, run, lumiblock):
        """Get existing bin group or create new one with resizable datasets."""
        if bin_key not in self.bin_datasets:
            # Create bin group
            bin_group = self.bins_group.create_group(bin_key)
            bin_group.attrs['run'] = run
            bin_group.attrs['lumiblock'] = lumiblock
            bin_group.attrs['n_events'] = 0
            
            # Create resizable datasets with initial size
            self.bin_datasets[bin_key] = {
                'events': bin_group.create_dataset(
                    'events',
                    shape=(150000,),
                    maxshape=(None,),
                    dtype=self.event_dtype,
                    chunks=(min(self.chunk_size, 50000),),
                    compression=self.compression,
                    compression_opts=self.compression_level
                ),

                'group': bin_group
            }
            self.bin_counts[bin_key] = 0
            
        return self.bin_datasets[bin_key]
    
    def _resize_bin_datasets(self, bin_key, final=False):
        """Resize datasets for a bin."""
        if bin_key not in self.bin_datasets:
            return
            
        current_size = self.bin_counts[bin_key]
        datasets = self.bin_datasets[bin_key]
        
        if final:
            # Final resize to exact size
            new_size = current_size
            datasets['group'].attrs['n_events'] = current_size
        else:
            new_size = datasets['events'].shape[0] + self.chunk_size

        # Resize datasets
        datasets['events'].resize((new_size,))
    
    def append_to_bin(self, bin_key, run, lumiblock, data):
        """Append data to a bin's datasets."""
        datasets = self._get_or_create_bin(bin_key, run, lumiblock)
        
        n_new = len(data)
        if n_new == 0:
            return
            
        current_size = self.bin_counts[bin_key]

        # Ensure we have enough space
        while current_size + n_new > datasets['events'].shape[0]:
            self._resize_bin_datasets(bin_key)
            
        # Write data
        datasets['events'][current_size:current_size + n_new] = data
        
        # Update count
        self.bin_counts[bin_key] = current_size + n_new
        datasets['group'].attrs['n_events'] = self.bin_counts[bin_key]
    
    def set_metadata(self, **kwargs):
        """Set metadata attributes."""
        for key, value in kwargs.items():
            if isinstance(value, dict):
                self.metadata_group.attrs[key] = json.dumps(value)
            else:
                self.metadata_group.attrs[key] = value


def process_root_to_hdf5_streaming(
    root_file_path,
    hdf5_output_path,
    cuts,
    bin_width=50,
    chunk_size="200 MB",
    compression='lzf',
    use_temp_file=True
):
    """
    Convert ROOT file to HDF5 with true streaming (minimal memory usage).
    
    Parameters
    ----------
    root_file_path : str
        Path to input ROOT file
    hdf5_output_path : str
        Path to output HDF5 file
    cuts : dict
        Dictionary of cuts: {branch_name: (low, high)}
    bin_width : int
        Width of lumiblock bins
    chunk_size : str or int
        Size for uproot.iterate
    compression : str
        HDF5 compression algorithm
    compression_level : int
        Compression level (1-9)
    use_temp_file : bool
        Write to temp file first (safer for large files)
    """
    
    print(f"Converting ROOT file to HDF5")
    print(f"  Input: {root_file_path}")
    print(f"  Output: {hdf5_output_path}")
    
    # Branches to read
    branches_to_read = ["Run", "LumiBlock", "B_J1_mass", "nPV"]
    
    # Add cut branches if not already included
    for cut_var in cuts.keys():
        if cut_var not in branches_to_read:
            branches_to_read.append(cut_var)
    
    # Get total entries for progress bar
    with up.open(root_file_path) as f:
        total_entries = f["ntuple"].num_entries
    print(f"Total entries in ROOT file: {total_entries:,}")
    
    # Use temporary file if requested (safer for crashes)
    if use_temp_file:
        temp_dir = os.path.dirname(hdf5_output_path)
        temp_fd, temp_path = tempfile.mkstemp(suffix='.h5', dir=temp_dir)
        os.close(temp_fd)
        write_path = temp_path
    else:
        write_path = hdf5_output_path
    
    processed = 0
    kept_total = 0
    unique_bins = set()

    print(f"Using chunk size: {chunk_size}")
    print(f"Using bin width: {bin_width}")
    print(f"Using compression: {compression}")
    print(f"Using temporary file: {use_temp_file}")
    print(f"Writing to: {write_path}")
    print("Starting processing...\n")
    
    try:
        with IncrementalHDF5Writer(write_path, bin_width, compression) as writer:
            print("Processing ROOT file in chunks...")
            # Process ROOT file in chunks
            with tqdm(total=total_entries, desc="Processing ROOT file") as pbar:
                for arrays in up.iterate(
                    f"{root_file_path}:ntuple",
                    branches_to_read,
                    step_size=chunk_size,
                    library="np"
                ):
                    chunk_size_entries = len(arrays["Run"])
                    
                    # Apply cuts
                    mask = np.ones(chunk_size_entries, dtype=bool)
                    for var, (low, high) in cuts.items():
                        if var in arrays:
                            if low is not None:
                                mask &= (arrays[var] > low)
                            if high is not None:
                                mask &= (arrays[var] < high)
                    
                    kept = int(mask.sum())
                    kept_total += kept
                    
                    if kept > 0:
                        # Extract filtered data
                        runs = arrays["Run"][mask]
                        lumis = arrays["LumiBlock"][mask]
                        masses = arrays["B_J1_mass"][mask]
                        npvs = arrays["nPV"][mask]

                        chunk_data = np.zeros(len(masses), dtype=writer.event_dtype)
                        chunk_data['mass'] = masses
                        chunk_data['npv'] = npvs
                        
                        # Group by bins and write immediately
                        # Create a structured array for efficient grouping
                        bin_keys = np.char.add(
                            np.char.add(runs.astype(str), '_'),
                            ((lumis // bin_width) * bin_width).astype(str)
                        )

                        # Process each unique bin in this chunk
                        for bin_key in np.unique(bin_keys):
                            bin_mask = bin_keys == bin_key
                            
                            # Parse bin info
                            run_str, lumibin_str = bin_key.split('_')
                            run = int(run_str)
                            lumiblock = int(lumibin_str)
                            
                            # Write this bin's data immediately
                            writer.append_to_bin(
                                bin_key,
                                run,
                                lumiblock,
                                chunk_data[bin_mask]
                            )
                            
                            unique_bins.add(bin_key)
                    
                    processed += chunk_size_entries
                    pbar.update(chunk_size_entries)
                    pbar.set_postfix({
                        'kept': f'{kept_total:,}',
                        'bins': len(unique_bins),
                        'mem_MB': get_memory_usage()
                    })
            
            # Set metadata
            writer.set_metadata(
                total_entries_original=total_entries,
                total_entries_filtered=kept_total,
                bin_width=bin_width,
                cuts=cuts,
                n_bins=len(unique_bins)
            )
        
        # Move temp file to final location
        if use_temp_file:
            shutil.move(temp_path, hdf5_output_path)
        
        # Get file size
        file_size_mb = os.path.getsize(hdf5_output_path) / (1024 * 1024)
        
        print(f"\n✓ HDF5 file created: {hdf5_output_path}")
        print(f"  File size: {file_size_mb:.1f} MB")
        print(f"  Original entries: {total_entries:,}")
        print(f"  Filtered entries: {kept_total:,} ({kept_total/total_entries*100:.1f}%)")
        print(f"  Number of bins: {len(unique_bins)}")
        print(f"  Average events/bin: {kept_total/len(unique_bins):.0f}")
        print(f"  Peak memory usage: ~{get_memory_usage():.1f} MB")
        
    except Exception as e:
        # Clean up temp file on error
        if use_temp_file and os.path.exists(temp_path):
            os.remove(temp_path)
        raise e


def get_memory_usage():
    """Get current memory usage in MB."""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        # If psutil not available, return 0
        return 0


def verify_hdf5_file(hdf5_path):
    """Verify the created HDF5 file structure and print statistics."""
    print(f"\nVerifying HDF5 file: {hdf5_path}")
    
    with h5py.File(hdf5_path, 'r') as hf:
        # Check metadata
        metadata = dict(hf['metadata'].attrs)
        print(f"Metadata:")
        for key, value in metadata.items():
            if key == 'cuts':
                value = json.loads(value)
            print(f"  {key}: {value}")
        
        # Check bins
        bins_group = hf['bins']
        n_bins = len(bins_group.keys())
        print(f"\nBins: {n_bins}")
        
        # Sample a few bins
        total_events = 0
        min_events = float('inf')
        max_events = 0
        
        for i, bin_key in enumerate(bins_group.keys()):
            bin_group = bins_group[bin_key]
            n_events = bin_group.attrs['n_events']
            total_events += n_events
            min_events = min(min_events, n_events)
            max_events = max(max_events, n_events)
            
            if i < 3:  # Show first 3 bins as examples
                print(f"  {bin_key}: {n_events:,} events")
        
        if n_bins > 3:
            print(f"  ... and {n_bins - 3} more bins")
        
        print(f"\nStatistics:")
        print(f"  Total events: {total_events:,}")
        print(f"  Min events/bin: {min_events:,}")
        print(f"  Max events/bin: {max_events:,}")
        print(f"  Avg events/bin: {total_events/n_bins:.0f}")


def main():
    parser = argparse.ArgumentParser(description='Convert ROOT file to HDF5 format (memory-optimized)')
    parser.add_argument('root_file', help='Input ROOT file path')
    parser.add_argument('hdf5_file', help='Output HDF5 file path')
    parser.add_argument('--bin-width', type=int, default=50, help='Lumiblock bin width')
    parser.add_argument('--chunk-size', default='200 MB', help='Chunk size for reading')
    parser.add_argument('--compression', default='lzf', choices=['gzip', 'lzf'],
                       help='HDF5 compression algorithm')
    parser.add_argument('--no-temp-file', action='store_true',
                       help='Write directly to output file (faster but less safe)')
    parser.add_argument('--verify', action='store_true',
                       help='Verify the output HDF5 file after creation')
    
    args = parser.parse_args()
    
    # Define cuts (you can modify this or load from config)
    cuts = {
        "B_J1_mass": (2.7, 3.5),
        "B_Mu1_pt": (5000, None),
        "B_Mu2_pt": (5000, None)
    }
    
    # Process the file
    process_root_to_hdf5_streaming(
        args.root_file,
        args.hdf5_file,
        cuts,
        bin_width=args.bin_width,
        chunk_size=args.chunk_size,
        compression=args.compression,
        use_temp_file=not args.no_temp_file
    )
    
    # Optionally verify the output
    if args.verify:
        verify_hdf5_file(args.hdf5_file)


if __name__ == "__main__":
    main()