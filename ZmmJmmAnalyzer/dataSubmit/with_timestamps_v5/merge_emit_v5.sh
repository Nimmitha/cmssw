#!/bin/bash
# Merge the with_timestamps_v5 CRAB outputs: one file per era dataset with all 8 PDs,
# PDMLM_mm_<year><era><version>_emit_v5.root (as PDMLM_mm_2024F1_emit_v4.root).
# Usage: ./merge_emit_v5.sh [output_dir]   (run after cmsenv, for hadd)

EOS_DIR=/eos/uscms/store/user/wkarunar/emit
OUTPUT_DIR="${1:-merged_root_files/emit_v5}"
mkdir -p "$OUTPUT_DIR"

# Task tags (e.g. 2024D1) from the PD0 task directories
for tag in $(ls -d $EOS_DIR/ParkingDoubleMuonLowMass0/crab_PDMLM0_mm_*_emit_v5 | sed 's/.*_mm_\(.*\)_emit_v5/\1/'); do
    output_file="$OUTPUT_DIR/PDMLM_mm_${tag}_emit_v5.root"
    input_files=$(find $EOS_DIR/ParkingDoubleMuonLowMass[0-7]/crab_PDMLM[0-7]_mm_${tag}_emit_v5 -mindepth 3 -maxdepth 3 -type f -name "*.root" | sort)
    echo ""
    echo "$tag: $(echo "$input_files" | wc -w) files from $(ls -d $EOS_DIR/ParkingDoubleMuonLowMass[0-7]/crab_PDMLM[0-7]_mm_${tag}_emit_v5 | wc -l) PDs -> $output_file"
    if [ -n "$input_files" ]; then
        hadd -f "$output_file" $input_files
    else
        echo "No root files found for $tag"
    fi
done

echo ""
echo "Script execution complete."
