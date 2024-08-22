#!/bin/bash
# Run this from the src directory

# Set up the environment
cmsenv

# Run the rootupler for each of the 10 input files
echo "Starting the loop to process files..."
for i in {1..10}; do
    echo "Processing file with index $i..."
    sed -i "s/zeejmm_2018\/zeejmm_ss\/MiniAOD\/MiniAOD_[0-9]\+\.root/zeejmm_2018\/zeejmm_ss\/MiniAOD\/MiniAOD_$i.root/g" miniAODeemm/python/miniAODmuonsRootupler_1.py
    sed -i "s/preselection\/zeejmm_mc_2018_[0-9]\+\.root/preselection\/zeejmm_mc_2018_$i.root/g" miniAODeemm/python/miniAODmuonsRootupler_1.py
    cmsRun miniAODeemm/python/miniAODmuonsRootupler_1.py
done

# combine the output files
echo "Combining the output files..."
hadd preselection/zeejmm_mc_2018.root preselection/zeejmm_mc_2018_*.root

echo "All done!"

# # Clean up the output directory
# echo "Cleaning up the output directory..."
# rm ../output/ZmmMCOut_*.root
