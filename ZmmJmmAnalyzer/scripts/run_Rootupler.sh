#!/bin/bash
# Run this from the src directory

# Set up the environment
cmsenv

# Run the rootupler for each of the 10 input files
echo "Starting the loop to process files..."
for i in {1..10}; do
    echo "Processing file with index $i..."
    sed -i "s/datasets\/mc\/run2_zmmJpsimm\/ZmmJpsimm_2018\/MiniAOD\/MiniAOD_[0-9]\+\.root/datasets\/mc\/run2_zmmJpsimm\/ZmmJpsimm_2018\/MiniAOD\/MiniAOD_$i.root/g" miniAODmmmm/python/miniAODmuonsRootupler_1.py
    sed -i "s/zmmjmm_mc_2018\/zmmjmm_mc_2018_[0-9]\+\.root/zmmjmm_mc_2018\/zmmjmm_mc_2018_$i.root/g" miniAODmmmm/python/miniAODmuonsRootupler_1.py
    cmsRun miniAODmmmm/python/miniAODmuonsRootupler_1.py
done

# combine the output files
echo "Combining the output files..."
hadd inputFiles/zmmjmm_mc_2018/zmmjmm_mc_2018.root inputFiles/zmmjmm_mc_2018/zmmjmm_mc_*.root

echo "All done!"

# # Clean up the output directory
# echo "Cleaning up the output directory..."
# rm ../output/ZmmMCOut_*.root