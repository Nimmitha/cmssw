#!/bin/bash
# Run this from the src/x directory

# Set up the environment
cmsenv

# Run the rootupler for each of the 10 input files
echo "Starting the loop to process files..."
for i in {1..10}; do
    echo "Processing file with index $i..."
    # sed -i "s/datasets\/mc\/zmmymm_run3\/mmmm_v1\/MiniAOD\/MiniAOD_[0-9]\+\.root/datasets\/mc\/zmmymm_run3\/mmmm_v1\/MiniAOD\/MiniAOD_$i.root/g" miniAODmmmm/python/miniAODmuonsRootupler_1.py
    # sed -i "s/mc_mmmm_v1_[0-9]\+\.root/mc_mmmm_v1_$i.root/g" miniAODmmmm/python/miniAODmuonsRootupler_1.py
    
    sed -i "s/datasets\/mc\/run3_zmmJpsimm\/ZmmJpsimm\/MiniAOD\/MiniAOD_[0-9]\+\.root/datasets\/mc\/run3_zmmJpsimm\/ZmmJpsimm\/MiniAOD\/MiniAOD_$i.root/g" miniAODmmmm/python/miniAODmuonsRootupler_1.py
    sed -i "s/mc_zmmjpsimm_v1_[0-9]\+\.root/mc_zmmjpsimm_v1_$i.root/g" miniAODmmmm/python/miniAODmuonsRootupler_1.py
    
    cmsRun miniAODmmmm/python/miniAODmuonsRootupler_1.py
done

# combine the output files
echo "Combining the output files..."
hadd inputFiles/mc_zmmjpsimm_v1/mc_zmmjpsimm_v1.root inputFiles/mc_zmmjpsimm_v1/mc_zmmjpsimm_v1_*.root

echo "All done!"

# Clean up the output directory
# echo "Cleaning up the output directory..."
# rm ../output/ZmmMCOut_*.root