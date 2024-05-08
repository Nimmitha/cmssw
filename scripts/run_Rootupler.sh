#!/bin/bash
# Run this from the src directory

# Set up the environment
cmsenv

# Run the rootupler for each of the 10 input files
echo "Starting the loop to process files..."
for i in {1..10}; do
    echo "Processing file with index $i..."
    sed -i "s/datasets\/ZmmYee\/MiniAOD\/MiniAOD_[0-9]\+\.root/datasets\/ZmmYee\/MiniAOD\/MiniAOD_$i.root/g" miniAODmuonsRootupler.py
    sed -i "s/output\/ZmmMCOut_[0-9]\+\.root/output\/ZmmMCOut_$i.root/g" miniAODmuonsRootupler.py
    cmsRun miniAODmuonsRootupler.py
done

# combine the output files
echo "Combining the output files..."
hadd output/ZmmMCOut.root output/ZmmMCOut_*.root

echo "All done!"

# # Clean up the output directory
# echo "Cleaning up the output directory..."
# rm ../output/ZmmMCOut_*.root