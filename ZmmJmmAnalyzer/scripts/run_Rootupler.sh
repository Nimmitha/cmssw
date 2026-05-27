#!/bin/bash
# Run this from the src directory

cmsenv

CFG="miniAODmmmm/python/miniAODmuonsRootupler_mc.py"

INPUT_DIR="/uscms/home/wkarunar/nobackup/datasets/mc/miniAOD/run2/zmmjmm/2018_ss/MiniAOD"
OUTPUT_DIR="selection/2018_ss"

mkdir -p "$OUTPUT_DIR"

echo "Starting the loop to process files..."

for i in {1..10}; do
    INPUT_FILE="file:${INPUT_DIR}/MiniAOD_${i}.root"
    OUTPUT_FILE="${OUTPUT_DIR}/zmmjmm_mc_2018_v1_${i}.root"

    echo "Processing file with index $i..."
    echo "Input : $INPUT_FILE"
    echo "Output: $OUTPUT_FILE"

    cmsRun "$CFG" \
        inputFiles="$INPUT_FILE" \
        outputFile="$OUTPUT_FILE"

    if [ $? -ne 0 ]; then
        echo "cmsRun failed for file index $i"
        exit 1
    fi
done

echo "All done!"