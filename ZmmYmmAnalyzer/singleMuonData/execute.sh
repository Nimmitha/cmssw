#!/bin/bash

# Check if directory argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <year>"
    exit 1
fi

year=$1

# Loop through each subdirectory in the specified directory
for dir in "$year"/*/; do
    echo "Processing directory: $dir"
    # Change into the directory and execute the command
    # (cd "$dir" && crab submit -c crabConfig_TTree.py)
    # (cd "$dir" && crab status -d crab_TTree_13TeV*)
    # (cd "$dir" && crab resubmit -d crab_TTree_13TeV*)
    (cd "$dir" && crab report -d crab_TTree_13TeV*)
done