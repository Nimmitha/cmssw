#!/bin/bash

# Function to print usage
usage() {
    echo "Usage: $0 <year> <action>"
    echo "Actions:"
    echo "  submit   - Run 'crab submit -c crabConfig_TTree.py'"
    echo "  status   - Run 'crab status -d crab_TTree_13TeV*'"
    echo "  resubmit - Run 'crab resubmit -d crab_TTree_13TeV*'"
    echo "  kill     - Run 'crab kill -d crab_TTree_13TeV*'"
    echo "  report   - Run 'crab report -d crab_TTree_13TeV*'"
    echo "  lumi     - Run 'brilcalc lumi -i processedLumis.json > onlinelumi.csv'"
    exit 1
}

# Default values for variables
CRAB_DIR_PATTERN="crab_TTree_13TeV*"
CRAB_CONFIG="crabConfig_TTree.py"

# Check if correct number of arguments is provided
if [ $# -ne 2 ]; then
    usage
fi

year=$1
action=$2

# Loop through each subdirectory in the specified directory
for dir in "$year"/*/; do
    echo "Processing directory: $dir"
    
    # Check if the directory exists and is indeed a directory
    if [ ! -d "$dir" ]; then
        echo "Directory $dir does not exist. Skipping..."
        continue
    fi

    # Execute the action based on input parameter
    case $action in
        submit)
            (cd "$dir" && crab submit -c $CRAB_CONFIG)
            ;;
        status)
            (cd "$dir" && crab status -d $CRAB_DIR_PATTERN)
            ;;
        resubmit)
            (cd "$dir" && crab resubmit -d $CRAB_DIR_PATTERN)
            ;;
        report)
            (cd "$dir" && crab report -d $CRAB_DIR_PATTERN)
            ;;
        kill)
            (cd "$dir" && crab kill -d $CRAB_DIR_PATTERN)
            ;;
        lumi)
            found_subdir=false
            for subdir in "$dir"/$CRAB_DIR_PATTERN/results; do
                # Check if the subdirectory exists
                if [ -d "$subdir" ]; then
                    found_subdir=true
                    echo "Found subdirectory: $subdir"
                    # Check if processedLumis.json exists
                    if [ -f "$subdir/processedLumis.json" ]; then
                        (cd "$subdir" && brilcalc lumi -i processedLumis.json > onlinelumi.csv)
                    else
                        echo "processedLumis.json not found in $subdir. Skipping..."
                    fi
                fi
            done

            if [ "$found_subdir" = false ]; then
                echo "No matching subdirectory found in $dir"
            fi
            ;;
        *)
            echo "Invalid action: $action"
            usage
            ;;
    esac
done
