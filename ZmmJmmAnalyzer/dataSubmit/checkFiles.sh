#!/bin/bash

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <year> <line_number>"
    exit 1
fi

year=$1
line_number=$2
filename_pattern="crabConfig_TTree.py"
# filename_pattern="miniAODmuonsRootupler_*.py"

echo "Searching for line $line_number in files matching pattern '$filename_pattern' for year $year..."

# Loop over files matching the pattern
for file in "$year"/*/$filename_pattern; do
    # Check if the file exists
    if [ -f "$file" ]; then
        echo "$file"
        # Print the specified line number from the file
        sed -n "${line_number}p" "$file"
        echo
    fi
done
