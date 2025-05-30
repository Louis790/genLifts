#!/bin/bash

# Check if exactly 7 arguments are passed
if [ "$#" -ne 7 ]; then
    echo "Usage: <partition> <modulo> <n> <k> <groupMin> <groupMax> <minGirth>"
    echo "Press Enter to exit..."
    read
    exit 1
fi

# Assign arguments to variables for clarity
partition=$1
modulo=$2
n=$3
k=$4
groupMin=$5
groupMax=$6
minGirth=$7

# Construct pipe paths based on parameters
pipe_prefix="/tmp/pipe_${partition}"
dreadnaut_in="${pipe_prefix}_in"
dreadnaut_out="${pipe_prefix}_out"

# Clean up any existing pipes
rm -f "${dreadnaut_in}" "${dreadnaut_out}"

# Create new named pipes
mkfifo "${dreadnaut_in}"
mkfifo "${dreadnaut_out}"

# Start dreadnaut in the background
cat "${dreadnaut_in}" | ./bin/dreadnaut > "${dreadnaut_out}" &
disown

# Run the main generator
./bin/genLifts "$partition" "$modulo" "$n" "$k" "$groupMin" "$groupMax" "$minGirth"


