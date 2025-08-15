#!/bin/bash

if [ "$#" -ne 8 ] && [ "$#" -ne 9 ]; then
    echo "Usage: <partition> <modulo> <n> <k> <groupMin> <groupMax> <minGirth> <pureGraph6> <Optional:groupfile>"
    echo "Press Enter to exit..."
    read
    exit 1
fi

partition=$1
modulo=$2
n=$3
k=$4
groupMin=$5
groupMax=$6
minGirth=$7
pureGraph6=$8
groupfile=${9:-""}

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
if [ -n "$groupfile" ]; then
    ./bin/genLifts "$partition" "$modulo" "$n" "$k" "$groupMin" "$groupMax" "$minGirth" "$pureGraph6" "$groupfile"
else
    ./bin/genLifts "$partition" "$modulo" "$n" "$k" "$groupMin" "$groupMax" "$minGirth" "$pureGraph6"
fi



