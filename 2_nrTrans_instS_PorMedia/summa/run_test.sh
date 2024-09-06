#!/bin/bash

# Get the directory of the script
script_dir=$(dirname "$(realpath "$0")")

# Expand the full path of the openWQ_master.json file
file_path="$script_dir/openWQ_master.json"

# Export the full path to the environment variable
export master_json="$file_path"

# Run the test
time /code/summa/bin/summa_openwq_sundials.exe -g 1 1 -m ./summa/SUMMA/summa_fileManager_OpenWQ_systheticTests_BGQ.txt

