#!/bin/bash

# ==============================================================================
# Driver Script for Sampling MERRA-21C Time-Invariant Constants
# along flight tracks in HSRL2/DIAL/HALO files
# ==============================================================================

# 1. Environment Setup
source /discover/nobackup/acollow/AeroApps/AeroApps/install/bin/g5_modules.sh

export PYTHONPATH="/discover/nobackup/acollow/GMAOpyobs/install/lib/Python:/discover/nobackup/acollow/GMAOpyobs/src:${PYTHONPATH}"

ulimit -s unlimited
export OMP_STACKSIZE=4G

# 2. Parse YAML configuration
CONFIG_FILE=$1

if [ -z "$CONFIG_FILE" ]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file $CONFIG_FILE not found!"
    exit 1
fi

# Only parsing the parameters needed for the constants script
CAMPAIGN=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['campaign'])")
OUTPUT_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['output_dir'])")
FLIGHT_DATES=$(python3 -c "import yaml; print(' '.join(yaml.safe_load(open('$CONFIG_FILE'))['flight_dates']))")

# Ensure output directory exists
mkdir -p "${OUTPUT_DIR}"

echo "=========================================================="
echo "Starting Constants Sampling for ${CAMPAIGN}..."
echo "Strategy: Flight-by-Flight processing"
echo "=========================================================="

for base_date in $FLIGHT_DATES; do
    
    echo ""
    echo "**********************************************************"
    echo "Processing Date: ${base_date}"
    echo "**********************************************************"

    # ==============================================================================
    # Sample MERRA-21C Constants
    # ==============================================================================
    echo "--- Sampling MERRA-21C Time-Invariant Constants ---"
    
    python3 samplecontantsfile.py "$CONFIG_FILE" "$base_date"

done # End of FLIGHT_DATES loop

echo "=========================================================="
echo "All constants processing complete!"
