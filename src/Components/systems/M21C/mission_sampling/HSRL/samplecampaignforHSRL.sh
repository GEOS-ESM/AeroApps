#!/bin/bash

# ==============================================================================
# Driver Script for Gridding, Sampling, and AOP of MERRA-2 / MERRA-21C
# Files at 60 second intervals according to flight tracks in HSRL2/DIAL/HALO files
# ==============================================================================

# 1. Environment Setup
source /discover/nobackup/acollow/AeroApps/AeroApps/install/bin/g5_modules.sh

export PYTHONPATH="/discover/nobackup/acollow/GMAOpyobs/install/lib/Python:/discover/nobackup/acollow/GMAOpyobs/src:${PYTHONPATH}"

ulimit -s unlimited
export OMP_STACKSIZE=4G

REGRID_EXEC="/discover/nobackup/acollow/AeroApps/AeroApps/install/bin/Regrid_Util.x"
AOP_EXEC="/discover/nobackup/acollow/GMAOpyobs/install/bin/aop"
GRID_SPEC="PC1152x721-DC"
FVINPUT_DIR="/discover/nobackup/projects/gmao/share/dasilva/fvInput/"

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

CAMPAIGN=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['campaign'])")
OUTPUT_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['output_dir'])")
TMP_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['tmp_dir'])")
EXP_ID=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['exp_id'])")
MERRA21C_ARCHIVE=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['merra21c_archive'])")
FLIGHT_DATES=$(python3 -c "import yaml; print(' '.join(yaml.safe_load(open('$CONFIG_FILE'))['flight_dates']))")

PLANE=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['plane'])")
WAVELENGTHS=$(python3 -c "import yaml; print(' '.join(map(str, yaml.safe_load(open('$CONFIG_FILE'))['wavelengths'])))")
AOP_CONF_M2=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['aop_config_m2'])")
AOP_CONF_M21C=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['aop_config_m21c'])")

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

echo "=========================================================="
echo "Starting Processing for ${CAMPAIGN}..."
echo "Strategy: Flight-by-Flight processing"
echo "=========================================================="

for base_date in $FLIGHT_DATES; do
    
    echo ""
    echo "**********************************************************"
    echo "Processing Date: ${base_date}"
    echo "**********************************************************"

    # ==============================================================================
    # STEP 1: Regrid MERRA-21C Files (For this specific date + 1 day)
    # ==============================================================================
    echo "--- STEP 1: Regridding MERRA-21C Files ---"
    
    for offset in 0 1; do
        target_date=$(date -d "${base_date} + ${offset} days" +%Y%m%d)
        YYYY=${target_date:0:4}
        MM=${target_date:4:2}
        DD=${target_date:6:2}
        
        for HH in "00" "03" "06" "09" "12" "15" "18" "21"; do
            input_file="${MERRA21C_ARCHIVE}/${EXP_ID}/chem/Y${YYYY}/M${MM}/${EXP_ID}.aer_inst_3hr_glo_C360x360x6_v72.${YYYY}-${MM}-${DD}T${HH}00Z.nc4"
            if [ ! -f "$input_file" ]; then continue; fi
            
            basename_file=$(basename "$input_file")
            regridded_file="${TMP_DIR}/regridded_${basename_file}"
            
            if [ -f "$regridded_file" ]; then 
                echo "  -> Regridded file already exists, moving on: $basename_file"
                continue
            fi

            echo "  -> Regridding: $basename_file"
            mpirun -np 6 ${REGRID_EXEC} \
                -i "${input_file}" \
                -o "${regridded_file}" \
                -ogrid "${GRID_SPEC}" \
                -file_weights 2>/dev/null || true
        done
    done
    
    # ==============================================================================
    # STEP 2: Sample MERRA-2 and MERRA-21C (For this specific date)
    # ==============================================================================
    echo "--- STEP 2: Sampling ---"
    
    python3 samplecampaign.py "$CONFIG_FILE" "$base_date"

    # ==============================================================================
    # STEP 3: Compute Aerosol Optical Properties (Extinction)
    # ==============================================================================
    echo "--- STEP 3: Computing Extinction (aop) ---"
    
    # Find all MERRA-2 sampled files for this date (handles L1, L2, etc. automatically)
    for m2_in in "${OUTPUT_DIR}/${CAMPAIGN}-MERRA2-aer-Nv-${PLANE}_Model_${base_date}"*.nc; do
        # If no files found, the glob string itself is returned, so we skip if it's not a real file
        [ -f "$m2_in" ] || continue 
        
        # Extract the suffix (everything after "Model_YYYYMMDD")
        prefix="${OUTPUT_DIR}/${CAMPAIGN}-MERRA2-aer-Nv-${PLANE}_Model_${base_date}"
        suffix="${m2_in#$prefix}"
        
        for wl in $WAVELENGTHS; do
            m2_out="${OUTPUT_DIR}/${CAMPAIGN}-MERRA2-ext${wl}-${PLANE}_Model_${base_date}${suffix}"
            
            if [ ! -f "$m2_out" ]; then
                echo "  -> Running aop for MERRA-2 | Suffix: ${suffix} | Wavelength: ${wl}"
                ${AOP_EXEC} -c ${AOP_CONF_M2} -r ${FVINPUT_DIR} --noaback -v -a ext -w ${wl} -o ${m2_out} ${m2_in}
            else
                echo "  -> aop output exists for MERRA-2 | Suffix: ${suffix} | WL: ${wl}"
            fi
        done
    done

    # Find all MERRA-21C sampled files for this date
    for m21c_in in "${OUTPUT_DIR}/${CAMPAIGN}-MERRA21C-aer-Nv-${PLANE}_Model_${base_date}"*.nc; do
        [ -f "$m21c_in" ] || continue 
        
        prefix="${OUTPUT_DIR}/${CAMPAIGN}-MERRA21C-aer-Nv-${PLANE}_Model_${base_date}"
        suffix="${m21c_in#$prefix}"
        
        for wl in $WAVELENGTHS; do
            m21c_out="${OUTPUT_DIR}/${CAMPAIGN}-MERRA21C-ext${wl}-${PLANE}_Model_${base_date}${suffix}"
            
            if [ ! -f "$m21c_out" ]; then
                echo "  -> Running aop for MERRA-21C | Suffix: ${suffix} | Wavelength: ${wl}"
                ${AOP_EXEC} -c ${AOP_CONF_M21C} -r ${FVINPUT_DIR} --noaback -v -a ext -w ${wl} -o ${m21c_out} ${m21c_in}
            else
                echo "  -> aop output exists for MERRA-21C | Suffix: ${suffix} | WL: ${wl}"
            fi
        done
    done

    # ==============================================================================
    # STEP 4: Clean up Temporary Files for this Date
    # ==============================================================================
    echo "--- STEP 4: Cleaning up temporary regridded files ---"
    
    for offset in 0 1; do
        target_date=$(date -d "${base_date} + ${offset} days" +%Y-%m-%d)
        #rm -f "${TMP_DIR}"/*"${target_date}"*.nc4
    done
    echo "  -> Cleanup for ${base_date} complete."

done # End of FLIGHT_DATES loop

echo "=========================================================="
echo "All processing complete!"
