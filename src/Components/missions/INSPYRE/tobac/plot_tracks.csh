#!/bin/csh

module load python/GEOSpyD/Min24.4.0-0_py3.11

# Define your directories
set input_dir = "./Save"
set output_dir = "./Plot"

mkdir -p "$output_dir"

echo "==================================================="
echo " Starting Processing for Track files"
echo " Reading from: $input_dir"
echo " Saving to:    $output_dir"
echo "==================================================="

foreach file ( "$input_dir"/Track_*.nc )

    # Safety check in case the folder is empty
    if ( ! -e "$file" ) then
        echo "Error: No Track_*.nc files found in $input_dir."
        exit 1
    endif

    # 3. Extract the event name from the filename
    set base = `basename "$file" .nc`
    set event_name = `echo "$base" | sed 's/^Track_//'`

    echo ""
    echo "---------------------------------------------------"
    echo " Processing Event: $event_name"
    echo "---------------------------------------------------"

    python -u plot_tracks_refined.py --event "$event_name"
    set py_status = $status
    # Check the specific exit status from Python
    if ( $py_status == 2 ) then
        echo " -> [Skipped] PyroCb cloud did not meet duration/distance criteria."
        continue
    else if ( $py_status != 0 ) then
        echo " [!] CRITICAL ERROR: Python script crashed on $event_name."
        continue
    endif
end

echo ""
echo "==================================================="
echo " Batch processing complete! All results are in $output_dir"
echo "==================================================="
