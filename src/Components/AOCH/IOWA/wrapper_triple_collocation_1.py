"""
1/2/2025
wrapper_triple_collocation.py: Wrapper to process multiple IOWA-CALIOP files through
collocated_MERRA-CALIOP-IOWA_3.2.py
"""
import subprocess
import sys
import os
import glob


def get_matching_files(input_directory, file_pattern):
   """Get list of files matching the pattern in the specified directory."""
   full_pattern = os.path.join(input_directory, file_pattern)
   file_list = glob.glob(full_pattern)
   return sorted(file_list)

def modify_and_run_script(filename):
   # Read the original script
   with open('collocated_MERRA-CALIOP-IOWA_3.2.py', 'r') as file:
      lines = file.readlines()

   # Find and modify the filename line
   for i, line in enumerate(lines):
      if 'CAL_IOWA_File =' in line:
         # Preserve the original indentation
         indentation = len(line) - len(line.lstrip())
         new_line = ' ' * indentation
         new_line += f"CAL_IOWA_File = '{os.path.basename(filename)}' #\n"
         lines[i] = new_line
         break

   # Write the modified script to a temporary file
   temp_script_name = 'temp_collocation.py'
   with open(temp_script_name, 'w') as file:
      file.writelines(lines)

   # Execute the modified script
   print(f"\nProcessing file: {os.path.basename(filename)}")
   try:
      subprocess.run([sys.executable, temp_script_name], check=True)
   except subprocess.CalledProcessError as e:
      print(f"Error processing file: {e}")

   # Clean up temporary file
   os.remove(temp_script_name)


###### -------- Main Code
if __name__ == "__main__":
   # Set input directory and file pattern
   INPUT_DIR = "/discover/nobackup/sgasso/Files/Satellite/SliceTROCAL/"
   FILE_PATTERN = "*_batch.nc"
#    FILE_PATTERN ="IOWA-CALIOP_2020-10-01T19-30-09ZD_test.nc"

   # Get list of all matching files
   all_files = get_matching_files(INPUT_DIR, FILE_PATTERN)
   
   if not all_files:
      print(f"No files found matching pattern: {FILE_PATTERN}")
      sys.exit()
   
   # Print all available files
   print(f"Found {len(all_files)} matching files:")
   for i, f in enumerate(all_files):
      print(f"{i}: {os.path.basename(f)}")
   
   # Option to process all files or specific indices
#    process_all = False  # Set to True to process all files
   process_all = True  # Set to True to process all files
   
   if process_all:
      files_to_process = all_files
   else:
      indices = [0,1]  # Specify which files to process
      try:
         files_to_process = [all_files[i] for i in indices]
      except (ValueError, IndexError) as e:
         print(f"Error in file selection: {e}")
         sys.exit()
   
   # Process selected files
   for filename in files_to_process:
      modify_and_run_script(filename)
      print("--------------------------------------------------------")
