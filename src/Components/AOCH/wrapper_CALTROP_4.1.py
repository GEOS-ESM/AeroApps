"""
1/1/2025
wrapper_CALTROP_4.1.py: this is a wrapper created with chat.gsfc that loops over a bunch of dates
and executes CALTROP_overlap_4.2.py for each of them

"""
import subprocess
import sys

# List of cases with their parameters
cases = [
   {'yyyy': '2020', 'mm': '08', 'dd': '21', 'cal_orb_time': [19,58], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '08', 'dd': '26', 'cal_orb_time': [13,16], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '07', 'cal_orb_time': [20,49], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '12', 'cal_orb_time': [20,40], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '14', 'cal_orb_time': [18,39], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '14', 'cal_orb_time': [20,17], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '15', 'cal_orb_time': [17,37], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '15', 'cal_orb_time': [19,17], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '16', 'cal_orb_time': [16,37], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '16', 'cal_orb_time': [18,15], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '16', 'cal_orb_time': [19,55], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '17', 'cal_orb_time': [20,33], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '18', 'cal_orb_time': [19,32], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '20', 'cal_orb_time': [19,9], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '09', 'dd': '27', 'cal_orb_time': [15,20], 'out_label': 'batch'},
   {'yyyy': '2020', 'mm': '10', 'dd': '01', 'cal_orb_time': [19,30], 'out_label': 'batch'}
]
# cases = [
#    {'yyyy': '2020', 'mm': '09', 'dd': '15', 'cal_orb_time': [17,39], 'out_label': 'batch'},
#    {'yyyy': '2020', 'mm': '10', 'dd': '01', 'cal_orb_time': [19,30], 'out_label': 'batch'}
#    
# ]
# cases = [
#    {'yyyy': '2020', 'mm': '08', 'dd': '21', 'cal_orb_time': [19,58], 'out_label': 'batch'},
#    {'yyyy': '2020', 'mm': '08', 'dd': '26', 'cal_orb_time': [13,16], 'out_label': 'batch'}
#    ]
# cases = [
#    {'yyyy': '2020', 'mm': '10', 'dd': '01', 'cal_orb_time': [19,30], 'out_label': 'test'},
#  ]

def modify_and_run_script(case):
   # Read the original script
   with open('CALTROP_overlap_4.25.py', 'r') as file:
      lines = file.readlines()

   # Find and modify the parameter line
   for i, line in enumerate(lines):
      if 'yyyy, mm, dd =' in line:
         # Preserve the original indentation
         indentation = len(line) - len(line.lstrip())
         new_line = ' ' * indentation  # Recreate the original indentation
         new_line += f"yyyy, mm, dd = '{case['yyyy']}', '{case['mm']}', '{case['dd']}' ; cal_orb_time=[{case['cal_orb_time'][0]},{case['cal_orb_time'][1]}] ;out_label='{case['out_label']}' #\n"
         lines[i] = new_line
#          print(lines[i])
         break

   # Write the modified script to a temporary file
   temp_script_name = 'temp_CALTROP.py'
   with open(temp_script_name, 'w') as file:
      file.writelines(lines)

   # Execute the modified script
   print(f"\nProcessing case: {case['yyyy']}-{case['mm']}-{case['dd']} "
         f"Time: {case['cal_orb_time'][0]}:{case['cal_orb_time'][1]}")
   try:
      subprocess.run([sys.executable, temp_script_name], check=True)
   except subprocess.CalledProcessError as e:
      print(f"Error running case: {e}")

   # Clean up temporary file
   import os
   os.remove(temp_script_name)

# Run all cases
for case in cases:
   modify_and_run_script(case)
   print("--------------------------------------------------------")