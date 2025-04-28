#
#
# Modules use by code test_get_mplAOCH_at_overpass.py
#
# created Apr28/2025 by @sgassoumd
#------------------------------



import os, sys 
from datetime import datetime

class InputData:
    def __init__(self):
        self.date_time = []
        self.orbit_number = []
        self.line = []
        self.column = []
        self.me_hw532 = []
    def __len__(self):
        return len(self.date_time)
    def __getitem__(self, key):
        #Allow dictionary-like access with records['column_name']"""
        if hasattr(self, key):
            return getattr(self, key)
        else:
            raise KeyError(f"Column '{key}' not found")

def f_get_records(filename):
    try:
        # Read the data
        records = read_file(filename)
        
        # Display results
        print(f"Successfully read {len(records)} records from {filename}")
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"Error processing file: {e}")
    return records 
        
def read_file(filename):
    # Read a MERRA data file into a list of InputData objects.
    # IN :filename: Path to the MERRA data file
    # OUT: List of InputData objects

    dataset = InputData()
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split()
            
            dataset.date_time.append(datetime.fromisoformat(parts[0]))
            dataset.orbit_number.append(int(parts[1]))
            dataset.line.append(int(parts[2]))
            dataset.column.append(int(parts[3]))
            dataset.me_hw532.append(float(parts[4]))
    
    return dataset
