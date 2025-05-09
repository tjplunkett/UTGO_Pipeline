"""
reductionBot.py

Author: Thomas Plunkett

Date: 5th Jan, 2024

Purpose:

A script that reads a text file of reduced folders and compares to current folder.
If a new folder is present, passes this to autoReduce.py for automatic reduction.

"""
#!/bin/bash

# Import necessary packages 
import pandas as pd
import os
import glob
import subprocess
from datetime import datetime

# Find current year and define important paths
dir_path = os.path.abspath('/home/obs/UTGO_Pipeline/directories.csv')
dir_df = pd.read_csv(dir_path)
year_str = str(datetime.now().year)
path2data = os.path.abspath(dir_df.data_dir[0])
path2yr = os.path.join(path2data,year_str)
path2txt = os.path.join(path2yr, 'ReducedFolders.txt')

# Read text file. Might fail if the file is empty/doesnt exist
try:
    red_df = pd.read_table(path2txt)
    red_df = red_df.astype(str)
    old_list = red_df['Folders'].values.tolist()
except:
    old_list = []
    
# Check the directory for new folders
dir_list = []
for item in os.listdir(path2yr):
    if os.path.isdir(os.path.join(path2yr, item)):
        dir_list += [item]

# Compare if entries are in the old list or not
dir_list = sorted(dir_list, key=lambda x: datetime.strptime(x, "%Y%m%d"))

for dirt in dir_list:
    if dirt not in old_list:
        print('Working on night: {:}'.format(dirt))
        
        # Re-name files here, so shit doesn't break
        file_list = glob.glob(os.path.join(os.path.join(path2yr, dirt), '*.fits'))
        for file in file_list:
            if ' ' in str(file):
                os.rename(file, str(file).replace(' ','_'))

        target_dir = str(os.path.join(path2yr, dirt))
        calib_dataset = 'Calibrations_'+dirt
        raw_dataset = 'Observations_'+dirt

        # Upload raw and calibration data
        cmd1 = 'python /home/obs/UTGO_Pipeline/uploader.py ' + target_dir + ' raw ' + raw_dataset + ' y y'
        process = subprocess.Popen([cmd1], shell=True)
        process.wait()

        cmd2 = 'python /home/obs/UTGO_Pipeline/uploader.py ' + target_dir + ' calib ' + calib_dataset + ' y y'
        process = subprocess.Popen([cmd2], shell=True)
        process.wait()
        
        # Extract any calibration frames 
        cmd3 = 'python /home/obs/UTGO_Pipeline/update_calibs.py ' + str(os.path.join(path2yr, dirt)) + ' n'
        process = subprocess.Popen([cmd3], shell=True)
        process.wait()
        
        # Auto reduction
        cmd4 = 'python /home/obs/UTGO_Pipeline/autoReduce.py ' + str(os.path.join(path2yr, dirt)) + ' n n'
        process = subprocess.Popen([cmd4], shell=True)
        process.wait()
        
# Write new folders into the text file 
with open(path2txt, 'w') as fp:
    fp.write("Folders\n")
    for fldr in dir_list:
        # write each item on a new line
        fp.write("%s\n" % fldr)

print('Reduction Bot ran successfully on: {:}'.format(datetime.now()))
