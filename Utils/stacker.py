"""
stacker.py

Author: Thomas Plunkett

Purpose: Stack a selection of images together, using either SWARP or a weighted pixel combination map.

"""
# Import necessary packages 
from prose import Image, Sequence, blocks, FitsManager, Observation
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import subprocess as sub
import os
import glob

# Get path for SWarp
cwd = os.getcwd()
swarp_file = os.path.join(cwd, 'default.swarp')

def splitter(fm, freq):
    """
    A function to split a set of images by no. of images, day or hour ready to pass to a stacker.
    
    params:
    fm - A prose FitsManager object 
    freq - The frequency to group images by
    """
    df_list = []
    
    # Access all files and extract JD and paths. Create datetime.
    total_df = fm.to_pandas(query='select path, jd, exposure from files')
    total_df['dt'] = pd.to_datetime(total_df.jd, unit='D', origin='julian')
    total_df = total_df.sort_values('jd').reset_index()
    
    # Only stack images of like exp. time
    for df in [total_df[total_df.exposure == e] for e in total_df.exposure.unique()]:
        df = df.reset_index()
        
        # Depending on the input to 'freq' param, create lists of data frames 
        if freq == 'day' or freq == 'Day' or freq == 'd' or freq == 'D':
            df = df.set_index('dt')
            df_list += [df[df.index.day == d] for d in df.index.day.unique()]

        elif freq == 'hour' or freq == 'Hour' or freq == 'h' or freq == 'H':
            df = df.set_index('dt')
            days_df = [df[df.index.day == d] for d in df.index.day.unique()]
            for day in days_df:
                df_list += [day[day.index.hour == h] for h in day.index.hour.unique()]
        else:
            df_groups = df.groupby(df.index // int(freq))

            for _, dfs in df_groups:
                df_list += [dfs]
        
    return df_list
    
def writer(lst, title):
    """
    A function to write a list of files to an ascii file to pass to SWarp.
    
    params:
    lst - The list of files to write to ascii
    title - The title of the output ascii file
    """
    file_name = title
    with open(file_name, 'w') as f:
        for path in lst:
            f.write(f"{path}\n")
    f.close()
    
def remover(path2stacks):
    """
    Function to remove previous stacks from folder to avoid confusion!
    
    param:
    path2stacks - The output folder where stacks may have been run before
    """
    stack_list = glob.glob(os.path.join(path2stacks, '*_Stack_*.fits'))
    ascii_list = glob.glob(os.path.join(path2stacks, '*_Stack_*.ascii')) 
    
    if len(stack_list) != 0 or len(ascii_list) != 0:
        print('Removing previous stacks and lists!')
        for stack in stack_list:
            os.remove(stack)
        for asci in ascii_list:
            os.remove(asci)

# Set up the parser
parser = argparse.ArgumentParser(description='Reduce and stack images')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
parser.add_argument('method', type=str, help='Use SWarp? (type swarp)')
parser.add_argument('period', type=str, help='Stack by hour or day? (hour/day)')
args = parser.parse_args()

# Safety and checks
target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Begin the action!
remover(target_dir)
fm = FitsManager(target_dir, depth=args.depth)
print(fm)
ref = Image(fm.all_images[0])
objct = ref.header['OBJECT']

# Define Calibration sequence from reduced
CalibFromRed = Sequence([
    blocks.Trim(),
    blocks.SegmentedPeaks(), # stars detection
    blocks.Cutouts(),                   # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
])
   
CalibFromRed.run(ref, show_progress=False)

default_stacker = Sequence([
    *CalibFromRed[0:-1],                   # apply the same calibration to all images
    blocks.psf.Moffat2D(reference=ref),   # providing a reference improves the PSF optimisation
    blocks.Twirl(ref.stars_coords),       # compute image transformation
    # Stack image
    blocks.AffineTransform(stars=False, data=True),
    blocks.MedianStack(name='stack'),
])
      
# Start action    
df_list = splitter(fm, str(args.period)) 
exp_list = []

# Try use Swarp
if args.method == 'swarp':
    try:
        for df in df_list:
            # Get exposure time and store for later
            exp = df.exposure.unique()[0]
            counter = exp_list.count(exp) + 1
            
            # Define paths to pass to swarp, write these to ascii files for posterity
            path_list = df.path.to_list()
            list_file = str(objct+'_Stack_'+str(exp)+'s_'+str(counter)+'.ascii')
            list_file = os.path.join(target_dir, list_file)
            writer(path_list, list_file)
            
            # Run swarp
            command = 'SWarp @' + list_file + ' -c ' + swarp_file + ' -IMAGEOUT_NAME ' + list_file.replace('.ascii', '.fits') 
            process = sub.Popen([command], shell=True)
            process.wait() 
            
            exp_list += [exp]
            
    except:
        print('Unable to stack with SWarp! Check installation, or use default stack method')
        raise SystemExit(1)

# Use default stacker from Prose
else:
    for df in df_list:
        # Get exposure time and store for later
        exp = df.exposure.unique()[0]
        counter = exp_list.count(exp) + 1
        
        # Define paths to pass to swarp, write these to header for posterity
        output_name = str(objct+'_Stack_'+str(exp)+'s_'+str(counter)+'.fits')
        output_path = os.path.join(target_dir, output_name)
        path_list = df.path.to_list()
        writer(path_list, output_path.replace('.fits', '.ascii'))
        default_stacker.run(path_list)
        default_stacker.stack.stack.header = Image(path_list[0]).header
        default_stacker.stack.stack.header['COMMENT'] = 'Created with stacker.py using default method.'
        default_stacker.stack.stack.writeto(output_path)
        
        exp_list += [exp]