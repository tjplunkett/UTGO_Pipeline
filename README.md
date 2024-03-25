# UTGO_Pipeline
## Scripts to aid in astronomical image manipulation, reduction and photometry for UTGO

#### Installation:
This code is mostly written in Python, leaning very heavily on the modular pipeline building package ‘Prose’ (Garcia et al. 2021). You will also need to install SExtractor, SWarp and Astrometry.net if you want full functionality. To use the scripts, you will need to install the required dependencies. I recommend the user to utilise Anaconda or Miniconda in order to manage python environments. To start, create a new environment (let’s call it ‘red’) by going to the terminal (or Anaconda Prompt on Windows) and typing:

conda create -n red python=3.9.16

Once activated, (red) should display on the far left of the terminal window. Download the project folder from github and then change into this directory, using:

cd [path2folder] 

Then, assuming PIP was installed successfully to this environment, type:

pip install -r requirements.txt 

For more assistance, please see the UTGOPipeline_Documentation file. 

Note: If you would like to run this code on your own system, you will need to add paths to the 'directories.csv'. This file tells the code where the calibration frames, data and output folders live. This is mainly implemented for security reasons.