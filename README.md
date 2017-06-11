# SpectralAnalysis


This program will calculate absorption coefficients via the Beer-Lambert law from transmission spectra.


### Installation
The first step will be to download this code onto your computer. You can do this by pasting the following into your terminal:

`git clone https://github.com/NauIceLab/SpectralAnalysis.git`


This will provide all of the necessary files to perform a test case to understand how the program operates step by step. 


### Running the code
Included are two CSV files with example data. BackgroundData.csv and RawData.csv are spectra collected of an empty cell chamber and pure methane (CH4) liquid. These spectra were taken at 92 Kelvin. First, you will run the code via the terminal with: `python GTanalysis.py`. This will start the program.The program will ask you for a directory name, this will be where all of the work will be stored. Various plots, csv files, and txt files will be generatred during the program and saved into this directory.

![first gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part1.gif) 


Next will be to import the data used for analysis. The program will first ask for your raw dataset, followed by your background spectrum. These data files do not need to be in the same directory as `GTanalysis.py`, but make sure to give a full path name if this is the case. 

![second gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part2.gif) 

After importing these data files, the program will save the temperature and concenntration of the mixture you are analyzing. The program will ask for both, and will store each answer in a separate txt file within the directory.

![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part3.gif)
