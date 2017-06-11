# SpectralAnalysis


This program will calculate absorption coefficients via the Beer-Lambert law from transmission spectra.


### Installation
The first step will be to download this code onto your computer. You can do this by pasting the following into your terminal:

`git clone https://github.com/NauIceLab/SpectralAnalysis.git`


This will provide all of the necessary files to perform a test case to understand how the program operates step by step. 


### Running the code
Included are two CSV files with example data. BackgroundData.csv and RawData.csv are spectra collected of an empty cell chamber and pure methane (CH4) liquid. These spectra were taken at 92 Kelvin. At the very end of this section is a full gif of the entire process of executing this program for reference.

1. First, you will run the code via the terminal with: `python GTanalysis.py`. This will start the program. The program will ask you for a directory name, this will be where all of the work will be stored. Various plots, csv files, and txt files will be generated during the program and saved into this directory. 

2. Next will be to import the data used for analysis. The program will first ask for your raw dataset, followed by your background spectrum. These data files do not need to be in the same directory as `GTanalysis.py`, but make sure to give a full path name if this is the case. Included in the download are `RawData.csv` and `RawBackground.csv`, and can be used to test the program. See the gif at the end of this section where these files are used.

3. After importing these data files, the program will save the temperature and concentration of the mixture you are analyzing. The program will ask for both, and will store each answer in a separate txt file within the directory.

4. Once the data have been imported, the program will plot both data sets and ask which filtering method the user desires. The options are utilizing a ![Savitzky-Golay filter](https://en.wikipedia.org/wiki/Savitzky–Golay_filter) (essentially a low-pass filter) or a Fast Fourier Transformation method (FFT) to manually remove unwanted frequencies. This filter step is used to remove any signal that might come form experimental methods. In the case of my project, there is a signal from relfection of sapphire windows in our laboratory setup. See ![Protopapa 2015](https://arxiv.org/pdf/1503.00703.pdf) and ![Grundy 2002](http://www.sciencedirect.com/science/article/pii/S0019103501967260) for more on this topic. I would recommend utilizing the FFT method, as this allows the user more control over the data being processed.

5. Once the user has determined a filtering method, the program will ask for further input for the user's choice. If the Savtizky-Golay filter was chosen, the user definse a window box size and a polynomial for the alogorithm. See ![the wikipedia article](https://en.wikipedia.org/wiki/Savitzky–Golay_filter) for more understanding of how it works. If the user chose the FFt option, a FFT is taken of the data, and the user hand selects the frequency to be cut out (see gif at bottom for example).

<!---
![first gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part1.gif) 
![second gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part2.gif) 
![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part3.gif)
-->
