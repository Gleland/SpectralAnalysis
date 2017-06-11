# Spectral Analysis


This program will calculate absorption coefficients via the Beer-Lambert law from transmission spectra.


### Installation
The first step will be to download this code onto your computer. You can do this by pasting the following into your terminal:

`git clone https://github.com/NauIceLab/SpectralAnalysis.git`


This will provide all of the necessary files to perform a test case to understand how the program operates step by step. 


### Running the code
Included are two CSV files with example data. BackgroundData.csv and RawData.csv are spectra collected of an empty cell chamber and pure methane (CH4) liquid. These spectra were taken at 92 Kelvin. At the very end of this section is a full gif of the entire process of executing this program for reference.

1. First, you will run the code via the terminal with: `python GTanalysis.py`. This will start the program. The program will ask you for a directory name, this will be where all of the work will be stored. Various plots, csv files, and txt files will be generated during the program and saved into this directory. 


![first gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part1.gif) 



2. Next will be to import the data used for analysis. The program will first ask for your raw dataset, followed by your background spectrum. These data files do not need to be in the same directory as `GTanalysis.py`, but make sure to give a full path name if this is the case. Included in the download are `RawData.csv` and `RawBackground.csv`, and can be used to test the program. See the gif at the end of this section where these files are used.


![second gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part2.gif) 


3. After importing these data files, the program will save the temperature and concentration of the mixture you are analyzing. The program will ask for both, and will store each answer in a separate txt file within the directory.

![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part3.gif)



4. Once the data have been imported, the program will plot both data sets and ask which filtering method the user desires. The options are utilizing a [Savitzky-Golay filter](https://en.wikipedia.org/wiki/Savitzky–Golay_filter) (essentially a low-pass filter) or a Fast Fourier Transformation method (FFT) to manually remove unwanted frequencies. This filter step is used to remove any signal that might come form experimental methods. In the case of my project, there is a signal from relfection of sapphire windows in our laboratory setup. See [Protopapa 2015](https://arxiv.org/pdf/1503.00703.pdf) and [Grundy 2002](http://www.sciencedirect.com/science/article/pii/S0019103501967260) for more on this topic. I would recommend utilizing the FFT method, as this allows the user more control over the data being processed.

![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/fft_choice.gif)



5. Once the user has determined a filtering method, the program will ask for further input for the user's choice. If the Savtizky-Golay filter was chosen, the user definse a window box size and a polynomial for the alogorithm. See [the wikipedia article](https://en.wikipedia.org/wiki/Savitzky–Golay_filter) for more understanding of how it works. If the user chose the FFT option, a FFT is taken of the data, and the user hand selects the frequency to be cut out. The program listens for mouse clicks on the plot that is displayed, a left click will add a vertical line, a middle click will remove the nearest line, and a right click will tell the program that the user is satisfied with the selection and wants to proceed.

![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/fft_cut.gif)


6. After the filtering process has been completed, the next step is to fit a continuum. This continuum defines where no absorption is occuring. The program will ask the user to zoom in on the region where the continuum wants to be fitted, and then will ask what order polynomial (linear, quadratic, etc) is to be used for creating the continuum. IN the gif below I utilized a third order (cubic) but any fit is available to the user. I used a cubic fit for my work as it gave the most consisent results for my project. While the user clicks on the graph, the program fits a polynomial to all the data points that have been added to the graph. A left click adds a point, a middle click remoes the nearest point, and a right click ends the collection.


![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/continuum.gif)
![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/polynomial_order.gif)



7. After the continuum has been fit, the program asks the user to input the thcikness of the cell used during spectra colelction. This is because the absorption of light depends on the thickness of the cell chamber, via the [Beer-Lambert law](http://life.nthu.edu.tw/~labcjw/BioPhyChem/Spectroscopy/beerslaw.htm). The cell chamber used in my experiemnts was 2 cm, and thus a value of `2` will be entered. Finally, the program will generate absorption coefficients and display them for the user to inspect.

8. When finished viewing the plot, the user enters `y`, and the program closes. All of the data and plots generated during the program will be saved into the direcotry that was input at the beginning.




## Demo of program, with a pure methane (Ch4) mixture at a temperature of 92.0 Kelvin.
![first gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/full_video.gif) 



# If you have any questions, feel free to contact me at garrett.leland.thompson@gmail.com

<!---
-->
