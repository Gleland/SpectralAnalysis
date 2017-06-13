# Spectral Analysis


This program will calculate absorption coefficients via the Beer-Lambert law from transmission spectra.


### Installation
The first step will be to download this code onto your computer. You can do this by pasting the following into your terminal:

`git clone https://github.com/NauIceLab/SpectralAnalysis.git`


This will include all of the necessary files to perform a test case to understand how the program operates step by step. 


### Running the code
Included are two CSV files with example data. RawBackground.csv and RawData.csv are spectra collected of an empty cell chamber and pure methane (CH4) liquid. These spectra were taken at 92 Kelvin. 

1. First, you will run the code via the terminal with: `python GTanalysis.py`. This will start the program. The program will ask you for a directory name, this will be where all of the work will be stored. Various plots, csv files, and txt files will be generated during the program and saved into this directory. 


![first gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part1.gif) 



2. Next, import the data used for analysis. The program will first ask for your raw dataset, followed by your background spectrum. These data files do not need to be in the same directory as `GTanalysis.py`, but make sure to give a full path name if this is the case. Included in the download are `RawData.csv` and `RawBackground.csv` and can be used to test the program.


![second gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part2.gif) 


3. After importing these data files, the program will save the temperature and concentration of the mixture you are analyzing. The program will ask for both, and will store each answer in a separate txt file within the directory.

![third gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/part3.gif)



4. Once the data have been imported, the program will plot both data sets for user inspection (to make sure the data look okay) and ask which filtering method the user desires. The options are utilizing a [Savitzky-Golay filter](https://en.wikipedia.org/wiki/Savitzky–Golay_filter) (essentially a low-pass filter) or a Fast Fourier Transformation method (FFT) to manually remove unwanted frequencies. This filter step is used to remove any signal that might come form experimental methods. In the case of my project, there is a signal from reflection of sapphire windows in our laboratory setup. See [Protopapa 2015](https://arxiv.org/pdf/1503.00703.pdf) and [Grundy 2002](http://www.sciencedirect.com/science/article/pii/S0019103501967260) for more on this topic. I would recommend utilizing the FFT method, as this allows the user more control over the data being processed.

![fourth gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/fft_choice.gif)



5. Once the user has determined a filtering method, the program will ask for further input for the user's choice. If the Savtizky-Golay filter was chosen, the user defines a window box size and a polynomial for the algorithm. See [the wikipedia article](https://en.wikipedia.org/wiki/Savitzky–Golay_filter) for more understanding of how it works. If the user chose the FFT option, a FFT is taken of the data, and the user hand selects the frequency to be cut out. The program listens for mouse clicks on the plot that is displayed, a left click will add a vertical line, a middle click will remove the nearest line, and a right click will tell the program that the user is satisfied with the selection and wants to proceed. In the gif provided one can see that the FFT is symmetric. The program will take the selection of frequencies (in the gif it is ~0.7 cm) and set them to zero, as well as their counterparts that are on the other side of the origin.

![fifth gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/fft_cut.gif)


6. After the filtering process has been completed, the next step is to fit a continuum. This continuum defines where no absorption is occurring. The program will ask the user to zoom in on the region where the continuum wants to be fitted, and then will ask what order polynomial (linear, quadratic, etc) is to be used for creating the continuum. In the gif below I utilized a third order (cubic) but any fit is available to the user. I used a cubic fit for my work as it gave the most consistent results for my project. While the user clicks on the graph, the program fits a polynomial to all the data points that have been added to the graph. A left click adds a point, a middle click removes the nearest point, and a right click ends the collection.


![sixth gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/continuum.gif)
![seventh gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/polynomial_order.gif)



7. After the continuum has been fit, the program asks the user to input the thickness of the cell used during spectra collection. This is because the absorption of light depends on the thickness of the cell chamber, via the [Beer-Lambert law](http://life.nthu.edu.tw/~labcjw/BioPhyChem/Spectroscopy/beerslaw.htm). The cell chamber used in my experiments was 2 cm, and thus a value of `2` will be entered. Finally, the program will generate absorption coefficients and display them for the user to inspect. Notice that not all of the spectrum can be relied upon for legitimate data, the end points (less than 90000 cm-1 and greater than 15000 cm-1) are not considered significant, so be sure not to include these values in calculations.

![eighth gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/thickness-and-finished.gif)
![ninth gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/abs_coeffs.gif)



8. When finished viewing the plot, the user enters `y`, and the program closes. All of the data and plots generated during the program will be saved into the directory that was input at the beginning.


![tenth gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/done.gif)



## Files saved during program execution

1. __10000_peak.csv__: To save computation later on, two peaks of interest (10,000 and 11,200 cm-1) are automatically saved as a csv for later analysis. This is helpful for methane and ethane mixtures, but can be commented out if undesired.
2. __11200_peak.csv__: Same as above, this time for the 11,200 cm-1 peak.
3. __10000area.txt__: The trapezoidal rule is invoked on the 10000_peak.csv file to calculate the integrated area.
4. __11200area.txt__: same as above, this time for the 11,200 cm-1 peak.
5. __alpha_coeffs.csv__: Absorption coefficients calculated during the program, saved as a csv file for later plotting and calculations as needed
6. __alpha_coeffs.pdf__: screenshot of the absorption coefficients taken during program execution for later reference
7. __continuum_chosen.pdf__: screenshot of the continuum chosen during program execution for later reference
8. __continuum_polynomial.txt__: The continuum polynomial is saved for future use if someone wants to plot the curve again
9. __dv_filt_specturm.pdf__: The spectrum after it has been filtered and a background divide (raw/background) has been performed
10. __continuum_polynomial.txt__: The continuum polynomial is saved for future use if someone wants to plot the curve again
11. __dv_filt_specturm.pdf__: The spectrum after it has been filtered and a background divide (raw/background) has been performed
12. __continuum_polynomial.txt__: The continuum polynomial is saved for future use if someone wants to plot the curve again
13. __dv_filt_specturm.pdf__: The spectrum after it has been filtered and a background divide (raw/background) has been performed
14. __fft_data.csv__: the raw data after the FFT filter has been applied for future reference if needed
15. __FFT_filter.pdf__: A screenshot of the filter window chosen during program execution
16. __freq_window.csv__: a csv file of the values picked for the FFT filter
17. __FFT_Raw_bg_data.csv__: same as above, but for the background data
18. __fft_background.pdf__: a screenshot of the FFT performed on the background data during program execution
19. __rawbackground.pdf__: __a screenshot of the raw background data being plotted at the beginning of the program for future reference
20. __rawspecturm.pdf__: same as above but for the raw data.
21. __temperature.txt__: a txt file with the user's input of the temperature of the sample
22. __concentration.txt__: a txt file with the user's input of the concentration of the sample
23. __data_files_used.txt__: a txt file with the user's input of the data files used during program execution
24. __time_created.txt__: a timestamp of when the program was executed

### If you have any questions, feel free to contact me at garrett.leland.thompson@gmail.com


<!---
## Demo of program, with a pure methane (Ch4) mixture at a temperature of 92.0 Kelvin.
![first gif](https://github.com/Gleland/SpectralAnalysis/blob/master/images/full_video.gif) 

-->
