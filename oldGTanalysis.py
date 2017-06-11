########################################################################################
#Authored by Garrett Thompson
#Graphical User Interface for Data Analysis
#Created at Northern Arizona University
#	for use at the Astrophysical Ice Laboratory
#Advisors: Jennifer Hanley, Will Grundy, Henry Roe
#garrett.leland.thompson@gmail.com
########################################################################################
"""
NOTE: two warnings are being silenced, a RankWarning, when creating a polynomial for continuum dividng,
      and a runtime invalid wraning, when doing -log()/thickness. Both are normal during operation, and
			shouldn't be mistaken as bad programming. See documentation for further explanation before removing
			lines of code that silence warnings
	
"""

#only external package dependencies are Numpy, Scipy, and Matplotlib
import os
import sys
import time
import warnings #for when creating polynomial curves
import numpy as np
import matplotlib
matplotlib.use("TkAgg")#allows matplotlib and Tkinter to get along

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf #for building continuums
from scipy.fftpack import fft, fftfreq, ifft #fft filtering
from scipy.signal import savgol_filter as sgf #savitkzy-golay filter

#for choosing files through Finder interface
from Tkinter import Tk
from tkFileDialog import askopenfilename
# writing out data files generated for later use
import csv

#warning when using sgf option
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
#https://github.com/scipy/scipy/issues/5998 if wanting to see more

#assumes a csv file, as all data stored from ice lab is in CSV format
def import_data(filename):
	raw_data = np.loadtxt(open(filename,"rb"),delimiter=",")
	xdat = raw_data[:,0]
	ydat = raw_data[:,1]
	return xdat,ydat


def freq_click(event):
	#plots data to cutout frequencies
	#global fig3,frq_x,fft_yb,fft_y,vert_lines, norm_fft
	plt.xlim(plt.gca().get_xlim())
	plt.ylim(plt.gca().get_ylim())
	#print(('button: ',event.button))
	if event.button==1: #left click, adds vertical line
		#if button_click = left: add left line
		#if button_click = middle: removes closest line
		#if button_lick = right: finish
		#		try:
		#			#only delete if a curve_fit line has already been drawn
		#			if len(ax.lines) !=1: ax.lines.remove(ax.lines[-1])
		#	except: UnboundLocalError
		#add clicked data points to list
		vert_lines.append(event.xdata)
		#print(vert_lines)
		#confirmation, and displays values in terminal for reference
		#print(('left click! ',(event.xdata, event.ydata)))
		
		frq_ax.plot(frq_x,np.log(np.abs(fft_ybg)),color='blue')
		plt.axvline(x=vert_lines[-1],color='black')
		plt.xlabel('Cycles/Wavenumber')
		plt.ylabel('Relative Intensity')
		#draws points as they are added
		plt.draw()
	#calc_coeffs(pvals)
	if event.button==2:#middle click, remove closest vertical line
		print ('pop!')
		#gets x,y limits of graph,saves them before destroying figure
		xlims = plt.gca().get_xlim()
		ylims = plt.gca().get_ylim()
		#clears axes, to get rid of old scatter points
		frq_ax.cla()
		#re-plots spectrum
		frq_ax.plot(frq_x,np.log(np.abs(fft_ybg)),color='blue')
		#sets axes limits to original values
		plt.xlim(xlims)
		plt.ylim(ylims)
		plt.xlabel('Cycles/Wavenumber')
		plt.ylabel('Relative Intensity')
		#deletes last recorded points (x,y coords)
		#xcoords.pop()
		#ycoords.pop()
		#deletes point closest to mouse click
		#takes the difference of every element with the coords, and finds the smallest value.
		xindx = np.abs(vert_lines-event.xdata).argmin()
		del vert_lines[xindx]
		print(vert_lines)
		for line in vert_lines:
			plt.axvline(x=line,color='black')
		#draws the new set of vertical lines
		#plt.axvline(x=vert_lines[-1],color='black')
		#tells the figure to update
		plt.draw()
	#calc_coeffs(pvals)
	if event.button==3:#right click, ends collection
		#ends clicking awareness
		frq_fig.canvas.mpl_disconnect(frq_cid)
		#title = raw_input('Enter title to save, or press enter to skip\n:')
		#if len(title) > 0:
		#plt.title(title)
		plt.savefig('FFT_filter.pdf')
		with open("freq_window.csv", "w") as f:
			writer = csv.writer(f)
			writer.writerow(["Xposition of vert. line"])
			writer.writerows(zip(vert_lines))
		
		window_filter(*vert_lines[:2])#first window
		if len(vert_lines) > 2:
			window_filter(*vert_lines[2:4])#second window if desired
		elif len(vert_lines) > 3:
			print('you might have added an extra window, only first four lines will be recorded')

		fft_calc()




def fft_calc():
	global norm_fft,frq_x
	#dividing filtered y data from filtered bg data
	norm_fft = ifft(filt_y)/ifft(filt_ybg)
	#cutout outliers if necessary
	#for i in range(len(raw_x)):
	#	if div_fft[i] > 3 or div_fft[i] < 0: div_fft[i]=0.
	#	if div_smooth[i] >3 or div_smooth[i] < 0: div_smooth[i]=0.
	#normailze if desired
	#norm_fft = norm_fft / max(norm_fft)
	rows = zip(raw_x,norm_fft.real)
	with open("fft_data.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["raw_x","fft_filt"])
		writer.writerows(rows)
	plot_data(raw_x,norm_fft.real)
#plot_data(raw_x,norm_smooth)



def sgf_calc(window_param, poly_param, raw_y, raw_ybg):
	global norm_smooth
	smoothed_y = sgf(raw_y,window_length=window_param,polyorder=poly_param,delta=(abs(raw_y)[1]-raw_y)[0])
	smoothed_ybg =sgf(raw_ybg,window_length=window_param,polyorder=poly_param,delta=(abs(raw_ybg)[1]-raw_ybg)[0])
	#dividing filtered y data from filtered bg data
	norm_smooth = smoothed_y/ smoothed_ybg
	#cutout outliers if necessary
	#for i in range(len(raw_x)):
	#	if div_fft[i] > 3 or div_fft[i] < 0: div_fft[i]=0.
	#	if div_smooth[i] >3 or div_smooth[i] < 0: div_smooth[i]=0.
	#normailze if desired
	#norm_smooth = norm_smooth / max(norm_smooth)
	rows = zip(raw_x,norm_smooth)
	with open("sgf_data.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["window","polynomail order"])
		writer.writerow([window_param,poly_param])
		writer.writerow(["raw_x","sgf_filt"])
		writer.writerows(rows)
	plot_data(raw_x,norm_smooth)




#range of frequenices to cut out
#visually determined to be around window_min,window_max = 0.2,0.8
def window_filter(window_min,window_max):
	for i in range(len(frq_x)):
		if (frq_x[i] >= window_min and frq_x[i] <=window_max) or (frq_x[i]>-1*window_max and frq_x[i]<-1*window_min):
			filt_y[i] = 0
			filt_ybg[i] = 0


def plot_data(a,b):
	#variables I will need in different scopes
	global fig,ax,cid,x,y,xcoords,ycoords,order
	#filenames so it doesn't have to be typed out
	x = a
	y = b
	#list of x and y values
	xcoords = []
	ycoords = []
	#make a figure for graphing
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#plot data with relevant settings
	ax.plot(x,y)
	plt.suptitle('Divide and Filtered Spectrum')
	plt.xlabel('Wavenumber cm-1')
	plt.ylabel('Relative Intensity')
	plt.show(block=False)
	plt.xlim(plt.gca().get_xlim())
	plt.ylim(plt.gca().get_ylim())
	plt.savefig('dv_filt_spectrum.pdf')
	#asks user what order polynomial for the continuum fit
	order = input('Zoom to liking and then enter what order polynomial for continuum fit\n:')
	#tells python to turn on awareness for button presses
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	print 'Left to add, middle to remove nearest, and right to finish'
	#displays plot
	plt.show()



#for creating continuum fit to divide out
def onclick(event):
	global fig,fig2,ax,ax2,cid,x,y,xcoords,ycoords,order,pvals
	plt.xlim(plt.gca().get_xlim())
	plt.ylim(plt.gca().get_ylim())
	#print(('button: ',event.button))
	if event.button==1: #left click
	#if button_click = left: add
	#if button_click = middle: remove
	#if button_lick = right: finished
		try:
			#only delete if a curve_fit line has already been drawn
			if len(ax.lines) !=1: ax.lines.remove(ax.lines[-1])
		except: UnboundLocalError
		#add clicked data points to list
		xcoords.append(event.xdata)
		ycoords.append(event.ydata)
		#confirmation, and displays values in terminal for reference
		#print(('left click! ',(event.xdata, event.ydata)))
		ax.scatter(xcoords,ycoords,color='black')
		plt.xlabel('Wavenumber cm-1')
		plt.ylabel('Relative Intensity')
		#draws points as they are added
		plt.draw()
		#want to draw curve as soon as enough data points present
		xvals = np.array(xcoords)
		yvals = np.array(ycoords)
		#fits values to polynomial, because we aren't using polynomial for stats, rankwarning is irrelevant
		warnings.simplefilter('ignore', np.RankWarning)
		p_fit = np.polyfit(xvals,yvals,order)
		pvals = np.poly1d(p_fit)
		ax.plot(x,pvals(x),color='black')
		plt.show(block=False)
		#calc_coeffs(pvals)
	if event.button==2:#middle click, remove closest point to click
		print ('pop!')
		#gets x,y limits of graph,saves them before destroying figure
		xlims = plt.gca().get_xlim()
		ylims = plt.gca().get_ylim()
		#clears axes, to get rid of old scatter points
		ax.cla()
		#re-plots spectrum
		ax.plot(x,y)
		#sets axes limits to original values
		plt.xlim(xlims)
		plt.ylim(ylims)
		plt.xlabel('Wavenumber cm-1')
		plt.ylabel('Relative Intensity')
		#deletes point closest to mouse click
		#takes the difference of every element with the coords, and finds the smallest value.
		xindx = np.abs(xcoords-event.xdata).argmin()
		del xcoords[xindx]
		yindx = np.abs(ycoords-event.ydata).argmin()
		del ycoords[yindx]
		#draws the new set of scatter points, and colors them
		ax.scatter(xcoords,ycoords,color='black')
		#tells the figure to update
		plt.draw()
		#making the lists into numpy arrays for curve fitting
		xvals = np.array(xcoords)
		yvals = np.array(ycoords)
		#fits values to polynomial, because we aren't using polynomial for stats, rankwarning is irrelevant
		warnings.simplefilter('ignore', np.RankWarning)
		p_fit = np.polyfit(xvals,yvals,order)
		pvals = np.poly1d(p_fit)
		ax.plot(x,pvals(x),color='black')
		plt.draw()
	if event.button==3:#right click, ends collection
		#ends clicking awareness
		#fig.canvas.mpl_disconnect(cid)
		#save plot image for furture inspection
		plt.savefig('continuum_chosen.pdf')
		#Saving polynomial eq'n used in continuum divide for reference

		with open("continuum_polynomial.txt", "w") as save_file:
			save_file.write("%s *x^ %d  " %(pvals[0],0))
			for i in (xrange(len(pvals))):
				save_file.write("+ %s *x^ %d  " %(pvals[i+1],i+1))
		calc_coeffs(pvals)


#calculates alpha coefficients
def calc_coeffs(pvals):
	global fig,fig2,ax,ax2,cid,x,y,xcoords,ycoords,order
	fit_y = pvals(x)
	#flattens the continuum
	new_continuum = y / fit_y
	thickness = int(raw_input('\nEnter thickness of cell in cm\n')) #2 cm for our work in 2016
	#remove runtime errors when taking negative log and dividing
	err_settings = np.seterr(invalid='ignore')
	alpha_coeffs = -np.log(new_continuum) / thickness
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	#plot data with relevant settings
	plt.suptitle('Alpha Coefficients')
	plt.show(block=False)
	ax2.plot(x,alpha_coeffs)
	plt.xlabel('Wavenumber cm-1')
	plt.ylabel('Alpha')
	plt.savefig('alpha_coeffs.pdf')
	plt.draw()
	alpha_rows = zip(x,alpha_coeffs)
	with open("alpha_coeffs.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["x","alpha"])
		writer.writerows(alpha_rows)
	#calculating without dividing by thickness for later processing
	#won't plot because user can manipulate data later and plot separately
	no_divide_alpha = -np.log(new_continuum)
	no_divide_rows = zip(x,no_divide_alpha)
	with open("no_divide_alpha_coeffs.csv", "w") as g:
		writer = csv.writer(g)
		writer.writerow(["x","no_divide_alpha"])
		writer.writerows(no_divide_rows)
	finish_prog = raw_input("Press 'y' when finished\n:")
	check = True
	while check:
		if (finish_prog =="y"): check = False
	plt.close('all')
	print "Finished!"
	quit() # end of program



#start of program, everything above are just functions being defined
########################################################################################




# Where all work to follow will be saved
folder_to_save = raw_input('Type name of directory to save all data being created\n:')

#saves terminal arguments as variables
#arg1 = sys.argv[1] # name of folder to create


#make and change to directory named by user
os.mkdir(folder_to_save)
os.chdir(folder_to_save)

#recording date and time that program is run, saving it to folder
with open("time_created.txt", "w") as text_file:
	text_file.write("Time this program was run: {} \n".format(time.strftime("%Y-%m-%d %H:%M")))


#Creating a window for user to pick files by hand, rather than typing out long strings
#instantiate a Tk window
root = Tk()
# prevents Tk window from opening needlessly
root.withdraw()

#dunno what this does, fixes askopenfilename if I use it.
root.update() #allows popup to close and open multiple times


#import raw data
print "choose a raw dataset for analysis in the pop up window"
raw_import = askopenfilename()
print  "\nGot it! Importing now... \n"
raw_x,raw_y = import_data(raw_import)

#import raw background
print "choose a raw background for analysis in the pop up window"
bg_import = askopenfilename()
print  "\nGot it! Importing now... \n"
raw_xbg,raw_ybg = import_data(bg_import)

#saving text file so user knows in future what datasets were used to create alpha coeffs
with open("data_files_used.txt", "w") as text_file:
	text_file.write("Raw data file used: {} \n".format(raw_import))
	text_file.write("Raw background data file used: {}".format(bg_import))


#plot raw spectum for inspections
print "plotting imported data..."
raw_plot = plt.figure()
raw_ax = raw_plot.add_subplot(111)
plt.suptitle('Raw Spectrum')
raw_ax.plot(raw_x,raw_y)
plt.xlabel('Wavenumber (cm-1)')
plt.ylabel('% Transmittance')
plt.show(block=False)
plt.savefig('rawspectrum.pdf')


#plot bg spectrum for inspection
bg_plot = plt.figure()
bg_ax = bg_plot.add_subplot(111)
plt.suptitle('Raw Background')
plt.xlabel('Wavenumber (cm-1)')
plt.ylabel('% Transmittance')
bg_ax.plot(raw_xbg,raw_ybg)
plt.show(block=False)
plt.savefig('rawbackground.pdf')

"""
#writing raw and bg files for use later
raw_rows = zip(raw_x,raw_y,raw_xbg,raw_ybg)
with open("raw_data.csv", "w") as f:
	writer = csv.writer(f)
	writer.writerow(["raw_x","raw_y","raw_xbg","raw_ybg"])
	writer.writerows(raw_rows)
"""

#based on the data, which method does the user want?
user_method = raw_input('Press "s" for savitsky-golay filter, or any other key for fft filter\n:')
if user_method.lower() == 's':
	#savitsky-golay option was chosen
	window_param = int(raw_input('Input window box size (must be odd number)\n:'))
	poly_param = int(raw_input('Input polynomial order for smoothing\n:'))
	#saving parameters chosen for future inspection
	with open("sgf_params.txt", "w") as sgf_file:
		sgf_file.write("Window parameter used: {} \n".format(window_param))
		sgf_file.write("Polynomial paramter used: {}".format(poly_param))
	
	sgf_calc(window_param,poly_param,raw_y,raw_ybg)
else:
	#fft option was chosen
	
	#finds FFT of ydata
	fft_y = fft(raw_y)
	fft_ybg = fft(raw_ybg)

	#gets corresponding frequencies for FFT of data by using length of array, and sample spacing
	frq_x = fftfreq(len(fft_y),((max(raw_x)-min(raw_x))/len(fft_y)))
	frq_xbg = fftfreq(len(fft_ybg),((max(raw_xbg)-min(raw_xbg))/len(fft_ybg)))

	#replacing values within window with zero to nip ringing from sapphire windows
	filt_y = fft_y.copy()
	filt_ybg = fft_ybg.copy()
	#plots graph of fft spectrum, for nipping freqs that are undesired
	vert_lines = []
	frq_fig = plt.figure()
	frq_ax = frq_fig.add_subplot(111)
	plt.suptitle('FFT of raw bg')
	frq_ax.plot(frq_x,np.log(np.abs(fft_ybg)))
	raw_fft_rows = zip(frq_x,np.log(abs(fft_ybg)))
	#saving raw fft data for later inspection
	with open("FFT_Raw_bg_data.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["frq_x","fft_bg"])
		writer.writerows(raw_fft_rows)
	#comparing if left and right are symmetric, but they are not
	#frq_ax.plot(frq_x,np.log(abs(fft_ybg-fft_ybg[::-1])))
	plt.show(block=False)
	raw_input('zoom to liking, then press enter to start')
	print 'Left to add, middle to remove nearest, and right to finish'
	frq_cid = frq_fig.canvas.mpl_connect('button_press_event', freq_click)
	plt.show()




