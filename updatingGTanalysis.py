##############################################################################
# Created by Garrett Thompson
# Graphical User Interface for Data Analysis
# Created at Northern Arizona University
#	for use in the Astrophysical Ice Laboratory
# Advisors: Jennifer Hanley, Will Grundy, Henry Roe
# garrett.leland.thompson@gmail.com
##############################################################################

import os
import time
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
from scipy.fftpack import fft, fftfreq, ifft
from scipy.signal import savgol_filter as sgf
from scipy.integrate import trapz
import csv

def main():
    print("main called")
    folder_to_save = choose_dir()
    #choose files for analysis
    raw_x,raw_y, raw_xbg,raw_ybg = choose_files(folder_to_save)

    print("plotting imported data...")
    plotting_data_for_inspection(raw_x,raw_y,'Raw Data','Wavenumber (cm-1)','% Transmittance','rawspectrum.pdf',False)
    plotting_data_for_inspection(raw_xbg,raw_ybg,'Raw Background','Wavenumber (cm-1)','% Transmittance','rawbackground.pdf',False)

    #user chooses method after inspecting plots
    user_method = str(raw_input('Press "s" for savitsky-golay filter, or "f" for fft filter\n:'))
    choosing = True
    while choosing:
        if user_method.lower() == 's':
            # savitsky-golay option was chosen
            choosing = False
            raw_x, norm_smooth = sgf_calc()
            plot_data(raw_x,norm_smooth)
        elif user_method.lower() == 'f':
            # fft option was chosen
            choosing = False
            frq_x,frq_xbg,fft_y,fft_ybg = fft_calculation(raw_x,raw_y,raw_xbg,raw_ybg)
            plot_figure, plot_axis = plotting_data_for_inspection(frq_x,np.log(abs(fft_ybg)),'FFT of raw bg','Cycles/Wavenumber (cm)','Log(Power/Frequency)','fft_background.pdf',False)
            filt_y = fft_y.copy()
            filt_ybg = fft_ybg.copy()
            raw_input('zoom to liking, then press enter to start')
            print 'Left to add, middle to remove nearest, and right to finish'
            global frq_cid 
            #frq_cid = frq_fig.canvas.mpl_connect('button_press_event',lambda event: freq_click(event, [fft_ybg,frq_fig]))
            vert_lines=[]
            frq_cid = plot_figure.canvas.mpl_connect('button_press_event',lambda event: freq_click(event, [frq_x,fft_ybg,plot_figure,plot_axis,vert_lines,filt_y,filt_ybg]))
            plt.show()





def fft_calculation(raw_x,raw_y,raw_xbg,raw_ybg): 
    """ calculates FFT of data for use in nipping unwanted frequencies"""
    # finds FFT of ydata
    fft_y = fft(raw_y)
    fft_ybg = fft(raw_ybg)

    # gets frequencies for FFT of data from array, and sample spacing
    frq_x = fftfreq(len(fft_y),((max(raw_x)-min(raw_x))/len(fft_y)))
    frq_xbg = fftfreq(len(fft_ybg),((max(raw_xbg)-min(raw_xbg))/len(fft_ybg)))

    raw_fft_rows = zip(frq_x,np.log(abs(fft_ybg)))
    # saving raw fft data for later inspection
    with open("FFT_Raw_bg_data.csv", "w") as f:
            writer = csv.writer(f)
            writer.writerow(["frq_x","fft_bg"])
            writer.writerows(raw_fft_rows)
    return frq_x, frq_xbg, fft_y, fft_ybg


def choose_dir():
    """
    User chooses where all work will be saved and 
    time stamp is created for future reference
    """
     
    # Where all work to follow will be saved
    folder_to_save = raw_input('Type name of directory to save all data being created\n:')
    # make and change to directory named by user
    os.mkdir(folder_to_save)
    os.chdir(folder_to_save)

    # recording date and time that program is run, saving it to folder
    with open(str(folder_to_save) + "time_created.txt", "w") as text_file: 
        text_file.write("Time this program was run: {} \n".format(time.strftime("%Y-%m-%d %H:%M")))
    os.chdir('..')
    return folder_to_save

def plotting_data_for_inspection(xdata,ydata,plot_title,plot_xlabel,plot_ylabel,filename_for_saving,block_boolean):
    """ 
    Plots data for user to look at within program

    parameters
    ----------
    xdata,ydata: x and y data to be plotted
    plot_xlabel,plot_ylabel: label x and y axes in plot
    file_name_for_saving: string given for saving file for later referece
    block_boolean: True or False, tells if program waits for figure to close
    """
    print 'plot_data_for_inspection was called'
    plot_figure, plot_axis = plt.subplots()
    plt.plot(xdata,ydata)
    plt.xlabel(plot_xlabel)
    plt.ylabel(plot_ylabel)
    plt.suptitle(plot_title)
    plt.show(block=block_boolean)
    plt.savefig(filename_for_saving)
    return plot_figure, plot_axis


def choose_files(folder_to_save):
    """ 
    Lets user determine which files will be imported for analysis
    and saves preferences for reference later on
    """
    #raw_import = str(raw_input('Enter a raw dataset for analysis\n:'))
    raw_import ='RawData.csv' 
    print  "\nGot it! Importing now... \n"
    raw_x,raw_y = import_data(raw_import)

    #bg_import = str(raw_input('Enter a raw background for analysis\n:'))
    bg_import ='BackgroundData.csv' 
    print  "\nGot it! Importing now... \n"
    raw_xbg,raw_ybg = import_data(bg_import)
    os.chdir(folder_to_save)
    with open("data_files_used.txt", "w") as text_file:
            text_file.write("Raw data file used: {} \n".format(raw_import))
            text_file.write("Raw background data file used: {}".format(bg_import))

    # saving text file of concentration for later use in plotting
    with open("concentration.txt","w") as f:
            f.write(folder_to_save[-12:-9])

    # saving text file of temperature for later use in plotting
    with open("temperature.txt","w") as f:
            f.write(folder_to_save[-5:-1])
    os.chdir('..')
    return raw_x, raw_y,raw_xbg,raw_ybg


# assumes a csv file, as all data stored from ice lab is in CSV format
def import_data(filename):
        print "import_data called"
	raw_data = np.loadtxt(open(filename,"rb"),delimiter=",")
	xdat = raw_data[:,0]
	ydat = raw_data[:,1]
	return xdat,ydat


def freq_click(event, args_list):
        print "freq_click called"
        #fft_ybg,frq_fig = args_list
        frq_x,fft_ybg,plot_figure,plot_axis,vert_lines, filt_y, filt_ybg = args_list

	plt.xlim(plt.gca().get_xlim())
	plt.ylim(plt.gca().get_ylim())
	if event.button==1: 
		# if button_click = left: add left line
		# if button_click = middle: removes closest line
		# if button_lick = right: finish
		# add clicked data points to list
		vert_lines.append(event.xdata)
		plot_axis.plot(frq_x,np.log(np.abs(fft_ybg)),color='blue')
		plt.axvline(x=vert_lines[-1],color='black')
		plt.xlabel('Cycles/Wavenumber')
		plt.ylabel('Relative Intensity')
		# draws points as they are added
		plt.draw()
	if event.button==2:
                # middle click, remove closest vertical line
		print ('pop!')
		# gets x,y limits of graph,saves them before destroying figure
		xlims = plt.gca().get_xlim()
		ylims = plt.gca().get_ylim()
		# clears axes, to get rid of old scatter points
		plot_axis.cla()
		# re-plots spectrum
		plot_axis.plot(frq_x,np.log(np.abs(fft_ybg)),color='blue')
		# sets axes limits to original values
		plt.xlim(xlims)
		plt.ylim(ylims)
		plt.xlabel('Cycles/Wavenumber')
		plt.ylabel('Relative Intensity')
		# deletes point closest to mouse click
		xindx = np.abs(vert_lines-event.xdata).argmin()
		del vert_lines[xindx]
		for line in vert_lines:
			plt.axvline(x=line,color='black')
		# draws the new set of vertical lines
		# tells the figure to update
		plt.draw()
	if event.button==3:
                # right click, ends collection
		# ends clicking awareness
		plot_figure.canvas.mpl_disconnect(frq_cid)
		plt.savefig('FFT_filter.eps')
		with open("freq_window.csv", "w") as f:
			writer = csv.writer(f)
			writer.writerow(["Xposition of vert. line"])
			writer.writerows(zip(vert_lines))
		# first window
                args_list =[vert_lines,frq_x,filt_y,filt_ybg] 
		filt_y,filt_ybg = window_filter(args_list)
		#if len(vert_lines) > 2:
                #        # second window if desired
		#	window_filter(vert_lines[2:4],frq_x)
		#elif len(vert_lines) > 3:
                #        print('you might have added an extra window, only first four lines will be recorded')
		fft_calc(filt_y, filt_ybg, frq_x)

def fft_calc(filt_y, filt_ybg, frq_x):
	# dividing filtered y data from filtered bg data
	norm_fft = ifft(filt_y)/ifft(filt_ybg)
	rows = zip(raw_x,norm_fft.real)
	with open("fft_data.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["raw_x","fft_filt"])
		writer.writerows(rows)
	plot_data(raw_x,norm_fft.real)



def sgf_calc():
        print "sgf_calc called"
        # warning when using sgf option
        warnings.filterwarnings(action="ignore", module="scipy",message="^internal gelsd")
        window_param = int(raw_input('Input window box size (must be odd number)\n:'))
        poly_param = int(raw_input('Input polynomial order for smoothing\n:'))
        # saving parameters chosen for future inspection
        with open("sgf_params.txt", "w") as sgf_file:
                sgf_file.write("Window parameter used: {} \n".format(window_param))
                sgf_file.write("Polynomial paramter used: {}".format(poly_param))
	global norm_smooth
	smoothed_y = sgf(raw_y,window_param,poly_param,delta=(abs(raw_y)[1]-raw_y)[0])
	smoothed_ybg =sgf(raw_ybg,window_param,poly_param,delta=(abs(raw_ybg)[1]-raw_ybg)[0])
	# dividing filtered y data from filtered bg data
	norm_smooth = smoothed_y / smoothed_ybg
	rows = zip(raw_x,norm_smooth)
	with open("sgf_data.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["window","polynomail order"])
		writer.writerow([window_param,poly_param])
		writer.writerow(["raw_x","sgf_filt"])
		writer.writerows(rows)
	return raw_x,norm_smooth

# range of frequenices to cut out
def window_filter(args_list):
        vert_lines, frq_x, filt_y, filt_ybg = args_list
        window_min, window_max= vert_lines[-2], vert_lines[-1]
        print "window_filter called"
	for i in range(len(frq_x)):
		if (frq_x[i] >= window_min and frq_x[i] <=window_max) or (frq_x[i]>-1*window_max and frq_x[i]<-1*window_min):
			filt_y[i] = 0
			filt_ybg[i] = 0
        return filt_y,filt_ybg


def plot_data(x,y):
        print "plot_data called"
	# variables I will need in different scopes
        global xcoords,ycoords
	# list of x and y values
	xcoords = []
	ycoords = []
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(x,y)
	plt.suptitle('Divide and Filtered Spectrum')
	plt.xlabel('Wavenumber cm-1')
	plt.ylabel('Relative Intensity')
	plt.show(block=False)
	plt.xlim(plt.gca().get_xlim())
	plt.ylim(plt.gca().get_ylim())
	plt.savefig('dv_filt_spectrum.pdf')
	order = input('Zoom to liking and then enter what order polynomial for continuum fit\n:')
	# tells python to turn on awareness for button presses
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	print 'Left to add, middle to remove nearest, and right to finish'
	plt.show()

# for creating continuum fit to divide out
def onclick(event):
        print "onclick called"
	global fig,fig2,ax,ax2,cid,x,y,xcoords,ycoords,order,pvals
	plt.xlim(plt.gca().get_xlim())
	plt.ylim(plt.gca().get_ylim())
	# print(('button: ',event.button))
	if event.button==1: 
                # left click
		try:
			# only delete if curve_fit line already drawn
			if len(ax.lines) !=1: ax.lines.remove(ax.lines[-1])
		except: UnboundLocalError
		# add clicked data points to list
		xcoords.append(event.xdata)
		ycoords.append(event.ydata)
		ax.scatter(xcoords,ycoords,color='black')
		plt.xlabel('Wavenumber cm-1')
		plt.ylabel('Relative Intensity')
		plt.draw()
		xvals = np.array(xcoords)
		yvals = np.array(ycoords)
		# fits values to polynomial, rankwarning is irrelevant
		warnings.simplefilter('ignore', np.RankWarning)
		p_fit = np.polyfit(xvals,yvals,order)
		pvals = np.poly1d(p_fit)
		ax.plot(x,pvals(x),color='black')
		plt.show(block=False)
	if event.button==2:
                # middle click, remove closest point to click
		print ('pop!')
		# gets x,y limits of graph,saves them before destroying figure
		xlims = plt.gca().get_xlim()
		ylims = plt.gca().get_ylim()
		# clears axes, to get rid of old scatter points
		ax.cla()
		# re-plots spectrum
		ax.plot(x,y)
		# sets axes limits to original values
		plt.xlim(xlims)
		plt.ylim(ylims)
		plt.xlabel('Wavenumber cm-1')
		plt.ylabel('Relative Intensity')
		# deletes point closest to mouse click
		xindx = np.abs(xcoords-event.xdata).argmin()
		del xcoords[xindx]
		yindx = np.abs(ycoords-event.ydata).argmin()
		del ycoords[yindx]
		# draws the new set of scatter points, and colors them
		ax.scatter(xcoords,ycoords,color='black')
		plt.draw()
		xvals = np.array(xcoords)
		yvals = np.array(ycoords)
		# fits values to polynomial, rankwarning is ignored
		warnings.simplefilter('ignore', np.RankWarning)
		p_fit = np.polyfit(xvals,yvals,order)
		pvals = np.poly1d(p_fit)
		ax.plot(x,pvals(x),color='black')
		plt.draw()
	if event.button==3:
		# right click,ends clicking awareness
		plt.savefig('continuum_chosen.pdf')
		# Saving polynomial eqn used in continuum divide for reference
		with open("continuum_polynomial.txt", "w") as save_file:
			save_file.write("%s *x^ %d  " %(pvals[0],0))
			for i in (xrange(len(pvals))):
				save_file.write("+ %s *x^ %d  " %(pvals[i+1],i+1))
		calc_coeffs(pvals)

def calc_coeffs(pvals):
        print "calc_coeffs called"
	global fig,fig2,ax,ax2,cid,x,y,xcoords,ycoords,order
	fit_y = pvals(x)
	# flattens the continuum
	new_continuum = y / fit_y
	thickness = int(raw_input('\nEnter thickness of cell in cm\n'))
        # 2 cm thickness for our work in 2016
	# remove runtime errors when taking negative log and dividing
	err_settings = np.seterr(invalid='ignore')
	alpha_coeffs = -np.log(new_continuum) / thickness
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
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

	# defining where the interesting peaks are
	c1= (x>10000) & (x<10500)
	c2 = (x>11200) & (x<12000)
	# creating masks, only around each peak
	xm1,ym1 = x[c1],alpha_coeffs[c1]
	xm2,ym2 = x[c2],alpha_coeffs[c2]

	# writing data for plotting later
	r1 = zip(xm1,ym1)
	r2 = zip(xm2,ym2)

	with open("10000_peak.csv","w") as f:
		writer = csv.writer(f)
		writer.writerow(["x","y"])
		writer.writerows(r1)

        with open("11200_peak.csv","w") as f:
		writer = csv.writer(f)
		writer.writerow(["x","y"])
		writer.writerows(r2)

	area10000=trapz(ym1,xm1)
	area11200=trapz(ym2,xm2)
	with open("10000area.txt","w") as f:
		f.write(str(area10000))
	with open("11200area.txt","w") as f:
		f.write(str(area11200))
	finish_prog = raw_input("Press 'y' when finished\n:")
	check = True
	while check:
		if (finish_prog =="y"): check = False
	plt.close('all')
	print "Finished!"
	quit() # end of program


if __name__ == '__main__':
    main()

