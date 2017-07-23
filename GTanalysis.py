##############################################################################
# Created by Garrett Thompson
# Graphical User Interface for Data Analysis
# Created at Northern Arizona University
#       for use in the Astrophysical Ice Laboratory
# Advisors: Jennifer Hanley, Will Grundy, Henry Roe
# garrett.leland.thompson@gmail.com
##############################################################################

import os
import csv
import time
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
from scipy.fftpack import fft, fftfreq, ifft
from scipy.signal import savgol_filter as sgf
from scipy.integrate import trapz

def main():
    folder_to_save = choose_dir()
    #choose files for analysis
    raw_x,raw_y, raw_xbg,raw_ybg = choose_files(folder_to_save)

    print("Plotting imported data...")
    plotting_data_for_inspection(raw_x,raw_y,'Raw Data','Wavenumber (cm-1)','% Transmittance','rawspectrum.pdf',folder_to_save, False)
    plotting_data_for_inspection(raw_xbg,raw_ybg,'Raw Background','Wavenumber (cm-1)','% Transmittance','rawbackground.pdf',folder_to_save, False)

    #user chooses method after inspecting plots
    user_method = str(input('Press "s" for savitsky-golay filter, or "f" for fft filter\n:'))
    choosing = True
    while choosing:
        if user_method.lower() == 's':
            # savitsky-golay option was chosen
            choosing = False
            args_list = [folder_to_save, raw_y, raw_ybg, raw_x]
            raw_x, norm_smooth = sgf_calc(args_list)
            plot_data(raw_x,norm_smooth,folder_to_save)
        elif user_method.lower() == 'f':
            # fft option was chosen
            choosing = False
            frq_x,frq_xbg,fft_y,fft_ybg = fft_calculation(raw_x,raw_y,raw_xbg,raw_ybg,folder_to_save)
            plot_figure, plot_axis = plotting_data_for_inspection(frq_x,np.log(abs(fft_ybg)),'FFT of raw bg','Cycles/Wavenumber (cm)','Log(Power/Frequency)','fft_background.pdf',folder_to_save, False)
            filt_y = fft_y.copy()
            filt_ybg = fft_ybg.copy()
            input('Zoom to liking, then press enter to start')
            print('Left to add, middle to remove nearest, and right to finish')
            # global frq_cid 
            vert_lines=[]
            frq_cid = plot_figure.canvas.mpl_connect('button_press_event',lambda event: freq_click(event, [frq_x,fft_ybg,plot_figure,plot_axis,vert_lines,filt_y,filt_ybg,folder_to_save,raw_x]))
            plt.show()
            plot_figure.canvas.mpl_disconnect(frq_cid)
            # vert_lines, frq_x, filt_y, filt_ybg = args_dict["vert_lines"],args_dict["frq_x"],args_dict["filt_y"],args_dict["filt_ybg"]



def save_as_csv(folder_to_save,title, column1_title,column2_title,column1_data,column2_data):
    os.chdir(folder_to_save)
    with open(title,"w") as f:
        writer = csv.writer(f)
        writer.writerow([column1_title,column2_title])
        writer.writerows(list(zip(column1_data,column2_data)))
    os.chdir('..')


def fft_calculation(raw_x,raw_y,raw_xbg,raw_ybg,folder_to_save): 
    """ calculates FFT of data for use in nipping unwanted frequencies"""
    # finds FFT of ydata
    fft_y = fft(raw_y)
    fft_ybg = fft(raw_ybg)
    # gets frequencies for FFT of data from array, and sample spacing
    frq_x = fftfreq(len(fft_y),((max(raw_x)-min(raw_x))/len(fft_y)))
    frq_xbg = fftfreq(len(fft_ybg),((max(raw_xbg)-min(raw_xbg))/len(fft_ybg)))
    save_as_csv(folder_to_save,"FFT_Raw_bg_data.csv","frq_x","log(abs(fft_bg))",frq_x,np.log(abs(fft_ybg)))
    return frq_x, frq_xbg, fft_y, fft_ybg


def choose_dir():
    """
    User chooses where all work will be saved and 
    time stamp is created for future reference
    """
    # Where all work to follow will be saved
    folder_to_save = input('Type name of directory to save all data being created\n:')
    # make and change to directory named by user
    os.mkdir(folder_to_save)
    os.chdir(folder_to_save)
    # recording date and time that program is run, saving it to folder
    with open("time_created.txt", "w") as text_file: 
        text_file.write("Time this program was run: {} \n".format(time.strftime("%Y-%m-%d %H:%M")))
    os.chdir('..')
    return folder_to_save

def plotting_data_for_inspection(xdata,ydata,plot_title,plot_xlabel,plot_ylabel,filename_for_saving,folder_to_save, block_boolean):
    """ 
    Plots data for user to look at within program

    parameters
    ----------
    xdata,ydata: x and y data to be plotted
    plot_xlabel,plot_ylabel: label x and y axes in plot
    file_name_for_saving: string given for saving file for later referece
    block_boolean: True or False, tells if program waits for figure to close
    """
    plot_figure, plot_axis = plt.subplots()
    plt.plot(xdata,ydata,color='blue')
    plt.xlabel(plot_xlabel)
    plt.ylabel(plot_ylabel)
    plt.suptitle(plot_title)
    plt.show(block=block_boolean)
    os.chdir(folder_to_save)
    plt.savefig(filename_for_saving)
    os.chdir('..')
    return plot_figure, plot_axis


def choose_files(folder_to_save):
    """ 
    Lets user determine which files will be imported for analysis
    and saves preferences for reference later on
    """
    raw_import = str(input('Enter a raw dataset for analysis\n:'))
    print("\nGot it! Importing now... \n")
    raw_x,raw_y = import_data(raw_import)

    bg_import = str(input('Enter a raw background for analysis\n:'))
    print("\nGot it! Importing now... \n")
    raw_xbg,raw_ybg = import_data(bg_import)
    os.chdir(folder_to_save)
    with open("data_files_used.txt", "w") as text_file:
            text_file.write("Raw data file used: {} \n".format(raw_import))
            text_file.write("Raw background data file used: {}".format(bg_import))

    concentration = str(input('Enter concentration of mixture\n:'))
    # saving text file of concentration for later use in plotting
    with open("concentration.txt","w") as f:
            f.write(concentration)

    temperature = str(input('Enter temperature of mixture\n:'))
    # saving text file of temperature for later use in plotting
    with open("temperature.txt","w") as f:
            f.write(temperature)
    os.chdir('..')
    return raw_x, raw_y,raw_xbg,raw_ybg


# assumes a csv file, as all data stored from ice lab is in CSV format
def import_data(filename):
        raw_data = np.loadtxt(open(filename,"rb"),delimiter=",")
        xdat = raw_data[:,0]
        ydat = raw_data[:,1]
        return xdat,ydat


def freq_click(event, args_list):
        # if button_click = left: add left line
        # if button_click = middle: removes closest line
        # if button_lick = right: finish
        # add clicked data points to list
        frq_x,fft_ybg,plot_figure,plot_axis,vert_lines, filt_y, filt_ybg,folder_to_save, raw_x = args_list

        plt.xlim(plt.gca().get_xlim())
        plt.ylim(plt.gca().get_ylim())
        if event.button==1: 
                vert_lines.append(event.xdata)
                plot_axis.plot(frq_x,np.log(np.abs(fft_ybg)),color='blue')
                #plt.axvline(x=vert_lines[-1],color='black')
                for val in vert_lines:
                    plt.axvline(x=val,color='black')
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
                plt.draw()
        if event.button==3:
                # right click, ends clicking awareness 
                # plot_figure.canvas.mpl_disconnect(frq_cid)
                os.chdir(folder_to_save)
                plt.savefig('FFT_filter.pdf')
                with open("freq_window.csv", "w") as f:
                        writer = csv.writer(f)
                        writer.writerow(["Xposition of vert. line"])
                        writer.writerows(list(zip(vert_lines)))
                os.chdir('..')
                # first window
                args_dict ={"vert_lines":vert_lines,"frq_x":frq_x,"filt_y":filt_y,"filt_ybg":filt_ybg}
                plt.close("all")
                argslist =  [vert_lines,frq_x,filt_y,filt_ybg] 
                filt_y,filt_ybg = window_filter(argslist)
                fft_calc(filt_y, filt_ybg, raw_x,folder_to_save)

def fft_calc(filt_y, filt_ybg, raw_x,folder_to_save):
        # dividing filtered y data from filtered bg data
        norm_fft = ifft(filt_y)/ifft(filt_ybg)
        save_as_csv(folder_to_save,"fft_data.csv","raw_x","fft_filt",raw_x,norm_fft.real)
        plot_data(raw_x,norm_fft.real,folder_to_save)



def sgf_calc(args_list):
        folder_to_save, raw_y, raw_ybg, raw_x = args_list
        # warning when using sgf option
        warnings.filterwarnings(action="ignore", module="scipy",message="^internal gelsd")
        window_param = int(input('Input window box size (must be odd number)\n:'))
        poly_param = int(input('Input polynomial order for smoothing\n:'))
        # saving parameters chosen for future inspection
        os.chdir(folder_to_save)
        with open("sgf_params.txt", "w") as sgf_file:
                sgf_file.write("Window parameter used: {} \n".format(window_param))
                sgf_file.write("Polynomial paramter used: {}".format(poly_param))
        #global norm_smooth
        smoothed_y = sgf(raw_y,window_param,poly_param,delta=(abs(raw_y)[1]-raw_y)[0])
        smoothed_ybg =sgf(raw_ybg,window_param,poly_param,delta=(abs(raw_ybg)[1]-raw_ybg)[0])
        # dividing filtered y data from filtered bg data
        norm_smooth = smoothed_y / smoothed_ybg
        rows = list(zip(raw_x,norm_smooth))
        with open("sgf_data.csv", "w") as f:
                writer = csv.writer(f)
                writer.writerow(["window","polynomail order"])
                writer.writerow([window_param,poly_param])
                writer.writerow(["raw_x","sgf_filt"])
                writer.writerows(rows)
        os.chdir('..')
        return raw_x,norm_smooth

# range of frequenices to cut out
def window_filter(args_list):
        vert_lines, frq_x, filt_y, filt_ybg = args_list
        window_min, window_max= vert_lines[-2], vert_lines[-1]
        for i in range(len(frq_x)):
                if (frq_x[i] >= window_min and frq_x[i] <=window_max) or (frq_x[i]>-1*window_max and frq_x[i]<-1*window_min):
                        filt_y[i] = 0
                        filt_ybg[i] = 0
        return filt_y,filt_ybg


def plot_data(x,y,folder_to_save):
        plot_figure,plot_axis = plotting_data_for_inspection(x,y,"Divide and Filtered Spectrum","Wavenumber cm-1","Relative Intensity","dv_filt_spectrum.pdf",folder_to_save, False)
        order = int(input('Zoom to liking and then enter what order polynomial for continuum fit\n:'))
        xcoords,ycoords = [],[]
        # tells python to turn on awareness for button presses
        global cid
        cid = plot_figure.canvas.mpl_connect('button_press_event', lambda event: onclick(event, [xcoords,ycoords,plot_figure,plot_axis,order,folder_to_save,x,y]))
        print('Left to add, middle to remove nearest, and right to finish')
        plt.show()

# for creating continuum fit to divide out
def onclick(event,argslist):
        xcoords,ycoords,plot_figure,plot_axis,order,folder_to_save,x,y = argslist
        global pvals
        if event.button==1: 
                # left click
                plt.xlim(plt.gca().get_xlim())
                plt.ylim(plt.gca().get_ylim())
                #plt.cla()
                try:
                        # only delete if curve_fit line already drawn
                        if len(plot_axis.lines) !=1: plot_axis.lines.remove(plot_axis.lines[-1])
                except: UnboundLocalError
                # add clicked data points to list
                xcoords.append(event.xdata)
                ycoords.append(event.ydata)
                plot_axis.scatter(xcoords,ycoords,color='black')
                plt.xlabel('Wavenumber cm-1')
                plt.ylabel('Relative Intensity')
                plt.draw()
                xvals = np.array(xcoords)
                yvals = np.array(ycoords)
                # fits values to polynomial, rankwarning is irrelevant
                warnings.simplefilter('ignore', np.RankWarning)
                p_fit = np.polyfit(xvals,yvals,order)
                pvals = np.poly1d(p_fit)
                plot_axis.plot(x,pvals(x),color='black')
                plt.draw()
                # plt.show(block=False)
        if event.button==2:
                # middle click, remove closest point to click
                print ('pop!')
                # gets x,y limits of graph,saves them before destroying figure
                xlims = plt.gca().get_xlim()
                ylims = plt.gca().get_ylim()
                # clears axes, to get rid of old scatter points
                plot_axis.cla()
                # re-plots spectrum
                plot_axis.plot(x,y)
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
                plot_axis.scatter(xcoords,ycoords,color='black')
                plt.draw()
                xvals = np.array(xcoords)
                yvals = np.array(ycoords)
                # fits values to polynomial, rankwarning is ignored
                warnings.simplefilter('ignore', np.RankWarning)
                p_fit = np.polyfit(xvals,yvals,order)
                pvals = np.poly1d(p_fit)
                plot_axis.plot(x,pvals(x),color='black')
                plt.draw()
        if event.button==3:
                # right click,ends clicking awareness
                plot_figure.canvas.mpl_disconnect(cid)
                os.chdir(folder_to_save)
                plt.savefig('continuum_chosen.pdf')
                # Saving polynomial eqn used in continuum divide for reference
                with open("continuum_polynomial.txt", "w") as save_file:
                        save_file.write("%s *x^ %d  " %(pvals[0],0))
                        for i in (range(len(pvals))):
                                save_file.write("+ %s *x^ %d  " %(pvals[i+1],i+1))
                os.chdir('..')
                calc_coeffs(pvals,x,y,folder_to_save)

def calc_coeffs(pvals,x,y,folder_to_save):
        fit_y = pvals(x)
        # flattens the continuum
        new_continuum = y / fit_y
        thickness = int(input('\nEnter thickness of cell in cm\n:'))
        # 2 cm thickness for our work in 2016
        # remove runtime errors when taking negative log and dividing
        err_settings = np.seterr(invalid='ignore')
        alpha_coeffs = -np.log(new_continuum) / thickness
        plotting_data_for_inspection(x,alpha_coeffs,"Alpha Coefficients","Wavenumber cm-1","Absorption cm-1","alpha_coeffs.pdf",folder_to_save,False)
        save_as_csv(folder_to_save,"alpha_coeffs.csv","x","alpha",x,alpha_coeffs)

        # creating masks around each peak
        x_mask1 = x[(x>10000) & (x<10500)]
        x_mask2 = x[(x>11200) & (x<12000)]
        y_mask1 = alpha_coeffs[(x>10000) & (x<10500)]
        y_mask2 = alpha_coeffs[(x>11200) & (x<12000)]

        # writing data for plotting later
        save_as_csv(folder_to_save,"10000_peak.csv","x","y",x_mask1,y_mask1)
        save_as_csv(folder_to_save,"11200_peak.csv","x","y",x_mask2,y_mask2)

        # integrated area calcs
        area10000=trapz(y_mask1,x_mask1)
        area11200=trapz(y_mask2,x_mask2)
        os.chdir(folder_to_save)
        with open("10000area.txt","w") as f:
                f.write(str(area10000))
        with open("11200area.txt","w") as f:
                f.write(str(area11200))
        os.chdir('..')
        finish_prog = input("Press 'y' when finished\n:")
        check = True
        while check:
                if (finish_prog =="y"): check = False
        plt.close('all')
        print("Finished!")
        quit() # end of program


if __name__ == '__main__':
    main()

