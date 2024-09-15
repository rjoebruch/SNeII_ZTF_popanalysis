

# from ztfquery import marshal
from astropy.coordinates import SkyCoord 
from astropy.table import Table, vstack, hstack
from astropy.io import ascii 
import numpy as np
import matplotlib.pylab as plt
from numpy import random
import pandas as pd

import random
import math
import scipy 
from scipy.optimize import curve_fit

from scipy.special import gamma
from scipy.integrate import quad

from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.stats import median_abs_deviation
from scipy.signal import savgol_filter



import sfdmap 
import extinction
from PyAstronomy import pyasl
from astropy.coordinates import Distance


from astropy import cosmology

from iminuit import Minuit




import pprint




def to_mag(f):
    '''
    this funtion converts flux to magnitude
    
    parameters
    ----------
    f.  column or array?
    
    returns
    -------
    '''
    f = abs(f)
    return -2.5*np.log10(f)

def error_on_conv_to_mag(eflux,flux):
    '''
    this function computes the error on the conversion from flux to mag
    
    parameters
    ----------
    eflux   error on flux
    flux.   flux values
    returns
    -------
    '''
    # eflux = abs(eflux)
    # flux  = abs(flux)
    return (2.5*eflux)/(np.log(10)*flux)



def bazin(t,A,tau_fall, tau_rise, B):
    '''
    This function is a phenomelogical function for SN lightcurves 

    parameters
    ----------


    returns
    -------
    
    '''
    
    fluxy =  np.where( t > 0,  ( A * np.exp(-(t)/tau_fall) ) / (1+np.exp((t)/tau_rise)) + B, 0 )
    
    return fluxy

  

def bazin_log(t,A,t_0,tau_fall, tau_rise, B):
    '''
    This function is a phenomelogical function for SN lightcurves 

    parameters
    ----------
    A
    t_0
    tau_fall
    tau_rise
    B

    returns
    -------
    
    '''

    return -2.5*np.log10(A*((  np.exp(-(t-t_0)/tau_fall) ) / (1+np.exp((t-t_0)/tau_rise))     ) + B)




def broken_PL(t,A, alpha1, t_break, alpha2, n):
    '''
    This function is a phenomelogical function used to describe therise to peak, with a change of regime at t_break. 
    See Eq. 2 in https://iopscience.iop.org/article/10.1086/503294/pdf

    parameters
    ----------
    t   [array or list or column?]
    A
    alpha1
    t_break
    alpha2
    n


    returns
    -------
    
    '''


    f = A * ( (t/t_break)**(-alpha1*n) + (t/t_break)**(-alpha2*n) )**(-1/n)
    

    return f


def bounded_broken_PL(t, A, alpha1, t_break, alpha2, n, t_peak, peak_val):
    
    '''
    This function is a phenomelogical function used to describe therise to peak, with a change of regime at t_break. 
    See Eq. 2 in https://iopscience.iop.org/article/10.1086/503294/pdf

    parameters
    ----------
    t   [array or list or column?]
    A
    alpha1
    t_break
    alpha2
    n
    t_peak   [float] determined time of peak 
    peak_val [float] flux value at peak calculated earlier
    
    returns
    -------
    
    '''

    # f = A * ( (t/t_break)**(-alpha1*n) + (t/t_break)**(-alpha2*n) )**(-1/n)

    # hiuv = np.where( t < t_peak, f , peak_val )

    f =  - A * ( (-(t-t_peak)/t_break)**(-alpha1*n) + (-(t-t_peak)/t_break)**(-alpha2*n) )**(-1/n) + peak_val

    return f



def simple_EXP_law(t,A, n):
    '''
    This function is a simple power law used to fit the early time data. 

    parameters
    ----------
    t
    A
    n


    returns
    -------
    
    '''



    f = A * t**n

    return f


def bounded_simple_EXP_law(t,A, n, t_peak, peak_val):
    '''
    This function is a simple power law used to fit the early time data. 

    parameters
    ----------
    t
    A
    n


    returns
    -------
    
    '''

    
    # f = A * t**n 

    f = - A * (-(t-t_peak))**n + peak_val

    # hiuv = np.where( t < t_peak, f , peak_val )


    return f

def bounded_linear_decline(t,a,b,t_peak):
    '''
    This function is used to fit a linear decline to the light curves of standard SNe II

    parameters
    ----------
    a
    b [float] peak value at peak time 
    t_peak

    
    '''

    return a*(t-t_peak) + b




def function_s1_s2(t,a1,t0,m0,a2 ):
    '''
    
    This function is broken line with one break point used to describe the decline from peak
    of SN II (from the silver sample)

    parameters
    ----------
    t : interval of time over which we are fitting
    a1: coefficient directeur de la premiere courbe
    t0 - m0 : breakpoint time and magnitude 

    # a1
    # b1 [float] peak value at peak time -FIXED
    # t_peak - FIXED 
    # a2
    # b2 - Il faut que b2 = a1*(t_break - t_peak) + b1


    

    >> Le meilleur outil c'est le piecewise de numpy je pense


    
    '''

    # a1*(t-t_peak) + b1 + a2*(t-t_peak)+b2

    # def piecewise_linear(x, x0, y0, k1, k2):
    # T = t-
    s1s2 = np.piecewise(t, [t < t0], [lambda x:a1*x + m0-a1*t0, lambda x:a2*x + m0-a2*t0])

    ## ça c'est pas anchored au peak... comment on force qu'a t = tpeak c'est egal a mpeak? 


    return s1s2





def bin_LC(table, day_min, day_max,sug = 'median', plot = False):
    '''
    This function smoothes the magnitude light curve (considering only the detections)
    Three types: mean binning, median binning and moving average
    
    parameters
    ----------
    table     [table] containing time and mag
    day_min   [int]   
    day_max   [int]
    sug       [string] -optional- 'mean', 'median'
    plot      [string] 'None' for no plot, 'mag' in magnitude, 'flux' in flux
    
    returns
    -------
    
    
    '''
    
    bin_phot_table_ = table[['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag', 'absmag','e_absmag']][ table['tfromexplo_zc'] < day_min ]
    if sug == 'mean':
        for i in range(day_min,day_max, 1):
            temp_phot_ =   table[['tfromexplo_zc', 'extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']][ (i <= table['tfromexplo_zc']) & (table['tfromexplo_zc'] < i+1) ]
            if len(temp_phot_) == 1:
                mu_time    = temp_phot_['tfromexplo_zc'][0]
                mu_flux    = temp_phot_['extcorrforcediffimflux'][0]
                sig_flux   = temp_phot_['extcorrforcediffimfluxunc'][0]
                mu_mag     = temp_phot_['mag'][0]
                sig_mag    = temp_phot_['emag'][0]
                mu_amag    = temp_phot_['absmag'][0]
                sig_amag   = temp_phot_['e_absmag'][0]

                bin_phot_table_.add_row([ mu_time ,mu_flux, sig_flux, mu_mag , sig_mag, mu_amag, sig_amag ] )
                
            if len(temp_phot_) > 1:
                mu_time    = np.average(temp_phot_['tfromexplo_zc'])
                mu_flux    = get_weighted_mean( np.array(temp_phot_['extcorrforcediffimflux']), np.array(temp_phot_['extcorrforcediffimfluxunc']) )
                sig_flux   = get_std_err_wmean(np.array(temp_phot_['extcorrforcediffimflux']), np.array(temp_phot_['extcorrforcediffimfluxunc']) )
                mu_mag     = get_weighted_mean( np.array(temp_phot_['mag']), np.array(temp_phot_['emag']) )
                sig_mag    = get_std_err_wmean(np.array(temp_phot_['mag']), np.array(temp_phot_['emag']) )
                mu_amag    = get_weighted_mean( np.array(temp_phot_['absmag']), np.array(temp_phot_['e_absmag']) )
                sig_amag   = get_std_err_wmean(np.array(temp_phot_['absmag']), np.array(temp_phot_['e_absmag']) )

                bin_phot_table_.add_row([ mu_time ,mu_flux, sig_flux, mu_mag , sig_mag, mu_amag, sig_amag ] )

    elif sug == 'median':
        for i in range(day_min,day_max, 1):
            temp_phot_ =   table[['tfromexplo_zc', 'extcorrforcediffimflux','extcorrforcediffimfluxunc','mag', 'emag', 'absmag','e_absmag']][ (i <= table['tfromexplo_zc']) & (table['tfromexplo_zc'] < i+1) ]
            if len(temp_phot_) == 1:
                mu_time    = temp_phot_['tfromexplo_zc'][0]
                mu_flux    = temp_phot_['extcorrforcediffimflux'][0]
                sig_flux   = temp_phot_['extcorrforcediffimfluxunc'][0]
                mu_mag     = temp_phot_['mag'][0]
                sig_mag    = temp_phot_['emag'][0]
                mu_amag    = temp_phot_['absmag'][0]
                sig_amag   = temp_phot_['e_absmag'][0]

                bin_phot_table_.add_row([ mu_time ,mu_flux, sig_flux, mu_mag , sig_mag, mu_amag, sig_amag ] )
               
            
            if len(temp_phot_) > 1:
                mu_time   = np.average(temp_phot_['tfromexplo_zc'])
                mu_flux   = np.median( np.array(temp_phot_['extcorrforcediffimflux']) )
                sig_flux  = get_std_err_wmean(np.array(temp_phot_['extcorrforcediffimflux']), np.array(temp_phot_['extcorrforcediffimfluxunc']) )
                mu_mag    = np.median( np.array(temp_phot_['mag']) )
                sig_mag   = get_std_err_wmean( np.array(temp_phot_['mag']), np.array(temp_phot_['emag']) )
                mu_amag    = np.median( np.array(temp_phot_['absmag']) )
                sig_amag   = get_std_err_wmean( np.array(temp_phot_['absmag']), np.array(temp_phot_['e_absmag']) )


                bin_phot_table_.add_row([ mu_time, mu_flux, sig_flux ,mu_mag , sig_mag, mu_amag, sig_amag ] )
                    
        

                
    if plot == 'mag':
        plt.figure()
        plt.plot(table['tfromexplo_zc'],table['mag'],'x', color = 'grey', alpha = 0.5)
        plt.plot(bin_phot_table_['tfromexplo_zc'],bin_phot_table_['mag'],'.', color = 'black')
        plt.gca().invert_yaxis()
        plt.xlabel('Time from FD')
        plt.ylabel('App. Mag')

    elif plot == 'flux':
        plt.figure()
        plt.plot(table['tfromexplo_zc'],table['extcorrforcediffimflux'],'x', color = 'grey', alpha = 0.5)
        plt.plot(bin_phot_table_['tfromexplo_zc'],bin_phot_table_['extcorrforcediffimflux'],'.', color = 'black')
        plt.xlabel('Time from FD')
        plt.ylabel('Flux')
    

    return bin_phot_table_




def moving_average(a, n=3) :
    '''
    this function computes the moving average over n index shift (by default 3)
    
    parameters
    ----------
    a [array/column] 
    n [int]
    
    returns
    -------
    smoothed data
    
    '''
    ret     = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n






def get_weighted_mean(val, e_val):
    '''
         this function computes the weighted mean of an esemble of points assuming that they are distributed following a normal distribution
    
        parameters
        ----------
        val      [numpy array] collection of values (e.g. magnitude, flux)
        e_val    [numpy array] collection of errors on the sus-mentioned values
        /!\ they need to be the same length 
    
        returns
        -------
        mu    [float] weighted mean
    '''
    
    num_   = np.sum(val/(e_val**2))
    denom_ = np.sum(1/(e_val**2))

    return num_/denom_




def get_std_err_wmean(vals,err):

    '''
    computes the standard error of the weighted mean 
    
    parameters
    ----------
    err : array of errors values
    
    returns
    -------
    the standard error of the weighted mean
    '''
    _muweigt = get_weighted_mean(vals,err)
    
    _weights = 1/(err**2)
    
    _var     = (np.sum(_weights*(vals**2))/np.sum(_weights) - _muweigt**2 )*(len(vals)/(len(vals)-1))
    
    return np.sqrt(_var/len(vals))




def interpolate_LC(t,mag,len_interpo = 1000, plot=False):
    '''
    this function interpolates the lightcurve prior to computing the autocorrelation
    Recommend to cut the lightcurve to <150 days 
    
    parameters
    ----------
    t           [column/array] original column/array of time
    mag         [column/array] original column/array of magnitudes
    len_interpo [int] -optional- len of the interpolated light curve. Can be changed according to number of points
    plot        [bool] -optional- plot the lightcurve to visualise the interpolation
    
    
    returns
    -------
    array of time, array of interpolated magnitudes 
    
    '''
    

    
    interp_mag = interp1d(t, mag,   bounds_error=False, fill_value='nan')

    interp_t   = np.linspace(min(t), max(t), len_interpo)

    if plot is True:

        plt.figure()
        plt.plot(t,mag,'.')
        plt.plot(interp_t, interp_mag(interp_t),'-')
        plt.gca().invert_yaxis()
        plt.xlabel('Time from FD')
        plt.ylabel('App. Mag')

    return interp_t, interp_mag(interp_t)

def interpolate_LC_with_errors(t,mag,emag,mininterp,maxinterp,peaktg,len_interpo = 1000, plot=False):
    '''
    this function interpolates the lightcurve prior to computing the autocorrelation
    Recommend to cut the lightcurve to <150 days 

    /!\
    THIS FUNCTION WAS WRITTEN TO INTERPOLATE THE INTERPOLATION (just so that we have the right datapoints foir color evolution)
    

    parameters
    ----------
    t           [column/array] original column/array of time
    mag         [column/array] original column/array of magnitudes
    emag        [column/array] original column/array of errors on magnitudes
    peaktg      [float] time of peak g. We're adding this specific point to calculate the color at peak g /!\\\\ 
    len_interpo [int] -optional- len of the interpolated light curve. Can be changed according to number of points
    plot        [bool] -optional- plot the lightcurve to visualise the interpolation
    
    
    returns
    -------
    array of time, array of interpolated magnitudes, array of interpolated errors on magnitude
    
    '''
    

    
    interp_mag  = interp1d(t, mag,   bounds_error=False, fill_value='nan')
    interp_emag = interp1d(t, emag,   bounds_error=False, fill_value='nan')

    interp_t   = np.linspace(mininterp, maxinterp, len_interpo)
    interp_t   = np.insert(interp_t,1,peaktg)
    interp_t.sort()

    if plot is True:
        
        plt.figure()

        h = plt.plot(t, mag, color = 'grey', alpha = 0.8 )
        plt.fill_between(t,mag + np.sqrt(emag), mag - np.sqrt(emag), 
                 alpha=0.5, color=h[0].get_color())

        i = plt.plot(interp_t, interp_mag(interp_t), color = 'gold', alpha = 0.8 )
        plt.fill_between(interp_t, interp_mag(interp_t) + np.sqrt(interp_emag(interp_t)), interp_mag(interp_t) - np.sqrt(interp_emag(interp_t)), 
                 alpha=0.4, color=i[0].get_color())


        plt.gca().invert_yaxis()
        plt.xlabel('Time from FD')
        plt.ylabel('App. Mag')

    return interp_t, interp_mag(interp_t), interp_emag(interp_t)
    
    

        




def interpolate_detec_LC(t,mag,len_interpo = 1000):
    '''
    this function interpolates the lightcurve prior to computing the autocorrelation
    Recommend to cut the lightcurve to <150 days 
    
    parameters
    ----------
    t           [column/array] original column/array of time
    mag         [column/array] original column/array of magnitudes
    len_interpo [int] -optional- len of the interpolated light curve. Can be changed according to number of points

    
    returns
    -------
    table of interpolation
    
    '''
    
    interp_mag = interp1d(t, mag,   bounds_error=False, fill_value='nan')
    interp_t   = np.linspace(min(t), max(t), len_interpo)
    
    interpopo  = Table(data=[interp_t, interp_mag(interp_t)], names = ('days', 'flux'))
    return interpopo


################### Functions for light curve smoothing and interpolation

def interpolate_fullLC_with_errors(t,mag,emag,mininterp,maxinterp,len_interpo = 1000, plot=False):
    '''
    this function interpolates the lightcurve prior to smoothing withSG filter
    Recommend to cut the lightcurve to <150 days 


    parameters
    ----------
    t           [column/array] original column/array of time
    mag         [column/array] original column/array of magnitudes
    emag        [column/array] original column/array of errors on magnitudes
    peaktg      [float] time of peak g. We're adding this specific point to calculate the color at peak g /!\\\\ 
    len_interpo [int] -optional- len of the interpolated light curve. Can be changed according to number of points
    plot        [bool] -optional- plot the lightcurve to visualise the interpolation
    
    
    returns
    -------
    Table containaing 
    array of time, array of interpolated magnitudes, array of interpolated errors on magnitude
    
    '''
    


    interp_mag  = interp1d(t, mag,   bounds_error=False, fill_value='extrapolate')
    interp_emag = interp1d(t, emag,   bounds_error=False, fill_value='extrapolate') ##interpolating errors 

    interp_t   = np.linspace(mininterp, maxinterp, len_interpo)
    interp_t.sort()

    if plot is True:
        
        plt.figure()

        h = plt.plot(t, mag, color = 'grey', alpha = 0.8 )
        plt.fill_between(t,mag + np.sqrt(emag), mag - np.sqrt(emag), 
                 alpha=0.5, color=h[0].get_color())

        i = plt.plot(interp_t, interp_mag(interp_t), color = 'gold', alpha = 0.8 )
        plt.fill_between(interp_t, interp_mag(interp_t) + np.sqrt(interp_emag(interp_t)), interp_mag(interp_t) - np.sqrt(interp_emag(interp_t)), 
                 alpha=0.4, color=i[0].get_color())


        plt.gca().invert_yaxis()
        plt.xlabel('Time from FD')
        plt.ylabel('App. Mag')
        
    table = Table(data = ([interp_t, interp_mag(interp_t), interp_emag(interp_t)]), names = ('t','magn','e_mag'))
    return table



def sg_perpart(interpolated_lc, peak_T):
    '''
    This function applies a smoothing with different windows to the light curve, ide est:
    for the rise part (0 to day 5): window of three data points
    for the "to peak" part (day 5 to 15): window of five data points 
    for the "to plateau" part (day 15 onwards): window of 9 datapoints
    
    parameters
    ----------
    tje grossly interpolated lc with "interpolate_LC_with_errors"
    
    returns
    ------
    two arrays of the smoothed absmag lc and error on abs mag
    
    '''
    # TODO: optimise the parameters: with a X2 test? or a grid of model > willl still be a least squared test 
    
    if len(peak_T)>0:
    
        rise      = interpolated_lc[interpolated_lc['t']<=peak_T]
        # print(len(rise))
        afterrise = interpolated_lc[interpolated_lc['t']>peak_T]
        # topeak    = interpolated_lc[(interpolated_lc['t']>3)&(interpolated_lc['t']<=15)]
        # toplateau = interpolated_lc[interpolated_lc['t']>15]
    else:
        rise      = interpolated_lc[interpolated_lc['t']<=15]
        afterrise = interpolated_lc[interpolated_lc['t']>15]
       
    
    if len(rise)>1:
        rise_sg      = savgol_filter(rise['magn'],3,2 )
        topeak_sg    = savgol_filter(afterrise['magn'],101,3 )
#         toplateau_sg = savgol_filter(toplateau['magn'],501,3 )

        e_rise_sg      = savgol_filter(rise['e_mag'],3,2 )
        e_topeak_sg    = savgol_filter(afterrise['e_mag'],101,3 )

        full_sg   = np.concatenate((rise_sg, topeak_sg))
        e_full_sg = np.concatenate((e_rise_sg, e_topeak_sg))

    elif len(rise)<=1:
        topeak_sg    = savgol_filter(afterrise['magn'],101,3 )

        e_topeak_sg    = savgol_filter(afterrise['e_mag'],101,3 )

        full_sg   = topeak_sg
        e_full_sg = e_topeak_sg

        full_sg   = np.concatenate((rise['magn'], topeak_sg))
        e_full_sg = np.concatenate((rise['e_mag'], e_topeak_sg))

    full_sg   = savgol_filter(full_sg,91,3)
    e_full_sg = savgol_filter(e_full_sg,91,3)
    
    full_sg   = savgol_filter(full_sg,201,3)
    e_full_sg = savgol_filter(e_full_sg,201,3)

    return full_sg,e_full_sg