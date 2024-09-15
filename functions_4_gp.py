
'''

This file contains the functions necessary to compute a gaussian processed interpolation of the ZTF P48 lightcurves

'''

import numpy as np
from astropy.io import ascii
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.interpolate import splev, splrep
from scipy.optimize import curve_fit
import george
import os
import math
from astropy.table import vstack, Table, hstack



###### FUNCTIONS #######

def sq_exp(x, l):
    '''
    function to fit autocorrelated function. Squared exponential. should be the on eused for gaussian processes
    (need to fit only half of it)
    parameters
    ----------
    x [array] time sequence (time lag )
    l [float] (parameter to optimize) typical size of
    
    
    returns
    -------
    
    
    '''
    
    return np.exp(-( (x**2) / (2*(l**2)) ) )


def bin_LC(table, day_min, day_max,sug = 'mean', plot = False):
    '''
    This function smoothes the magnitude light curve (considering only the detections)
    Three types: mean binning, median binning and moving average
    
    parameters
    ----------
    table     [table] containing time and mag
    day_min   [int]   
    day_max   [int]
    sug       [string] -optional- 'mean', 'median'
    plot
    
    returns
    -------
    
    
    '''
    
    bin_phot_table_ = table[['tfromexplo_zc','mag','emag']][ table['tfromexplo_zc'] < day_min ]
    if sug == 'mean':
        for i in range(day_min,day_max, 1):
            temp_phot_ =   table[['tfromexplo_zc', 'mag','emag']][ (i <= table['tfromexplo_zc']) & (table['tfromexplo_zc'] < i+1) ]
            if len(temp_phot_) == 1:
                mu_time  = temp_phot_['tfromexplo_zc'][0]
                mu_mag   = temp_phot_['mag'][0]
                sig_mag  = temp_phot_['emag'][0]

                bin_phot_table_.add_row([ mu_time , mu_mag , sig_mag ] )
                
            if len(temp_phot_) > 1:
                mu_time  = np.average(temp_phot_['tfromexplo_zc'])
                mu_mag   = get_weighted_mean( np.array(temp_phot_['mag']), np.array(temp_phot_['emag']) )
                sig_mag  = get_std_err_wmean(np.array(temp_phot_['mag']), np.array(temp_phot_['emag']) )

                bin_phot_table_.add_row([ mu_time , mu_mag , sig_mag ] )

    elif sug == 'median':
        for i in range(day_min,day_max, 1):
            temp_phot_ =   table[['tfromexplo_zc', 'mag', 'emag']][ (i <= table['tfromexplo_zc']) & (table['tfromexplo_zc'] < i+1) ]
            if len(temp_phot_) == 1:
                mu_time  = temp_phot_['tfromexplo_zc'][0]
                mu_mag   = temp_phot_['mag'][0]
                sig_mag  = temp_phot_['emag'][0]

                bin_phot_table_.add_row([ mu_time , mu_mag , sig_mag ] )
               
            
            if len(temp_phot_) > 1:
                mu_time_     = np.average(temp_phot_['tfromexplo_zc'])
                mu_flux_     = np.median( np.array(temp_phot_['mag']) )
                sig_mag      = get_std_err_wmean( np.array(temp_phot_['mag']), np.array(temp_phot_['emag']) )


                bin_phot_table_.add_row([ mu_time , mu_mag , sig_mag ] )
                    
        

                
    if plot is True:
        plt.figure()
        plt.plot(table['tfromexplo_zc'],table['mag'],'x', color = 'grey', alpha = 0.5)
        plt.plot(bin_phot_table_['tfromexplo_zc'],bin_phot_table_['mag'],'.', color = 'black')
        plt.gca().invert_yaxis()
        plt.xlabel('Time from FD')
        plt.ylabel('App. Mag')

    return bin_phot_table_


def interpolate_LC(t,mag,len_interpo = 1000, plot='False'):
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



def autocorrelate_LC(t,mag, plot = False, table = False):
    '''
    this function computes the autocorrelation of the magnitude 
    The code goes through an array of  time shifts (of one day) and
    computes the Pearson product-moment correlation coefficients (see numpy doc)
    
    > should take as input both the time and the magnitude
    > time should be a "linear space" (equally spaced) array and magnitudes should be the interpolation of the LC 
    
    parameters
    ----------
    t     [array] interpolated time SHOULD HAVE EQUALLY SPACED T 
    mag   [array] interpolated magnitude 
    plot  [bool] -optional- plot the autocorrelation 
    table [bool] -optional- if you want to return it in a table format
    
    returns
    -------
    two arrays: lags and the autocorrelation
    
    omelette du fromage
    
    '''

    _autoco = []
    _laggy  = []
    
    for shift in range(len(mag) - 1):
        if shift == 0 :
            autoco = np.corrcoef( np.array([mag , mag]) )[0,1]
            _autoco.append(autoco)
            _laggy.append(0)
            
        else:
            autoco = np.corrcoef( np.array([mag[:-shift] , mag[shift:]]) )[0,1]
            diff_T = t[shift]-t[0]
            if diff_T <= 20:
                _autoco.append(autoco)
                _laggy.append(diff_T)
            else:
                pass
    

    
    if plot is True:
        plt.figure()
        plt.plot(_laggy,_autoco)
        plt.xlabel("Lags [days]")
        plt.xlabel("N-P correlation coef.")
    
    if table is True: 
        _taby = Table([_laggy, _autoco], names = ('lags', 'autocorrelation'))
        return _taby
    elif table is False:
        return  _laggy, _autoco


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


def moving_average(a, n=20) :
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



def gp_interp_lc(lc, len_kernel , filtre, length_interp = 1000, mode = 'compute'):
    '''
    This function performs the guassian process interpolaiton based on the kernel size fitted before
    
    parameters
    ----------
    lc            [astropy table] magnitude light curve (tfromexplo_zc, mag, e_mag)
    len_kernel    [array/list]    The length of the kernel (squared exponential). Should be popt[0] if you fit the autocorrelation of the lightcurve prior 
    filtre        [string]        filter of the LC interp
    length_interp [int] how many data points to interpolate between min and max of th light curve. 
    mode          [string] 'compute' when you want to have the interpolation at time which were not measured and  'residual' when testing for chi2
    
    returns
    -------
    gp_fit  [astropy table] the GP interpolated fit
    
    '''

    if mode == 'compute':
        # time scale over which we're interpolating
        x_fit  = np.linspace(min(lc['tfromexplo_zc']),max(lc['tfromexplo_zc']),length_interp)

        # introducing the parameters of the kernel we fit earlier
        kernel =  np.var(lc['mag']) * george.kernels.ExpSquaredKernel(len_kernel) 

        # generating and computing the GP
        gp     = george.GP(kernel)
        gp.compute(lc['tfromexplo_zc'], yerr=lc['emag']*2)

        # applying the GP to the LC
        # NOTE: the GP wants to go to "0" (defined mean)
        # , hence we need to put the data arounbd a mean of zero. We add back the mean of the LC afterwards
        y_fit, y_var = gp.predict(lc['mag'] - np.mean(lc['mag']), x_fit, return_var=True) #added here the "back to 0"
        y_fit       += np.mean(lc['mag'])

        gp_fit = Table([x_fit, y_fit, y_var, [filtre]*len(x_fit)], names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))

        return gp_fit

    elif mode == 'residual':
         # time scale over which we're interpolating
        x_fit  = lc['tfromexplo_zc']

        # introducing the parameters of the kernel we fit earlier
        kernel =  np.var(lc['mag']) * george.kernels.ExpSquaredKernel(len_kernel) 

        # generating and computing the GP
        gp     = george.GP(kernel)
        gp.compute(lc['tfromexplo_zc'], yerr=lc['emag']*2)

        # applying the GP to the LC
        # NOTE: the GP wants to go to "0" (defined mean)
        # , hence we need to put the data arounbd a mean of zero. We add back the mean of the LC afterwards
        y_fit, y_var = gp.predict(lc['mag'] - np.mean(lc['mag']), x_fit, return_var=True) #added here the "back to 0"
        y_fit       += np.mean(lc['mag'])

        gp_fit = Table([x_fit, y_fit, y_var, [filtre]*len(x_fit)], names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))

        return gp_fit




def gp_interp_lc_linear_mean(lc, len_kernel , filtre,  mean_param, length_interp = 1000 ):
    '''
    This function performs the guassian process interpolaiton based on the kernel size fitted before
    
    parameters
    ----------
    lc            [astropy table] magnitude light curve (tfromexplo_zc, mag, e_mag)
    len_kernel    [array/list]    The length of the kernel (squared exponential). Should be popt[0] if you fit the autocorrelation of the lightcurve prior 
    filtre        [string]        filter of the LC interp
    length_interp [int] how many data points to interpolate between min and max of th light curve. 
    min_param     [list] list of parameters to include the mean in the function [a,b,peakt]
    
    returns
    -------
    gp_fit  [astropy table] the GP interpolated fit
    
    '''
    class linear_mean_decline(george.modeling.Model):

        parameter_names = ("a", "b", "peakt")

        def get_value(self, t):
            return self.b + self.a * (t-self.peakt)


    
    mean_lin = linear_mean_decline(*mean_param)

    # time scale over which we're interpolating
    x_fit  = np.linspace(min(lc['tfromexplo_zc']),max(lc['tfromexplo_zc']),length_interp)

    # introducing the parameters of the kernel we fit earlier
    kernel =  np.var(lc['mag']) * george.kernels.ExpSquaredKernel(len_kernel) # not happy yet with this expression... 

    # generating and computing the GP
    gp     = george.GP(kernel , mean=mean_lin)
    gp.compute(lc['tfromexplo_zc'], yerr=lc['emag']*2)

    # applying the GP to the LC
    
    # , hence we need to put the data arounbd a mean of zero. We add back the mean of the LC afterwards
    y_fit, y_var = gp.predict(lc['mag'] , x_fit, return_var=True) 

    gp_fit = Table([x_fit, y_fit, y_var, [filtre]*len(x_fit)], names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))
    return gp_fit



#########################

### FUNCTION FOR THE GLOBAL COLOR EVOLUTION 

def error_on_col_ev(e_g, e_r):
    '''
    this function calculates the error on the color evolution based on the errors extracted from the gaussian process
    
    paramaters
    ----------
    e_r 
    e_g 
    
    
    returns
    -------
    error on the color ev 
    '''
    
    return np.sqrt(e_r**2 + e_g**2)

def col_ev(g,r):
    '''
    this function returns the color difference g-r
    
    parameters
    ----------
    g    
    r   
    
    
    returns
    ------- 
    color evolution 
    '''
    
    return g-r



############### CHI_2 test of the LC

def chi_2_lc(lc, gp_lc, e_lc ):
    '''
    This function computes the Chi_2 of the light curve vs gaussian processes
    Also returns the redisuals of the function

    parameters
    ----------
    lc      magnitude lightcurve
    gp_lc   gaussian process magnitude interpolaiton
    e_lc    error on the magnitude lightcurve

    returns
    -------
    residuals [array] 
    Chi_2     [float]
    
    '''

    residuals = lc - gp_lc

    pull_     = ((residuals)/ e_lc )**2
    chi_sq    = np.sum(pull_)
    

    return residuals, chi_sq



def perform_gp_plot_full(lc, filter):
    '''
    This function was made to facilitate the use of Jupyter notebooks. 

    parameters
    ----------
    lc        [table] magnitude light curve with  time from rxplosion
    filter    [string] band to treat, 'g' or  'r'




    returns
    -------
    Figure of the GP with details on the autocorrelation
    
    '''
    lc_b = lc[lc['filter']== filter] 
    bined_lc_b = bin_LC(lc_b,4,150,plot=True)

    if len(lc_b) > 5: 

        t_b, mag_b     = interpolate_LC(bined_lc_b['tfromexplo_zc'], bined_lc_b['mag'])
        autoco_b       = autocorrelate_LC(t_b,mag_b, plot=False, table=True)
        autoco_b       = autoco_b[autoco_b['autocorrelation']>=0]
        popt_b, pcov_b = curve_fit(sq_exp, autoco_b['lags'], autoco_b['autocorrelation'])

        # The gaussian proc
        gp_lc_b       = gp_interp_lc(lc_b, popt_b[0], filter)


        # TESTING THE INTERPOLATION CHI_2

        gp_lc_b_interp       = gp_interp_lc(lc_b, popt_b[0], filter, mode = 'residual')
        residual_b, chi2_b   = chi_2_lc(lc_b['mag'], gp_lc_b_interp['gp_mag'], lc_b['emag'])

        gp_lc_b_interp['residuals'] = residual_b



        ###### PLOTTING
        fig_b   = plt.figure(figsize=(6*np.sqrt(2),6))
        grid    = plt.GridSpec(3, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)
        interp_b     = fig_b.add_subplot(grid[0,  0])
        autocop_b    = fig_b.add_subplot(grid[0,  1])
        residuals_b  = fig_b.add_subplot(grid[1, :2]) 
        gp_interp_b  = fig_b.add_subplot(grid[2, :2])

        # Binning an interpolation

        interp_b.plot(lc_b['tfromexplo_zc'], lc_b['mag'],'x', color = 'grey', alpha = 0.4 )
        interp_b.plot(bined_lc_b['tfromexplo_zc'], bined_lc_b['mag'],'.', color = 'black' )
        if filter == 'g':
            interp_b.plot(t_b, mag_b,'g-') 
        elif filter == 'r':
            interp_b.plot(t_b, mag_b,'r-') 

        interp_b.invert_yaxis()
        interp_b.set_xlabel('Time from EED [days]')
        interp_b.set_ylabel('App. Mag')

        # Autocorrelation and fir to the autocorrelation

        autocop_b.plot(autoco_b['lags'], autoco_b['autocorrelation'], 'grey', alpha = 0.35)
        if filter == 'g':
            autocop_b.plot(autoco_b['lags'], sq_exp(autoco_b['lags'], *popt_b), 'g-',  label = f'l={popt_b[0]:.2f}')
        elif filter == 'r':
            autocop_b.plot(autoco_b['lags'], sq_exp(autoco_b['lags'], *popt_b), 'r-',  label = f'l={popt_b[0]:.2f}')
        autocop_b.set_xlabel('Lags [days]')
        autocop_b.set_ylabel('Autocorr coef')
        autocop_b.legend()


        # Residual of GP vs data

        residuals_b.plot(gp_lc_b_interp['gp_time'], gp_lc_b_interp['residuals'], '.', color = 'grey')
        residuals_b.set_ylabel('Residuals')
        residuals_b.text(5,0.2,f'X_2 = {chi2_b:.3f}')


        # GP interp 

        gp_interp_b.errorbar(lc_b['tfromexplo_zc'], lc_b['mag'], lc_b['emag'], fmt='o', ms = 3, color = 'grey')
        if filter == 'g':    
            h = gp_interp_b.plot(gp_lc_b['gp_time'], gp_lc_b['gp_mag'], color = 'darkolivegreen' )
        elif filter == 'r':
            h = gp_interp_b.plot(gp_lc_b['gp_time'], gp_lc_b['gp_mag'], color = 'orangered' )

        gp_interp_b.fill_between(gp_lc_b['gp_time'], gp_lc_b['gp_mag'] + np.sqrt(gp_lc_b['gp_e_mag']), gp_lc_b['gp_mag'] - np.sqrt(gp_lc_b['gp_e_mag']), 
                                alpha=0.5, color=h[0].get_color())

        gp_interp_b.set_xlabel('Days from EED ')
        gp_interp_b.set_ylabel('App. Mag')
        gp_interp_b.invert_yaxis()

        plt.show()
        return gp_lc_b

    else:
        print(f'Not enough data in {filter} band to perform the gp interpolation')
        gp_lc_b = Table(names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))
        return gp_lc_b




def perform_gp_plot_by_parts(lc, filter, peak_time, part ):
    '''
    This function was made to facilitate the use of Jupyter notebooks. 
    Here, we perform the gaussian processes for a specific part of the light curve : the rise or the decline. The separation is based on the peak time. 


    parameters
    ----------
    lc        [table] magnitude light curve with  time from rxplosion
    filter    [string] band to treat, 'g' or  'r'
    part      [string] 'rise' or 'decline'




    returns
    -------
    Figure of the GP with details on the autocorrelation
    
    '''
    lc_b = lc[lc['filter']== filter] 


    if part == 'rise':
        lc_b       = lc_b[lc_b['tfromexplo_zc'] <= peak_time]
        bined_lc_b = lc_b
        _limit_fit = 3

    elif part == 'decline':
        lc_b       = lc_b[lc_b['tfromexplo_zc']> peak_time]
        lwbound    = int(math.floor(peak_time)+1)
        bined_lc_b = bin_LC(lc_b,lwbound,180,plot=True)
        _limit_fit = 5



    if len(lc_b) > _limit_fit: 

        t_b, mag_b     = interpolate_LC(bined_lc_b['tfromexplo_zc'], bined_lc_b['mag'])
        autoco_b       = autocorrelate_LC(t_b,mag_b, plot=False, table=True)
        autoco_b       = autoco_b[autoco_b['autocorrelation']>=0]
        popt_b, pcov_b = curve_fit(sq_exp, autoco_b['lags'], autoco_b['autocorrelation'])

        # The gaussian proc
        gp_lc_b       = gp_interp_lc(lc_b, popt_b[0], filter)


        # TESTING THE INTERPOLATION CHI_2

        gp_lc_b_interp       = gp_interp_lc(lc_b, popt_b[0], 'g', mode = 'residual')
        residual_b, chi2_b   = chi_2_lc(lc_b['mag'], gp_lc_b_interp['gp_mag'], lc_b['emag'])

        gp_lc_b_interp['residuals'] = residual_b



        ###### PLOTTING
        fig_b   = plt.figure(figsize=(6*np.sqrt(2),6))
        grid    = plt.GridSpec(3, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)
        interp_b     = fig_b.add_subplot(grid[0,  0])
        autocop_b    = fig_b.add_subplot(grid[0,  1])
        residuals_b  = fig_b.add_subplot(grid[1, :2]) 
        gp_interp_b  = fig_b.add_subplot(grid[2, :2])

        # Binning an interpolation

        interp_b.plot(lc_b['tfromexplo_zc'], lc_b['mag'],'x', color = 'grey', alpha = 0.4 )
        interp_b.plot(bined_lc_b['tfromexplo_zc'], bined_lc_b['mag'],'.', color = 'black' )
        if filter == 'g':
            interp_b.plot(t_b, mag_b,'g-') 
        elif filter == 'r':
            interp_b.plot(t_b, mag_b,'r-') 

        interp_b.invert_yaxis()
        interp_b.set_xlabel('Time from EED [days]')
        interp_b.set_ylabel('App. Mag')

        # Autocorrelation and fir to the autocorrelation

        autocop_b.plot(autoco_b['lags'], autoco_b['autocorrelation'], 'grey', alpha = 0.35)
        if filter == 'g':
            autocop_b.plot(autoco_b['lags'], sq_exp(autoco_b['lags'], *popt_b), 'g-',  label = f'l={popt_b[0]:.2f}')
        elif filter == 'r':
            autocop_b.plot(autoco_b['lags'], sq_exp(autoco_b['lags'], *popt_b), 'r-',  label = f'l={popt_b[0]:.2f}')
        autocop_b.set_xlabel('Lags [days]')
        autocop_b.set_ylabel('Autocorr coef')
        autocop_b.legend()


        # Residual of GP vs data

        residuals_b.plot(gp_lc_b_interp['gp_time'], gp_lc_b_interp['residuals'], '.', color = 'grey')
        residuals_b.set_ylabel('Residuals')
        residuals_b.text(5,0.2,f'X_2 = {chi2_b:.3f}')


        # GP interp 

        gp_interp_b.errorbar(lc_b['tfromexplo_zc'], lc_b['mag'], lc_b['emag'], fmt='o', ms = 3, color = 'grey')
        if filter == 'g':    
            h = gp_interp_b.plot(gp_lc_b['gp_time'], gp_lc_b['gp_mag'], color = 'darkolivegreen' )
        elif filter == 'r':
            h = gp_interp_b.plot(gp_lc_b['gp_time'], gp_lc_b['gp_mag'], color = 'orangered' )

        gp_interp_b.fill_between(gp_lc_b['gp_time'], gp_lc_b['gp_mag'] + np.sqrt(gp_lc_b['gp_e_mag']), gp_lc_b['gp_mag'] - np.sqrt(gp_lc_b['gp_e_mag']), 
                                alpha=0.5, color=h[0].get_color())

        gp_interp_b.set_xlabel('Days from EED ')
        gp_interp_b.set_ylabel('App. Mag')
        gp_interp_b.invert_yaxis()

        plt.show()
        return gp_lc_b

    else:
        print(f'Not enough data in {filter} band to perform the gp interpolation')
        gp_lc_b = Table(names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))
        return gp_lc_b














