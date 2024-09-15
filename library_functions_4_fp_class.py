# Functions used for the FP class


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

import sfdmap 
import extinction
from PyAstronomy import pyasl
from astropy.coordinates import Distance


from astropy import cosmology

from iminuit import Minuit



import pprint



###############################################

#               FUNCTIONS                     #

###############################################


def peak_mag_finder(table,dmin,dmax,plot=False,expon = 3,**kwargs):
    '''
    This function interpolates close to the peak and estimates the peak flux and peak magnitude.
    The interpolation is made with a third order polynomial
    
    
    parameters
    ----------
    table    [astropytable] reduced forced photometry table
    dmin     [float]        lower day bound to perform the fit
    dmax     [float]        lower day bound to perform the fit
    plot     [bool]         plot the interpolation
    
    returns
    -------
    
    maximax['leflux'][0]    maximum flux
    maximax['lesjours'][0]  day of maximum flux
    magximax                magnitude at max. flux
    
    '''
    
    if len(table[(table['tfromexplo_zc']<=dmax) & (table['tfromexplo_zc']>=dmin) ]) > 0:

        days  = table[(table['tfromexplo_zc']<=dmax) & (table['tfromexplo_zc']>=dmin) ]['tfromexplo_zc']
        flux  = table[(table['tfromexplo_zc']<=dmax) & (table['tfromexplo_zc']>=dmin)]['extcorrforcediffimflux']
        eflux = table[(table['tfromexplo_zc']<=dmax) & (table['tfromexplo_zc']>=dmin)]['extcorrforcediffimfluxunc']

        flux_interpol = np.polyfit(days, flux, expon)
        flux_fit      = np.poly1d(flux_interpol)


        interays      = np.linspace(dmin,dmax,150)



        #### STOCK THE INTERPOLATION IN A TABLE
        linterpo_r = Table([interays,flux_fit(interays)], names = ['lesjours','leflux'])

        #### FIND MAX FLUX
        # This solution comes from stackoverflow, whose idea is to compute some sort of difference between each points, then second derivative estimate the sign and difference again 
        # (some sort of second derivative), pad the array with zeros (r/l) and then extract the index of the maximum 
        seconddev  = np.pad(np.diff(np.sign( np.diff ( linterpo_r['leflux'] ) ) ), (1,1), 'constant')
        
        try:
            index_max  = np.where(seconddev == -2)[0][0] # <= -2 is where you reach a maximum, taking the first maximum
            fluxmax    = linterpo_r['leflux'][index_max]

            maximax    = linterpo_r[linterpo_r['leflux'] == fluxmax ]

            # old way of doing it: 
            # maximax  = linterpo_r[linterpo_r['leflux'] == max(linterpo_r['leflux'])]
            #print("The maximum flux is " , maximax['leflux'][0]  )
            #print("reached at day" , maximax['lesjours'][0] )

            #### CONVERT TO APPARENT MAGNITUDE
            magximax = -2.5*np.log10(maximax['leflux'][0])
            #print("The maximum apparent magnitude is" , magximax  )

            ##### PLOT THE FIT 
            if plot is True:
                #plt.figure()
                plt.errorbar(days, flux, eflux, fmt='o', alpha = 0.015,color = 'grey',**kwargs)
                plt.plot(interays, flux_fit(interays) ,lw = 1.2, label = f'{dmin}-{dmax} d')
                plt.axvline(maximax['lesjours'][0], alpha = 0.03, color = 'grey')
                plt.axhline(maximax['leflux'][0], alpha = 0.01, color = 'blue')

                plt.xlabel('Time from EED')
                plt.ylabel('Flux')

            return maximax['leflux'][0] ,  maximax['lesjours'][0], magximax


        except IndexError:
            print('No maximum was found, here is the extrema matrix: ')
            print(seconddev)
            print(f'The code failed for the following boundaries : {dmin} - {dmax}')

            return None, None, None
            
  
    else:
        return None, None, None 
    
    
    
    
    
    
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
    eflux = abs(eflux)
    flux  = abs(flux)
    return (2.5*eflux)/(np.log(10)*flux)



def read_table_forced_phot(path_name):
        '''
        This function reads the table and sets the column names
        
        
        '''



        colname = ['index', 'field', 'ccdid', 'qid', 'filter', 'pid', 'infobitssci', 'sciinpseeing', 'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms', 
                    'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2', 'diffmaglim', 'zpdiff', 'programid', 'jd', 'rfid', 'forcediffimflux', 'forcediffimfluxunc',
                    'forcediffimsnr', 'forcediffimchisq', 'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr', 'dnearestrefsrc', 'nearestrefmag', 
                    'nearestrefmagunc', 'nearestrefchi', 'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatus']

        # table            = ascii.read(path_name, delimiter = ',')
        table            = ascii.read(path_name, names = colname)

        table.remove_row(0)
        # #print(table)
        # table.meta = {}
        # lenomdescols     = list(table.colnames)
        # #print(table[1])
        # # For some reason, the update of astropy made it so that it would not interpret the first row of the table as a collection of string. so I had to modify this
        # lenewnomsdescols = [x.split(',')[0] for x in table[0] ] 
        # # lenewnomsdescols = [str(x).split(',')[0] for x in table[0] ]
        # table.rename_columns(lenomdescols, lenewnomsdescols)
        
        
        return table
        
def remove_procstatus_61(table):
    '''
    This procstatus causes to have "null" measurements which are preventing the creation of the photometric table
    
    '''

    table = table[table['procstatus'] != '61']
        
    return table 


def procstatus0_only(table):
    '''
    This procstatus causes to have "null" measurements which are preventing the creation of the photometric table
    
    '''

    table = table[table['procstatus'] == '0']
        
    return table 

# def remove_procstatus(table, procnumber = '56'):
#     '''
#     remove the procstatus specified by the user
#     Procstatus are : 56, 57, 58,59,60,61,62,63,64,65,255 (HAS TO BE A STRINF)

#     In the process of creating the light curve we always renove procstatus 57 and 61
#     There is however a chance you might need to remove more. 


#     '''

#     table = table[table['procstatus'] != procnumber]
#     return table

def filter_bad_seeing(table, seeing = 4):
    '''
    This function filters the FP table on the science image seeing value. 
    By default we filter out any measurement with a seeing > 4"

    parameters
    ----------
    table [astropy.table]
    seeing [float] seeing value in as 
    '''

    table = table[table['sciinpseeing']<=4]
    return table

    
def remove_procstatus_57(table):
    '''
    This procstatus causes to have "null" measurements which are preventing the creation of the photometric table
    
    '''
    zaple = []
    zaple = [x for x in table if x['procstatus'][:2] != '57' and x['procstatus'][3:5] != '57'  ]
    if len(zaple) > 0 :
                  
        pable = vstack(zaple)
    
        return pable
    else :
        pable = table[table['procstatus']=='0']
        return pable


def remove_null_measurements(table):
    '''
    This function removes null measurement, independantly of procstatus... 
    (note: in the future, I'd just rely on this and forget about these stupid procstatus which are inconsistent....! )
    '''
    #print(type(table['forcediffimflux'][3]))
    table = table[table['forcediffimflux'] != 'null']
    #table = table[table['forcediffimfluxunc'] != 'null']
    #table = table[table['forcediffimchisq'] != 'null']
    #print(table)
    
    return table

def remove_infobits_nonzero(table):
    '''
    Removes all the entries with infobits different from 0 (which are supposed to be baddly processed?)
    '''
    table = table[table['infobitssci'] == '0.0']
    return table


def filter_StN_ratio(table, threshold):
    '''
    This function filters out on the S/N provided for the difference flux estimation

    parameters
    ----------
    threshold [float] S/N ratio
    '''


    table = table[table['forcediffimsnr'] >= threshold]
    return table
        
def convert_table_float(table):
        '''
        All the entries in the forced photometry table are strings, so we need to convert the relevant ones into floats
        
        '''
            
        table['procstatus'] = [ x[0:2] for x in table['procstatus'] ]
            
        
        list_keys = list(table.colnames)


        #NOTE: sometimes the the measurements are not null, but the nearest reference source within 5" returns nothing so it still gives a null. 
        # Make sure to remove the null flux measurements before converting to nans and flaots


        for lacle in list_keys:
            for _ in range(len(table[lacle])):
                if table[lacle][_] == 'null':
                    table[lacle][_] = 'nan'

        for lacle in list_keys:
            # table[lacle] = ['nan' for x in table[lacle] if x == 'null'] 
            table[lacle] = [float(x) for x in table[lacle]]
        table['index'] = [int(x) for x in table['index']]
        #print(table)
        return table
    
def rise_typeII(jd, a, t_exp, n):
    '''
    This function defines the type of rise we want to fit to the early light curve.

    f(x) = | a(t-T_exp)**n if t >= T_exp
           | 0             if t <= T_exp
    
    parameters
    ----------
    jd     [array] time interval
    a      [float] amplitude/normalisation scale
    T_exp  [float] estimated time of zero flux
    n      [float] fit exponent
    
    
    returns
    -------
    array
    
    '''
  
    hiuv = np.where( jd-t_exp > 0, a*((jd-t_exp)**n) , 0 )

    # rise   = a*((jd+T_exp)**n)
    # #a*np.sign(jd-T_exp)*(np.abs(jd-T_exp)**n)
    # # a*((jd-T_exp)**n)
    # condi0 = ( (jd + T_exp) <= 0 )

    # #condi0 = ((rise) >= 0 )
    # rise   = condi0*rise

    # return rise
    return hiuv



def proba_chi_sq(chi_2_val, N):
    '''
    This function computes the probability to get a value of chi_2 given a set of data 

    parameters
    ----------
    chi_2 [float] computed chi_2 (from function chi_sq) of the fit vs.data
    N     [float] degrees of freedom of the distribution n - m. n is the number of data points and m is the number of parameters 
    returns
    -------
    P(chi_2,N) [float] the probability for chi_2 to have this value (an estimation of the goodness of fit)

    '''
    chi2_pdf_ = lambda chi2_ : (1/gamma(N/2))*(2**(-N/2))*(chi2_**(N/2-1))*np.exp(-(chi2_/2))

    return quad(chi2_pdf_, chi_2_val, np.inf)



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



def interp_spline_baseline(table,t_b):
    '''
    this function interpolates the baseline
    '''
    table = table[table['tfrommarshalfdec'] <= t_b ]
    spl   = UnivariateSpline(table['tfrommarshalfdec'],table['extcorrforcediffimflux'])
    return spl
    
    
def interp_spline_rise(table,t_b, t_r):
    '''
    this function interpolates the rise
    '''
    table = table[(table['tfrommarshalfdec'] >= t_b) & (table['tfrommarshalfdec'] <= t_r) ]
    spl   = UnivariateSpline(table['tfrommarshalfdec'],table['extcorrforcediffimflux'])
    return spl


def interp_spline_rest(table, t_r):
    '''
    this function interpolates the rest of the lightcurve
    '''
    table = table[table['tfrommarshalfdec'] >= t_r ]
    spl   = UnivariateSpline(table['tfrommarshalfdec'],table['extcorrforcediffimflux'])
    return spl




def remove_big_err_point(table, threshold = 5):
    '''
    Removes data points with large error bars

    should implement a sigma clipping method? 
    '''

    thres  = threshold * np.median(table['extcorrforcediffimfluxunc'])

    _table = table[table['extcorrforcediffimfluxunc'] <= thres ]   

    return _table


def create_new_id(table):
    '''
    This function creates a unique identification key for the field/quadrant...etc 
    '''

    new_id = []
    for _ in range(len(table['field'])):
        bli = float(str(math.floor(table['field'][_]))+str(math.floor(table['ccdid'][_]))+str(math.floor(table['qid'][_])))
        new_id.append(bli)
    table['obs_id'] = new_id
    return table
    



def clean_baseline_phot(table):
    '''
    This function is used to clean the baseline of outliers. It looks for the baseline median and std and kicks out points which are 5*std away 
    from the median
    '''
    _keep       = table[table['tfrommarshalfdec'] >= -2.5 ] # want to conserve the information about the rise if it happens earlier 

    _tempbase   = table[table['tfrommarshalfdec'] < -2.5  ]

    if len(_tempbase) != 0:
        _med, _mad  = np.median(_tempbase['forcediffimflux']) , median_abs_deviation(_tempbase['forcediffimflux'])

        _tempbase   = _tempbase[(_tempbase['forcediffimflux'] < _med + 5 * _mad ) & (_tempbase['forcediffimflux'] > _med - 5 * _mad )]
    #print(_tempbase)

        return vstack([_tempbase,_keep])

    elif len(_tempbase) == 0:
        return table



def correct_4_baseline_offsets(table):
    '''
    This function corrects the baseline for any offset for each 

    NOTE: you need to create the unique ID before using this function. see "create_new_id" function

    parameter
    ---------
    table [astrppytable]

    '''
    
    obids       = list(np.unique(table['obs_id']))

    # a very disgusting way to create a new table with the same keys
    
    lacorrectab = []
    # lacorrectab.remove_rows(slice(0,len(table)))
    

    for obs_id in obids: 
        
        _namebase           = 'temp_'+str(math.floor(obs_id)) 
        locals()[_namebase] = table[table['obs_id'] == obs_id] 

        _temp       = locals().get(_namebase)

        _tempbase   = _temp[_temp['tfrommarshalfdec'] < -2.5  ]

        if len(_tempbase) != 0: 
            _med                            = np.median(_tempbase['forcediffimflux']) 
        
            _temp['forcediffimflux'] = _temp['forcediffimflux'] - _med
            
            lacorrectab = vstack([lacorrectab,_temp])


        else: 
            lacorrectab = vstack([lacorrectab,_temp])

    if len(lacorrectab) !=0:
        # print(lacorrectab)
        lacorrectab.sort('jd')
        return lacorrectab

    elif len(lacorrectab) == 0: 
        print('I could not compute a baseline or correct it')
        return table


   

def correct_4_baseline_offset_withoutsep(table):
    '''
    This function corrects the baseline for any offset for each 

    NOTE: you need to create the unique ID before using this function. see "create_new_id" function

    parameter
    ---------
    table [astrppytable]

    '''
    
    
    
    _tempbase   = table[table['tfrommarshalfdec'] < -2.5  ]


        

    if len(_tempbase) != 0: 
        _med                            = np.median(_tempbase['forcediffimflux']) 
    
        table['forcediffimflux'] = table['forcediffimflux'] - _med

        return table
        
        
    else: 
        print('Not enough data to compute baseline shift')
        return table
            
    
    




    
def bin_phot_perday(phot_table,day_min, day_max, sug = 'median' ):
    
    #save_phot = False, dirout = None, returns_binphot = True, candi
    '''
        this functions bins photometry per day (if several measurement per day, like for P48). It is assuming that each measurement ~ G(mu,sigma). For each bin of day, the weighted mean and error on the weighted mean are computed
        The binning is applied on both the flux and the mag measurements

        /!\ WARNING
        -> instrument specific
        -> filter specific 
        -> table need the columns "mag,flux,emag,eflux,real_time" 

        parameters
        ----------
        phot_table [astropy table] photometric table. WARNING: photometric measurement have to be taken with the same filter and same instrument /!\ They shhould also have the column 'real_time'
        day_min    [int]         when you want to start the binning
        day_max    [int]         when you want to stop the binning 
        sug       [string]      'median' or 'mean', this is for the type of binning you want. In the case of outliers, I recommand using 'median'


        /!\ Provide round numbers
        /!\ In the case of a median, the errors are still computed as weighted mean.
       
        returns
        -------
        bin_phot_table_ [astropy table] table containeing the binned photometry
        hThe code also saves the lightcurve in the folder specified by cdirout
        /!\ the returned table contains less keys, essentially the redshift corrected JD and corresponding marshal time as well as the flux and the 
        uncertainty on the flux (corrected for redshift, ISM extinction and all.)
    '''

    
#     filte       = phot_table['filter'][0] 
#     inst        = phot_table['instrument'][0]

  
    # assuming you want to bin only what comes after the rise ... have to rethink the strategy if want to append 
    bin_phot_table_ = phot_table[['jd_zcorr','tfrommarshalfdec', 'extcorrforcediffimflux', 'extcorrforcediffimfluxunc']][ phot_table['tfrommarshalfdec'] < day_min ]
    if sug == 'mean':
        for i in range(day_min,day_max, 1):
            temp_phot_ =   phot_table[['jd_zcorr','tfrommarshalfdec', 'extcorrforcediffimflux', 'extcorrforcediffimfluxunc']][ (i <= phot_table['tfrommarshalfdec']) & (phot_table['tfrommarshalfdec'] < i+1) ]
            if len(temp_phot_) > 0:
                mu_time_     = np.average(temp_phot_['tfrommarshalfdec'])
                mu_jd_       = np.average(temp_phot_['jd_zcorr'])
                mu_flux_     = get_weighted_mean( np.array(temp_phot_['extcorrforcediffimflux']), np.array(temp_phot_['extcorrforcediffimfluxunc']) )
                sigma_flux_  = get_error_on_weighted_mean( np.array(temp_phot_['extcorrforcediffimfluxunc']) )


                bin_phot_table_.add_row([ mu_jd_ ,mu_time_,mu_flux_ , sigma_flux_ ] )

    elif sug == 'median':
        for i in range(day_min,day_max, 1):
            temp_phot_ =   phot_table[['jd_zcorr','tfrommarshalfdec', 'extcorrforcediffimflux', 'extcorrforcediffimfluxunc']][ (i <= phot_table['tfrommarshalfdec']) & (phot_table['tfrommarshalfdec'] < i+1) ]
            if len(temp_phot_) > 0:
                mu_time_     = np.average(temp_phot_['tfrommarshalfdec'])
                mu_jd_       = np.average(temp_phot_['jd_zcorr'])
                mu_flux_     = np.median( np.array(temp_phot_['extcorrforcediffimflux']) )
                sigma_flux_  = get_error_on_weighted_mean( np.array(temp_phot_['extcorrforcediffimfluxunc']) )


                bin_phot_table_.add_row([ mu_jd_ ,mu_time_,mu_flux_ , sigma_flux_ ] )


    return bin_phot_table_





