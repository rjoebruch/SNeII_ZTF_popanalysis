'''
This script is to generate the color evolution curves


'''
from class_fit_bazin_function import * 
from class_fit_brokenPL_iminuit2 import * 
from class_fit_simpleEXPlaw_iminuit2 import *
from class_separated_fp_explotimin import *
from scipy import integrate
import pandas as pd
from scipy.stats import pearsonr
import os
import glob


################################################## LINKS ######################################################################

save_forcedphot  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'

# BROKENPL = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Broken_PL_alpha1_0-2/results/'
# SIMPLEPL = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Simple_PL/results/'

SAVEFIG  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/color_evolution/early_time/'

##################################################################################################################################

################################################## TABLES ########################################################################

peaks_andall  = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_24012022.ascii', delimiter=',')
table_infants = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflash_22122021.csv', delimiter = ',')
standard_SNII = table_infants[(table_infants['rach-classification']=='SN II')|(table_infants['rach-classification']=='SN IIP')]

interp_early_time_info = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/summary_full/interpolation_both_bands_24012022.ascii'
                                   , delimiter = ',')


###################################################################################################################################

#################################################### FUNCTIONS ####################################################################
def generate_LC(candiname, filt, peak):
    '''
    This function generates the LC 

    parameters
    ----------
    candiname [string] 
    filt      [string]
    '''

    
    _filt = f'ZTF_{filt}'

    band         = ForcedPhot(la_force_ohoto, _filt , candiname )
    band.correct_lc(correct_unc = True, correct_z = True, correct_ext = True, add_jd_f_texp= True)
    band.add_magnitudes_detvndet_meas()
    band.table   = band.table[(band.table['mag'] <= 21.)|(band.table['mag'] >= 99.)] 
    # band.plot_maglc()
    photable = band.table['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']
    photable.add_row([peak['pday'], peak['pflux'], peak['e_pflux'],peak['papmag'],peak['e_papmag'],peak['pabmag'], peak['e_pabmag']])
    photable.sort('tfromexplo_zc')

    detec    = photable[photable['mag']!=99.0]
    detec    = detec[detec['tfromexplo_zc']>=-0.5]
    return detec


def generate_BPL_fit(detec, candiname, filt, peak_time, peak_flux, zez):
    '''
    This function generates the fit with the broken power law 

    parameters
    ----------
    detec
    candiname
    filt
    peak_time 
    peak_flux 
    zez
    
    returns
    -------
    fitty [class]
    '''

    BROKENPL  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Broken_PL_alpha1_0-2/results/'
    param_BPL = ascii.read(BROKENPL+f'{candiname}_table_earlylc_parameters_{filt}_band.csv', delimiter = ',')

    guess = ( param_BPL['A'], param_BPL['alpha1'] , param_BPL['t_break'], param_BPL['alpha2'], 
             1, peak_time, peak_flux)
    bound = {'A': (0,np.inf), 'alpha1':(0,2), 't_break':(0.1,20), 'alpha2': (0,2) } 
    step  = {'alpha1':0.01, 't_break':0.5, 'alpha2': 0.01} 
    fixed = ['n', 't_peak','peak_val']
    fitty = Fit_broken_PL(detec, 0 , peak_time, peak_flux , filter = filt, zez=zez)
    fitty.fit_minuit(guess, boundaries = bound, fixed = fixed, step_sizes = step)



def generate_SPL_fit(detec, candiname, filt, peak_time, peak_flux, zez):
    '''
    This function generates the fit with the simple power law

    parameters
    ----------
    detec
    candiname
    filt
    peak_time 
    peak_flux 
    zez
    
    returns
    -------
    fitty [class]
    '''
    SIMPLEPL  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Simple_PL/results/'
    param_SPL = ascii.read(SIMPLEPL+f'{candiname}_table_earlylc_parameters_{filt}_band.csv', delimiter = ',')


    guess = ( param_SPL['A'], param_SPL['n'], peak_time, peak_flux)
    bound = {'A': [0,np.inf], 'n':[-2,2] } 
    fixed = ['t_peak','peak_val']
    fitty = Fit_simple_EXP_law(detec, 0 , peak_time, peak_flux,  filter = filt, zez=zez )
    fitty.fit_minuit(guess, bound, fixed = fixed)

    return fitty

def generate_g_r_LC(fittyr, fittyg, peak_time_g):
    '''
    This function generates the color light curve g-r based on the interpolation caculated

    parameters
    ----------
    fittyr [class] 
    fittyg [class]
    peak_time_g [float]

    returns
    -------
    g_r       [table]
    _interp_r []
    _interp_g []
    

    '''

    _min_interp = max(min(fittyg.t_fit), min(fittyr.t_fit))
    # _max_interp = min(max(fittyg.t_fit), max(fittyr.t_fit))

    _interp_r = fittyr.get_fit_table_result(_min_interp,peak_time_g)
    _interp_g = fittyg.get_fit_table_result(_min_interp,peak_time_g)

    t_from_peakg = _interp_g['t']   - peak_time_g
    g_r          = _interp_g['mag'] - _interp_r['mag']
    e_g_r        = np.sqrt(_interp_g['emag']**2 + _interp_r['emag']**2)

    g_r = Table([t_from_peakg, g_r, e_g_r], names = ('t_from_gpeak', 'g-r', 'e_g-r'))

    return g_r, _interp_r, _interp_g

    




###################################################################################################################################

####################################################### MAIN #######################################################################

the_weirdos = []

for candiname in standard_SNII['name']:
    
    canditab = interp_early_time_info[interp_early_time_info['name']==candiname]

    if len(canditab) == 2: 
        meta      = standard_SNII[standard_SNII['name']==candiname]
        zez       = [meta['redshift'], meta['e_redshift']]

        la_force_ohoto = save_forcedphot  +  f'{candiname}_fp.ascii' 

        ## R BAND 
        #### LIGHT CURVE COMPUTATION
        peak_r   = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]
        detec_r  = generate_LC(candiname, 'r', peak_r)
        fitype_r = canditab[canditab['filter']=='r']['fit_type'][0]

        #### FIT COMPUTATION
        if fitype_r == 'brokenPLa2':

            fittyr     = generate_BPL_fit(detec_r, candiname, 'r', peak_r['pday'][0], peak_r['pflux'][0], zez)

        elif fitype_r == 'simplePL':

            fittyr     = generate_SPL_fit(detec_r, candiname, 'r', peak_r['pday'][0], peak_r['pflux'][0], zez)
        

        ## G BAND
        #### LIGHT CURVE COMPUTATION 
        peak_g   = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]
        detec_g  = generate_LC(candiname, 'g', peak_g)
        fitype_g = canditab[canditab['filter']=='g']['fit_type'][0]
        #### FIT COMPUTATION
        if fitype_g == 'brokenPLa2':

            fittyg     = generate_BPL_fit(detec_g, candiname, 'g', peak_g['pday'][0],peak_g['pflux'][0], zez)

        elif fitype_g == 'simplePL':

            fittyg     = generate_SPL_fit(detec_g, candiname, 'g', peak_g['pday'][0],peak_g['pflux'][0], zez)

        try:
            ### COMPUTE THE COLOR CURVE

            peak_time_g = peak_g['pday'][0]

            color, interp_r, interp_g = generate_g_r_LC(fittyr, fittyg, peak_time_g)

            ### PLOT RESULTS 

            fig     = plt.figure(figsize=(6*np.sqrt(2),6))
            grid    = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)

            r_b_interpo    = fig.add_subplot(grid[0,  0])
            g_b_interpo    = fig.add_subplot(grid[0,  1])
            color_interpo  = fig.add_subplot(grid[1, :2])



            r_b_interpo.errorbar(fittyr.t_fit, fittyr.mag, fittyr.e_mag, fmt = 'o', ms=3.5 , color = 'red' )
            r_b_interpo.plot(interp_r['t'] , interp_r['mag'], ls = '--' , color = 'black' )
            r_b_interpo.fill_between(interp_r['t'], interp_r['mag'] - interp_r['emag'], interp_r['mag'] + interp_r['emag'],
                                     color = 'grey', alpha=0.5)
            r_b_interpo.set_xlabel('Time from EED [days]')
            r_b_interpo.set_ylabel('Apparent Magnitude [mag]')
            r_b_interpo.invert_yaxis()


            g_b_interpo.errorbar(fittyg.t_fit, fittyg.mag, fittyg.e_mag, fmt = 'o', ms=3.5 , color = 'green' )
            g_b_interpo.plot(interp_g['t'] , interp_g['mag'], ls = '--' , color = 'black' )
            g_b_interpo.fill_between(interp_g['t'], interp_g['mag'] - interp_g['emag'], interp_g['mag'] + interp_g['emag'],
                                     color = 'grey', alpha=0.5)
            g_b_interpo.set_xlabel('Time from EED [days]')
            g_b_interpo.set_ylabel('Apparent Magnitude [mag]')
            g_b_interpo.invert_yaxis()



            color_interpo.plot(color['t_from_gpeak'] , color['g-r'], ls = '--' , color = 'orange' )
            color_interpo.fill_between(color['t_from_gpeak'], color['g-r'] - color['e_g-r'], color['g-r'] + color['e_g-r'], color = 'orange', alpha=0.5)
            color_interpo.set_xlabel('Time from Peak g ')
            color_interpo.set_ylabel('g-r [mag]')

            plt.savefig(SAVEFIG+f'{candiname}_early_lc_g-r.pdf')

            plt.close()
        except AttributeError:
            print(f'{candiname} was weirdly computed, check it out')
            the_weirdos.append(candiname)

with open('check_early_time_color_bug.txt', 'w') as f:
    for _ in range(len(the_weirdos)):
        f.write(the_weirdos[_])
        f.write('\n')

