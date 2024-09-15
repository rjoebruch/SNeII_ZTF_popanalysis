'''

This script generates the decline fit of standard SN II LC
This version systematically fits until x days, there's no option for interactiveness 

Note that below 30 days, some candidates fail 
ZTF20av..moak should not be included in these fit, cause the amount of data is too poor 
TODO: identify those with poor data and that should be kicked out 
> see the TODO LIST in the class file: to implement before next run
> we need a consistent definition... the data points should all fit the line for a set amount of time. 
> maybe use a binned light curve... 


'''



from class_Linear_decline import *
from class_separated_fp_explotimin import *


from scipy import integrate
import pandas as pd
from scipy.stats import pearsonr
import os
import glob


####################### PATHS AND TABLE FILES ######################################################################################################

# save_forcedphot  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'

# peaks_andall   = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_01022022.ascii', delimiter=',')

# table_infants  = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeak_24012022.csv', delimiter = ',')

# standard_SNII  = table_infants[(table_infants['rach-classification']=='SN II')|(table_infants['rach-classification']=='SN IIP')]



# SAVE           = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/decline_linear_fit/'



save_forcedphot  = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'

peaks_andall   = ascii.read('/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_01022022.ascii', delimiter=',')

table_infants  = ascii.read('/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeakclimb_15022022.csv', delimiter = ',')

standard_SNII  = table_infants[(table_infants['rach-classification']=='SN II')|(table_infants['rach-classification']=='SN IIP')]



SAVE           = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/decline_linear_fit/'

SAVENAME       = 'result_fit_declinelinear_fixed50.ascii'

####################################################################################################################################################



####################### FUNCTIONS ##################################################################################################################

def generate_LC(candiname, filt, peak):
    '''
    This function generates the LC 

    parameters
    ----------
    candiname [string] 
    filt      [string]
    peak      [table]

    '''

    la_force_ohoto = save_forcedphot  +  f'{candiname}_fp.ascii' 


    _filt = f'ZTF_{filt}'

    band         = ForcedPhot(la_force_ohoto, _filt , candiname )
    band.correct_lc(correct_unc = True, correct_z = True, correct_ext = True, add_jd_f_texp= True)
    band.add_magnitudes_detvndet_meas()
    band.table   = band.table[(band.table['mag'] <= 21.)|(band.table['mag'] >= 99.)] 


    photable = band.table['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']

    photable.add_row([peak['pday'], peak['pflux'], peak['e_pflux'], peak['papmag'],peak['e_papmag'] ,peak['pabmag'], peak['e_pabmag']])
    photable.sort('tfromexplo_zc')

    detec    = photable[photable['mag']!=99.0]

    detec    = detec[(detec['tfromexplo_zc']>=-0.5)&(detec['tfromexplo_zc']<=100)]

    # if filt == 'r':
    #     plt.figure()
    #     plt.errorbar(detec['tfromexplo_zc'],detec['mag'],detec['emag'],
    #                  fmt='o', alpha = 0.8 ,ms = 3.5,  color = 'red' )
    #     plt.errorbar(peak['pday'], peak['papmag'], peak['e_papmag'], peak['e_pday'], 
    #                 fmt= '*', alpha = 0.9,ms = 7 ,elinewidth=4,color = 'black')
    #     plt.gca().invert_yaxis()
    #     plt.xlabel('Time from EED [days]')
    #     plt.ylabel('Apparent Magnitude')
    #     plt.title(f'{candiname}')
    #     plt.show()

    # elif filt == 'g':
    #     plt.figure()
    #     plt.errorbar(detec['tfromexplo_zc'],detec['mag'],detec['emag'],
    #                  fmt='o', alpha = 0.8 ,ms = 3.5,  color = 'green' )
    #     plt.errorbar(peak['pday'], peak['papmag'], peak['e_papmag'], peak['e_pday'], 
    #                 fmt= '*', alpha = 0.9,ms = 7 ,elinewidth=4,color = 'black')
    #     plt.xlabel('Time from EED [days]')
    #     plt.ylabel('Apparent Magnitude')
    #     plt.title(f'{candiname}')
    #     plt.gca().invert_yaxis()
    #     plt.show()

    return detec


def generate_linear_fit(late_detec, peak_mag, peak_time, filtre, max_fit = 40 ):
    '''
    Generates the linear fit

    parameters
    ----------
    late_detec [table]
    peak_mag
    peak_time
    filtre
    max_fit
    
    '''
    FIG_SAVE       = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/decline_linear_fit/figures_fixed50/'

    guess = [0.03 , peak_mag, peak_time ]
    bound = {'a': [-0.05,np.inf] } 
    fixed = ['b','t_peak']

    if filtre == 'r' : 
        try:
            fitty = Fit_Linear_Decline(late_detec, peak_time , max_fit ,  filter = filtre )
            fitty.fit_minuit(guess, bound, fixed = fixed)
            fitty.plot_fit(add_plot=None)
            # plt.show()
            plt.savefig(FIG_SAVE+f'{candiname}_decline_fit_{filtre}.pdf')
            
            return fitty
        except TypeError:
            print('Did not go through')
            return None

    elif filtre == 'g':
        try:
            fitty = Fit_Linear_Decline(late_detec, peak_time , max_fit ,  filter = filtre )
            fitty.fit_minuit(guess, bound, fixed = fixed)
            fitty.plot_fit(add_plot=None)
            # plt.show()
            plt.savefig(FIG_SAVE+f'{candiname}_decline_fit_{filtre}.pdf')
            
            return fitty
        except TypeError:
            print('Did not go through...')
            return None

    


####################################################################################################################################################


check_later        = Table(names=('name', 'filter'), dtype = ('S20', 'S20'))

## TODO: change to complete to add the camdidates I had removed 

### WARNING : in a new case, please use the extension "_incomplete" !!!! 

# if os.path.exists(SAVE+ SAVENAME) is True :
    

#     result_fit_decline = ascii.read(SAVE+SAVENAME)

#     print(' Opening incomplete file from folder  ')

# elif os.path.exists(SAVE+f'result_fit_declinelinear_fixed30_incomplete.ascii') is False :

result_fit_decline = Table(names=('name', 'filter', 'a', 'e_a', 'max_fit'), dtype = ('S20', 'S20', 'f8','f8','f8'))

print('Creating new table')

# else :
#     print('No file was found and could not create it')
#     print('Stopping execution of the script')
#     exit()



    

for candiname in standard_SNII['name']:

    if len( result_fit_decline[result_fit_decline['name'] == candiname] ) == 0 :

        print(f'################### {candiname} ################### ')

        peak_r = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]

        if len(peak_r)>0:
            peak_timer = peak_r['pday'][0]
            peak_magr  = peak_r['papmag'][0]

            ### LIGHT CURVE GENERATION 
            print('## R BAND ##')
            print(' PLOTTING LIGHT CURVE, please select zone of interest')

            detec_r  = generate_LC(candiname, 'r', peak_r)
            ### LOOKING ONLY AT THE LC AFTER PEAK 
            # late_detec_r = detec_r[detec_r['tfromexplo_zc'] >= peak_timer]
            print(' Until when would you like to perform the decline fit? ')


            # satisfaction_r_check = 'n'

            # while satisfaction_r_check != 'y': 
            #     max_fit_r = input('Max fit for decline? ')
            fittyr    = generate_linear_fit(detec_r,peak_magr,peak_timer,'r',50)


                # satisfaction_r_check = input('Are you satisfied with the fit? (y/n) ')

                # if satisfaction_r_check == 'y':
                #     print('Moving on...')
                # elif satisfaction_r_check == 'n':
                #     print('Input again the fit parameter')
                # else : 
                #     print('Not understandable, please enter y or n')
                #     satisfaction_r_check = input('Are you satisfied with the fits? (y/n) ')

            if fittyr is not None:
                result_fit_decline.add_row([candiname, 'r', fittyr.minuit_output.params[0].value, fittyr.minuit_output.params[0].error, 40])
            else:
                check_later.add_row([candiname, 'r'])



        else:
            print(' NO PEAK MAG IN R BAND ')

        peak_g = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]

        if len(peak_g)>0:
            peak_timeg = peak_g['pday'][0]
            peak_magg  = peak_g['papmag'][0]

            ### LIGHT CURVE GENERATION 
            detec_g  = generate_LC(candiname, 'g', peak_g)
            ### LOOKING ONLY AT THE LC AFTER PEAK 
            # late_detec_g = detec_g[detec_g['tfromexplo_zc'] >= peak_timeg]

            # satisfaction_g_check = 'n'

            # while satisfaction_g_check != 'y': 
            #     max_fit_g    = input('Max fit for decline? ')
            fittyg       = generate_linear_fit(detec_g,peak_magg,peak_timeg,'g',50)


                # satisfaction_g_check = input('Are you satisfied with the fit? (y/n) ')

                # if satisfaction_g_check == 'y':
                #     print('Moving on...')
                # elif satisfaction_g_check == 'n':
                #     print('Input again the fit parameter')
                # else : 
                #     print('Not understandable, please enter y or n')
                #     satisfaction_g_check = input('Are you satisfied with the fits? (y/n) ')


            if fittyg is not None:
                result_fit_decline.add_row([candiname, 'g', fittyg.minuit_output.params[0].value, fittyg.minuit_output.params[0].error, 40])
            else:
                check_later.add_row([candiname, 'g'])

        else:
            print(' NO PEAK MAG IN R BAND ')




    #     saveloop = True
    #     while saveloop is True:
        
    #         save_peakmag = input('Would you like to save progress? (y/n) ')
    #         if save_peakmag == 'y': 
    #             ascii.write(result_fit_decline, SAVE+'result_fit_declinelinear_04062022_incomplete.ascii', delimiter = ',', overwrite = True)
    #             print('Saved!')
    #             saveloop = False

    #         elif save_peakmag == 'n':
    #             print(' The table was not saved! ') 
    #             saveloop = False

    #         else :
    #             print('Wrong input, save? y/n ')

    #     print(' Moving on! ')

    # elif len( result_fit_decline[result_fit_decline['name'] == candiname] ) > 0 :
    #     print(f'{candiname} was already treated, moving on')
    
    



print(result_fit_decline)


ascii.write(result_fit_decline, SAVE+ SAVENAME, delimiter = ',')
ascii.write(check_later, SAVE+'rchecklater_fixed50.ascii', delimiter = ',')