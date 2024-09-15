ls# Python 3.8
'''
This script is using the forced photometry class (class_separated_fp_explotimin.py) and the class to fit the broken power law (class_fit_brokenPL.py) 
to fit the early light curve of SN II. The fit is from EED until EPM (estimated peak magnitude) which was estiamted using the interactive technique 
(see script_mag_finder_zextcorr_interactive.py)

'''
from class_fit_brokenPL_iminuit2 import * 
from class_separated_fp_explotimin import *



######################################################### PATHS ########################################################################################

# table_sample   = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeak_24012022.csv', delimiter = ',')
# peaks_andall   = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_01022022.ascii', delimiter=',')


# save_forcedphot  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'
# TAB_SAVE         = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Broken_PL_alpha1_0-2/results/'
# FIG_SAVE         = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Broken_PL_alpha1_0-2/figures/'


table_sample   = ascii.read('/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeak_24012022.csv', delimiter = ',')
peaks_andall   = ascii.read('/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_01022022.ascii', delimiter=',')


save_forcedphot  = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'
TAB_SAVE         = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Broken_PL_alpha1_0-2/results/'
FIG_SAVE         = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/Broken_PL_alpha1_0-2/figures/'



######################################################### FUNCTIONS ########################################################################################

def get_fit_results(detection_table, peak_time, peak_flux, filter , candiname):
    '''
    This function calls the fitting class and returns the fit results. 
    It also plots the fit and saves it 

    parameters
    ----------
    returns
    -------
    
    '''
    _reduced_to_peak = detection_table[detection_table['tfromexplo_zc']<=peak_time]

    if len(_reduced_to_peak['tfromexplo_zc']) > 4  : # there are 4 parameters over which we are optimising, so there cannot be less datapoints than parameters

        if min(detection_table['tfromexplo_zc'])<=2.5:

            if candiname != 'ZTF18aatlfus': ## ZTF18aatlfus has weird behaviour at early time and cannot be fit in g band either... 

                inter       = interpolate_detec_LC(detection_table['tfromexplo_zc'], detection_table['extcorrforcediffimflux'])
                # guess_break = inter[(inter['flux']<=0.85*peak_flux)&(inter['days']<=peak_time)]['days'][-1] 
                # the gussing of the break is located at 75% of the max_flux 

                # guess = [peak_flux/2, 0.8 , 5, 0.2, 1, peak_time, peak_flux ]
                # guess = [peak_flux/10, 0.8 , 5, 0.2, 1, peak_time, peak_flux ]
                
                guess = ( peak_flux/10, 5 , peak_time-2 , 2, 1, peak_time, peak_flux)

                # bound = {'A': [0,np.inf], 'alpha1':[0,10], 't_break':[0,50], 'alpha2': [0,10] } 

                bound = {'A': (0,np.inf), 'alpha1':(0,5), 't_break':(0.1,25), 'alpha2': (0,5) } 

                fixed = ['n', 't_peak','peak_val']
                step  = {'alpha1': 0.01, 't_break': 0.5, 'alpha2': 0.01} 



                fitty = Fit_broken_PL(detection_table, 0 , peak_time, peak_flux, filter = filter )
                fitty.fit_minuit(guess, boundaries = bound, fixed = fixed, step_sizes = step)
                fitty.plot_fit(add_plot=None)
                

                plt.savefig(FIG_SAVE+f'{candiname}_early_lc_fit_{filter}.pdf')
                plt.close()

                fitty.plot_fit_conv_mag()

                plt.savefig(FIG_SAVE+f'{candiname}_early_lc_fit_mag_{filter}.pdf')
                plt.close()


                _params = [fitty.minuit_output.params[x].value for x in range(7)]
                _errors = [fitty.minuit_output.params[x].error for x in range(7)]

                kishta = fitty.get_redchi_2_fit(_params)

                fit_early_rise = {
                    'name'       : [candiname],
                    'filter'     : [filter],
                    'max_fit'    : [peak_time],
                    'A'          : [_params[0]],
                    'e_A'        : [_errors[0]],
                    'alpha1'     : [_params[1]],
                    'e_alpha_1'  : [_errors[1]],
                    't_break'    : [_params[2]],
                    'e_t_break'  : [_errors[2]],
                    'alpha2'     : [_params[3]],
                    'e_alpha_2'  : [_errors[3]],
                    'n'          : [_params[4]],
                    'e_n'        : [_errors[4]],
                    'chi2fit'    : [kishta]
                 }

                return fit_early_rise
            elif candiname == 'ZTF18aatlfus':
                print('######### This is ZTF18aatlfus, which is unfittable with this law ##########')

                fit_early_rise = {
                            'name'       : [candiname],
                            'filter'     : [filter],
                            'max_fit'    : [9999],
                            'A'          : [9999],
                            'e_A'        : [9999],
                            'alpha1'     : [9999],
                            'e_alpha_1'  : [9999],
                            't_break'    : [9999],
                            'e_t_break'  : [9999],
                            'alpha2'     : [9999],
                            'e_alpha_2'  : [9999],
                            'n'          : [9999],
                            'e_n'        : [9999],
                            'chi2fit'    : [9999],
                                    }
                return fit_early_rise
        else:
            print('Data starts more than 2 days from the EED')
            fit_early_rise = {
                'name'       : [candiname],
                'filter'     : [filter],
                'max_fit'    : [9999],
                'A'          : [9999],
                'e_A'        : [9999],
                'alpha1'     : [9999],
                'e_alpha_1'  : [9999],
                't_break'    : [9999],
                'e_t_break'  : [9999],
                'alpha2'     : [9999],
                'e_alpha_2'  : [9999],
                'n'          : [9999],
                'e_n'        : [9999],
                'chi2fit'    : [9999],
            }
        return fit_early_rise
    else: 
        print('Not enough data at early time to fit')
        fit_early_rise = {
            'name'       : [candiname],
            'filter'     : [filter],
            'max_fit'    : [9999],
            'A'          : [9999],
            'e_A'        : [9999],
            'alpha1'     : [9999],
            'e_alpha_1'  : [9999],
            't_break'    : [9999],
            'e_t_break'  : [9999],
            'alpha2'     : [9999],
            'e_alpha_2'  : [9999],
            'n'          : [9999],
            'e_n'        : [9999],
            'chi2fit'    : [9999],
        }
        return fit_early_rise

##################################################### MAIN ###############################################################



# keys = ['max_fit', 'A','e_A','alpha1','e_alpha_1','t_break','e_t_break','alpha2','e_alpha_2','n','e_n','chi2fit']

revisit_later = Table(names= ('name', 'filter'), dtype = ('S20', 'S20'))

standard_SNII = table_sample[(table_sample['rach-classification']=='SN II')|(table_sample['rach-classification']=='SN IIP')]



for candiname in standard_SNII['name']:

    # if candiname == 'ZTF18aatlfus':
    #     pass
    # else: 

    print(f'##########################{candiname}##########################')
    la_force_ohoto = save_forcedphot  +  f'{candiname}_fp.ascii'

    if len(peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]) != 0 :

        print('####################### R BAND ###############################')

        peak_r = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]

        r_band         = ForcedPhot(la_force_ohoto,'ZTF_r', candiname )

        if len(r_band.table)>0:

            r_band.correct_lc(correct_unc = True, correct_z = True, correct_ext = True, add_jd_f_texp= True)
            r_band.add_magnitudes_detvndet_meas()

            r_band.table = r_band.table[(r_band.table['mag'] <= 21.)|(r_band.table['mag'] >= 99.)] 
            photable_r   = r_band.table['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']

            ## Adding artificially the peak flux in the table with very low errors to force it to be bound at peak in the end
            # photable_r.add_row([peak_r['pday'], peak_r['pflux'], 1e-18, peak_r['papmag'], 1e-12 ,peak_r['pabmag'], peak_r['e_pabmag']])
            photable_r.sort('tfromexplo_zc')
            
            #keeping only the detections for the fit

            detec_r    = photable_r[photable_r['mag']!=99.0]
            detec_r    = detec_r[detec_r['tfromexplo_zc']>=0]

            peak_fluxr = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]['pflux'][0]
            peak_timer = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]['pday'][0]


            try:
                d_r = get_fit_results(detec_r, peak_timer, peak_fluxr, 'r', candiname)

                table_earlyfit_r = pd.DataFrame.from_dict(d_r)
                table_earlyfit_r.to_csv(TAB_SAVE+f'{candiname}_table_earlylc_parameters_r_band.csv',
                                        sep = ",", index=False, 
                                        columns=('name','filter','max_fit','A','e_A','alpha1','e_alpha_1','t_break','e_t_break','alpha2','e_alpha_2','n','e_n','chi2fit'), 
                                        encoding= 'ascii')

                del d_r, table_earlyfit_r
            except TypeError:
                print('############################# ERROR ####################################')
                print('The covariance did not compute. Try this candidate again later ')
                revisit_later.add_row([candiname,'r'])


        else: 

            print(f'There is no peak magnitude for {candiname} in r band')


    if len(peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]) != 0 :

        print('####################### G BAND ###############################')
        
        peak_g = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]

        g_band         = ForcedPhot(la_force_ohoto,'ZTF_g', candiname )

        if len(g_band.table)>0:

            g_band.correct_lc(correct_unc = True, correct_z = True, correct_ext = True, add_jd_f_texp= True)
            g_band.add_magnitudes_detvndet_meas()
            
            g_band.table = g_band.table[(g_band.table['mag'] <= 21.)|(g_band.table['mag'] >= 99.)] 
            photable_g   = g_band.table['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']

            ## Adding artificially the peak flux in the table with very low errors to force it to be bound at peak in the end
            # photable_g.add_row([peak_g['pday'], peak_g['pflux'], 1e-18, peak_g['papmag'], 1e-12 ,peak_g['pabmag'], peak_g['e_pabmag']])
            photable_g.sort('tfromexplo_zc')
            
            #keeping only the detections for the fit

            detec_g    = photable_g[photable_g['mag']!=99.0]
            detec_g    = detec_g[detec_g['tfromexplo_zc']>=-1]

            peak_fluxg = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]['pflux'][0]
            peak_timeg = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]['pday'][0]

            try:
                d_g = get_fit_results(detec_g, peak_timeg, peak_fluxg, 'g', candiname)


                table_earlyfit_g = pd.DataFrame.from_dict(d_g)
                table_earlyfit_g.to_csv(TAB_SAVE+f'{candiname}_table_earlylc_parameters_g_band.csv', 
                                        sep = ",", index=False, 
                                        columns=('name','filter','max_fit','A','e_A','alpha1','e_alpha_1','t_break','e_t_break','alpha2','e_alpha_2','n','e_n','chi2fit'), 
                                        encoding= 'ascii')

                del d_g, table_earlyfit_g
            except TypeError:
                print('############################# ERROR ####################################')
                print('The covariance did not compute. Try this candidate again later ')
                revisit_later.add_row([candiname,'g'])
        else:

            print(f'There is no peak magnitude for {candiname} in g band')

print(f'############################### THERE ARE {len(revisit_later)} TO INVESTIGATE #####################################')
print(revisit_later)
ascii.write(revisit_later, '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/early_lc_fit/revisit_candidate_with_BPL_bug.ascii', delimiter = ',')