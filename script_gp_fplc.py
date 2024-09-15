# python 3.8

'''
This script is automated to perform guassian processes on the ext and z corrected light curves 
Update on the 22/11/2021: I am re-renning this script with the corrected magnitude lightcurves (there was an extinction and redshift). 
                          I'm also increasing the number of data points to consider. 
                          The bining of the lightcurve is also operated one day later. 


Update on the 25/01/2022: We are using this code to only interpolate the light curve after peak. We assume that afterpeak the timescale of change
                          of the light curve are on a slower and "fixed" timescale compared to the early time. 

'''




from functions_4_gp import *
from class_separated_fp_explotimin import *


# tablepath        = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/'
# # forcedphot_lc    = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/case_by_case/lc/'
# save_forcedphot  = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'
# # forcedphot_maglc = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/mag_lc_ex_z_corr_2/'
# save_path        = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/interp_gp/decline_interp-wmean/'


# peaks_andall     = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_24012022.ascii', delimiter=',')
# sample_table     = ascii.read(tablepath + 'RISNeIIfull_radeczeztexpflashpeak_24012022.csv', delimiter = ',')
# declinaison      = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/decline_linear_fit/result_fit_declinelinear_26012022.ascii', delimiter = ',')




tablepath        = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/tables/'
# forcedphot_lc    = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/case_by_case/lc/'
save_forcedphot  = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'
# forcedphot_maglc = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/mag_lc_ex_z_corr_2/'
save_path        = '/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/interp_gp/decline_interp-wmean/'


peaks_andall     = ascii.read('/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_24012022.ascii', delimiter=',')
sample_table     = ascii.read(tablepath + 'RISNeIIfull_radeczeztexpflashpeak_24012022.csv', delimiter = ',')
declinaison      = ascii.read('/Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/decline_linear_fit/result_fit_declinelinear_26012022.ascii', delimiter = ',')






standard_SNII  = sample_table[(sample_table['rach-classification']=='SN II')|(sample_table['rach-classification']=='SN IIP')]

######################################### FUNCTIONS 
def generate_LC(candiname, filt, peak):
    '''
    This function generates the LC 

    parameters
    ----------
    candiname [string] 
    filt      [string]
    '''

    la_force_ohoto = save_forcedphot  +  f'{candiname}_fp.ascii' 


    _filt = f'ZTF_{filt}'

    band         = ForcedPhot(la_force_ohoto, _filt , candiname )
    band.correct_lc(correct_unc = True, correct_z = True, correct_ext = True, add_jd_f_texp= True)
    band.add_magnitudes_detvndet_meas()
    band.table   = band.table[(band.table['mag'] <= 21.)|(band.table['mag'] >= 99.)] 

#     band.plot_maglc()
    photable = band.table['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']

    photable.add_row([peak['pday'], peak['pflux'], 1e-18 ,peak['papmag'],1e-9 ,peak['pabmag'], peak['e_pabmag']])
    photable.sort('tfromexplo_zc')

    detec    = photable[photable['mag']!=99.0]

    detec    = detec[(detec['tfromexplo_zc']>=-0.5)&(detec['tfromexplo_zc']<=100)]

    return detec


################################# MAIN 

kernel_size = Table(names = ('name', 'filter', 'ks', 'e_ks'), dtype = ('S20', 'S20', 'f8','f8'))


for candiname in standard_SNII['name']:

    print(candiname)

    if os.path.isdir(save_path+candiname) is True:
        pass
    else:
        os.mkdir(save_path+candiname)

    peak_r = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]
    peak_g = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]


    ############## RED BAND

    if len(peak_r)>0:
        peaktime_r = peak_r['pday'][0]
        peakmag_r  = peak_r['papmag'][0]

        ### LIGHT CURVE GENERATION 
        detec_r  = generate_LC(candiname, 'r', peak_r)

        ### LOOKING ONLY AT THE LC AFTER PEAK 
        late_detec_r = detec_r[detec_r['tfromexplo_zc'] >= peaktime_r]
        bined_lc_r   = bin_LC(late_detec_r,math.floor(peaktime_r),100,plot=True) # binning after rise
        

        if len(detec_r) > 5  :
            t_r, mag_r     = interpolate_LC(bined_lc_r['tfromexplo_zc'], bined_lc_r['mag'])
            autoco_r       = autocorrelate_LC(t_r,mag_r, plot=False, table=True)    
            autoco_r       = autoco_r[autoco_r['autocorrelation']>=0]

            popt_r, pcov_r = curve_fit(sq_exp, autoco_r['lags'], autoco_r['autocorrelation'])

            kernel_size.add_row([candiname, 'r', popt_r[0], np.sqrt(pcov_r[0]) ])

            # THE GAUSSIAN PROC 
            dec_param_r    = declinaison[(declinaison['name']==candiname)&(declinaison['filter']=='r')]['a'][0]

            mean_param_r   = [dec_param_r, peakmag_r, peaktime_r]
            gp_lc_r        = gp_interp_lc_linear_mean(late_detec_r, popt_r[0], 'r', mean_param_r)

            ##### PLOTTING ##### 
            fig_r   = plt.figure(figsize=(6*np.sqrt(2), 6))

            grid    = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)

            interp_r     = fig_r.add_subplot(grid[0,  0])
            autocop_r    = fig_r.add_subplot(grid[0,  1])
            gp_interp_r  = fig_r.add_subplot(grid[1, :2])

            # Binning an interpolation
            interp_r.plot(detec_r['tfromexplo_zc'], detec_r['mag'],'x', color = 'grey', alpha = 0.4 )
            interp_r.plot(peak_r['pday'][0], peak_r['papmag'][0], '*', ms = 5 ,color = 'mediumpurple', alpha = 0.9 )
            interp_r.plot(bined_lc_r['tfromexplo_zc'], bined_lc_r['mag'],'.', color = 'black' )
            interp_r.plot(t_r, mag_r,'r-') 
            interp_r.invert_yaxis()
            interp_r.set_xlabel('Time from EED (rest frame) [d]')
            interp_r.set_ylabel('App. Mag')

            # Autocorrelation and fir to the autocorrelation
            autocop_r.plot(autoco_r['lags'], autoco_r['autocorrelation'],color = 'grey', alpha = 0.35)
            autocop_r.plot(autoco_r['lags'], sq_exp(autoco_r['lags'], *popt_r), 'r-', label = f'l={popt_r[0]:.2f}')
            autocop_r.set_xlabel('Lags [days]')
            autocop_r.set_ylabel('Autocorr coef')
            autocop_r.legend()


            # GP interp 

            gp_interp_r.errorbar(detec_r['tfromexplo_zc'], detec_r['mag'], detec_r['emag'], fmt='o', ms = 3, color = 'grey')

            h = gp_interp_r.plot(gp_lc_r['gp_time'], gp_lc_r['gp_mag'], color = 'orangered' )
            gp_interp_r.fill_between(gp_lc_r['gp_time'], gp_lc_r['gp_mag'] + np.sqrt(gp_lc_r['gp_e_mag']), gp_lc_r['gp_mag'] - np.sqrt(gp_lc_r['gp_e_mag']), 
                                    alpha=0.5, color=h[0].get_color())

            gp_interp_r.set_xlabel('Days from EED')
            gp_interp_r.set_ylabel('App. Mag')
            gp_interp_r.invert_yaxis()

            # gp_interp_r.legend()

            plt.savefig(save_path+candiname+f'/{candiname}_gp_interp_results_r.pdf')

            plt.close(fig_r)


            # plt.plot(lagt_r,autoco_r, 'red', alpha = 0.15)
        else: 
            print(f'Not enough data in r band to perform the gp interpolation for {candiname}')
            gp_lc_r = Table(names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))
    else:
        print(f'No peak magnitude in r band for {candiname}')
        gp_lc_r = Table(names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))


    ############## GREEN BAND


    if len(peak_g)>0:
        peaktime_g = peak_g['pday'][0]
        peakmag_g  = peak_g['papmag'][0]

        detec_g    = generate_LC(candiname, 'g', peak_g)

        late_detec_g = detec_g[detec_g['tfromexplo_zc'] >= peaktime_g]
        bined_lc_g   = bin_LC(late_detec_g,math.floor(peaktime_g),100,plot=True) # binning after rise

        

        if len(detec_g) >5:
            t_g, mag_g     = interpolate_LC(bined_lc_g['tfromexplo_zc'], bined_lc_g['mag'])
            autoco_g       = autocorrelate_LC(t_g,mag_g, plot=False, table=True)
            autoco_g       = autoco_g[autoco_g['autocorrelation']>=0]


            popt_g, pcov_g = curve_fit(sq_exp, autoco_g['lags'], autoco_g['autocorrelation'])

            kernel_size.add_row([candiname, 'g', popt_g[0], np.sqrt(pcov_g[0]) ])

            # The gaussian proc
            dec_param_g   = declinaison[(declinaison['name']==candiname)&(declinaison['filter']=='g')]['a'][0]

            mean_param_g  = [dec_param_g, peakmag_g, peaktime_g]
            gp_lc_g       = gp_interp_lc_linear_mean(late_detec_g, popt_g[0], 'g', mean_param_g)

            ###### PLOTTING


            fig_g   = plt.figure(figsize=(6*np.sqrt(2),6))

            grid    = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)

            interp_g     = fig_g.add_subplot(grid[0,  0])
            autocop_g    = fig_g.add_subplot(grid[0,  1])
            gp_interp_g  = fig_g.add_subplot(grid[1, :2])

            # Binning an interpolation
            interp_g.plot(detec_g['tfromexplo_zc'], detec_g['mag'],'x', color = 'grey', alpha = 0.4 )
            interp_g.plot(peak_g['pday'][0], peak_g['papmag'][0],'*', ms = 5, color = 'mediumpurple', alpha = 0.9 )
            interp_g.plot(bined_lc_g['tfromexplo_zc'], bined_lc_g['mag'],'.', color = 'black' )
            interp_g.plot(t_g, mag_g,'g-') 
            interp_g.invert_yaxis()
            interp_g.set_xlabel('Time from EED [days]')
            interp_g.set_ylabel('App. Mag')

            # Autocorrelation and fir to the autocorrelation
            autocop_g.plot(autoco_g['lags'], autoco_g['autocorrelation'], 'grey', alpha = 0.35)
            autocop_g.plot(autoco_g['lags'], sq_exp(autoco_g['lags'], *popt_g), 'g-',  label = f'l={popt_g[0]:.2f}')
            autocop_g.set_xlabel('Lags [days]')
            autocop_g.set_ylabel('Autocorr coef')
            autocop_g.legend()

            # GP interp 

            gp_interp_g.errorbar(detec_g['tfromexplo_zc'], detec_g['mag'], detec_g['emag'], fmt='o', ms = 3, color = 'grey')

            h = gp_interp_g.plot(gp_lc_g['gp_time'], gp_lc_g['gp_mag'], color = 'darkolivegreen' )
            gp_interp_g.fill_between(gp_lc_g['gp_time'], gp_lc_g['gp_mag'] + np.sqrt(gp_lc_g['gp_e_mag']), gp_lc_g['gp_mag'] - np.sqrt(gp_lc_g['gp_e_mag']), 
                                    alpha=0.5, color=h[0].get_color())

            gp_interp_g.set_xlabel('Days from EED ')
            gp_interp_g.set_ylabel('App. Mag')
            gp_interp_g.invert_yaxis()

            # gp_interp_g.legend()

            plt.savefig(save_path+candiname+f'/{candiname}_gp_interp_results_g.pdf')

            plt.close(fig_g)

        else:
            print(f'Not enough data in g band to perform the gp interpolation for {candiname}')
            gp_lc_g = Table(names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))
    else:
        print(f'No peak mag in g band for{candiname}')
        gp_lc_g = Table(names = ('gp_time', 'gp_mag', 'gp_e_mag', 'filter'), dtype = ('f8','f8','f8','S20'))

        

    
    
    gp_lc = vstack(gp_lc_r,gp_lc_g)

    ascii.write(gp_lc, save_path+candiname+f'/{candiname}_gp_interpolation_mag.ascii', delimiter = ',')

    
    plt.close('all')

ascii.write(kernel_size, save_path + 'kernel_size_summary.csv', delimiter  = ',')






