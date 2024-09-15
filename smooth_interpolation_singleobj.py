from library_4_bazin_functions import *
from class_separated_fp_explotimin import *

'''
We want to generate the smoothed interpolated magnitude light curve of the SN II sample 
Provide the name of transient that sparks your interest in the command line after the python filename 
SINGLE OBJECT - NO SAVE - DISPLAY 
'''


import sys 


print(sys.argv[1])
#######################################################################################################################


save_forcedphot  = '/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'

SAVE_INTERP      = '/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/forced_photometry/SNII_LC_analysis/maglc_snii/'
SAVE_PLOT        = '/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/forced_photometry/SNII_LC_analysis/plot_maglc/'

peaks_andall  = ascii.read('/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/peak_mag/peak_mag_fulinfantsample_01022022.ascii', delimiter=',')
table_infants = ascii.read('/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeakclimb_15022022.csv', delimiter = ',')

standard_SNII = table_infants[(table_infants['rach-classification']=='SN II')|(table_infants['rach-classification']=='SN IIP')]



############################################################################ FUNCTIONS ######################################

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


    if len(band.table)>0:
        
        band.correct_lc(correct_unc = True, correct_z = True, correct_ext = True, add_jd_f_texp= True)
        band.add_magnitudes_detvndet_meas()
        band.table   = band.table[(band.table['mag'] <= 21.)|(band.table['mag'] >= 99.)] 
        # band.plot_maglc()

        photable = band.table['tfromexplo_zc','extcorrforcediffimflux','extcorrforcediffimfluxunc','mag','emag','absmag','e_absmag']
        if len(peak)>0:
            photable.add_row([peak['pday'], peak['pflux'], 1e-18 ,peak['papmag'],1e-9 ,peak['pabmag'], peak['e_pabmag']])
            photable.sort('tfromexplo_zc')

        detec    = photable[photable['mag']!=99.0]
        detec    = detec[(detec['tfromexplo_zc']>=0)&(detec['tfromexplo_zc']<=200)]

        return detec
    else : 
        print(f'No Data in {filt} band')
        return None


############# MAIN 


candiname = sys.argv[1]
#ZTF18aadsuxd



print(f'########### {candiname} ###########')
meta      = table_infants[table_infants['name']==candiname]
zez       = [meta['redshift'], meta['e_redshift']]


print('## R BAND ##')
#### LIGHT CURVE COMPUTATION ###

print(' GENERATING LIGHT CURVE IN R BAND ')
peak_r   = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='r')]
detec_r  = generate_LC(candiname, 'r', peak_r)
    
if detec_r is not None:
    print(detec_r)
    low_bg       = min(detec_r['tfromexplo_zc']) 
    upp_br       = max(detec_r['tfromexplo_zc']) 


    interp_amag_r = interpolate_fullLC_with_errors(detec_r['tfromexplo_zc'], detec_r['absmag'],
                                                    detec_r['e_absmag'],low_bg, upp_br, plot=False)
    smoothed_lc_Test = sg_perpart(interp_amag_r,  peak_T=peak_r['pday'])

    plt.figure()
    plt.errorbar(detec_r['tfromexplo_zc'], detec_r['absmag'], detec_r['e_absmag'],
                 fmt = 'o', ms = 3,alpha = 0.7, color = 'grey' )
    plt.plot(interp_amag_r['t'], interp_amag_r['magn'], color = 'blue', alpha = 0.5)
    i = plt.plot(interp_amag_r['t'], smoothed_lc_Test[0], color = 'red')
    plt.fill_between(interp_amag_r['t'], smoothed_lc_Test[0] + np.sqrt(smoothed_lc_Test[1]), 
                     smoothed_lc_Test[0] - np.sqrt(smoothed_lc_Test[1]), 
                     alpha=0.4, color=i[0].get_color())
    
    plt.gca().invert_yaxis()
    # plt.savefig(f'/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/forced_photometry/SNII_LC_analysis/interpolation/{candiname}_rband_interp.pdf' )
    # plt.close()   
    



elif detec_r is None:
    print(f'No table, no plot for {candiname} in r band ')




print('## G BAND ##')


print(' GENERATING LIGHT CURVE IN G BAND ')
peak_g   = peaks_andall[(peaks_andall['name']==candiname)&(peaks_andall['filter']=='g')]
print(peak_g)

detec_g  = generate_LC(candiname, 'g', peak_g)
    
if detec_g is not None:
    print(detec_g)
    low_bg       = min(detec_g['tfromexplo_zc']) 
    upp_bg       = max(detec_g['tfromexplo_zc']) 

    interp_amag_g = interpolate_fullLC_with_errors(detec_g['tfromexplo_zc'], detec_g['absmag'], detec_g['e_absmag'], low_bg, upp_bg, plot=False)
    smoothed_lc_Test = sg_perpart(interp_amag_g, peak_T=peak_g['pday'])
    plt.figure()
    plt.errorbar(detec_g['tfromexplo_zc'], detec_g['absmag'], detec_g['e_absmag'],
                 fmt = 'o', ms = 3,alpha = 0.7, color = 'grey' )
    plt.plot(interp_amag_g['t'], interp_amag_g['magn'], color = 'blue', alpha = 0.5)
    i =plt.plot(interp_amag_g['t'], smoothed_lc_Test[0], color = 'green' )
    plt.fill_between(interp_amag_g['t'], smoothed_lc_Test[0] + np.sqrt(smoothed_lc_Test[1]), 
                     smoothed_lc_Test[0] - np.sqrt(smoothed_lc_Test[1]), 
                     alpha=0.4, color=i[0].get_color())
    plt.gca().invert_yaxis()
    # plt.savefig(f'/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/forced_photometry/SNII_LC_analysis/interpolation/{candiname}_gband_interp.pdf' )
    # plt.close()
    plt.show()
        
elif detec_g is None:
    print(f'No table, no plot for {candiname} in g band ')
