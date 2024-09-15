'''
This script is for computing the g-r color evolution given the gaussian processes interpolations

'''





from functions_4_gp import *


tablepath        = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/'
# forcedphot_lc    = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/case_by_case/lc/'
forcedphot_maglc = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/mag_lc_ex_z_corr_2/'
save_path        = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/interp_gp/results/'

sample_table = ascii.read(tablepath + 'RISNeII_radeczeztexp_30062021.csv', delimiter = ',')




for candiname in sample_table['name']:

    print(candiname)

    gp_lc   = ascii.read(save_path+candiname+f'/{candiname}_gp_interpolation_mag.ascii', delimiter = ',')

    gp_lc_r = gp_lc[gp_lc['filter']=='r'] 
    gp_lc_g = gp_lc[gp_lc['filter']=='g'] 

    if len(gp_lc_g['gp_time'])>0 and len(gp_lc_r)>0:
        t_min = min(gp_lc_r['gp_time'][0],gp_lc_g['gp_time'][0])
        t_max = max(gp_lc_r['gp_time'][-1],gp_lc_g['gp_time'][-1])

        