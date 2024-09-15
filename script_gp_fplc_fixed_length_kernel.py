# python 3.8

'''
This script is automated to perform guassian processes on the ext and z corrected light curves but this time without fitting the autocorrelation length. 
(23/11/2021): We are choosing a kernel of length 100 days for the autocorroleation length and leave it fixed as it is. 
The variance remains the standard deviation... 

Question: is it worth to perform the gaussian process in flux space and only then convert to magnitude or to perform it directly on the magnitude light curve? 

'''




from functions_4_gp import *


tablepath        = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/'
# forcedphot_lc    = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/case_by_case/lc/'
forcedphot_maglc = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/mag_lc_ex_z_corr_2/'
save_path        = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/interp_gp/results/sq_exp_fixed100_kernel/'

sample_table = ascii.read(tablepath + 'RISNeII_radeczeztexp_30062021.csv', delimiter = ',')



#####################  PARAMETERS ###################


length_kernel = 100

#####################################################


##### PLOTTING ##### 

fig   = plt.figure(figsize=(6*np.sqrt(2), 6))
grid    = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)


######################


for candiname in sample_table['name']:

    print(candiname)

    if os.path.isdir(save_path+candiname) is True:
        pass
    else:
        os.mkdir(save_path+candiname)

    lc   = ascii.read(forcedphot_maglc+f'{candiname}_mag_fp_lc_exzcorr.ascii')

    # remove non detections 
    lc   = lc[lc['mag']!=99.]

    # restrain light curve in time
    lc   = lc[lc['tfromexplo_zc']<=180]
    lc   = lc[lc['tfromexplo_zc']>=-1]

    # attention aux conditions de s'il existe observations dans le filtre, etc... 


    ############## RED BAND
    lc_r = lc[lc['filter']== 'r'] 
    bined_lc_r = bin_LC(lc_r,5,150,plot=True) # binning after rise


    if len(lc_r) >5:
        # t_r, mag_r     = interpolate_LC(bined_lc_r['tfromexplo_zc'], bined_lc_r['mag'])
        # autoco_r       = autocorrelate_LC(t_r,mag_r, plot=False, table=True)    
        # autoco_r       = autoco_r[autoco_r['autocorrelation']>=0]

        # popt_r, pcov_r = curve_fit(sq_exp, autoco_r['lags'], autoco_r['autocorrelation'])


        # THE GAUSSIAN PROC 

        gp_lc_r       = gp_interp_lc(lc_r, length_kernel , 'r')
        


        

        interp_r     = fig_r.add_subplot(grid[0,  0])
        autocop_r    = fig_r.add_subplot(grid[0,  1])
        gp_interp_r  = fig_r.add_subplot(grid[1, :2])

        # # Binning an interpolation
        # interp_r.plot(lc_r['tfromexplo_zc'], lc_r['mag'],'x', color = 'grey', alpha = 0.4 )
        # interp_r.plot(bined_lc_r['tfromexplo_zc'], bined_lc_r['mag'],'.', color = 'black' )
        # interp_r.plot(t_r, mag_r,'r-') 
        # interp_r.invert_yaxis()
        # interp_r.set_xlabel('Time from EED (rest frame) [d]')
        # interp_r.set_ylabel('App. Mag')

        # # Autocorrelation and fir to the autocorrelation
        # autocop_r.plot(autoco_r['lags'], autoco_r['autocorrelation'],color = 'grey', alpha = 0.35)
        # autocop_r.plot(autoco_r['lags'], sq_exp(autoco_r['lags'], *popt_r), 'r-', label = f'l={popt_r[0]:.2f}')
        # autocop_r.set_xlabel('Lags [days]')
        # autocop_r.set_ylabel('Autocorr coef')
        # autocop_r.legend()
        

        # GP interp 

        gp_interp_r.errorbar(lc_r['tfromexplo_zc'], lc_r['mag'], lc_r['emag'], fmt='o', ms = 3, color = 'grey')

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
    

    ############## GREEN BAND
    lc_g = lc[lc['filter']== 'g'] 
    bined_lc_g = bin_LC(lc_g,4,150,plot=True) # binning after rise 

    if len(lc_g) >5:
        t_g, mag_g     = interpolate_LC(bined_lc_g['tfromexplo_zc'], bined_lc_g['mag'])
        autoco_g       = autocorrelate_LC(t_g,mag_g, plot=False, table=True)
        autoco_g       = autoco_g[autoco_g['autocorrelation']>=0]


        popt_g, pcov_g = curve_fit(sq_exp, autoco_g['lags'], autoco_g['autocorrelation'])

        # The gaussian proc
        gp_lc_g       = gp_interp_lc(lc_g, popt_g, 'g')

        ###### PLOTTING


        fig_g   = plt.figure(figsize=(6*np.sqrt(2),6))

        grid    = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3, bottom = 0.1, top = 0.9)

        interp_g     = fig_g.add_subplot(grid[0,  0])
        autocop_g    = fig_g.add_subplot(grid[0,  1])
        gp_interp_g  = fig_g.add_subplot(grid[1, :2])

        # Binning an interpolation
        interp_g.plot(lc_g['tfromexplo_zc'], lc_g['mag'],'x', color = 'grey', alpha = 0.4 )
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

        gp_interp_g.errorbar(lc_g['tfromexplo_zc'], lc_g['mag'], lc_g['emag'], fmt='o', ms = 3, color = 'grey')

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
        

    
    
    gp_lc = hstack(gp_lc_r,gp_lc_g)

    ascii.write(gp_lc, save_path+candiname+f'/{candiname}_gp_interpolation_mag.ascii', delimiter = ',')

    
    plt.close('all')


    # plt.plot(lagt_g,autoco_g, 'green', alpha = 0.15)





