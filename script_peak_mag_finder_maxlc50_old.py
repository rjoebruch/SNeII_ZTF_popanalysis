from class_separated_fp_explotimin import *


# realinfants  = ascii.read('/Users/r.saturn/PhD/WIP/Infant_Sne/tables/real_infant_fromexcell.ascii')
save_forcedphot = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/sample422/lc/'


peak_mag_summ = Table(names = ('name', 'redshift', 'e_redshift','filter',
                               'pflux','e_pflux', 'papmag', 'e_papmag', 'pday','e_pday'),
                      dtype = ('S20',  'f8' ,'f8','S20', 'f8',  'f8',  'f8',  'f8',  'f8', 'f8'))


finris  = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeII_radeczeztexp_28062021.csv')

# lesub7s = finris[finris['ttofirstspec']<=7]
# lesub7s.info.format = '7.3f'

# for pipi in range(len(lesub7s['iauname'])):
#     #print(lesub7s[pipi]['iauname'][0:2])
#     temp_ = lesub7s[pipi]['iauname']
#     if lesub7s[pipi]['iauname'][0:2] == 'AT':
#         lesub7s[pipi]['iauname'] = 'SN' + temp_[2:]
#         #print(lesub7s[pipi]['iauname'])


PATH_SAVE = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/peak_mag/p3_0-2p5_15-20/'




for candiname in finris['name']:
    print(candiname)
    la_force_ohoto = save_forcedphot  +  f'{candiname}_fp.ascii'

    _tempinfo      = finris[finris['name']==candiname]

    _jd_t_exp_fin  = _tempinfo['jd_t_exp'][0]
    _e_t_exp_fin   = _tempinfo['e_t_exp'][0]



    ### R BAND
    r_band       = ForcedPhot(la_force_ohoto,'ZTF_r', candiname )
    if len(r_band.table)>0: 
        r_band.correct_lc() # this is corrected for MW extinction but not for z 
        r_band.add_magnitudes_detvndet_meas()


        
        # get only the relevant information
        photable_r = r_band.table['jd','tfrommarshalfdec','extcorrforcediffimflux','extcorrforcediffimfluxunc']
        photable_r['jd_fromexplo'] = photable_r['jd'] - _jd_t_exp_fin

        # trim the spec table
        photable_r = photable_r[photable_r['jd_fromexplo']<=50]


        # initialise the collection or peak flux day and magnitudes

        listrf = []
        listrm = []
        listrd = []
        # randomise boundaries
        plt.figure()
        for _ in range(100):

            # _dmin = random.uniform(0,2)
            # _dmax = random.uniform(10,25)

            _dmin = random.uniform(0,2.5)
            _dmax = random.uniform(15, 20)

            #print(f'Lower bound is {_dmin}')
            #print(f'Upper bound is {_dmax}')
            rf, rd, rm = peak_mag_finder(photable_r,_dmin,_dmax, plot=True)
            
            rbpeak     = [rf, rd, rm]

            if rbpeak == [ None, None, None ]: 
                pass
            else :

                listrf.append(rf)
                listrd.append(rd)
                listrm.append(rm)

        plt.savefig(PATH_SAVE+f'{candiname}_rband.pdf')
        plt.close()
        #calculate the average and std of the peak values -> should take the median, not the
        rmf = np.median(listrf)
        rmm = np.median(listrm)
        rmd = np.median(listrd)

        remm = np.std(listrm)
        remf = np.std(listrf)
        remd = np.std(listrd)

        print(f'Maximum flux is {rmf} pm {remf}')
        print(f'Maximum mag is {rmm} pm {remm}')
        print(f'Peak reached at {rmd} pm {remd}')

        peak_mag_summ.add_row([ _tempinfo['name'][0],  _tempinfo['redshift'][0], _tempinfo['e_redshift'][0],'r', rmf, remf,
                       rmm, remm, rmd, np.sqrt(remd**2+_e_t_exp_fin**2)
                              ])        
    

    ####  G BAND 
    g_band = ForcedPhot(la_force_ohoto,'ZTF_g', candiname )
    if len(g_band.table)>0:
        g_band.correct_lc()

        photable_g = g_band.table['jd','tfrommarshalfdec','extcorrforcediffimflux','extcorrforcediffimfluxunc']
        photable_g['jd_fromexplo'] = photable_g['jd'] - _jd_t_exp_fin


        photable_g = photable_g[photable_g['tfrommarshalfdec']<=50]

        listgf = []
        listgm = []
        listgd = []

        plt.figure()
        for _ in range(100):
            # _dmin = random.uniform(0,2)
            # _dmax = random.uniform(10,25)

            _dmin = random.uniform(0,2.5)
            _dmax = random.uniform(15, 20)
            #print(f'Lower bound is {_dmin}')
            #print(f'Upper bound is {_dmax}')

            gf, gd, gm = peak_mag_finder(photable_g,_dmin,_dmax, plot=True)

            gbpeak     = [gf, gd, gm]
            if gbpeak == [ None, None, None ]: 
                pass
            else :
                listgf.append(gf)
                listgd.append(gd)
                listgm.append(gm)


        plt.savefig(PATH_SAVE+f'{candiname}_gband.pdf')
        plt.close()

        gmf = np.median(listgf)
        gmm = np.median(listgm)
        gmd = np.median(listgd)

        gemm = np.std(listgm)
        gemf = np.std(listgf)
        gemd = np.std(listgd)

        print(f'Maximum flux is {gmf} pm {gemf}')
        print(f'Maximum mag is {gmm} pm {gemm}')
        print(f'Peak reached at {gmd} pm {gemd}')

        peak_mag_summ.add_row([ _tempinfo['name'][0],  _tempinfo['redshift'][0], _tempinfo['e_redshift'][0],'g', gmf, gemf,
                       gmm, gemm, gmd,  np.sqrt(gemd**2+_e_t_exp_fin**2)
                      ])




print(peak_mag_summ)
ascii.write(peak_mag_summ, PATH_SAVE+'peak_magnitude_estimation_30062021.ascii')