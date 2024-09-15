from class_forced_phot_lvlup_110421 import *

# from matplotlib import rc
# rc('font', **{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

# plt.rcParams['xtick.labelsize']=15
# plt.rcParams['ytick.labelsize']=15

#------------------------Functions-------



# def to_abs_mag(x,dm = r_band._DM):
#     return x - dm

# def to_app_mag(x,dm = r_band._DM):
#     return x + dm


#----------------------/functions-------


# papertab = ascii.read('/Users/r.saturn/PhD/WIP/Infant_Sne/tables/table_final2june_paper.ascii') 
papertab = ascii.read('/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/pythoncode/sample_construction/table/test_subsample_2019_RIinfant_rate.ascii') 
# realinfants  = ascii.read('/Users/r.saturn/PhD/WIP/Infant_Sne/tables/final_feb27_rinfantsne2_fromexcel.ascii')
# phath = '/Users/r.saturn/PhD/WIP/Infant_Sne/forcedphot_infsne2018/forced_photometry_RIsample_final/resub_feb20/named/' 
phath = '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/subsamp_testRI_fplc/' 

discomagtab = Table(names=('name','appmagdisco','jddisco','filter'), dtype = ('S20','f8','f8','S20'))


for candiname in papertab['name']:

    

    la_force_ohoto = phath  +  f'{candiname}_forcedphotometry.ascii'

    # plt.figure(figsize=(6*np.sqrt(2),5.4))
    # ax = plt.subplot(111)

# R BAND
    r_band       = ForcedPhot(la_force_ohoto,'ZTF_r', candiname )


    if len(r_band.table) != 0 :
        r_band.correct_lc(only_meas = False)
        r_band.clean_baseline()
        r_band.correct_baseline_offset()
        mag_rband = r_band.compute_nd_vs_det()
        det_rband = mag_rband[mag_rband['mag']!=99.]
        det_rband['filter'] = 'r'
        #_filtyr   = ['r']*len(det_rband['mag'])
        #det_rband.add_column('r', name = 'filter')
        # ndet_rband = mag_rband[mag_rband['mag']==99.]


        # these functions are useful for plotting the secondary axis. 
        # def to_abs_mag(x,dm = r_band._DM):
        #     return x - dm

        # def to_app_mag(x,dm = r_band._DM):
        #     return x + dm

        # PLOT and save the lightcurve
        
        # ax.errorbar(det_rband['tfrommarshalfdec'], det_rband['mag'], det_rband['emag'], fmt = 'o', ms = 5,color = 'orangered', label = 'ZTF-r')
        # ax.plot(ndet_rband['tfrommarshalfdec'],ndet_rband['limmag'], lw = 0, marker = 'v', alpha = 0.2, color = 'orangered')
        # ax2 = ax.secondary_yaxis('right', functions=(to_abs_mag, to_app_mag))
        

    
# G BAND
    g_band = ForcedPhot(la_force_ohoto,'ZTF_g', candiname )

    if len(g_band.table) != 0 : 
        g_band.correct_lc(only_meas = False)
        g_band.clean_baseline()
        g_band.correct_baseline_offset()
        mag_gband = g_band.compute_nd_vs_det()
        det_gband = mag_gband[mag_gband['mag']!=99.]

        det_gband['filter'] = 'g'
        #_filtyg   = ['g']*len(det_gband['mag'])
        #det_gband.add_column('g', name = 'filter')
        # ndet_gband = mag_gband[mag_gband['mag']==99.]

        # def to_abs_mag2(x,dm = g_band._DM):
        #     return x - dm

        # def to_app_mag2(x,dm = g_band._DM):
        #     return x + dm

        # PLOT and save the lightcurve
        # plt.figure(figsize=(6*np.sqrt(2),5.4))
        # ax = plt.subplot(111)
        # ax.errorbar(det_gband['tfrommarshalfdec'], det_gband['mag'], det_gband['emag'], fmt = 'o', ms = 5,color = 'darkolivegreen', label = 'ZTF-g')
        # ax.plot(ndet_gband['tfrommarshalfdec'],ndet_gband['limmag'], lw = 0, marker = 'v', alpha = 0.2, color = 'darkolivegreen')
        # ax2 = ax.secondary_yaxis('right', functions=(to_abs_mag2, to_app_mag2))

    if len(det_rband['mag'])!=0 :
        if  len(det_gband['mag'])!=0 : 
            _fullmag = vstack([det_rband,det_gband])
            print(_fullmag)
        elif len(det_gband['mag'])==0 : 
            _fullmag = det_rband
            print(_fullmag)
    elif len(det_rband['mag'])==0 : 
        if  len(det_gband['mag'])!=0 : 
            _fullmag = det_gband
            print(_fullmag)
        else :
            print(f'There is no data for {candiname}')
            pass
    
    if len(_fullmag)!= 0:
        _fullmag.sort('jd') 

        discomagtab.add_row([candiname, _fullmag['mag'][0], _fullmag['jd'][0], _fullmag['filter'][0]])

print(discomagtab)

    
ascii.write(discomagtab, '/Users/r.saturn/PhD/WIP/Infant_Sne/tables/discovery_mag_fullsample_fp.ascii' )
    # ax.set_ylabel('Apparent Magnitude', size = 15)
    # ax2.set_ylabel('Absolute Magnitude', size = 15)
    # plt.axvline(0,alpha = 0.3, ls = ':', lw = 0.75)
    # plt.xlabel('Days from A.S. first detection', size = 15)
    # plt.xlim([-10,50])
    # plt.ylim([16.5,21.5])
    # plt.text(0.5, 21, candiname,  fontsize = 15)

    # plt.legend(fontsize = 13)
    # plt.gca().invert_yaxis()
    # plt.savefig(f'/Users/r.saturn/PhD/WIP/Infant_Sne/figures/LC_mag_FP_4mid/{candiname}_COMBINED_zoomedx2.pdf')
    # plt.close()
    











