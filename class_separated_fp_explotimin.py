from library_functions_4_fp_class import *


#TODO: make this class ONLY a provider of correct photometry. let the epxlosion time and the peak magnitude be satellite functions or classes ? 
#TODO: adapt to the new filescheme




###############################################

#               CLASS                         #

###############################################



class ForcedPhot( object ): 
    '''
    This class defines the object "photometric table" from the forced photometry table. It takes as input the forced
    photometry file and the filter in which you want to obtain the "final" product of the photometric table. 
    For more information, check the pdf file with al the information on forced photometry
    
    '''
    

    ###############################################

    #             GLOBAL VARIABLES                #

    ###############################################

    #  Global variable include the general table which contains all the necessary meta data for the FP lightcurve treatment
    # sample_table =  '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/SNeII_032018-122020_ndlim2p5_zclasscentroidradec.ascii'


    '''
    WARNING: this table has all the metadata needed. It's defined locally in the class... 
    '''
    
    sample_table =  '/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeakclimb_15022022.csv'
    # sample_table =  '/Users/r.saturn/PhD/WIP/Full_RISNeII_ZTF1/tables/RISNeIIfull_radeczeztexpflashpeak_24012022.csv'
    


    # /Users/r.saturn/Dropbox (Weizmann Institute)/PhD/WIP/Full_RISNeII_ZTF1/forced_photometry/case_by_case



    #--------#
    # SETTER #
    #--------#

    
    
    def __init__(self, forced_phot_filename , filtre , name , plots = False ):
        '''
        This function initisalises the class
        
        Parameters
        ----------
        self
        forced_phot_file [path+string] path + name of the file which we have to be treated
        filtre [string] : each object will be discriminated based on the ztf filter of the forced photometry
        name   [string] : object name
        plots  [bool]-optional- Whether you want to visualize some plots like the raw LC plot and the flux vs chi2 for the uncertainties. 
        other things can be added afterwards
        -optional- codes of processing status you want removed form the table
                                          ['255','64','56','61'] are the default status filtered out  
        
        Returns
        -------
        The actual things you're going to be working with goddamit
        
        '''
        self.cosmo = cosmology.Planck13 # We are using the Planck 13 cosmology for our calculations 
        
        self._filter_name   = {'ZTF_g':'darkolivegreen','ZTF_r':'orangered', 'ZTF_i':'gold'}
        
        self.set_name(name)
        self.set_metadata()
        
        self.set_rawtable(forced_phot_filename)
        self.set_filtre(filtre)
        
        self.set_table()
    
        self._corrected_LC     = False 
        self._corrected_4_z    = False
        self._corrected_4_ext  = False
        self._converted_to_mag = False
        self._normalised_LC    = False
        self._corrected_unc    = False
        self._interpolatedLC   = False
        self._binned_LC        = False
        
        if plots is  True :
            self.plot_check_unc_chisq()
            self.plot_lc()
        
    
        
        
    def set_rawtable(self, forced_phot_filename):
        '''
        
        '''
        self._rawtable = read_table_forced_phot(forced_phot_filename)
    
    
    
    def set_name(self, name):
        '''
        sets the name of the supernova (important to query other parameters such as redshift and coordinates):
        
        Parameters
        ----------
        
        Returns
        -------
        '''
        
        self.objectname = name
        
       
        
    def set_metadata(self, table = sample_table ):
       

        '''
        Sets the metadata such as redshift and coordinates which will be use afterwards 
        parameters
        ----------
        table [string] path to the table of objexts of list of objects or just a line where there's all the metadata information

        /!\ This table is defined as the global variable of the class. One can change it by setting "self.sample_table" BEFORE working with the LC

        '''
        
        # NOTE: changed the coordinates to the centroid coordinates
        
        _table_infoobj     = ascii.read(self.sample_table)
        _table_infoobj     = _table_infoobj[_table_infoobj['name'] == self.objectname]
        self._redshift     = float(_table_infoobj['redshift'])

        # Which cosmology should we assume? Plack 13.  Why? Cause 18 is fucked up

        self._DM           = self.cosmo.distmod(z = self._redshift).value
        # self._e_DM         = float(_table_infoobj['e_redshift'])*self._DM/self._redshift
        self._e_DM         = (5*_table_infoobj['e_redshift'])/(np.log(10)*_table_infoobj['redshift']) # for small z!! 
        
        self._coordinates  = SkyCoord(np.array(_table_infoobj['ra_med']), np.array(_table_infoobj['dec_med']),unit = 'deg', frame = 'fk5')
        
        #TODO: can change the first detection to the first detection with the forced photometry light curve. 
        self._marshalfdet  = _table_infoobj['First Detection'] # This is given by the full table. 
        #  self._fpfirstdet   = _table_infoobj['First Detection'] 
        try:
            self._jd_texp = _table_infoobj['jd_t_exp'][0]
        except KeyError:
            print('There is no explosion date in this table')
            
                
                           
    
    def set_filtre(self, filtre):
        '''
        sets the filter in which you want to obtain the photometry. the variable is only an entry variable. 
        parameter
        ---------
        filtre [string] name of the filter/ 'ZTF_{g,i,r}'
        
        '''
        self.filtre = filtre
        
        
    def set_table(self):
        '''
        
        This function prepares the rawtable to be worked with in a single band (filter). 
        It removes all the 'null' measurements. 

        #TODO: In a future version, give the procstatus on which to filter. 
        In this version, we implement quality checks on the seeing, bad pixels 
        
        Parameters 
        ----------
        self 
        
        Returns
        -------
        self.table
        '''
        # Selects only measurements in the relevant filters
        self.table = self._rawtable[self._rawtable['filter'] == self.filtre]
        
        # removing the "null" measurements : this is when the fit failed
        self.table = remove_null_measurements(self.table)


        #NOTE: sometimes the the measurements are not null, but the nearest reference source within 5" returns nothing so it still gives a null. 
        # Make sure to remove the null flux measurements before converting to nans and flaots
        if len(self.table['jd']) > 0 : 
            #removing the actual strings from the column which will be of no use for this object
            self.table.remove_columns(['filter']) # string removal, each light curbe is per filter

            self.table = convert_table_float(self.table)

            # Giving each observation a unique ID based on the CCD ID, the Quadrant and the Field IDs. 
            self.table = create_new_id(self.table)

            self.table['tfrommarshalfdec'] = self.table['jd'] - self._marshalfdet
            
            
        elif len(self.table['jd']) == 0 : 
            
            print(f'There does not seem to be any data for {self.objectname} in {self.filtre}')
        

          
    #--------#
    # FILTER #
    #--------#

    def remove_big_errors(self, thres = 7):
        '''
        #IMPORTANT: this function removes a significant amount of data points? THis should be better investigated... 
       
        This function looks for the median of the error of the flux measurements and removes those with errors bigger than 3*median
        
        parameters
        ----------
        thres [float] the threshold of the MAD we remove 


        returns
        -------
        table without the faulty data points

        '''

        self.table = remove_big_err_point(self.table, threshold  = thres)
        # return remove_big_err_point(self.table, threshold  = thres)

    
    
    def reset_table(self):
        '''
        Useful function to reinitialise the table after applying cuts. 

        Warning
        -------
        Remember to correct the lightcurve and normalise it before using the fitter
        
        '''
       
        del self.table
        self.set_table()
        
        self._corrected_LC      = False 
        self._corrected_4_z     = False
        self._corrected_4_ext   = False
        self._normalised_LC     = False
        self._converted_to_mag  = False
        self._corrected_unc     = False
        self._binned_LC         = False
        self._clean_LC          = False
        self._added_jd_from_exp = False

    
    # def filter_procstatus(self, procnumber ):
    #     '''
    #     This function filters the procstatus

    #     parameters
    #     ----------

    #     procnumber [float] 56,57,58 ...etc, see documentation

    #     returns
    #     -------
    #     table filtered from chosen procstatus 
    #     not saved to class
    #     '''

    #     #self.table = remove_procstatus(self.table , procnumber=pn)

    #     return self.table[self.table['procstatus'] != procnumber]
        

    
    # def filter_ccdid(self,ccdid_code):
    #     '''
    #     this function filters out the specified filter 

    #     Parameters
    #     ----------
        
    #     Returns
    #     -------
    #     a filtered table 
    #     Not saved into the class
        
    #     '''

    #     # self.table = self.table[self.table['ccdid'] != ccdid_code]
    #     return self.table[self.table['ccdid'] != ccdid_code]
          
        
        
        
    #--------#
    # GETTER #
    #--------#
    
    def get_baseline(self, binned = False):
        '''
        this functions aims at estimating the baseline of the non detections 
        Baseline correction is important vefore conversion to magnitude / upper limit. 
        Parameters
        ----------
        indiex [int]
        
        
        Returns
        -------
        Value of the median of the baseline and its standard deviation
        '''
        if self._corrected_LC  is True: 
            if binned is False:
                _temp = self.table['extcorrforcediffimflux'][self.table['tfrommarshalfdec'] < -2.5]  
                return np.median(_temp ), np.std( _temp )

            elif binned is True:
                if self._binned_LC is True:
                    _temp = self.binned_lc['extcorrforcediffimflux'][self.binned_lc['tfrommarshalfdec'] <= -2.5] 
                    return np.median( _temp ), np.std( _temp )
                else : 
                    print('You need to bin the LC first')
        else:
            print('You need to correct the LC first')

    
    
    #TODO: get the  band explosion date 
    def get_band_t_exp(self, weights = 'error'):

        '''
        This function gets the weighted mean explosion date. The fits are performed from1.5 to 5.5 days after first detection, with a step of 0.5 days 
        
        parameters
        ----------
        weights [string] 'errors' to weight with the errors, 'chi2' to weight with the chi2 of the fit

        returns
        -------
        band_tx  - weighted mean on the the estimated explosion dates
        band_etx - error on the wieghted mean 
        '''

        texp    = []
        etexp   = []
        chi2fit = []

        guess = [2e-7, -1 , 0.5]
        bound = {'a': [0,np.inf], 't_exp':[-5,0.], 'n':[0,5]} 
        interv    = np.linspace(1.5,5.5,8) 

        kishta = 0.

        for _ in interv:
            if kishta <= 10.:
            
                _temp_fit = Fit_t_exp(self.table, -10 , _ )
                _temp_fit.fit_minuit(guess, bound)

                params = [_temp_fit.minuit_output.params[x].value for x in range(3)]

                kishta      = _temp_fit.get_redchi_2_fit(params)
                _tx, _e_tx  = _temp_fit.get_explosion_t() 
                
                texp.append(_tx)
                etexp.append(_e_tx)
                chi2fit.append(kishta)

        texp    = np.asarray(texp)
        etexp   = np.asarray(etexp)
        chi2fit = np.asarray(chi2fit)


        if weights == 'error':
            
            band_tx  = get_weighted_mean(texp,etexp)
            band_etx = get_std_err_wmean(texp,etexp)
        elif weights == 'chi2':
            band_tx  = get_weighted_mean(texp,chi2fit)
            band_etx = get_std_err_wmean(texp,chi2fit)
        else:
            print('You need to specify the type of weights to estimate the weighted mean and errors. Either error or chi2')

        print(f'tx in this band is {band_tx} ± {band_etx}')

        return band_tx, band_etx

    def get_t_exp_(self, max_fit = 2.5, min_fit = -10, plot_fit = True, add_plot=None, 
                    special_check = False):

        '''
        This function computes the explosion date for a specified time zone 

        parameters
        ---------
        max_fit   [float] upperlimit of the data to fit
        min_fit   [float] lower limit of the data to fit
        plot_fit  [bool] whether you want to plot the fit
        add_plot  [None or ax.] whether you want to verlay several fits
        special_check [bool] bypass the limit on the number of point after the first marshal detection

        returns
        -------
        fit_t_exp [dictionnary] compose of :
                                d       time zon
                                t_exp   explosion time
                                e_t_exp error of the fit on the explosion date
                                n       index of the fit
                                e_n     error on the index of the fit
                                chi2fit chisquare on the zone specifief of the fit

        Returns None if the timezone to fit has less than 1 point

        '''

        guess = [2e-7, -1 , 1]
        bound = {'a': [0,np.inf], 't_exp':[-10,0.], 'n':[0,5]} 
        

    
        _temp_fit = Fit_t_exp(self.table, min_fit , max_fit , filter=self.filtre)
        lenfifi = _temp_fit.get_fit_roi_length()
        # print(lenfifi)
        
        if special_check is False: 
            if lenfifi > 1 :

                try : 
                    _temp_fit.fit_minuit(guess, bound)

                    _params = [_temp_fit.minuit_output.params[x].value for x in range(3)]
                    _errors = [_temp_fit.minuit_output.params[x].error for x in range(3)]

                    kishta = _temp_fit.get_redchi_2_fit(_params)

                    fit_t_exp = {
                        'd_fit'   : max_fit,
                        't_exp'   : _params[1],
                        'e_t_exp' : _errors[1],
                        'n'       : _params[2],
                        'e_n'     : _errors[2],
                        'chi2fit' : kishta
                    }

                    if plot_fit == True:
                        _temp_fit.plot_fit(add_plot=add_plot)

                    return fit_t_exp

                except ValueError :
                    print('Oops, something happened... Better check !')
                    return None
                

            elif lenfifi <=1 :
                print('Not enough datapoint to fit')
                return None

        elif special_check is True: 
            _temp_fit.fit_minuit(guess, bound)

            _params = [_temp_fit.minuit_output.params[x].value for x in range(3)]
            _errors = [_temp_fit.minuit_output.params[x].error for x in range(3)]
            kishta = _temp_fit.get_redchi_2_fit(_params)
            fit_t_exp = {
                'd_fit'   : max_fit,
                't_exp'   : _params[1],
                'e_t_exp' : _errors[1],
                'n'       : _params[2],
                'e_n'     : _errors[2],
                'chi2fit' : kishta
            }

            if plot_fit == True:
                _temp_fit.plot_fit(add_plot=add_plot)

            return fit_t_exp

    def get_ealylc_n_behaviour(self, guess_n = 0.5, RoC = 'concave', plot_fit = True, add_plot=None, 
                    special_check = False):

        '''
        This function computes the explosion date for a specified time zone 

        parameters
        ---------
        max_fit   [float] upperlimit of the data to fit -> fixed to 3 d
        min_fit   [float] lower limit of the data to fit -> fixed to -10 days
        guess_n   [float] guess for the index of fit law
        RoC       [string] Curvature of the fit law 
        plot_fit  [bool] whether you want to plot the fit
        add_plot  [None or ax.] whether you want to verlay several fits
        special_check [bool] bypass the limit on the number of point after the first marshal detection

        returns
        -------
        fit_t_exp [dictionnary] compose of :
                                RoC     curvature of fit index
                                t_exp   explosion time
                                e_t_exp error of the fit on the explosion date
                                n       index of the fit
                                e_n     error on the index of the fit
                                chi2fit chisquare on the zone specifief of the fit

        Returns None if the timezone to fit has less than 1 point

        '''     

        max_fit = 3 
        min_fit = -10

        if RoC == 'concave':
            lowbound_n = 0
            upbound_n  = 1
        elif RoC == 'convex':
            lowbound_n = 1
            upbound_n  = np.inf
        else:
            raise ValueError('You need to specify concave or convex')



        guess = [2e-7, -1 , guess_n]
        bound = {'a': [0,np.inf], 't_exp':[-10,0.], 'n':[lowbound_n,upbound_n]} 
        

    
        _temp_fit = Fit_t_exp(self.table, min_fit , max_fit , filter=self.filtre)
        lenfifi = _temp_fit.get_fit_roi_length()
        # print(lenfifi)
        
        if special_check is False: 
            if lenfifi > 1 :

                try : 
                    _temp_fit.fit_minuit(guess, bound)

                    _params = [_temp_fit.minuit_output.params[x].value for x in range(3)]
                    _errors = [_temp_fit.minuit_output.params[x].error for x in range(3)]

                    kishta = _temp_fit.get_redchi_2_fit(_params)

                    fit_t_exp = {
                        'RoC'     : RoC,
                        't_exp'   : _params[1],
                        'e_t_exp' : _errors[1],
                        'n'       : _params[2],
                        'e_n'     : _errors[2],
                        'chi2fit' : kishta
                    }

                    if plot_fit == True:
                        _temp_fit.plot_fit(add_plot=add_plot)

                    return fit_t_exp

                except ValueError :
                    print('Oops, something happened... Better check !')
                    return None
                

            elif lenfifi <=1 :
                print('Not enough datapoint to fit')
                return None

        elif special_check is True: 
            _temp_fit.fit_minuit(guess, bound)

            _params = [_temp_fit.minuit_output.params[x].value for x in range(3)]
            _errors = [_temp_fit.minuit_output.params[x].error for x in range(3)]
            kishta = _temp_fit.get_redchi_2_fit(_params)
            fit_t_exp = {
                'RoC'     : RoC,
                't_exp'   : _params[1],
                'e_t_exp' : _errors[1],
                'n'       : _params[2],
                'e_n'     : _errors[2],
                'chi2fit' : kishta
            }

            if plot_fit == True:
                _temp_fit.plot_fit(add_plot=add_plot)

            return fit_t_exp

   
    
    #-----------#
    # Corrector #
    #-----------#
    
    def correct_lc(self, correct_unc = True , correct_z = False , correct_ext = True, correct_baseline = True, SNR_thres = 3 ,
                        correct_baseline_method = 'separated', add_jd_f_texp = False):
        '''
        This function corrects the lightcurve for:
        0) remove known bad measurements
        1) recalibrating the uncertainties on the PSF flux measurement as suggested in the FP guide 
        2) correct for extinciton (by default)
        3) can correct for redshift (but this is turned off by default. ) #NOTE: this will be relevant when working with GP for LC parameters

        It can also "correct" the uncertainty which by default is true here. This only means multiplying the raw errors by the sqrt of the chi2
        Correcting the light curve is physical and is based on the measurements and on the redshift, zero point and etc 


        (il faut faire les quality cuts bien avant les corrections)
        
        Parameters
        ----------
        self
        correct_unc       [bool] whether you want to also correct the uncertainties. By default it will be true.  
        correct_z         [bool] whether to correct for redshift
        correct_ext       [bool] whether to correct for extinction
        correct_baseline  [bool] whether to correct for the baseline offset. By default True. 
        SNR_thres         [int] the S/N threshold to filter the SN LC from noisy points (not the baseline)
        
        
        returns
        -------
        corrected light curve (self table)
        
        
        '''
        ############# STEP 1: filter known bad measurement: 
        # Remove measurements with a seeing > 4" 
        self.table = filter_bad_seeing(self.table, seeing = 4)


        # remove high spatial nosie sigma per pixel , by default threshold is 25. Can be stricter.
        self.cuts_scisigpix()


        ############ STEP 2 : correct light curve of physical parameters

        # baseline correction (should be done before the zero point correction, no?)
        # important: should correct for the offset BEFORE cleaning the baseline... !!! 
        
        if correct_baseline is True:
            self.correct_baseline_offset(correction_method=correct_baseline_method)
            self.clean_baseline()
            
            

        # correction for zeropoint

        self.correct_zeropoint()


        # recalibration of uncertainties of diff flux
        
        if correct_unc is True:
            self.correct_uncertainties()
        self.compute_errors_on_final_diffflux()

        
        # Remove low S/N (3) in the Supernova LC (outside of the baseline)
    
        self.filter_on_SNR(SNR_thres=SNR_thres)


        # remove big big errors , there's a bug in this thing... 

        self.remove_big_errors()


        # correction for extinction and redshift

        self.correct_redshift(correct_z=correct_z)
        self.correct_extinction(correct_ext=correct_ext)

        if add_jd_f_texp is True: 
            self.add_t_from_explo()



        # toggles this variable to True to let the user know that corrections were made.
        self._corrected_LC = True

    # def get_peak(self):
    #     '''
    #     This function determines the peak value (in flux) and time of the light curve 
    #     '''

    #     'tfrommarshalfdec'


    def add_t_from_explo(self):
        '''
        If the table we provide has the explosion time already computed, this function adds a column of time from explosion

        parameters
        ----------

        returns
        -------

        
        '''

        if self._corrected_4_z is True: 
            self._jd_texp_zcor          = self._jd_texp/(1+self._redshift)
            self.table['tfromexplo_zc'] = self.table['jd_zcorr'] - self._jd_texp_zcor
            # print(self.table)
            self._added_jd_from_exp     = True
        else:
            self.table['tfromexplo_znc'] = self.table['jd'] - self._jd_texp
            self._added_jd_from_exp = True
        


    ##################
    #    MAG LC      #
    ##################
    

    def add_magnitudes_detvndet_meas(self, return_magtable = False):
        '''
        This function computes the detections vs Non detections. 
        According to the Forced photometry guide, we choose SNT = 3 and SNU = 5. 
        If a the flux divided by the uncertainty is bigger/equal than 3, it is a dectection. If flux/unc < 3, it's a non-detection.
        Limiting magnitudeds are computed as 5 times sigma.

        parameters
        ----------
        return_magtable [bool] whether you want to return only the magnitude table. it will still be added to the main table

        returns
        -------

        note
        ----

        When computing the flux, we already account for the zero point so we don't need to take it into account again.

        '''

        if self._converted_to_mag == False:

            SNT_ = 3
            SNU_ = 5

            _magtable = Table(names = ('mag','emag','limmag', 'absmag', 'e_absmag'))

            for indy in range(len(self.table['jd'])):

                if self.table['extcorrforcediffimflux'][indy]/self.table['extcorrforcediffimfluxunc'][indy] >= SNT_: 

                    mag_          = to_mag(self.table['extcorrforcediffimflux'][indy])
                    absmag_       =  mag_ - self._DM
                    emag_         = 1.0857 * self.table['extcorrforcediffimfluxunc'][indy] / self.table['extcorrforcediffimflux'][indy]
                    e_absmag_     = np.sqrt(emag_ **2 + self._e_DM**2)
                    limmag_       = -2.5   * np.log10(SNU_*abs(self.table['extcorrforcediffimfluxunc'][indy]))
                else:
                    mag_          = 99.
                    absmag_       = 99.
                    emag_         = 99.
                    e_absmag_     = 99.
                    limmag_       = -2.5*np.log10(SNU_*abs(self.table['extcorrforcediffimfluxunc'][indy]))



                _magtable.add_row([  mag_ , emag_ , limmag_ , absmag_ , e_absmag_  ]  )

            self.table = hstack([self.table, _magtable ])      

            self._converted_to_mag = True
            

            if return_magtable is True:
                if self._added_jd_from_exp is True:
                    if self._corrected_4_z is True:
                        _magtable = hstack([self.table['tfromexplo_zc'], _magtable ])    
                          
                        return _magtable
                    else:
                        _magtable = hstack([self.table['tfromexplo_znc','jd'], _magtable ])      
                        return _magtable

        else:
            print('The table was already converted to magnitudes. No magtable returned')
            return None




    def clean_baseline(self):
        '''
        this function computes the median value of the flux and the median absolute deviation in the baseline (i.e. <-2.5d from FD). 
        It removes any point more than 3 MAD away from the median flux

        NOTE: Clean the baseline BEFORE you correct for the baseline

        '''
        self.table = clean_baseline_phot(self.table)
    
    
    def correct_baseline_offset(self, correction_method = 'separated'):
        '''
        This function computes the median of the baseline and substracts this value to the full flux light curve.


        parameters
        ----------
        correction_method : choose 'separated' to correct the baseline according to the unique OBS ID (field/quadrant/ccd) to avoid any bias. Otherwise, choose 'all' to not distinguish between
        different field/quadrant/ccds. 


        '''
        if correction_method == 'separated':
            self.table = correct_4_baseline_offsets(self.table)

        elif correction_method == 'all':
            self.table = correct_4_baseline_offset_withoutsep(self.table)

        




    def normalise_flux(self):
        '''
        this functions normalises the flux to 1 in order to avoir very small numbers. 
        
        We normalise to the maximal flux point

        /!\ The flux which is being normalised is the extinction AND zero point corrected flux.
        parameters
        ----------
        photfluxe

        returns
        -------
        normalised photflux
        '''
        _normabite                                    = max(self.table['extcorrforcediffimflux'])
        self.table['extcorrforcediffimflux_norm']     = self.table['extcorrforcediffimflux']/_normabite
        self.table['extcorrforcediffimfluxunc_norm']  = self.table['extcorrforcediffimfluxunc']/_normabite

        self._normalised_LC = True
        


        
    def correct_uncertainties(self):
        '''
        This function re estimates the uncertainty on the forced flux measurement fits according to a chi_2 dependency
        As recommended by the forced photometry guide. 
        # Cprrect it before the zero point conversion? 
        
        parameters
        ----------
        self
        
        returns
        -------
        self.table['forcediffimfluxunc_chcorr'] corrected
        
        '''
        # for _ in range(len(self.table['forcediffimfluxunc'])):
        #     if self.table['forcediffimchisq'][_] >= 1.15: # completly arbitrary, let's be honest
        #         self.table['forcediffimfluxunc'][_] = self.table['forcediffimfluxunc'][_] * np.sqrt(self.table['forcediffimchisq'][_])
        
        self.table['forcediffimfluxunc'] = self.table['forcediffimfluxunc'] * np.sqrt(self.table['forcediffimchisq'])

        self._corrected_unc = True
    

    def compute_errors_on_final_diffflux(self):
        '''
      This function computes the errors on the computed Flux of the difference image (z_p corrected). Should also be 
      corrected for extinction, redshift and zeropoint
      Final means that these are the errors on the zp corrected flux... 
      
      parameters
      ----------
      the entire table, in the future this funciton will be in the class so it will be easier to handle
       
      returns
      -------
      
    
    
      '''
        _f_0 = 10**((-0.4) * self.table['zpdiff'])
        
        
        self.table['extcorrforcediffimfluxunc'] = np.sqrt( ( _f_0 *  self.table['forcediffimfluxunc'] )**2 + 
                                                            ( (-0.4) * np.log(10)* self.table['forcediffimflux'] * _f_0 * self.table['zpmaginpsciunc'] )**2  )  
            
    
    def filter_on_SNR(self,SNR_thres): 
        '''
        This function filters on the S/N ration of the diff flux

        parameters
        ----------
        S/R threshold
        '''

        _temp_SN       = self.table[self.table['tfrommarshalfdec'] >= 0.] #to not kick out important first few points
        _temp_Baseline = self.table[self.table['tfrommarshalfdec'] < 0.]

        _SN       = filter_StN_ratio(_temp_SN,SNR_thres)
        # _Baseline = filter_StN_ratio(_temp_Baseline,0.5)

        self.table = vstack([_temp_Baseline,_SN])

    
    def correct_zeropoint(self):
        '''
        This functions corrects for the zero point

        The zero point of an instrument, by definition, is the magnitude of an object that produces one count (or data number, DN) per second. 
        The magnitude of an arbitrary object producing DN counts in an observation of length EXPTIME is therefore:
        m = -2.5 x log10(DN / EXPTIME) + ZEROPOINT 

        For some obscur reason, the zero point converted to flux is already in DN ... So no need to divide by the exposure time. 
        
        parameters
        ----------
        table
    
        returns
        -------
        self.table[''] new column
    
        '''
    
        #return 10**((-0.4) * (-2.5*np.log10((table['forcediffimflux'])/table['exptime']) + table['zpdiff'] ) ) -> doesn't work???

        self.table['zpcorrfordifflux'] = (self.table['forcediffimflux']) * 10**((-0.4) * self.table['zpdiff'] ) 


    def correct_redshift(self, correct_z = False):
        '''
        This function creates a new column with the time corrected for redshift. We here do not consider correction on the flux because of redshift 
        
        parameters
        ----------
        table
        correct_z : bool. meaning we are considering only measured parameters and are not correcting for redshift. 
        
        NOTE: there's not real interest in correcting for redshift here? also not sustainable is no error on the redshift?? 
        
        returns
        -------
        
        '''
        
        
        # about correcting the flux with redsdhift, eran says we should, but ?
        if correct_z is True:
            # self.table['zp_zcorrfordifflux'] = self.table['zpcorrfordifflux'] * (1+self._redshift)
            self.table['zp_zcorrfordifflux'] = self.table['zpcorrfordifflux'] * (1+self._redshift) # correcting the flux for redshift effect ?

            self.table['jd_zcorr'] = self.table['jd']/ (1+self._redshift)

            self.table['tfrommarshalfdec'] = (self.table['tfrommarshalfdec']) / (1+self._redshift) 

            self._corrected_4_z = True
        else:
            self.table['zp_zcorrfordifflux'] = self.table['zpcorrfordifflux']
            self.table['jd_zcorr'] = self.table['jd']/ (1+self._redshift)
            

        
        # time from marshal first detection corrected for redshift
        
        
        
        
    def correct_extinction(self,correct_ext = True):
        '''
        ISM correction... 
        WARNING: make sure you have the dust maps and that they are stored in the right folder
        YOU NEED TO CORRECT FOR ZEROPOINT BEFORE USING THIS FUNCTION
        '''
        
        
        _p48   = { 'ZTF_r': 6339.6,'ZTF_g': 4722.7,'ZTF_i': 7886.1 }
        
        
        _gal_reddening = sfdmap.SFDMap('/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/code_library/python/sfddata-master/', scaling=0.86)
        _ebv           = _gal_reddening.ebv(self._coordinates)
        _alam          = extinction.ccm89(np.array([_p48.get(self.filtre)]) , _ebv * 3.1, 3.1, unit = 'aa')
        
        if correct_ext is True:
            self.table['extcorrforcediffimflux'] = self.table['zp_zcorrfordifflux']*(10**(0.4*_alam))
            self._corrected_4_ext = True

        else:
            self.table['extcorrforcediffimflux'] = self.table['zp_zcorrfordifflux']
        
        
        
        
        
        
        
    

        
        
# TODO: instead of a value change to how many MAD way from median to perform the cut

    #NOTE: can't find refsigpix in the new tables
    # def cuts_refsigpix(self, thres_ref = 5 ):
    #     '''
    #     This function permorms cuts on the light curve based on refsigpix 
    #     paramters
    #     ---------
    #     thres_ref [float/int] how many MAD away from median. by default = 5 
        
    #     returns
    #     -------
    #     '''
        
    #     _thres_ref = np.median(self.table['refsigpix']) + thres_ref * median_absolute_deviation(self.table['refsigpix'])

    #     self.table = self.table[self.table['refsigpix'] <= _thres_ref ]
        

        
    
    def cuts_scisigpix(self, thres_sci = None):
        '''
        We defined an abcolute value of 25 (Yuhan), twice what 

        Robust sigma per pixel in sci image [DN]

        parameters
        ---------
        thres_sci  [float] How many MAD away from median, by default = 5
        
        returns
        -------
        Does not return anything, just filters the table based on the scisigpix threshold 
        '''

        # faut checked en plus ou moins non? ...
        if type(thres_sci)==float:
            _thres_sci = np.median(self.table['scisigpix']) + thres_sci * median_absolute_deviation(self.table['scisigpix']) #actually more strict than what yuhan does... 
            self.table = self.table[self.table['scisigpix'] <= _thres_sci ]
        else:
            self.table = self.table[self.table['scisigpix'] <= 25 ]
       






    #-----------------# 
    #    BIN LC       #
    #-----------------# 


    def bin_lightcurve(self,tmin, tmax, sug = 'median'):
        '''


        note that you have to use ".binned_lc" afterwards 

        '''
        self.binned_lc  = bin_phot_perday(self.table, tmin, tmax, sug = sug)
        self._binned_LC = True




 
    #---------#
    # PLOTTER #
    #---------#
    
    def plot_lc(self, save_fig = False, path_save = None ,
                lc_corr = False, norma = False, add_plot = False, **kwargs):
        '''
        this function plots the light curve
        
        Parameters
        ----------
        add_plot [bool] if you want to plot several LC on the same figure, toggle to "True". Otherwise per default it opens a new figure per LC.
        
        Returns
        -------
        '''
        if add_plot is False:
            plt.figure()
        
        if lc_corr is False :
            plt.errorbar(self.table['jd'], self.table['forcediffimflux'], self.table['forcediffimfluxunc'], fmt='o',  ms = 2.5, color = self._filter_name.get(self.filtre),**kwargs)
            plt.xlabel('JD')
            plt.ylabel('#Counts')

        elif lc_corr is True and norma is False: 
            if self._corrected_LC is True :
                plt.errorbar(self.table['tfrommarshalfdec'], self.table['extcorrforcediffimflux'], self.table['extcorrforcediffimfluxunc'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre),**kwargs)
                plt.xlabel('Time from first ASFD', size = 15)
                plt.ylabel('Flux [Jy]', size = 15)

            else :
                raise Exception('You need to correct the LC before plotting it')
        elif norma is True: 

            if self._corrected_LC is True and self._normalised_LC is True : 
                plt.errorbar(self.table['jd_zcorr'], self.table['extcorrforcediffimflux_norm'], self.table['extcorrforcediffimfluxunc_norm'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre),**kwargs)
                plt.xlabel('JD, z corrected')
                plt.ylabel('Jy, norm')
            else :
                raise Exception('You need to normalise the light curve first')
        


        
        
        if save_fig == True :
            if path_save != None:
                plt.savefig(path_save)
                
            else :
                print('You need to provide a path to save the figure')




    def plot_maglc(self, save_fig = False, path_save = None ,
                add_plot = False, **kwargs):
        '''
        this function plots the magnitude light curve
        
        Parameters
        ----------
        add_plot [bool] if you want to plot several LC on the same figure, toggle to "True". Otherwise per default it opens a new figure per LC.
        
        Returns
        -------
        '''
         ##### FUNCTION FOR THE SECONDARY AXIS 
        def to_abs_mag(x,dm = self._DM):
            '''
            converts to absolute magnitude knowing the distance modulus
            '''
            return x - dm

        def to_app_mag(x,dm = self._DM):
            '''
            converts to apparent magnitude knowing the distance modulus
            '''
            return x + dm



        if add_plot is False:
            plt.figure()
            ax = plt.subplot(111)
 
        if self._converted_to_mag==True:

            _det  = self.table[self.table['mag']!=99.]
            _ndet = self.table[self.table['mag']==99.]


            ax.errorbar(_det['tfrommarshalfdec'], _det['mag'], _det['e_absmag'], fmt='o', alpha = 0.2 ,ms = 3.5, elinewidth=4, color = self._filter_name.get(self.filtre),**kwargs)
            ax.errorbar(_det['tfrommarshalfdec'], _det['mag'], _det['emag'], fmt='o',  ms = 3.5, color = self._filter_name.get(self.filtre),**kwargs)
            
            ax.plot(_ndet['tfrommarshalfdec'], _ndet['limmag'], lw=0 ,marker = 'v', ms = 3.5, color = self._filter_name.get(self.filtre),**kwargs)
            

            ax2 = ax.secondary_yaxis('right', functions=(to_abs_mag, to_app_mag))

        
            ax.set_ylabel('Apparent Magnitude', size = 15)
            ax2.set_ylabel('Absolute Magnitude', size = 15)
            plt.axvline(0,alpha = 0.3, ls = ':', lw = 0.75)
            plt.xlabel('Days from A.S. first detection', size = 15)
            # plt.xlim([-50,])
            plt.xlim(left=-3)
            # plt.ylim([16.5,21.5])

            plt.gca().invert_yaxis() 

            if save_fig == True :
                if path_save != None:
                    plt.savefig(path_save)
                else :
                    print('You need to provide a path to save the figure')

        else: 
            print('You need to add the magnitudes column before plotting it. See the function add_magnitudes_detvndet_meas. ')



    



#TODO: use a minimiser on the number of days to find the best fit? 
#NOTE: there's an interesting check to do at early time: what type of behaviour could be relevant? n=2 or n=0.5 ?






################################################### 

#    CLASS FOR FITTING EXPLOSION DATE 

##################################################





class Fit_t_exp( object ):

    #TODO: fix the chi2on the full data

    '''
        This class is defined to fit the early LC of Infant SNe II in order to estimate the explosion parameters

    '''

    def __init__(self, table, min_fit, max_fit, mini = -10, maxi = 7, filter = 'ZTF_r'):
        '''
        Parameters
        ----------
        table contraining the time from first marshal detection, the flux and the error on the flux 

        min_fit [float] minimal bound of the fit
        max_fit [float] maximal bound of the fit
        mini    [flaot] minimal bound of data considered
        maxi    [flaot] maximal bound of data considered for the final chi_2 fit 




        '''
        self._set_max_fit_val(max_fit)
        self._time_zone_selected = False
        self._filter_name   = {'ZTF_g':'darkolivegreen','ZTF_r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        self.set_filter(filter)

        # full data to comapre the chi_2 to
        self.set_time_zone_interest(mini, maxi)
        self.set_time_from_asfd()
        self.set_flux()
        self.set_e_flux()


        # reduced data to perform the fit onto
        self.set_time_zone_fit(min_fit, max_fit)
        # if len(self.fittable) <

        self.set_time_from_asfd_fit()
        self.set_flux_fit()
        self.set_e_flux_fit()

        

    #-----------#
    #  SETTERS  #
    #-----------#

    def _set_max_fit_val(self, max_fit):
        '''
        '''
        self._max_fit = max_fit

    def set_filter(self,filter):
        '''
        '''
        self.filter = filter

    def set_phot_table(self, table):
        '''

        '''
        self.table = table

    def set_time_zone_interest(self, mini = -10, maxi = 7):
        '''
        the early light curve we're interested in:
        taking from 10 days prior to FD until a week from FD
        '''
        self.table = self.table[(self.table['tfrommarshalfdec']>=mini)&(self.table['tfrommarshalfdec']<=maxi)]

    def set_time_from_asfd(self):
        '''
        '''
        self.t = self.table['tfrommarshalfdec']

    def set_flux(self):
        '''
        '''

        self.f = self.table['extcorrforcediffimflux']

    def set_e_flux(self):
        '''
        '''

        self.e_f = self.table['extcorrforcediffimfluxunc']




    def set_time_zone_fit(self, min_fit , max_fit ):
        '''

        '''

        self.fittable = self.table[(self.table['tfrommarshalfdec']>=min_fit)&(self.table['tfrommarshalfdec']<=max_fit)]

        self._time_zone_selected = True

    def set_time_from_asfd_fit(self):
        '''
        '''
        self.t_fit = self.fittable['tfrommarshalfdec']

    def set_flux_fit(self):
        '''
        '''

        self.f_fit = self.fittable['extcorrforcediffimflux']

    def set_e_flux_fit(self):
        '''
        '''
        self.e_f_fit = self.fittable['extcorrforcediffimfluxunc']


    ############
    # GETTERS  #
    ############


    def get_fit_roi_length(self):
        '''
        This function gets the length of the region of ionterest of fit (early rise)
        '''
        if self._time_zone_selected == True:
            len_roi = len(self.fittable[self.fittable['tfrommarshalfdec']>=0])
            return len_roi
        else:
            print('You need to select a region of interest got the fit')

    
    def get_chi_2(self, param):
        '''
        This function returns the Chi squared of the function fitted to the data. 

        parameters
        ----------
        self

        returns
        -------

        ''' 
        #a , T_exp , n  = param 
        model_ = rise_typeII(self.t_fit,*param)
        pull_  = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq = np.sum(pull_)
        dof    = len(self.f_fit)-len(param)
        if dof == 0: #avoid division by zero? 
            dof = 1

        return  chi_sq / dof
    
    def _get_chi2minuit(self, a, t_exp, n):
        '''

        '''

        param = [a,t_exp,n]

        return self.get_chi_2(param)



    def get_explosion_t(self):
        '''
        '''
        self.t_exp  = self.minuit_output.params[1].value
        self.dt_exp = self.minuit_output.params[1].error
        return self.t_exp, self.dt_exp
        

    def get_fit_index(self):
        '''
        '''
        self.n      = self.minuit_output.params[2].value
        self.dt_n   = self.minuit_output.params[2].error
        return self.n, self.dt_n




    #######################
    #        FITTER       #
    #######################
    
    def set_minuit(self, guess, boundaries=None, fixed=None, errordef=1 ,print_level=0):
        '''
        
        '''

        a, t_exp, n     = guess
        self._paramname = "a,t_exp,n".split(",")
        
        local_ = locals() # dict with the local variables
        self.minuit_kwargs = {}
        for p in self._paramname:
            self.minuit_kwargs[p] = local_[p]
            
        if boundaries is not None:
            for k,v in boundaries.items():
                self.minuit_kwargs["limit_"+k] = v
                
        if fixed is not None:
            for k in fixed:
                self.minuit_kwargs["fix_"+k] = True
        
        self.minuit = Minuit(self._get_chi2minuit, errordef=errordef, error_a = 5e-9 , error_t_exp = 0.1, error_n = 0.05,
                                    print_level=print_level, **self.minuit_kwargs)
        

    def fit_minuit(self, guess, boundaries=None, fixed=None):
        
        """ 
    
        """
        
        self.set_minuit(guess, boundaries, fixed)
        self.minuit_output = self.minuit.migrad()
        print(self.minuit_output)

    
    
    ##################
    #   PLOTTING     #
    ##################



    def plot_fit(self, add_plot = False):
        '''
        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_ = np.linspace(-10, max(self.t), 1000)
        
        a_      = self.minuit_output.params[0].value
        t_exp_  = self.minuit_output.params[1].value
        n_      = self.minuit_output.params[2].value
        
        
        dt_exp_ = self.minuit_output.params[1].error
        dn_     = self.minuit_output.params[2].error

        
        
        tx_string = f't_exp = {t_exp_:.3f}±{dt_exp_:.3f} d '
        n_string  = f'n = {n_:.3f}±{dn_:.3f} '

        chi2_      = self.get_redchi_2_gen([a_, t_exp_, n_])
        chi_string = f'Chi_2 early LC = {chi2_:.3f}'

        chi2_fit       = self.get_redchi_2_fit([a_, t_exp_, n_])
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        fit_interval_string = f'd = {self._max_fit}' 




        if add_plot == None:
            plt.figure()
            plt.errorbar(self.t, self.f, self.e_f, fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter) )
            plt.plot(x_ , rise_typeII(x_, a_, t_exp_, n_), ls = '--' , label= n_string , color = _fitcolor )

            plt.plot(x_[0] , rise_typeII(x_, a_, t_exp_, n_)[0], color = 'white',label= chi_string )
            plt.plot(x_[1] , rise_typeII(x_, a_, t_exp_, n_)[1], color = 'white',label= chi_fit_string )


            plt.axvline(t_exp_, color = _fitcolor ,ls=':',alpha = 0.5, label = tx_string )
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from ASFD [days]', size = 15)
            plt.ylabel('Flux [Jy]', size = 15)

            plt.legend()
        
        else: 
            add_plot.plot(0,0, ls = '--' , label= fit_interval_string , color = 'white' )
            add_plot.errorbar(self.t, self.f, self.e_f, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
            add_plot.plot(x_ , rise_typeII(x_, a_, t_exp_, n_), ls = '--' , label= n_string , color = _fitcolor )

            add_plot.plot(x_[0] , rise_typeII(x_, a_, t_exp_, n_)[0], color = 'white',label= chi_string )
            add_plot.plot(x_[1] , rise_typeII(x_, a_, t_exp_, n_)[1], color = 'white',label= chi_fit_string )


            add_plot.axvline(t_exp_, color = _fitcolor ,ls=':',alpha = 0.5, label = tx_string )
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from ASFD [days]', size = 15)
            plt.ylabel('Flux [Jy]', size = 15)

            plt.legend(fontsize = 10)

        
        
        

        
        
 


    def get_redchi_2_gen(self, param):
        '''
    
        parameters
        ----------
        self

        returns
        -------

        ''' 
    
        #a , T_exp , n  = param 
        model_ = rise_typeII(self.t,*param)
        pull_  = (( self.f - model_ )/ self.e_f )**2
        chi_sq = np.sum(pull_)
        dof    = len(self.f) - len(param)

        return  chi_sq/dof 


    def get_redchi_2_fit(self, param):
        '''
        parameters
        ----------
        self

        returns
        -------
        ''' 
        #a , T_exp , n  = param 
        model_ = rise_typeII(self.t_fit,*param)
        pull_  = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq = np.sum(pull_)
        dof    = len(self.f_fit)-len(param)

        return  chi_sq/dof


        






     
        # if zp_corr is False: 
            
        #     _temp_tab           = self.table[self.table['index']<=indiex]
        #     _temp_forcedflux    = np.array(_temp_tab['forcediffimflux'])
        
        
        #     return np.median(_temp_forcedflux)
        
        # elif zp_corr is True: 

        #     if 'zpcorrfordifflux' in self.table.colnames :
        #         _temp_tab           = self.table[self.table['index']<=indiex]
        #         _temp_forcedflux    = np.array(_temp_tab['zpcorrfordifflux'])
        #         return np.median(_temp_forcedflux)
        #     else :
        #         raise Exception('You need to correct for zeropoint prior to the baseline')
        
  
        
    # def get_std_dev(self,indiex):
    #     '''
    #     this is essentially for the baseline
    #     '''


    #     if self._corrected_LC  is True: 
    #         if binned is False:
    #             return np.median( self.table['extcorrforcediffimflux'][self.table['tfrommarshalfdec'] <= -0.5]  )
    #         elif binned is True:
    #             if self._binned_LC is True:
    #                 return np.median( self.binned_lc['extcorrforcediffimflux'][self.binned_lc['tfrommarshalfdec'] <= -0.5]  )
    #             else : 
    #                 print('You need to bin the LC first')
    #     else:
    #         print('You need to correct the LC first')
        
    #     self.table['index'] = [int(x) for x in self.table['index']]
    #     _temp_tab           = self.table[self.table['index']<=indiex]
    #     _temp_forcedflux    = [float(x['forcediffimfluxunc']) for x in _temp_tab  if x['forcediffimflux'][-1].isdigit()]
        
    #     return np.std(_temp_forcedflux)
    
    # #def get_explosion_time():
        
    

    # def compute_magnitudes_from_corrected_flux(self):
    #     '''
    #     This function computes the magnitude from the corrected flux (ie: extcorrforcediffimflux)

    #     parameters
    #     ----------

    #     returns
    #     -------
    #     '''

    #     self.table['magcorr']  = to_mag(self.table['extcorrforcediffimflux'])
    #     self.table['emagcorr'] = error_on_conv_to_mag(self.table['extcorrforcediffimfluxunc'], self.table['extcorrforcediffimflux'])
        

    # def compute_nd_vs_det(self):
    #     '''
    #     This function computes the detections vs Non detections. 
    #     According to the Forced photometry guide, we choose SNT = 3 and SNU = 5. 
    #     If a the flux divided by the uncertainty is bigger/equal than 3, it is a dectection. If flux/unc < 3, it's a non-detection.
    #     Limiting magnitudeds are computed as 5 times sigma.

    #     parameters
    #     ----------

    #     returns
    #     -------

    #     note
    #     ----

    #     When computing the flux, we already account for the zero point so we don't need to take it into account again.

    #     '''

    #     SNT_ = 3
    #     SNU_ = 5

    #     self.magtable = Table(names = ('jd','tfrommarshalfdec','mag','emag','limmag'))

    #     for indy in range(len(self.table['jd'])):
    #         if self.table['extcorrforcediffimflux'][indy]/self.table['extcorrforcediffimfluxunc'][indy] >= SNT_: 
    #             mag_     = to_mag(self.table['extcorrforcediffimflux'][indy])
    #             emag_    = 1.0857*self.table['extcorrforcediffimfluxunc'][indy]/self.table['extcorrforcediffimflux'][indy]
    #             limmag_  = -2.5*np.log10(SNU_*abs(self.table['extcorrforcediffimfluxunc'][indy]))
    #         else:
    #             mag_      = 99.
    #             emag_     = 99.
    #             limmag_   = -2.5*np.log10(SNU_*abs(self.table['extcorrforcediffimfluxunc'][indy]))



    #         self.magtable.add_row([ self.table['jd'][indy], self.table['tfrommarshalfdec'][indy],  mag_ , emag_ , limmag_   ]  )

    #     return self.magtable





    ##################
    #    CLEANERS
    ###############

    # def clean_lc(self,  remove_big_errors = True, clean_baseline_fo = True ,bin_late_time_lc = False) :
    #     '''
    #     This function cleans the lightcurve from potential outliers. 
    #     Cleaning the light curve is STATIUSTICAL and is part of further analysis.
    #     One option is also to remove "big errors" meaning data points with large uncertainties. These are kicked out by looking at most
    #     of the errors of the data and kicking out those which have an error 3*std bigger than the median error.

    #     parameters
    #     ----------
    #     remove_big_errors [bool] whether or not you want to also correct the uncertainties. By default it will be true 

    #     '''


    #     # _temp = remove_big_err_point(self.table)
        

    #     # self.cleantable = 


    #     # if remove_big_errors is True:
    #     #     #self.remove_big_err_point()
        
    #     if clean_baseline_fo is True:
    #         self.clean_baseline()

    #     if bin_late_time_lc is True :
    #         self.bin_lightcurve(tmin = 7, tmax = 100 , sug = 'median')
        
    #     self._clean_LC = True













# def cuts_scibckgnd(self, thres_bkg = 5):
    #     '''
    #     This function permorms cuts on the light curve based on scisigpix. By default it will kick out any points whos scisigpix is higher than the median threshold 
    #     +5*the medidan absolute deviation. This ensures the quality of the difference image in order to perform forced photometry. 

    #     Background level in sci image [DN]

    #     parameters
    #     ---------
    #     thres_bkg  [float] How many MAD away from median, by default = 5
        
    #     returns
    #     -------
    #     Does not return anything, just filters the table based on the scisigpix threshold 
    #     '''
        
    #     _thres_bkg = np.median(self.table['scibckgnd']) + thres_bkg * median_absolute_deviation(self.table['scibckgnd'])
        
    #     self.table = self.table[self.table['scibckgnd'] <= _thres_bkg ]
        
        
    # def cuts_sciseeing(self, thres_seeing = 4 ):
    #     '''
    #     This function performs cuts on the seeing of the science image
        
        
    #     '''
    #     self.table = self.table[self.table['sciinpseeing'] <= thres_seeing]
    

    

                
    
    # def plot_check_unc_chisq(self):
    #     '''
    #     This functions checks the PSF fit chi_2 dependancy on the flux. If a correlation is seen, I strongly recommenc
    #     to correct the uncertainties on the flux accordingly.
    #     Parameters
    #     ----------
    #     Self
        
    #     Returns
    #     -------
    #     A plot 
        
    #     '''
        
    #     _field_list = list(np.unique(self.table['field']))
        
    #     plt.figure()
    #     for _fieldx in _field_list:
    #         _temptab = self.table[self.table['field'] == _fieldx ]
    #         _mu      = np.mean(_temptab['forcediffimchisq'])
    #         plt.plot(_temptab['forcediffimflux'], _temptab['forcediffimchisq'], marker='x', ms = 4,lw = 0,label = f'Field {_fieldx:.0f}, m = {_mu:.2f}' )
    #         del _temptab
        
        
    #     plt.xlabel('forcediffimflux')
    #     plt.ylabel('forcediffimchisq')
    #     plt.legend()
                


        
    # def plot_interpolation(self, savefig = False, pathsavefig = None):
    #     '''
    #     this function plots the interpolation of the lightcurve 

    #     parameters
    #     ---------

    #     returns
    #     -------

    #     '''
    #     if self._interpolatedLC is True:
    #         if self._binned_LC is True:
    #             plt.figure()
    #             plt.errorbar(self.binned_lc['tfrommarshalfdec'], self.binned_lc['extcorrforcediffimflux'], self.binned_lc['extcorrforcediffimfluxunc'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre), alpha = 0.1, label = None)
    #             plt.plot(self.binned_lc['tfrommarshalfdec'], self.binned_lc['intrasplinflux'], ls = '--', color = 'dimgray', label = 'Full LC interpolation')
    #             plt.xlabel('Time from 1st detection [days]')
    #             plt.ylabel('Flux [Jy]')
    #             plt.legend()

    #             if savefig is True : 
    #                 if pathsavefig is not None:
    #                     plt.savefig(pathsavefig + f'fullLCinterpolation_{self.objectname}_{self.filtre}.pdf')

    #         elif self._binned_LC is False: 
    #             print('This plot does not use the binned LC ! ')
    #             plt.figure()
    #             plt.errorbar(self.table['tfrommarshalfdec'], self.table['extcorrforcediffimflux'], self.table['extcorrforcediffimfluxunc'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre), alpha = 0.1, label = None)
    #             plt.plot(self.table['tfrommarshalfdec'], self.table['intrasplinflux'], ls = '--', color = 'dimgray', label = 'Full LC interpolation')
    #             plt.xlabel('Time from 1st detection [days]')
    #             plt.ylabel('Flux [Jy]')
    #             plt.legend()

    #             if savefig is True : 
    #                 if pathsavefig is not None:
    #                     plt.savefig(pathsavefig + f'fullLCinterpolation_{self.objectname}_{self.filtre}.pdf')

    #     else:
    #         print('You need to interpolate the data first! ')

    # def plot_binned_lc(self):
    #     '''
    #     This function plots the binned light curve. Note: you need to correct the lightcurve prior to bin and so on...?
    #     '''
    #     if self._binned_LC is True:
    #         plt.figure()
    #         plt.errorbar(self.binned_lc['tfrommarshalfdec'], self.binned_lc['extcorrforcediffimflux'], self.binned_lc['extcorrforcediffimfluxunc'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre), alpha = 0.8, label = None)
    #         #plt.plot(self.binned_lc['tfrommarshalfdec'], self.binned_lc['intrasplinflux'], ls = '--', color = 'dimgray', label = 'Full LC interpolation')
    #         plt.xlabel('Time from 1st detection [days]')
    #         plt.ylabel('Flux [Jy]')
            
    #     else: 
    #         print('You need to bin the light curve prior to plot it! ')
       
    







    # #TODO: remove these functions and change it for a class itself which we will call explosion fitter.
    
    # def fit__early_lc(self,minx,maxx, special_condiation = False):
    #     '''
    #     this function uses curve_fit to fit the early rise function
    
    #     Parameters
    #     ----------
    #     dataphot
    
    #     minx
    #     maxx
    
    
    #     Returns
    #     -------
    #     popt [array]
    #     pcov [array]
    
    #     '''

    #     _phot = self.table[ (self.table['tfrommarshalfdec']>=minx)  & (self.table['tfrommarshalfdec']<=maxx) ]


    #     if special_condiation is True:
    #         # if self.objectname == 'ZTF20ablygyy': 
    #         _phot = _phot[(_phot['tfrommarshalfdec']<0) | (_phot['tfrommarshalfdec']>2.5) ]
    #     # print(_phot)

    #     _days  = _phot['tfrommarshalfdec']
    #     _flux  = _phot['extcorrforcediffimflux']
    #     _error = _phot['extcorrforcediffimfluxunc']

        
    #     if self.objectname == 'ZTF18abcptmt':
    #         popt , pcov = curve_fit(rise_typeII, _days, _flux, sigma=_error, p0 = [2e-8,-1.2,0.5],bounds = ( [0,-5,0] , [np.inf, 0 ,10] )   )
    #     else:
    #         popt , pcov = curve_fit(rise_typeII, _days, _flux, sigma=_error,bounds = ( [0,-5,0] , [np.inf, 0 ,10] )   )

    #     return popt,pcov



    
    
    # def fit_and_plot_phot_mod( self, peak_time,fit_upp_bound , bin_data = False, return_paramchi2 = False, figpath = None, add_plot = None , special_condiation = False):
    #     '''
    #         This function ... 

    #         parameters
    #         ----------

    #         candi        [string] name of the object
    #         filt         [string] which filter
    #         peak_time    [float]  Upper bound for fitting and plotting (should be t_stop, but for simplicity...)
    #         figpath      [string] where do you want to save the figure
    #         fit_upp_bound 
    #         bin_data     [bool] whether you want to use the binned light curve 
    #         add_plot     [plt]  If you want to overlay several plots of the fitting (say red and g band). This variable should be something like ax = ax = plt.subplot(111)
    #                             after a plt.figure(). If you want to make several comparison, please make sure you define different ax. 


    #         returns
    #         -------


    #     '''


        

    #     #if self._normalised_LC is True : 

    #     #_temptab         = self.table[self.table['tfrommarshalfdec']<=peak_time]

    #     param , covar   = self.fit__early_lc(-10, fit_upp_bound, special_condiation = special_condiation)

    #     print(param)
    #     print(covar)
    #     print(f'The explosion date is estimated at {param[1]} pm {covar[1][1]}')
    #     print(f'The exponent of the fit  is {param[2]} pm {covar[2][2]}')

    #     days            = np.linspace(-10,peak_time,2000)
    #     flux_flit       = rise_typeII(days, *param)
    #     chi2 , dof      = self.chi_sq( fit_type ='LCrise' , minx = -10 , maxx = fit_upp_bound, params = param )

    #     print(f'The Chi_2 is {chi2/dof}')
      
      
    #     if add_plot is None: 
    #         plt.figure()
    #         plt.errorbar(self.table['tfrommarshalfdec'], self.table['extcorrforcediffimflux'], self.table['extcorrforcediffimfluxunc'], fmt='o',  ms = 6 , alpha = 0.5 ,color = self._filter_name.get(self.filtre), label = None)
    #         plt.plot(days,flux_flit , lw = 1 , ls = '--', label = f'[0-{fit_upp_bound}] ')
    #         plt.xlabel('Days from 1st Marshal Detection')
    #         plt.ylabel('Flux')
    #         plt.title( f'Chi2 = {chi2/dof:.3f}; Texp = {param[1]:.3f} pm {covar[1][1]:.3f}', fontsize=15 )
    #         #P(chi2) = {    list(proba_chi2)[0]:.6f}
    #         plt.suptitle(f'{self.objectname}, {self.filtre} band')

    #         if figpath is not None: 
    #             plt.savefig(figpath+f'rise_fit_{self.objectname}_from0to{fit_upp_bound}days_{self.filtre}band.pdf')
            

    #     else: 
    #         add_plot.errorbar(self.table['tfrommarshalfdec'], self.table['extcorrforcediffimflux'], self.table['extcorrforcediffimfluxunc'], fmt='o',  ms = 6 , color = self._filter_name.get(self.filtre), alpha = 0.5, label = None)
    #         add_plot.plot(days,flux_flit , lw = 1 , ls = '--', color = self._filter_name.get(self.filtre) ,label = f'[0-{fit_upp_bound}] ')
    #         add_plot.set_xlabel('Days from 1st Marshal Detection', size =20)
    #         add_plot.set_ylabel('Flux [Jy]', size = 20)
    #         #plt.title( f'Chi2 = {chi2/dof:.3f}; Texp = {param[1]:.3f} pm {covar[1][1]:.3f}', fontsize=15 )
    #         #plt.suptitle(f'{self.objectname}, {self.filtre} ')

    #         if figpath is not None: 
    #             plt.savefig(figpath+f'rise_fit_{self.objectname}_from0to{fit_upp_bound}days_{self.filtre}band.pdf')
            
        

    #         if return_paramchi2 is True: 
    #             return param , covar, chi2/dof
        # else :
        #     raise Exception('You need to normalise the LC before fitting it ')

        
        # def final_expl_time_plot(self):
        #     '''
        #     this function plots the fit of the explosion time ( final version ? )
        #     '''



    #-------------------#
    #  SPLINE INTERPO   #
    #-------------------#
    # def nopt_spline_interpolation_LC(self, t_base, t_apeak , save_figure = False, pathsavefig = None):
    #     '''
    #     This is the non optimisez spline interpolation. Here, you specify the days at which you wish to perform the spline inteprolation.
    #     The spline interpolation in the case of these lightcurves has three regimes : the baseline, the rise and then the decline/plateau phase. These regimes
    #     are delimited by the time "t_base" and "t_apeak" which are supposed to delimit the region of the rise. However, I've noticed that giving a time slightly after
    #     peak helps providing a better interpolation. 
    #     This function was made in order to help the sigma clipping. 

    #     parameters
    #     ----------
    #     t_base  [float] 
    #     t_apeak [float]
    #     save_figure       [bool]   -optional-
    #     pathsavefig       [string] -optional-

    #     returns
    #     -------

    #     '''
    # #tparam, tcovar = self.optimised_spline_interp_fullLC()
        
    #     self.add_spline_interpolation(t_base, t_apeak)
    #     self.plot_interpolation(savefig = save_figure, pathsavefig = pathsavefig)
        
       
        
    # def spline_interpolation_LC(self, return_parameters = False, save_figure = False, pathsavefig = None):
    #     '''
    #     This is the full function to perform an "optimised" spline interpolation. 
    #     This function calls the "optimised_pline_interp_fullLC" which will look for the a priori base time to delimit the differet regimes of the lightcurve.
    #     It then calculates the chi_2 of tge interpolation vs the data, adds the spline interpolation to the data as a column and plots the spline interpoaltion
    #     for visualisaiton

    #     parameters
    #     ----------
    #     return_parameters [bool]   -optional-
    #     save_figure       [bool]   -optional-
    #     pathsavefig       [string] -optional-



    #     returns
    #     -------
    #     '''

    #     tparam, tcovar = self.optimised_spline_interp_fullLC()
    #     chi2           = self.chi_2_spline_interpo(tparam)
    #     self.add_spline_interpolation(tparam[0], tparam[1])
    #     self.plot_interpolation(savefig = save_figure, pathsavefig = pathsavefig)
    #     print(f'Chi 2 is {chi2}')

    #     if return_parameters is True :
    #         return tparam, tcovar, chi2





    # def fun_spline_full(self, _x, t_br , t_rd):
    #     '''
    #     The full spline interpolation function the lightcurve with the three different regimes delimited though boolean filters

    #     parameters
    #     ----------

    #     returns
    #     -------
    #     '''

    #     # _x              = np.array(table['tfrommarshalfdec'])
    #     self.spline_baseline = interp_spline_baseline(self.table,t_br)
    #     self.spline_rise     = interp_spline_rise(self.table,t_br,t_rd)
    #     self.spline_decl     = interp_spline_rest(self.table,t_rd)

    #     t_baseline      = ( _x < t_br  ) 
    #     _baseline       = np.array(self.spline_baseline(_x) * t_baseline)

    #     t_rise_c0       = ( _x >= t_br )
    #     t_rise_c1       = ( _x  < t_rd )
    #     _rise           = np.array(self.spline_rise(_x) * t_rise_c0 * t_rise_c1)

    #     t_decl          = ( _x >= t_rd )
    #     _decl           = np.array(self.spline_decl(_x) * t_decl )

    #     interpo         = _baseline + _rise + _decl
    #     condition_a_0   = ( interpo >= 0 )
    #     interpo         = interpo * condition_a_0

    #     return interpo
        

    # def add_spline_interpolation(self, t_br, t_rd):
    #     '''
    #    this functions interpolates the lightcurve according to phase, see the function called 

    #     parameters
    #     ----------
    #     t_br     [float] time at which you estimate that the baseline ends and the rise begins (e.g. : -1 days)
    #     t_rd     [float] time at which you estimate that the rise ends and that the decline or plateau phase begins (e.g.: 15 days)

    #     Note that these time have to be given in the timeframe of the marshal first detection. i.e.: day 0 is the day of first detection

    #     returns
    #     -------
    #     adds the column "intrasplinflux" which is the interpolated lightcurve to the table
    #     '''


    #     self.table['intrasplinflux'] = self.fun_spline_full( self.table['tfrommarshalfdec'], t_br, t_rd)
        
    #     self._interpolatedLC = True 



    # def optimised_spline_interp_fullLC(self):
    #     '''

    #     '''

    #     _phot = self.table

    #     _days  = _phot['tfrommarshalfdec']
    #     _flux  = _phot['extcorrforcediffimflux']
    #     _error = _phot['extcorrforcediffimfluxunc']


    #     popt , pcov = curve_fit(self.fun_spline_full, _days, _flux,  sigma=_error, bounds = ( [-2.5,0] , [0, 40] )   )
    #     return popt, pcov

    # def chi_2_spline_interpo(self, params = None):
    #     '''
    #     chi2 of the spline intepolation
    #     '''
    #     flux_           = self.table['extcorrforcediffimflux']  
    #     days_           = self.table['tfrommarshalfdec']
    #     eflux_          = self.table['extcorrforcediffimfluxunc']
   
    #     pull_  = (( flux_ - self.fun_spline_full(days_,*params) )/ eflux_ )**2
    #     chi_sq = np.sum(pull_)
    #     dof    = len(flux_)-len(params)
    #     print(f'Reduced Chi_sq is  {chi_sq/dof}')
    #     return chi_sq/dof


    # '''
    # Spline interpolation but only on the binned light curve
    # '''

    # def spline_interpolation_LC(self, return_parameters = False, save_figure = False, pathsavefig = None):
    #     '''
    #     '''

    #     tparam, tcovar = self.optimised_spline_interp_fullLC()
    #     chi2           = self.chi_2_spline_interpo(tparam)
    #     self.add_spline_interpolation(tparam[0], tparam[1])
    #     self.plot_interpolation(savefig = save_figure, pathsavefig = pathsavefig)
    #     print(f'Chi 2 is {chi2}')

    #     if return_parameters is True :
    #         return tparam, tcovar, chi2





    # def fun_spline_full(self, _x, t_br , t_rd):
    #     '''
    #     '''

    #     # _x              = np.array(table['tfrommarshalfdec'])
    #     self.spline_baseline = interp_spline_baseline(self.binned_lc,t_br)
    #     self.spline_rise     = interp_spline_rise(self.binned_lc,t_br,t_rd)
    #     self.spline_decl     = interp_spline_rest(self.binned_lc,t_rd)

    #     t_baseline      = ( _x < t_br  ) 
    #     _baseline       = np.array(self.spline_baseline(_x) * t_baseline)

    #     t_rise_c0       = ( _x >= t_br )
    #     t_rise_c1       = ( _x  < t_rd )
    #     _rise           = np.array(self.spline_rise(_x) * t_rise_c0 * t_rise_c1)

    #     t_decl          = ( _x >= t_rd )
    #     _decl           = np.array(self.spline_decl(_x) * t_decl )

    #     interpo         = _baseline + _rise + _decl
    #     condition_a_0   = ( interpo >= 0 )
    #     interpo         = interpo * condition_a_0

    #     return interpo
        

    # def add_spline_interpolation(self, t_br, t_rd):
    #     '''
    #    this functions interpolates the lightcurve according to phase, see the function called 

    #     parameters
    #     ----------
    #     t_br     [float] time at which you estimate that the baseline ends and the rise begins (e.g. : -1 days)
    #     t_rd     [float] time at which you estimate that the rise ends and that the decline or plateau phase begins (e.g.: 15 days)

    #     Note that these time have to be given in the timeframe of the marshal first detection. i.e.: day 0 is the day of first detection

    #     returns
    #     -------
    #     adds the column "intrasplinflux" which is the interpolated lightcurve to the table
    #     '''


    #     self.binned_lc['intrasplinflux'] = self.fun_spline_full( self.binned_lc['tfrommarshalfdec'], t_br, t_rd)
        
    #     self._interpolatedLC = True 



    # def optimised_spline_interp_fullLC(self):
    #     '''

    #     '''
    #     if self._binned_LC is True: 
    #         _phot = self.binned_lc
    #     else:
    #         print('You need to bin the LC before applying the interpolation')

    #     _days  = _phot['tfrommarshalfdec']
    #     _flux  = _phot['extcorrforcediffimflux']
    #     _error = _phot['extcorrforcediffimfluxunc']


    #     popt , pcov = curve_fit(self.fun_spline_full, _days, _flux,  sigma=_error, bounds = ( [-2.5,0] , [0, 40] )   )
    #     return popt, pcov

    # def chi_2_spline_interpo(self, params = None):
    #     '''
    #     chi2 of the spline intepolation
    #     '''
    #     flux_           = self.binned_lc['extcorrforcediffimflux']  
    #     days_           = self.binned_lc['tfrommarshalfdec']
    #     eflux_          = self.binned_lc['extcorrforcediffimfluxunc']
   
    #     pull_  = (( flux_ - self.fun_spline_full(days_,*params) )/ eflux_ )**2
    #     chi_sq = np.sum(pull_)
    #     dof    = len(flux_)-len(params)
    #     print(f'Reduced Chi_sq is  {chi_sq/dof}')
    #     return chi_sq/dof

    
    #-------------------#
    #  SIGMA CLIPPING   #
    #-------------------#

    # def sigma_clip_LC(self, sigma = 3, visualisation = False, save_figure = False, pathsavefig = None):
    #     '''
    #     this function clips the data which looks faulty based on a spline interpolation

    #     parameters
    #     ----------

    #     returns
    #     -------

    #     warning
    #     -------
    #     The table returned is reduced. (in the next version, there is a chance that I will remove straight away all "irrelevant" columns)

    #     '''

    #     _chloe = [x['jd','forcediffimflux','forcediffimfluxunc','zpcorrfordifflux','tfrommarshalfdec','jd_zcorr','extcorrforcediffimflux','extcorrforcediffimfluxunc','intrasplinflux'] 
    #                 for x in self.table 
    #                 if
    #                 (x['extcorrforcediffimflux'] < x['intrasplinflux'] + sigma * x['extcorrforcediffimfluxunc']) and 
    #                 (x['extcorrforcediffimflux'] > x['intrasplinflux'] - sigma * x['extcorrforcediffimfluxunc']) 
    #                 ]


    #     if visualisation is False:
    #         self.clipped_table = vstack(_chloe)
    #     elif visualisation is True:
    #         self.spline_interpolation_LC(save_figure=False,pathsavefig=None)
    #         self.table = vstack(_chloe)
    #         if self.filtre == 'ZTF_r':
    #             plt.errorbar(self.table['tfrommarshalfdec'], self.table['extcorrforcediffimflux'],self.table['extcorrforcediffimfluxunc'] ,fmt = 'o', ms = 2.5,color = 'darkred', alpha = 0.5)
    #         elif self.filtre == 'ZTF_g':
    #             plt.errorbar(self.table['tfrommarshalfdec'], self.table['extcorrforcediffimflux'],self.table['extcorrforcediffimfluxunc'] ,fmt = 'o', ms = 2.5,color = 'darkgreen', alpha = 0.5)
    #         if save_figure is True :
    #             if pathsavefig is not None : 
    #                 plt.savefig(pathsavefig + f'sigmaclippedLC_{self.objectname}_{self.filtre}.pdf')
    #             elif pathsavefig is None:
    #                 print('You need to specify a path where to save the sigma clipped LC')


    # def sigma_clip_binned_LC(self, sigma = 3, visualisation = False, save_figure = False, pathsavefig = None):
    #     '''
    #     this function clips the data which looks faulty based on a spline interpolation

    #     parameters
    #     ----------

    #     returns
    #     -------

    #     warning
    #     -------
    #     The table returned is reduced. (in the next version, there is a chance that I will remove straight away all "irrelevant" columns)

    #     '''
    #     # _chloe = [x['jd','forcediffimflux','forcediffimfluxunc','zpcorrfordifflux','tfrommarshalfdec','jd_zcorr','extcorrforcediffimflux','extcorrforcediffimfluxunc','intrasplinflux'] 
    #     #             for x in self.table 
    #     #             if
    #     #             (x['extcorrforcediffimflux'] < x['intrasplinflux'] + sigma * x['extcorrforcediffimfluxunc']) and 
    #     #             (x['extcorrforcediffimflux'] > x['intrasplinflux'] - sigma * x['extcorrforcediffimfluxunc']) 
    #     #             ]


    #     _chloe = [x['tfrommarshalfdec','extcorrforcediffimflux','extcorrforcediffimfluxunc','intrasplinflux'] 
    #                 for x in self.binned_lc 
    #                 if
    #                 (x['extcorrforcediffimflux'] < x['intrasplinflux'] + sigma * x['extcorrforcediffimfluxunc']) and 
    #                 (x['extcorrforcediffimflux'] > x['intrasplinflux'] - sigma * x['extcorrforcediffimfluxunc']) 
    #                 ]


    #     if visualisation is False:
    #         self.binned_lc = vstack(_chloe)
    #     elif visualisation is True:
    #         self.spline_interpolation_LC(save_figure=False,pathsavefig=None)
    #         self.binned_lc = vstack(_chloe)
    #         if self.filtre == 'ZTF_r':
    #             plt.errorbar(self.binned_lc['tfrommarshalfdec'], self.binned_lc['extcorrforcediffimflux'],self.binned_lc['extcorrforcediffimfluxunc'] ,fmt = 'o', ms = 2.5,color = 'darkred', alpha = 0.5)
    #         elif self.filtre == 'ZTF_g':
    #             plt.errorbar(self.binned_lc['tfrommarshalfdec'], self.binned_lc['extcorrforcediffimflux'],self.binned_lc['extcorrforcediffimfluxunc'] ,fmt = 'o', ms = 2.5,color = 'darkgreen', alpha = 0.5)
    #         if save_figure is True :
    #             if pathsavefig is not None : 
    #                 plt.savefig(pathsavefig + f'sigmaclippedLC_{self.objectname}_{self.filtre}.pdf')
    #             elif pathsavefig is None:
    #                 print('You need to specify a path where to save the sigma clipped LC')

    # def sigma_clip_LC(self):
    #     '''
    #     This function clips the light curve 

    #     Parameters
    #     ----------

    #     Returns
    #     -------
    #     '''


