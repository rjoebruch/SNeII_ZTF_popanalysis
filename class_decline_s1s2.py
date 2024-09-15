'''

This class fits a linear decline after peak to standard SN II 
Update of 16/01/2022: we impose a bound at the end of the fit so taht f(t_peak) = F_peak 


IMPORTANT NOTE: this verison of the class is using the new version of iminuit (>2) which has the "propagrate" option

TODO: - change the error computation type from propagation to bootstrap, see the docs, using a multivariate gaussian
- add the residual from the fit. if a data point is 0.1 mag * error away, mark it in red on the residual plot. 



'''
from library_4_bazin_functions import *
from iminuit.util import propagate
from astropy import cosmology


class Fit_s1s2_Decline( object ):

    

    '''

        This class is defined to fit the decline from peak of standard SNe II

        This class is used in MAGNITUDES

    '''

    def __init__(self, table, min_fit, max_fit, filter = 'r', zez = None, peakmag = None):
        '''
        Parameters
        ----------
        table [table] contraining the time from first marshal detection, the flux and the error on the flux 

        min_fit [float] minimal bound of the fit - peak time
        max_fit [float] maximal bound of the fit
        mini    [flaot] minimal bound of data considered
        maxi    [flaot] maximal bound of data considered for the final chi_2 fit 
        zez     [list] values of REDSHIFT and ERROR ON REDSHIFT [z, e_z]

        '''

        self.cosmo = cosmology.Planck13
        

        self.set_peak_mag_value(peakmag)
        # self._time_zone_selected = False
        self._filter_name   = {'g':'darkolivegreen','r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        self.set_filter(filter)

        # reduced data to perform the fit onto
        self.set_time_zone_fit(min_fit, max_fit)
        
        self.set_time_from_eed_fit()
        self.set_mag_fit()
        self.set_e_mag_fit()

        # self.set_mag()
        # self.set_e_mag()

        self.set_redshift_DM(zez=zez)
        

    #-----------#
    #  SETTERS  #
    #-----------#

    def set_redshift_DM(self, zez):
        '''
        This function sets the redshift and DM based on input
        '''
        if zez is not None:
            if type(zez) ==list:
                self._z   = zez[0]
                self._e_z = zez[1]

                self._DM           = self.cosmo.distmod(z = self._z).value
                self._e_DM         = (5*self._e_z)/(np.log(10)* self._z)
            else:
                print('ZEZ needs to be a list to set the DM and redshifst')
        else:
            print('Not setting redshift or DM, no value provided ')


    def set_peak_mag_value(self, peakmag):
        '''
        This functions sets the boundary of the fit so that f(t_peak = max_fit) == peak_flux
        '''
        
        self._peak_mag = peakmag

    # def _set_max_fit_val(self, max_fit):
    #     '''
    #     This function sets the temporal boundary of the fit where max_fit == t_peak 
    #     t_peak should be determined prior to the fit using methods in "script_mag_finders" 
    #     '''
    #     self._max_fit = max_fit

    def set_filter(self,filter):
        '''
        This function sets the filter we are working in. 'r' or 'g' 
        '''
        self.filter = filter

    def set_phot_table(self, table):
        '''
        This function sets the photometric table that we will use throughout the class.
        '''
        self.table = table


    # def set_mag(self):
    #     '''
    #     '''

    #     self.mag    = self.fittable['mag']
    #     self.absmag = self.fittable['absmag']

    # def set_e_mag(self):
    #     '''
    #     '''

    #     self.e_mag    = self.fittable['emag']
    #     self.e_absmag = self.fittable['e_absmag']




    def set_time_zone_fit(self, min_fit , max_fit ):
        '''
        This function cuts the table to the temporal zone relevant for the  fit
        '''

        self.fittable = self.table[(self.table['tfromexplo_zc']>=min_fit)&(self.table['tfromexplo_zc']<=max_fit)]

        self._time_zone_selected = True

    def set_time_from_eed_fit(self):
        '''
        This function sets the time array in the relevant temporal zone from the estimated explosion time
        '''
        self.t_fit = self.fittable['tfromexplo_zc']

    def set_mag_fit(self):
        '''
        This function sets the flux array relevant for the fit
        '''

        self.f_fit = self.fittable['mag']

    def set_e_mag_fit(self):
        '''
        This function sets the errors on the flux relevant for the fit.
        '''
        self.e_f_fit = self.fittable['emag']


    ############
    # GETTERS  #
    ############

    
    def _get_chi2minuit(self, a1,t0,m0,a2):
        '''
        a1
        b1

        t_peak [fixed]

        a2
        b2

        '''

        # param = [a,b,t_peak]

        # return self.get_chi_2(param)

        model_    = function_s1_s2(self.t_fit,a1,t0,m0,a2)
        pull_     = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq    = np.sum(pull_)

        num_param = 1
        dof       = len(self.f_fit) - num_param
        if dof == 0: #avoid division by zero? 
            dof = 1

        return  chi_sq / dof




    #######################
    #        FITTER       #
    #######################
    
    def set_minuit(self, guess, boundaries=None, fixed=None, step_sizes = None, print_level=0):
        '''
        This function sets the optimiser Minuit before running migrad 

        parameters
        ----------
        guess (tuple) starting point of the optimisation
        boundaries None or dictionnary 
        fixed    list of names of the variables you want to fix
        errordef 
        print_level
        step_size   None or dictionnary

        
        '''

        ######################### OLD VERSION OF CODING IMINUIT  ##################################################

        # A, alpha1, t_break, alpha2, n, t_peak, peak_val   = guess
        # self._paramname = "A,alpha1,t_break,alpha2,n,t_peak,peak_val".split(",")
        
        # local_ = locals() # dict with the local variables

        # self.minuit_kwargs = {}
        # for p in self._paramname:
        #     self.minuit_kwargs[p] = local_[p]
            
        # if boundaries is not None:
        #     for k,v in boundaries.items():
        #         self.minuit_kwargs["limit_"+k] = v
                
        # if fixed is not None:
        #     for k in fixed:
        #         self.minuit_kwargs["fix_"+k] = True
        
        # if step_sizes is not None: 
        #     for k,v in step_sizes.items():
        #         self.minuit_kwargs['error_'+k] = v
            

        # self.minuit = Minuit(self._get_chi2minuit, errordef=errordef, 
        #                             print_level=print_level, **self.minuit_kwargs)
        #############################################################################################################################

        
        a1,t0,m0,a2   = guess
        self._paramname = "a1,t0,m0,a2".split(",")

        local_ = locals() # dict with the local variables

        self.minuit_kwargs = {}
        for p in self._paramname:
            self.minuit_kwargs[p] = local_[p]


        # initialise the optimiser: we are using our own cost function defined as a reduced chi_2 and we feed it the guesses
        self.minuit = Minuit(self._get_chi2minuit, **self.minuit_kwargs)

        # since we are using a least square type of cost funtion, we need to define the error level definition here

        self.minuit.errordef = Minuit.LEAST_SQUARES

        self.minuit.print_level = print_level

        # we then set all the parameters supposedly fixed, their limits and step sizes 
        if boundaries is not None:
            for _ in boundaries.keys():
                self.minuit.limits[_] = boundaries[_]

        if fixed is not None:
            for _ in fixed:
                self.minuit.fixed[_] = True

        if step_sizes is not None: 
            for _ in step_sizes.keys():
                self.minuit.errors[_] = step_sizes[_]

    

    def fit_minuit(self, guess, boundaries=None, fixed=None, step_sizes=None):
        
        """ 
    
        """
        
        self.set_minuit(guess=guess, boundaries=boundaries, fixed=fixed, step_sizes = step_sizes )
        self.minuit_output = self.minuit.migrad()
        print(self.minuit_output)

    
    def error_propag(self, min_interp = None, max_interp = None, nump = 500 , returnop = False ) :
        '''
        this function propagates the error on the minuit fit and returns the estimated errors 

        parameters
        ----------
        min_interp [float]
        max_interp [float]
        nump       [float] number of points to interpolate, pd 500
        
        '''
        if (min_interp, max_interp) == (None, None):
            min_interp = min(self.t_fit)
            max_interp = max(self.t_fit)

        self.t_res = np.linspace(min_interp, max_interp, nump)

        self.fit_mag_res , self.fit_mag_cov_res  = propagate(lambda param: function_s1_s2(self.t_res, *param), self.minuit.values, self.minuit.covariance)
        
        self.fit_e_mag_res = np.diag(self.fit_mag_cov_res)**0.5

        if returnop is True:
            return self.t_res, self.fit_mag_res, self.fit_e_mag_res  


    def error_via_parametric_bootstrap(self, min_interp = None, max_interp = None, nump = 1000 , returnop = False ) :
        '''
        this function propagates the error on the minuit fit and returns the estimated errors 

        parameters
        ----------
        min_interp [float]
        max_interp [float]
        nump       [float] number of points to interpolate, pd 500
        
        '''
        if (min_interp, max_interp) == (None, None):
            min_interp = min(self.t_fit)
            max_interp = max(self.t_fit)

        self.t_res       = np.linspace(min_interp, max_interp, nump)

        self.fit_mag_res = function_s1_s2(self.t_res, *self.minuit.values[:])

        rng   = np.random.default_rng(1)
        par_b = rng.multivariate_normal(self.minuit.values, self.minuit.covariance, size=1000)

        y_b       = [function_s1_s2(self.t_res, *param) for param in par_b]
        
        self.fit_e_mag_res = np.std(y_b,axis=0)

        # self.fit_mag_res , self.fit_mag_cov_res  = propagate(lambda param: function_s1_s2(self.t_res, *param), self.minuit.values, self.minuit.covariance)
        
        # self.fit_e_mag_res = np.diag(self.fit_mag_cov_res)**0.5

        if returnop is True:
            return self.t_res, self.fit_mag_res, self.fit_e_mag_res 
    
    # def conv_mag_fit(self, returnop = False):

    #     '''
        
    #     '''
    #     self.fit_mag_res          = to_mag(self.fit_flux_res)
        
    #     self.fit_e_mag_res        = 1.0857 * self.fit_e_flux_res / self.fit_flux_res

    #     if returnop is True:
    #         return self.fit_mag_res, self.fit_e_mag_res

    def conv_abs_mag_fit(self):
        '''
        this function converts to absolute magnitudes

        '''

        self.fit_absmag_res    =  self.fit_mag_res - self._DM
                    
        self.fit_e_absmag_res  = np.sqrt(self.fit_e_mag_res **2 + self._e_DM**2)

    def get_fit_table_result(self, min_interp, max_interp):
        '''
        THis function returns the table with the interpolated light curve in flux, absolute magnitude and apparent magnitude

        parameters
        ----------
        min_interp [float] minimuma t which you want to perform the interpolation
        max_interp [float] maximum at which you want to perform the interpolation

        note
        ----
        /!\ min_interp cannot be smaller than min(self.t_fit)
        /!\ max_interp cannot be smaller than max(sel.t_fit)

        '''

        if min_interp < min(self.t_fit):
            print('/!\ min_interp cannot be smaller than min(self.t_fit)')
            return None
        else:
            if max_interp > max(self.t_fit):
                print('/!\ max_interp cannot be smaller than max(sel.t_fit)')
                return None
            else:
                self.error_propag(min_interp=min_interp, max_interp=max_interp)
                self.conv_mag_fit()
                self.conv_abs_mag_fit()

                interp_table = Table([self.t_res,self.fit_flux_res,self.fit_e_flux_res,self.fit_mag_res,self.fit_e_mag_res,self.fit_absmag_res,self.fit_e_absmag_res],
                                      names=('t', 'flux', 'e_flux', 'mag', 'emag', 'absmag', 'eabsmag'))
                return interp_table

    
    

    ##################
    #   PLOTTING     #
    ##################



    def plot_fit(self ,add_plot = False, **kwargs):
        '''

        This function plots the results of the optimisation with migrad

        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_         = np.linspace(min(self.t_fit), max(self.t_fit), 1000)
        

        _params = [self.minuit_output.params[x].value for x in range(len(self.minuit_output.params))]
        # _errors = [self.minuit_output.params[x].error for x in range(3)]
        
       

        chi2_fit       = self.get_redchi_2_fit(_params)
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        
        t_ , fit, e_fit = self.error_via_parametric_bootstrap(returnop=True)



        if add_plot == None:
            plt.figure()


            plt.plot(t_ , fit, ls = '--' , color = 'black' )
            plt.fill_between(t_, fit - e_fit, fit + e_fit,color = 'grey', alpha=0.5)
            
            plt.errorbar(self.table['tfromexplo_zc'], self.table['mag'], self.table['emag'], fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter), alpha=0.2,**kwargs )
            # plt.plot(x_ , function_s1_s2(x_, *_params ), ls = '--' , color = _fitcolor )

            

            
            # plt.plot(x_[1] , function_s1_s2(x_, *_params)[1], color = 'white',label= chi_fit_string )


            
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)
            # plt.axvline(_params[2], lw = 0.5, color = 'red')

            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Apparent Magnitude', size = 15)

            plt.gca().invert_yaxis()
            # plt.legend()
        
        else: 
            add_plot.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
            add_plot.plot(x_ ,function_s1_s2(x_, *_params), ls = '--'  , color = _fitcolor )

            
            add_plot.plot(x_[1] , function_s1_s2(x_,*_params)[1], color = 'white',label= chi_fit_string )



            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Apparent Magnitude', size = 15)

            plt.gca().invert_yaxis()

            plt.legend(fontsize = 10)


    def plot_fit_mag_normalised(self, add_plot = False):
        '''

        This function plots the results of the optimisation with migrad

        '''

        if self._peak_mag is not None: 

            _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
            _fitcolor = random.choice(_colors)


            x_         = np.linspace(min(self.t_fit), max(self.t_fit), 1000)


            _params = [self.minuit_output.params[x].value for x in range(3)]
            _errors = [self.minuit_output.params[x].error for x in range(3)]



            chi2_fit       = self.get_redchi_2_fit(_params)
            chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'


            t_ , fit, e_fit = self.error_propag(returnop=True)



            if add_plot == None:
                plt.figure()
                plt.errorbar(self.table['tfromexplo_zc'], self.table['mag']-self._peak_mag, self.table['emag'], fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter) )
                plt.plot(x_ , function_s1_s2(x_, *_params )-self._peak_mag, ls = '--' , color = _fitcolor )

                plt.plot(t_ , fit-self._peak_mag, ls = '--' , color = 'black' )
                plt.fill_between(t_, fit-self._peak_mag - e_fit, fit-self._peak_mag + e_fit, color = 'grey', alpha=0.5)


                plt.plot(x_[1] , function_s1_s2(x_, *_params)[1]-self._peak_mag, color = 'white',label= chi_fit_string )



                # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)
                # plt.axvline(_params[2], lw = 0.5, color = 'red')

                plt.xlabel('Time from EED [days]', size = 15)
                plt.ylabel('Apparent Magnitude', size = 15)

                plt.gca().invert_yaxis()
                plt.legend()

            else: 
                add_plot.errorbar(self.table['tfromexplo_zc'], self.table['mag']-self._peak_mag, self.table['emag'], fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
                
                add_plot.plot(x_ , function_s1_s2(x_, *_params )-self._peak_mag, ls = '--'  , color = _fitcolor )
                add_plot.fill_between(t_, fit-self._peak_mag - e_fit, fit-self._peak_mag + e_fit,color = 'grey', alpha=0.5)


                # add_plot.plot(x_[1] , function_s1_s2(x_,*_params)[1], color = 'white',label= chi_fit_string )



                # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


                plt.xlabel('Time from EED [days]', size = 15)
                plt.ylabel('Apparent Magnitude', size = 15)

                plt.gca().invert_yaxis()

                plt.legend()
        else:
            print('No peak mag provided when you initiated the class!')
        
        
    def get_redchi_2_fit(self, param):
        '''
        parameters
        ----------
        self

        returns
        -------
        ''' 
        #(a,b,t_peak)= param 
        model_  = function_s1_s2(self.t_fit,*param)
        pull_   = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq  = np.sum(pull_)
        n_param = 4
        dof     = len(self.f_fit)-n_param


        return  chi_sq/dof



    def get_chi2_fit(self, param):
        '''
        parameters
        ----------
        self
        param

        returns
        -------
        ''' 
        #(a,b,t_peak )= param 
        model_ = function_s1_s2(self.t_fit,*param)
        pull_  = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq = np.sum(pull_)
        # dof    = len(self.f_fit)-len(param)

        return  chi_sq

    def proba_chi_2(chi2, N):
        '''
        this function computes the probability for a chi_2 value 
        
        '''

        #### the PDF of Chi2 
        def pdf_chi2(N):
            return lambda chizwei : (1 /( 2**(N/2) * gamma(N/2) )) * (chizwei**(N/2 - 1)) * np.exp(-chizwei/2)

        return quad(pdf_chi2, chi2, np.inf)


