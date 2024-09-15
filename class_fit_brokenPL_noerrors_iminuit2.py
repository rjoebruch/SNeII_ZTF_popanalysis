'''

This class fits a broken power law for SN II?
Update of 16/01/2022: we impose a bound at the end of the fit so taht f(t_peak) = F_peak 
                      see function "bounded_Broken_PL" in the function file

Important note: here we fix n, t_peak and peak_val. The number of parameters is 4 

IMPORTANT NOTE: this verison of the class is using the new version of iminuit (>2) which has the "propagrate" option


VERSION OF 26/05/2022: 
The error propagation method using the Jacobian (mapping between input and output) fails often and does not reflect the reality. 
Moreover, the fit sometimes fail (cf ZTF21abkygyy)

We are updating the script and archituecture of the code to opt for a boot strapping method. 

HOWEVER /!!\ There exist in iminuit a bootstrapping method to evalutate the errors of the fit, but the bootstraping relies on varying the parameters of the fit within their own errors
It is not clear that the errors on the parameters as output of the fit arerealistic. This should be checked... 

We are hence opting for a manual bootstrapping method which involves feeding to iminuit random light curves and estimating the median and std light curve from this experiment. 




'''
from library_4_bazin_functions import *
from iminuit.util import propagate
from astropy import cosmology


class Fit_broken_PL( object ):

    

    '''

        This class is defined to fit the early LC of Infant SNe II in order to estimate the explosion parameters

    '''

    def __init__(self, table, min_fit, max_fit, peak_val, filter = 'r', zez = None):
        '''
        Parameters
        ----------
        table [table] contraining the time from first marshal detection, the flux and the error on the flux 

        min_fit [float] minimal bound of the fit
        max_fit [float] maximal bound of the fit
        mini    [flaot] minimal bound of data considered
        maxi    [flaot] maximal bound of data considered for the final chi_2 fit 
        zez     [list] values of REDSHIFT and ERROR ON REDSHIFT [z, e_z]

        '''

        self.cosmo = cosmology.Planck13
        

        self.set_peak_flux_value(peak_val=peak_val)
        self._filter_name   = {'g':'darkolivegreen','r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        self.set_filter(filter)

        # reduced data to perform the fit onto
        # self.set_time_zone_fit(min_fit, max_fit)
        
        self.set_time_from_eed_fit()
        self.set_flux_fit()
        # self.set_e_flux_fit()
        self.set_mag()
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


    def set_peak_flux_value(self, peak_val):
        '''
        This functions sets the boundary of the fit so that f(t_peak = max_fit) == peak_flux
        '''

        self._peak_val = peak_val

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


    # def set_time_zone_fit(self, min_fit , max_fit ):
    #     '''

    #     '''

    #     self.fittable = self.table[(self.table['tfromexplo_zc']>=min_fit)&(self.table['tfromexplo_zc']<=max_fit)]

    #     self._time_zone_selected = True

    def set_time_from_eed_fit(self):
        '''
        '''
        self.t_fit = self.table[0]

    def set_flux_fit(self):
        '''
        '''

        self.f_fit = self.table[1]


    def set_mag(self):
        '''
        WARNING: since we are generating random datapoints in FLUX space, we need to convert this flux into magnitude and not ntake the "pre" mag

        
        '''

        self.mag    = to_mag(self.f_fit)

    ############
    # GETTERS  #
    ############

    
    def _get_chi2minuit(self, A, alpha1, t_break, alpha2, n, t_peak, peak_val):
        '''
        This cost function does not involve the errors!! This is only for the bootstrapping method! 


        A
        alpha1
        t_break
        alpha2
        n fixed = 1
        t_peak, peak_val [fixed]

        '''

        # param = [A, alpha1, t_break, alpha2, n, t_peak, peak_val]


        model_    = bounded_broken_PL(self.t_fit,A, alpha1, t_break, alpha2, n, t_peak, peak_val)
        # pull_     = (( self.f_fit - model_ )/ self.e_f_fit )**2
        pull_     = ( self.f_fit - model_ )**2
        chi_sq    = np.sum(pull_)

        num_param = 4 # the actually really varying parameters. n, tpeak and peak val are fixed 
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

        
        A, alpha1, t_break, alpha2, n, t_peak, peak_val   = guess
        self._paramname = "A,alpha1,t_break,alpha2,n,t_peak,peak_val".split(",")

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

    

    ''' 
    Removing the error propagation from the fit. Over the full sample, it does not turn out to be very useful... 
    '''

    # def error_propag(self, min_interp = None, max_interp = None, nump = 500 , returnop = False ) :
    #     '''
    #     this function propagates the error on the minuit fit and returns the estimated errors 

    #     parameters
    #     ----------
    #     min_interp [float]
    #     max_interp [float]
    #     nump       [float] number of points to interpolate, pd 500
        
    #     '''
    #     if (min_interp, max_interp) == (None, None):
    #         min_interp = min(self.t_fit)
    #         max_interp = max(self.t_fit)

    #     self.t_res = np.linspace(min_interp, max_interp, nump)

    #     self.fit_flux_res , self.fit_flux_cov_res  = propagate(lambda param: bounded_broken_PL(self.t_res, *param), self.minuit.values, self.minuit.covariance)
        
    #     self.fit_e_flux_res = np.diag(self.fit_flux_cov_res)**0.5

    #     if returnop is True:
    #         return self.t_res, self.fit_flux_res, self.fit_e_flux_res  

    def generate_fit(self, min_interp, max_interp, nump = 500 ,returnop = False):

        '''
        This function generates the function form the fit parameter AFTER optimisation
         for nump = 500 (pd) datapoints

        parameters
        ----------

        min_interp [float]  lower bound of the fit
        max interp [float]  upper bound of the fit
        nump                number of data points 
        returnop   [bool]   whether to return the magnitude 
        '''

        self.t_res                = np.linspace(min_interp, max_interp, nump)

        self.fit_flux_res         = bounded_broken_PL(self.t_res , *self.minuit.values )

        # self.fit_e_mag_res        = 1.0857 * self.fit_e_flux_res / self.fit_flux_res

        if returnop is True:
            return self.t_res, self.fit_flux_res

    def conv_mag_fit(self,returnop = False):
        '''
        this function converts to magnitude and absolute magnitude absolute magnitudes

        Removing the error... this will be somethign external 

        '''

        self.fit_mag_res          = to_mag(self.fit_flux_res)


        self.fit_absmag_res    =  self.fit_mag_res - self._DM

        if returnop is True:
            return self.fit_mag_res, self.fit_absmag_res
                    
        # self.fit_e_absmag_res  = np.sqrt(self.fit_e_mag_res **2 + self._e_DM**2)

    # def get_fit_table_result(self, min_interp, max_interp):
    #     '''
    #     THis function returns the table with the interpolated light curve in flux, absolute magnitude and apparent magnitude

    #     parameters
    #     ----------
    #     min_interp [float] minimuma t which you want to perform the interpolation
    #     max_interp [float] maximum at which you want to perform the interpolation

    #     note
    #     ----
    #     /!\ min_interp cannot be smaller than min(self.t_fit)
    #     /!\ max_interp cannot be smaller than max(sel.t_fit)

    #     '''

    #     if min_interp < min(self.t_fit):
    #         print('/!\ min_interp cannot be smaller than min(self.t_fit)')
    #         return None
    #     else:
    #         if max_interp > max(self.t_fit):
    #             print('/!\ max_interp cannot be smaller than max(sel.t_fit)')
    #             return None
    #         else:
    #             self.error_propag(min_interp=min_interp, max_interp=max_interp)
    #             self.conv_mag_fit()
    #             self.conv_abs_mag_fit()

    #             interp_table = Table([self.t_res,self.fit_flux_res,self.fit_e_flux_res,self.fit_mag_res,self.fit_e_mag_res,self.fit_absmag_res,self.fit_e_absmag_res],
    #                                   names=('t', 'flux', 'e_flux', 'mag', 'emag', 'absmag', 'eabsmag'))
    #             return interp_table

    
    

    ##################
    #   PLOTTING     #
    ##################



    def plot_fit(self, add_plot = False):
        '''

        This function plots the results of the optimisation with migrad, version without the errors since thsi class is for only generating a fit 

        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)

        x_         = np.linspace(0, max(self.t_fit), 1000)
        _params = [self.minuit_output.params[x].value for x in range(7)]
        
        chi2_fit       = self.get_redchi_2_fit(_params)
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        # t_ , fit_  = self.generate_fit_and_conv_mag(returnop=True)

        if add_plot == None:
            plt.figure()
            plt.plot(self.t_fit, self.f_fit, '.', color = self._filter_name.get(self.filter) )
            plt.plot(x_ , bounded_broken_PL(x_, *_params ), ls = '--' , color = _fitcolor )


            plt.plot(x_[1] , bounded_broken_PL(x_, *_params)[1], color = 'white',label= chi_fit_string )
            plt.axvline(_params[2], lw = 0.5, color = 'red')

            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Flux', size = 15)

            plt.legend()
        
        else: 
            add_plot.plot(self.t_fit, self.f_fit, '.' , color = self._filter_name.get(self.filter) )
            add_plot.plot(x_ ,bounded_broken_PL(x_, *_params), ls = '--'  , color = _fitcolor )
            add_plot.plot(x_[1] , bounded_broken_PL(x_,*_params)[1], color = 'white',label= chi_fit_string )

            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Flux', size = 15)


            plt.legend(fontsize = 10)
        
        



    def plot_fit_conv_mag(self):
        '''
        This function converts the fit found to magnitude and overlays it on the actual magnitude measurements

        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)

        x_         = np.linspace(0, max(self.t_fit), 1000)
        
        _params = [self.minuit_output.params[x].value for x in range(7)]

        chi2_fit       = self.get_redchi_2_fit(_params)
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        
        plt.figure()
        plt.plot(self.t_fit, self.mag,'o', ms=3.5 , color = self._filter_name.get(self.filter) )
        plt.plot(x_ , to_mag(bounded_broken_PL(x_,*_params)), ls = '--' , color = _fitcolor )
        

        plt.xlabel('Time from EED [days]', size = 15)
        plt.ylabel('Apparent Magnitude', size = 15)
        # plt.ylim([min(mfit)-0.3,21])

        plt.gca().invert_yaxis()
        # plt.legend()
        
       


    
    

    def get_redchi_2_fit(self, param):
        '''
        parameters
        ----------
        self

        returns
        -------
        ''' 
        #(A_, alpha1_, t_break_, alpha2_, n_,t_peak_, peak_val_)= param 
        model_  = bounded_broken_PL(self.t_fit,*param)
        pull_   = (( self.f_fit - model_ ))**2
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
        #(A_, alpha1_, t_break_, alpha2_, n_,t_peak_, peak_val_ )= param 
        model_ = bounded_broken_PL(self.t_fit,*param)
        pull_  = (( self.f_fit - model_ ) )**2
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


