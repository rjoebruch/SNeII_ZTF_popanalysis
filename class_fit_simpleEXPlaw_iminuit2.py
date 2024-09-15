'''

This class fits a  simple power law for SN II ar early time
Update 16/01/2022 : we bound the early time fit so that f(t_peak) = peak flux


Important note   : here we fix t_peak and peak_val. The number of degree of freedom is 2
other important
note             : on the 17/01/2022 we updated iminuit to version 2 and had to re-adapt the whole code...

'''
from library_4_bazin_functions import *
from iminuit.util import propagate
from astropy import cosmology


class Fit_simple_EXP_law( object ):

    

    '''

        This class is defined to fit the early LC of Infant SNe II in order to estimate the explosion parameters

    '''

    def __init__(self, table, min_fit, max_fit, peak_val,  filter = 'r', zez = None ):
        '''
        Parameters
        ----------
        table contraining the time from first marshal detection, the flux and the error on the flux 

        min_fit [float] minimal bound of the fit
        max_fit [float] maximal bound of the fit
        mini    [flaot] minimal bound of data considered
        maxi    [flaot] maximal bound of data considered for the final chi_2 fit 




        '''

        self.cosmo = cosmology.Planck13
    


        self.set_peak_flux_value(peak_val)
        # self._time_zone_selected = False
        self._filter_name   = {'g':'darkolivegreen','r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        self.set_filter(filter)

        # full data to comapre the chi_2 to
        # self.set_time_zone_interest(mini, maxi)
        # self.set_time_from_eed()
        


        # reduced data to perform the fit onto
        self.set_time_zone_fit(min_fit, max_fit)
        # if len(self.fittable) <

        self.set_time_from_eed_fit()
        self.set_flux_fit()
        self.set_e_flux_fit()

        self.set_mag()
        self.set_e_mag()

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


    def set_filter(self,filter):
        '''
        '''
        self.filter = filter

    def set_phot_table(self, table):
        '''

        '''
        self.table = table


    def set_time_zone_fit(self, min_fit , max_fit ):
        '''

        '''

        self.fittable = self.table[(self.table['tfromexplo_zc']>=min_fit)&(self.table['tfromexplo_zc']<=max_fit)]

        self._time_zone_selected = True

    def set_time_from_eed_fit(self):
        '''
        '''
        self.t_fit = self.fittable['tfromexplo_zc']

    def set_flux_fit(self):
        '''
        '''

        self.f_fit = self.fittable['extcorrforcediffimflux']

    def set_e_flux_fit(self):
        '''
        '''
        self.e_f_fit = self.fittable['extcorrforcediffimfluxunc']


    def set_mag(self):
        '''
        '''

        self.mag    = self.fittable['mag']
        self.absmag = self.fittable['absmag']

    def set_e_mag(self):
        '''
        '''

        self.e_mag    = self.fittable['emag']
        self.e_absmag = self.fittable['e_absmag']

    ############
    # GETTERS  #
    ############

    
    def _get_chi2minuit(self, A, n,t_peak, peak_val):
        '''
        A
        n

        '''

        # param = [A, n,t_peak, peak_val]

        # return self.get_chi_2(param)

        model_    = bounded_simple_EXP_law(self.t_fit,A, n,t_peak, peak_val)
        pull_     = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq    = np.sum(pull_)

        num_param = 2

        dof    = len(self.f_fit) - num_param
        if dof == 0: #avoid division by zero? 
            dof = 1

        return  chi_sq / dof



    #######################
    #        FITTER       #
    #######################
    
    def set_minuit(self, guess, boundaries=None, fixed=None, errordef=1 ,print_level=0):
        '''
        
        '''

        A, n, t_peak, peak_val   = guess
        self._paramname = "A,n,t_peak,peak_val".split(",")
        
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


        

    def fit_minuit(self, guess, boundaries=None, fixed=None):
        
        """ 
    
        """
        
        self.set_minuit(guess, boundaries, fixed)
        self.minuit_output = self.minuit.migrad()
        print(self.minuit_output)

    


    def error_propag(self, min_interp = None, max_interp = None, nump = 500, returnop = False) :
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

        self.fit_flux_res , self.fit_flux_cov_res  = propagate(lambda param: bounded_simple_EXP_law(self.t_res, *param), self.minuit.values, self.minuit.covariance)
        
        self.fit_e_flux_res = np.diag(self.fit_flux_cov_res)**0.5
        
        if returnop == True:
            return self.t_res, self.fit_flux_res, self.fit_e_flux_res 


    def conv_mag_fit(self, returnop = False):

        '''
        this function converts to magnitudes the resulting fit with error propafation 
        
        '''
        self.fit_mag_res          = to_mag(self.fit_flux_res)
        
        self.fit_e_mag_res        = 1.0857 * self.fit_e_flux_res / self.fit_flux_res

        if returnop == True:
            return self.fit_mag_res, self.fit_e_mag_res

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



    def plot_fit(self, add_plot = False):
        '''
        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_         = np.linspace(0, max(self.t_fit), 1000)
        
        A_         = self.minuit_output.params[0].value
        n_         = self.minuit_output.params[1].value
        t_peak_    = self.minuit_output.params[2].value
        peak_val_  = self.minuit_output.params[3].value
        

        chi2_fit       = self.get_redchi_2_fit([A_, n_ , t_peak_, peak_val_])
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        
        t_, fit, e_fit = self.error_propag(returnop=True)



        if add_plot == None:
            plt.figure()
            plt.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter) )
            plt.plot(x_ , bounded_simple_EXP_law(x_, A_,  n_ ,t_peak_, peak_val_), ls = '--' , color = _fitcolor )

            
            plt.plot(t_ , fit, ls = '--' , color = 'black' )
            plt.fill_between(t_, fit - e_fit, fit + e_fit,color = 'grey', alpha=0.5)

            
            plt.plot(x_[1] , bounded_simple_EXP_law(x_,A_, n_,t_peak_, peak_val_ )[1], color = 'white',label= chi_fit_string )


            
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)
            # plt.axvline(t_break_, lw = 0.5, color = 'red')

            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Flux', size = 15)

            # plt.gca().invert_yaxis()
            plt.legend()
        
        else: 
            add_plot.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
            add_plot.plot(x_ ,bounded_simple_EXP_law(x_, A_, n_,t_peak_, peak_val_), ls = '--'  , color = _fitcolor )

            
            add_plot.plot(x_[1] , bounded_simple_EXP_law(x_,A_, n_,t_peak_, peak_val_)[1], color = 'white',label= chi_fit_string )



            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Flux', size = 15)

            # plt.gca().invert_yaxis()

            plt.legend(fontsize = 10)


    def plot_current_parameters(self, add_plot = False):
        '''
        This function plots the current parameters, and not the results of the migrad optimisation
        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_         = np.linspace(0, max(self.t_fit), 1000)
        
        A_         = self.minuit.params[0].value
        n_         = self.minuit.params[1].value
        t_peak_    = self.minuit.params[2].value
        peak_val_  = self.minuit.params[3].value
        

        chi2_fit       = self.get_redchi_2_fit([A_, n_ , t_peak_, peak_val_])
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        
        t_, fit, e_fit = self.error_propag(returnop=True)



        if add_plot == None:
            plt.figure()
            plt.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter) )
            plt.plot(x_ , bounded_simple_EXP_law(x_, A_,  n_ ,t_peak_, peak_val_), ls = '--' , color = _fitcolor )

            
            plt.plot(t_ , fit, ls = '--' , color = 'black' )
            plt.fill_between(t_, fit - e_fit, fit + e_fit,color = 'grey', alpha=0.5)

            
            plt.plot(x_[1] , bounded_simple_EXP_law(x_,A_, n_,t_peak_, peak_val_ )[1], color = 'white',label= chi_fit_string )


            
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)
            # plt.axvline(t_break_, lw = 0.5, color = 'red')

            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Flux', size = 15)

            # plt.gca().invert_yaxis()
            plt.legend()
        
        else: 
            add_plot.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
            add_plot.plot(x_ ,bounded_simple_EXP_law(x_, A_, n_,t_peak_, peak_val_), ls = '--'  , color = _fitcolor )

            
            add_plot.plot(x_[1] , bounded_simple_EXP_law(x_,A_, n_,t_peak_, peak_val_)[1], color = 'white',label= chi_fit_string )



            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Flux', size = 15)

            # plt.gca().invert_yaxis()

            plt.legend(fontsize = 10)


 


    def plot_fit_conv_mag(self):
        '''
        This function converts the fit found to magnitude and overlays it on the actual magnitude measurements

        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_         = np.linspace(0, max(self.t_fit), 1000)
        

        _params = [self.minuit_output.params[x].value for x in range(4)]
        _errors = [self.minuit_output.params[x].error for x in range(4)]

        chi2_fit       = self.get_redchi_2_fit(_params)
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        mfit, e_mfit = self.conv_mag_fit(returnop=True)

        
        plt.figure()
        plt.errorbar(self.t_fit, self.mag, self.e_mag, fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter) )

        plt.plot(x_ , to_mag(bounded_simple_EXP_law(x_,*_params)), ls = '--' , color = _fitcolor )

        plt.plot(self.t_res , mfit, ls = '--' , color = 'black' )
        plt.fill_between(self.t_res, mfit - e_mfit, mfit + e_mfit,color = 'grey', alpha=0.5)
        
        # plt.plot(x_[1] , bounded_simple_EXP_law(x_,*_params )[1], color = 'white',label= chi_fit_string )
        
        # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)
        # plt.axvline(_params[2], lw = 0.5, color = 'red')

        plt.xlabel('Time from EED [days]', size = 15)
        plt.ylabel('Apparent Magnitude', size = 15)
        plt.ylim([min(mfit)-0.3,21])
        
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
        #(A_, n_)= param 
        model_  = bounded_simple_EXP_law(self.t_fit,*param)
        pull_   = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq  = np.sum(pull_)
        n_param = 2
        dof     = len(self.f_fit)-n_param

        return  chi_sq/dof