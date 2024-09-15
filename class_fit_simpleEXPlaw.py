'''

This class fits a  simple power law for SN II ar early time
Update 16/01/2022 : we bound the early time fit so that f(t_peak) = peak flux


Important note   : here we fix t_peak and peak_val. The number of degree of freedom is 2


'''
from library_4_bazin_functions import *



class Fit_simple_EXP_law( object ):

    

    '''

        This class is defined to fit the early LC of Infant SNe II in order to estimate the explosion parameters

    '''

    def __init__(self, table, min_fit, max_fit, peak_val,  filter = 'r'):
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

        self.set_peak_flux_value(peak_val)
        # self._time_zone_selected = False
        self._filter_name   = {'g':'darkolivegreen','r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        self.set_filter(filter)

        # full data to comapre the chi_2 to
        # self.set_time_zone_interest(mini, maxi)
        # self.set_time_from_eed()
        # self.set_mag()
        # self.set_e_mag()


        # reduced data to perform the fit onto
        self.set_time_zone_fit(min_fit, max_fit)
        # if len(self.fittable) <

        self.set_time_from_eed_fit()
        self.set_mag_fit()
        self.set_e_mag_fit()

        

    #-----------#
    #  SETTERS  #
    #-----------#
    def set_peak_flux_value(self, peak_val):
        '''
        This functions sets the boundary of the fit so that f(t_peak = max_fit) == peak_flux
        '''

        self._peak_val = peak_val

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

    # def set_time_zone_interest(self, mini = 0, maxi = 50):
    #     '''
    #     the full LC here
    #     '''
    #     self.table = self.table[(self.table['tfromexplo_zc']>=mini)&(self.table['tfromexplo_zc']<=maxi)]

    # def set_time_from_eed(self):
    #     '''
    #     '''
    #     self.t = self.table['tfromexplo_zc']

    def set_mag(self):
        '''
        '''

        self.f = self.fittable['mag']

    def set_e_mag(self):
        '''
        '''

        self.e_f = self.fittable['emag']




    def set_time_zone_fit(self, min_fit , max_fit ):
        '''

        '''

        self.fittable = self.table[(self.table['tfromexplo_zc']>=min_fit)&(self.table['tfromexplo_zc']<=max_fit)]

        self._time_zone_selected = True

    def set_time_from_eed_fit(self):
        '''
        '''
        self.t_fit = self.fittable['tfromexplo_zc']

    def set_mag_fit(self):
        '''
        '''

        self.f_fit = self.fittable['extcorrforcediffimflux']

    def set_e_mag_fit(self):
        '''
        '''
        self.e_f_fit = self.fittable['extcorrforcediffimfluxunc']


    ############
    # GETTERS  #
    ############


    # def get_fit_roi_length(self):
    #     '''
    #     This function gets the length of the region of ionterest of fit (early rise)
    #     '''
    #     if self._time_zone_selected == True:
    #         len_roi = len(self.fittable[self.fittable['tfromexplo_zc']>=0])
    #         return len_roi
    #     else:
    #         print('You need to select a region of interest got the fit')

    
    def get_chi_2(self, param):
        '''
        This function returns the Chi squared of the function fitted to the data. 

        parameters
        ----------
        self

        returns
        -------

        ''' 
        
        #A,n = param 
        model_    = bounded_simple_EXP_law(self.t_fit,*param)
        pull_     = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq    = np.sum(pull_)

        num_param = 2

        dof    = len(self.f_fit) - num_param
        if dof == 0: #avoid division by zero? 
            dof = 1

        return  chi_sq / dof
    
    def _get_chi2minuit(self, A, n,t_peak, peak_val):
        '''
        A
        n

        '''

        param = [A, n,t_peak, peak_val]

        return self.get_chi_2(param)



    # def get_explosion_t(self):
    #     '''
    #     '''
    #     self.t_exp  = self.minuit_output.params[1].value
    #     self.dt_exp = self.minuit_output.params[1].error
    #     return self.t_exp, self.dt_exp
        

    # def get_fit_index(self):
    #     '''
    #     '''
    #     self.n      = self.minuit_output.params[2].value
    #     self.dt_n   = self.minuit_output.params[2].error
    #     return self.n, self.dt_n




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
            
        if boundaries is not None:
            for k,v in boundaries.items():
                self.minuit_kwargs["limit_"+k] = v
                
        if fixed is not None:
            for k in fixed:
                self.minuit_kwargs["fix_"+k] = True
        
        self.minuit = Minuit(self._get_chi2minuit, errordef=errordef, 
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


        x_         = np.linspace(0, max(self.t_fit), 1000)
        
        A_         = self.minuit_output.params[0].value
        n_         = self.minuit_output.params[1].value
        t_peak_    = self.minuit_output.params[2].value
        peak_val_  = self.minuit_output.params[3].value
        
        dA_         = self.minuit_output.params[0].error
        dn_         = self.minuit_output.params[1].error

        
        
        # A1_string       = f't_exp = {A_:.3f}±{dA_:.3f} d '
        # n1_string       = f'n = {alpha1_:.3f}±{dalpha1_:.3f} '
        # t_break_string  = f'n = {t_break_:.3f}±{dt_break_:.3f} '
        # A2_string       = f'n = {alpha2_:.3f}±{dalpha2_:.3f} '
        # n2_string       = f'n = {n_:.3f}±{dn_:.3f} '

    

        chi2_fit       = self.get_redchi_2_fit([A_, n_ , t_peak_, peak_val_])
        chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        




        if add_plot == None:
            plt.figure()
            plt.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=3.5 , color = self._filter_name.get(self.filter) )
            plt.plot(x_ , bounded_simple_EXP_law(x_, A_,  n_ ,t_peak_, peak_val_), ls = '--' , color = _fitcolor )

            
            plt.plot(x_[1] , bounded_simple_EXP_law(x_,A_, n_,t_peak_, peak_val_ )[1], color = 'white',label= chi_fit_string )


            
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)
            # plt.axvline(t_break_, lw = 0.5, color = 'red')

            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Apparent Mag', size = 15)

            # plt.gca().invert_yaxis()
            plt.legend()
        
        else: 
            add_plot.errorbar(self.t_fit, self.f_fit, self.e_f_fit, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
            add_plot.plot(x_ ,bounded_simple_EXP_law(x_, A_, n_,t_peak_, peak_val_), ls = '--'  , color = _fitcolor )

            
            add_plot.plot(x_[1] , bounded_simple_EXP_law(x_,A_, n_,t_peak_, peak_val_)[1], color = 'white',label= chi_fit_string )



            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from EED [days]', size = 15)
            plt.ylabel('Apparent Mag', size = 15)

            # plt.gca().invert_yaxis()

            plt.legend(fontsize = 10)
        
        




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