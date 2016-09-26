# Class Eisenman conceptual model
# BY: Mark Dekker

#Preambule
import sys,os
import time as T

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use('ggplot')
import csv
import matplotlib.cm as cm

class Eisenman(object):

    def __init__(self,parameters,method='Upwind',initialvalue=None):
        '''
        Initialization method. 

        INPUTS:
            -parameters: dictionary with model parameters: Requires. ['xmax','dx','tmax','dt','u0','k','E']
                xmax    : length of x-dimension, float
                dx      : spatial resolution, float
                tmax    : length of integration, float
                dt      : timestep, float
                u0      : windspeed, float
                k       : decay constant. float or array of the correct size (nx = xmax / dx)
                E       : sources. array of the correct size (nx = xmax / dx)
            -method: integration rule. Supports 'Upwind','Leapfrog', 'LaxWendroff'. by default 'Upwind'
            -initialvalue: initial condition. array of the correct size (nx = xmax / dx)
        OUTPUT: no output variables

        Does the following:
            -Assign parameters from the dictionary parameters to self.P
            -Set self.initialvalue to its value, or to 0, if not specified
            -Set integration method (by default Upwind)
            -Check if all essential parameters have been passed from parameters
            -call the method self.initialize
                -computes derived parameters
                -initializes arrays for time, space and results
                -assigns initial condition to first entry of results
                -checks whether initialvalue is either a single float or an array of size nx. If not, breaks.
        '''
        self.P = dict.fromkeys(['tmax','Dt','T0','S0','t_r','alpha_t','alpha_s','H','t_d','q/v','I','k0','M','alpha','alpha_p','alpha_e','Fs','V','StefBoltz','Emissivity','G'])
        self.P.update(parameters)
  
        # Set initial value to zero, if not specified
        self.initialvalue = 0
        if initialvalue is not None:
            self.initialvalue = initialvalue

        # Set integration rule to Upwind + Euler, if not otherwise specified
        self.method = method

        if None in self.P.values():
            print('WARNING: SOME NECESSARY PARAMETERS HAVE NOT BEEN SPECIFIED. FIX BEFORE RUN')
            print('These parameters are undefined:')
            for key in self.P:
                if self.P[key] == None:
                    print(key)

        self.initialize()
        print('Model initialized')

    def initialize(self):
        '''
        Initialization helper routine

        INPUTS: none (except for self)
        OUTPUTS: none

        computes derived parameters and creates arrays for time, space, and results
        assigns initial condition to first entry of results array
            -computes derived parameters
            -initializes arrays for time, space and results
            -assigns initial condition to first entry of results
            -checks whether initialvalue is either a single float or an array of size nx. If not, breaks.
        '''

        # initialize arrays
        self.time = np.arange(0,self.P['tmax'],self.P['Dt'])
        
        # compute derived parameters
        
        # initialize results array

        # set initial condition, if correct dimension
        #try:
            #self.results[0,:] = self.initialvalue
        #except ValueError:
            #sys.exit('Initial value has wrong dimension: %s instead of %s' % (str(self.initialvalue.shape),str(self.results[0,:])))

        print('Done with initialize')


    def updateParameters(self,newParams):
        '''
        INPUTS: newParams, dictionary

        Updates the model by loading new values into the parameter dictionary
        NOTE: this does not check if the parameters in newParams are the ones used by the model and/or in the correct data format
        '''
        self.P.update(newParams)
        self.initialize() # recompute derived parameters

    def reportParameters(self):
        print('Reporting Model parameters')
        paramList = list(self.P.keys())
        paramList.sort()
        for key in paramList:
            print('{:.12s} , {:}'.format(key,self.P[key]))

    def setInitial(self,initialvalue):
        '''
        Reinitialize the model with a new initial condition
        then calls self.initialize
        '''
        self.initialvalue = initialvalue
        self.initialize()


    def reportInitialState(self):   
        print('The initial state of the model is: later Ill do this')


    def set_method(self,method):
        '''
        set integration rule.
        available:
            -Upwind : Upwind in space, Euler forward in time
            -Leapfrog : Leapfrog in space, Euler forward in time
            -LaxWendroff: Lax-Wendroff
        '''
        possibleMethods = ['Upwind','Leapfrog','LaxWendroff']
        if method not in possibleMethods:
            print('Integration method incorrectly specified. Method "%s" does not exist' % method)
            print('Possible methods are %s' % (str(possibleMethods)))
            sys.exit('Exiting ...')
        self.method = method
           
    
    def integrateModel(self):
        
        tmax=self.P['tmax']     
        tlen=np.int(self.P['tmax']/self.P['Dt'])   
        
        T_eq_ocean=np.zeros(tlen)
        T_pole_ocean=np.zeros(tlen)
        T_eq_atm=np.zeros(tlen)
        T_pole_atm=np.zeros(tlen)
        S_eq_ocean=np.zeros(tlen)
        S_pole_ocean=np.zeros(tlen)
        A=np.zeros(tlen)
        P=np.zeros(tlen)
        Co=np.zeros(tlen)
        Ca=np.zeros(tlen)
        c_t=np.zeros(tlen)
        CO3=np.zeros(tlen)
        
        
        T_eq_ocean[0]=self.P['T_eq_ocean_0']
        T_pole_ocean[0]=self.P['T_pole_ocean_0']
        T_eq_atm[0]=self.P['T_eq_atm_0']
        T_pole_atm[0]=self.P['T_pole_atm_0']
        S_eq_ocean[0]=self.P['S_eq_ocean_0']
        S_pole_ocean[0]=self.P['S_pole_ocean_0']
        A[0]=self.P['A_0']
        P[0]=self.P['P_0']
        Ca[0]=self.P['Ca_0']
        Co[0]=self.P['Co_0']
        CO3[0]=self.P['CO3_0']
        

        def Qvector(self, T_nu,S_nu):
            gradientrho=1.-self.P['alpha_t']*(T_nu-self.P['T0'])+self.P['alpha_s']*(S_nu-self.P['S0'])
            return 1/self.P['t_d']+self.P['q/v']*gradientrho**2
            
            
        def kfunction(self, T_nu):
            return self.P['k0']*np.exp(self.P['alpha']*T_nu)
        
        def kgrowth(self,time):
            return self.P['k0']*(1.+self.P['beta']*np.cos(time*2.*np.pi/self.P['Period']+self.P['phase']))
        
        '''
        Integrate the model based on
            -specified parameters, including source array
            -specified initial condition
            -specified integration method
        Results are saved in self.results
        '''

        # test if E has correct dimension
        #if not self.P['E'].shape == self.results[0,:].shape:
        #    sys.exit('Source array does not have correct shape: %s . Needs same shape as x-dimemsion: %s')
        

        # =====================================
        # METHOD SELECTION LOOP
        # =====================================
        start_time = T.time()
        if self.method == 'Upwind':
            # MODEL RUN
            print('================================================')
            print('Starting model run with method %s' % self.method)
            
            #Ocean box
            self.time = np.linspace(0,self.P['tmax'],tlen)
                		                
            for t in range(1,tlen):
                Qvec_eq_ocean=Qvector(self,T_eq_ocean[t-1],S_eq_ocean[t-1])
                Qvec_pole_ocean=Qvector(self,T_pole_ocean[t-1],S_pole_ocean[t-1])
                kvec=kfunction(self,T_eq_ocean[t-1])
                kvec2=kgrowth(self,t)
                
                T_eq_ocean[t]=T_eq_ocean[t-1]+(-1/self.P['t_r']*(T_eq_ocean[t-1]-T_eq_atm[t-1])-1./2.*Qvec_eq_ocean*(T_eq_ocean[t-1]-T_pole_ocean[t-1]))*self.P['Dt']
                T_pole_ocean[t]=T_pole_ocean[t-1]+(-1/self.P['t_r']*(T_pole_ocean[t-1]-T_pole_atm[t-1])-1./2.*Qvec_pole_ocean*(T_pole_ocean[t-1]-T_eq_ocean[t-1]))*self.P['Dt']
                
                S_eq_ocean[t]=S_eq_ocean[t-1]+(self.P['Fs']*self.P['S0']/2./self.P['H']-1./2.*Qvec_eq_ocean*(S_eq_ocean[t-1]-S_pole_ocean[t-1]))*self.P['Dt']
                S_pole_ocean[t]=S_pole_ocean[t-1]+(0-self.P['Fs']*self.P['S0']/2./self.P['H']-1./2.*Qvec_pole_ocean*(S_pole_ocean[t-1]-S_eq_ocean[t-1]))*self.P['Dt']
                
    
                #Atmospheric box                
                T_eq_atm[t]=T_eq_atm[t-1]+(T_eq_atm[t-1]*(1-self.P['alpha_e'])+self.P['G']-self.P['StefBoltz']*self.P['Emissivity']*T_eq_atm[t-1]**4)*self.P['Dt']+A[t-1]*np.log(Ca[t-1]/Co[t-1])*self.P['Dt']
                T_pole_atm[t]=T_pole_atm[t-1]+(T_pole_atm[t-1]*(1-self.P['alpha_e'])+self.P['G']-self.P['StefBoltz']*self.P['Emissivity']*T_pole_atm[t-1]**4)*self.P['Dt']+A[t-1]*np.log(Ca[t-1]/Co[t-1])*self.P['Dt']
                

                A[t]=A[t-1]+(self.P['I']-kvec*A[t-1]*P[t-1])*self.P['Dt']
                P[t]=P[t-1]+(kvec*A[t]*P[t-1]-self.P['SedRate']*P[t-1])*self.P['Dt']                
                CO3[t]=A[t]-np.exp(P[t])
                Ca[t]=(A[t]-2.*CO3[t-1])**2./(CO3[t-1])*kfunction(self,(T_eq_atm[t-1]+T_pole_atm[t-1])/2.)
                Co[t]=A[t]-CO3[t-1]
                
                #A[t]=A[t-1]+(self.P['I']-kvec*CO3[t-1]*P[t-1])*self.P['Dt']
                #P[t]=P[t-1]+(kvec*CO3[t-1]*P[t-1]-self.P['SedRate']*P[t-1])*self.P['Dt']
               
                #print((T_eq_atm[t-1]+T_pole_atm[t-1])/2.)
                #print(kfunction(self,(T_eq_atm[t-1]+T_pole_atm[t-1])/2.))
                #Ca[t]=(A[t]-2.*CO3[t-1])**2./(CO3[t-1])*kfunction(self,(T_eq_atm[t-1]+T_pole_atm[t-1])/2.)
                #Co[t]=A[t]-CO3[t-1]
                #c_t[t]=self.P['V']/2.*A[t]
                #CO3[t]=CO3[t-1]*1.0001
               from sympy.solvers import solve
               from sympy import Symbol
               x = Symbol('x')
               CO3[t]=solve(self.P['M']*(A[t]-2*x)**2./(x*kfunction(self,(T_eq_atm[t-1]+T_pole_atm[t-1])/2.))+self.P['V']*(A[t]-x)==c_t[t], x)[0]
                
                    
                if np.mod(t,np.int(tlen/10.))==0:
                    print('Progress is at ', np.float(t)/(tlen)*100., 'percent')
    
    
            self.T_eq_ocean = T_eq_ocean
            self.T_pole_ocean = T_pole_ocean
            self.S_eq_ocean = S_eq_ocean
            self.S_pole_ocean = S_pole_ocean
            self.T_eq_atm = T_eq_atm
            self.T_pole_atm = T_pole_atm
            self.A = A
            self.P = P
            self.Ca = Ca
            self.Co = Co
            self.c_t = c_t
            self.CO3 = CO3
                     
            
            print('Total time required: %.2f seconds' % (T.time() - start_time))
            print('Model run finished')
            print('================================================')
            # END OF MODEL RUN

        else:
            sys.exit('Integration method incorrectly specified. Method "%s" does not exist' % self.method)

        # =====================================


    def saveModel(self,savename):
        '''
        INPUTS: savename: data is saved in the file 'savename.npy'

        saves model
            -parameters self.P
            -method, initialvalue
            -x,time, results
        '''
        container = self.P.copy()

        container['method'] = self.method
        container['initialvalue'] = self.initialvalue
        container['Mtot'] = self.Mtot
        container['M2'] = self.M2
        container['x'] = self.x
        container['time'] = self.time
        container['results'] = self.results

        np.save(savename + '.npy', container)
        print('Saving done')