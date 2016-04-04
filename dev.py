import data_handling
import nComp
import pure
import plot
import logging
import os
import numpy
import scipy
import scipy.interpolate
import Van_der_Waals
VdW = Van_der_Waals.VdW()
import argparse
import matplotlib.pyplot as plot

class Args:
    def __init__(self):
        self.compounds = ['acetone']
        # Options: 'acetone' 'benzene' 'hexane' 'cyclohexane' 'water' etc.
        self.phases = ['x', 'y']
        self.eos = 'DWPM'
        self.model = 'Adachi-Lu' # Options: 'Adachi-Lu'
                                 #          'Soave'
        self.r = None
        self.s = None
        self.temperature = None
        self.pressure = None
        self.z = None

        self.lle_only = False
        self.vle_only = False

        self.plot_iso = False
        self.plot_gibbs = False
        self.plot_pure = False
        # Saves
        self.optimise = False
        self.save = False
        self.save_pure = False
        self.force_pure_update = False

if __name__ == '__main__':
    args = Args()
    data = data_handling.ImportData()
    data.run_options(args)
    # Load all pure dictionaries data.c[i]
    data.load_pure_data()

    s, p = pure.pure_sim(data, i=0)
    
    s['b'] = p['b_c'] # (b is not temperature dependant)

    ## START HERE:
    # V_root
    V_l_old = []
    V_v_old = []

    V_l_new = []
    V_v_new = []
    
    errorVold = []
    errorLold = []
    
    errorVnew = []
    errorLnew = []
    

    for T, P in zip(p['T'], p['P']):
        s['T'] = T
        s['P'] = P
        s = VdW.a_T(s,p)
        # Old V_root functions
        s = VdW.V_root(s, p)
        V_l_old.append(s['V_l'])
        V_v_old.append(s['V_v'])
        #print s['V_l']
        #error Vapour
        errorV1= s['P'] - (p['R']*s['T'])/(s['V_v']-s['b']) + s['a']/s['V_v']**2
        #print errorV
        #error Liquid
        errorL1= s['P'] - (p['R']*s['T'])/(s['V_l']-s['b']) + s['a']/s['V_l']**2
        
        errorVold.append(errorV1)
        errorLold.append(errorL1)
        
        #######################################################################
        # Add your new function with same outputs as VdW.V_roots here
        s = VdW.V_root_new(s, p)
        #print s['V_l']
        #######################################################################
        V_l_new.append(s['V_l'])
        V_v_new.append(s['V_v'])
        #error Vapour
        errorV2= s['P'] - (p['R']*s['T'])/(s['V_v']-s['b']) + s['a']/s['V_v']**2
        #print errorV
        #error Liquid
        errorL2= s['P'] - (p['R']*s['T'])/(s['V_l']-s['b']) + s['a']/s['V_l']**2
        
        errorVnew.append(errorV2)
        errorLnew.append(errorL2)
       
     
    
   
    
    # Plot Volume root results on data
    plot.figure(1)
    plot.plot(V_l_old, p['P'], 'b*', label='Old liquid root')
    plot.plot(V_v_old, p['P'],  'r*', label='Old vapour root')
    plot.plot(V_l_new, p['P'], 'bx', label='New liquid root')
    plot.plot(V_v_new, p['P'],  'rx', label='New vapour root')
    plot.ylabel("Pressure / Pa")
    plot.title("Volume root differences on data P, "
               "T for {}".format(p['name'][0]))
    plot.legend()
    plot.show()
#
    plot.figure(2)
    plot.plot(range(len(V_l_old)),
              numpy.array(V_l_old) - numpy.array(V_l_new),
              'b*', label='$\Delta V_l$')
    plot.plot(range(len(V_v_old)),
              numpy.array(V_v_old) - numpy.array(V_v_new),
              'r*', label='$\Delta V_l$')
    plot.ylabel("$\Delta V^{old} - V^{new}$")
    plot.title("Volume root differences on data P, "
               "T for {}".format(p['name'][0]))
    plot.legend()
    plot.figure(3)
    plot.plot(p['P'],errorVold,'r',p['P'],errorLold,'b',p['P'],errorVnew,'r*',p['P'],errorLnew,'b*')
    plot.show()
#
#    ###########################################################################
#    # Psat_V_roots
#    # Phase envelopes below the crtical point
#    P_range = numpy.linspace(p['P'][0],                      # Lower limit
#                             0.99 * p['P'][len(p['P']) - 1], # Upper limit
#                             200)                           # Sampling points
#    T_range = numpy.linspace(p['T'][0],                      # Lower limit
#                             0.99 * p['T'][len(p['T']) - 1], # Upper limit
#                             200)                           # Sampling points
#    V_l_sat_old = []
#    V_v_sat_old = []
#    P_sat_old = []
#    V_l_sat_new = []
#    V_v_sat_new = []
#    P_sat_new = []

#    for T, P in zip(T_range, P_range):
#        #print T
#        #print P
#        s['T'] = T
#        s['P'] = P * 1e-2
#        # Old Volume roots
#
#
#        # Old P_sat functions
#        s = VdW.Psat_V_roots(s, p)
#        V_l_sat_old.append(s['V_l'])
#        V_v_sat_old.append(s['V_v'])
#        P_sat_old.append(s['P_sat'])
#        # New P_sat functions
#        s['T'] = T
#        s['P'] = P * 1e-2 # Reset starting point to range
#        #######################################################################
#        # Add your new function with same outputs as VdW.Psat_V_roots here
#        s = VdW.Psat_V_roots(s, p)
#        #######################################################################
#        V_l_sat_new.append(s['V_l'])
#        V_v_sat_new.append(s['V_v'])
#        P_sat_new.append(s['P_sat'])
#
#    #Plot phase envelopes
#    plot.figure()
#    plot.plot(V_l_sat_old, P_sat_old, 'b--', label='Old liquid phase root')
#    plot.plot(V_v_sat_old, P_sat_old,  'r--', label='Old vapour phase root')
#    plot.plot(V_l_sat_new, P_sat_new, 'b', label='New liquid phase root')
#    plot.plot(V_v_sat_new, P_sat_new,  'r', label='New vapour phase root')
#    plot.ylabel("Pressure$^{sat}$ / Pa" ) #,rotation=-1)
#    plot.xlabel("Volume / L")
#    plot.title("Phase envelope for {}".format(p['name'][0]))
#    plot.legend()
#    plot.yscale('log')
#    plot.xscale('log')
#    plot.figure(3)
#    plot.plot(V_l_sat_new,P_sat_new)
#    plot.show()
#    #print V_l_sat_new
