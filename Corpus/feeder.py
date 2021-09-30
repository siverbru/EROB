# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:39:35 2014

@author: Ruben
"""

import residential
import _pickle as cPickle
import numpy as np
import os


class IDEAS_Feeder(object):
    '''
    The Community class defines a set of hosueholds.
    '''
    
    def __init__(self, name, nBui, path, OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt):
        '''
        Create the community based on number of households and simulate for
        output towards IDEAS model.
        '''
        self.name = name
        self.nBui = nBui
        self.nbedr = nbedr
        self.nday = nday
        # we create, simulate and pickle all 'nBui' buildings
        self.simulate(path, OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt)
        # then we loop through all variables and output as single file
        # for reading in IDEAS.mo
        os.chdir(path)
        variables =['Window','DHW', 'P', 'sh_day','sh_night','sh_bath', 'QRad','QCon']
        for var in variables:
            self.output(var)
        # and conclude
        print ('\n')
        print (' - Feeder %s outputted.' % str(self.name))

    def simulate(self, path, OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt):
        '''
        Simulate all households in the depicted feeder
        '''
        #######################################################################
        # we loop through all households for creation, simulation and pickling.
        # whereas the output is done later-on.
        cwd = os.getcwd()
        for i in range(self.nBui):
            hou = residential.Household(str(self.name)+'_'+str(i), OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, i, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt)
            hou.simulate()
            hou.roundUp()
            os.chdir(path)
            hou.pickle()
            os.chdir(cwd)

    def output(self, variable):
        '''
        Output the variable for the dwellings in the feeder as a *.txt readable
        for Modelica.
        '''
        #######################################################################
        # we loop through all households for loading the depicted variable data
        # which is stored in the object pickle.
        nday = self.nday
        dat = np.zeros(0)
        for i in range(self.nBui):
            hou = cPickle.load(open(str(self.name)+'_'+str(i)+'.p','rb'))
            var = eval('hou.'+variable)
            if len(dat) != 0:
                dat = np.vstack((dat,var))
            else:
                dat = var
        #######################################################################
        # and output the array to txt
        sec = nday*24*3600
        step = nday*144
        tim = np.linspace(0,sec,step+1)
        if variable in ['P','DHW','sh_day','sh_bath','sh_night']:
            dat = np.vstack((tim,dat))
            hea ='#1 \n double data('+str(int(len(tim)))+','+str(self.nBui+1)+')'
        else:
            dat = np.vstack((tim,dat))
            hea ='#1 \n double data('+str(int(len(tim)))+','+str((self.nBui*11)+1)+')'

        directory = os.getcwd()
        os.chdir(directory + '/Output')
        np.savetxt(fname=variable+'.txt', X=dat.T, header=hea, comments='')
        os.chdir(directory)

