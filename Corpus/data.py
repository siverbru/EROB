# -*- coding: utf-8 -*-
"""
Created on Thu Feb 06 09:48:07 2014

@author: Ruben Baetens
"""

import os
import numpy as np
import stats

def get_clusters(employment, **kwargs):
    '''
    Find the clusters for each day of the week for a household member
    of the given employment type based on Belgian time use data of 2013.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    PATH = '../Data/Occupancy/Crosstables/'
    os.chdir(PATH)
    
    #create an empty dictionary
    keys = ['mon', 'tue', 'wed', 'thu', 'fri', 'sat', 'sun']
    cluDict = dict()

    # determine the cluster for each of the daytypes for the given employment
    for key in keys:
        if employment in ['FTE','Student']:
            if key in ['tue','wed', 'thu']:
                cluster = cluDict['mon']
            else:
                order = ['FTE','PTE','Unemployed','Retired','School','Student']
                emp_i = order.index(employment)
                data = np.loadtxt(key+'.txt', float).T[emp_i]
                rnd = np.random.random()
                cluster = stats.get_probability(rnd, data[0:], p_type='prob')
        elif employment in ['School']:
            if key in ['tue', 'thu', 'fri']:
                cluster = cluDict['mon']
            else:
                order = ['FTE','PTE','Unemployed','Retired','School','Student']
                emp_i = order.index(employment)
                data = np.loadtxt(key+'.txt', float).T[emp_i]
                rnd = np.random.random()
                cluster = stats.get_probability(rnd, data[0:], p_type='prob')
        else:
            order = ['FTE','PTE','Unemployed','Retired','School','Student']
            emp_i = order.index(employment)
            data = np.loadtxt(key+'.txt', float).T[emp_i]
            rnd = np.random.random()
            cluster = stats.get_probability(rnd, data[0:], p_type='prob')
        cluDict.update({key:cluster})

    # and return the final clusters
    os.chdir(cdir)
    return cluDict

def get_occDict(cluster, **kwargs):
    '''
    Create the dictionary with occupancy data based on the 
    cluster assignment. These occupancy data are based
    on Belgian time use statistics of 2013.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    DATA_PATH = '../Data/Occupancy'
    PATH = DATA_PATH + '/Pattern' + str(cluster)
    os.chdir(PATH)
    
    occDict = dict()

    # first we load the occupancy start states 'oss' from StartStates.txt
    ss = dict()
    data = np.loadtxt('StartStates.txt', float)
    for i in range(len(data)):
        ss.update({str(i+1):data[i]})
    occDict.update({'ss':ss})

    # Second we load the occupancy transitions state probabilities 'osn'
    # from TransitionProbability.txt
    data = np.loadtxt('TransitionProbability.txt', float)
    for i in range(3):
        os_i = dict()
        for j in range(48):
            os_i.update({str(j+1):data[i*48+j]})
        occDict.update({'os_'+str(i+1):os_i})

    # Third we load the Markov time density 'ol' from DurationProbability.txt
    data = np.loadtxt('DurationProbability.txt', float)
    for i in range(3):
        ol_i = dict()
        for j in range(48):
            ol_i.update({str(j+1):data[i*48+j]})
        occDict.update({'ol_'+str(i+1):ol_i})

    # and return the final occDict
    os.chdir(cdir)
    return occDict

def get_taskDict(weekday, HHsize, **kwargs):
    '''
    Create the dictionary with household task data based on Belgian time use statistics from 2013.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    DATA_PATH = '../Data/HouseholdTasks'
    os.chdir(DATA_PATH)
    
    # create an empty dictionary
    taskDict = dict()
    
    # first we define the dictionary used as legend for the load file
    task = {0:'cook', 1:'dishes', 2:'vacuum', 3:'iron', 4:'wash', 5:'dry'}

    
    # Second we load the start probability from the corresponding file
    # based on the day of the week and the household size
    if weekday == 0:
        WKD = 'Weekend'
    else:
        WKD = 'Week'
    FILNAM = str(HHsize) +'Pers'+ WKD +'.txt'
    data = np.loadtxt(FILNAM, float)
    for i in range(6):
        taskDict.update({task[i]:data.T[i]})

    # and return the final actDict
    taskDict.update({'period':600, 'steps':144})
    os.chdir(cdir)
    return taskDict

def get_taskdurDict(HHsize, **kwargs):
    '''
    Create the dictionary with probability duration of the household tasks.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    DATA_PATH = '../Data/HouseholdTasks'
    os.chdir(DATA_PATH)
    
    taskdurDict = dict()
    
    if HHsize == 1:
        SIZE = 1
    else:
        SIZE = 2
    # We load the Markov time density 'ol' from DurationProbability.txt
    data = np.loadtxt('DurationProbability'+str(SIZE)+'.txt', float)
    for i in range(6):
        dur_i = dict()
        for j in range(48):
            dur_i.update({str(j):data[i*48+j]})
        taskdurDict.update({'dur_'+str(i):dur_i})

    # and return the final taskdurDict
    os.chdir(cdir)
    return taskdurDict

def get_actDict(weekday, employment, **kwargs):
    '''
    Create the dictionary with startprobabilities of the individual activities.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    DATA_PATH = '../Data/IndividualActivities'
    os.chdir(DATA_PATH)
    
    # create an empty dictionary
    actDict = dict()

    # first we define the dictionary used as legend for the load file
    act = {0:'pc', 1:'tv', 2:'audio', 3:'adm'}

    # Second we load the start probability from the corresponding file
    # based on the day of the week and the employment type
    FILNAM = employment + str(weekday) +'.txt'
    data = np.loadtxt(FILNAM, float)
    for i in range(4):
        actDict.update({act[i]:data.T[i]})
        
    # and return the final actDict
    actDict.update({'period':600, 'steps':144})
    os.chdir(cdir)
    return actDict

def get_actdurDict(employment, **kwargs):
    '''
    Create the dictionary with probability duration of the household tasks.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    DATA_PATH = '../Data/IndividualActivities'
    os.chdir(DATA_PATH)
    
    actdurDict = dict()
    
    # We load the Markov time density 'ol' from DurationProbability.txt
    if employment in ['Student','School']:
        data = np.loadtxt('DurationProbabilitySTUDENT.txt', float)
    else:
        data = np.loadtxt('DurationProbabilityOTHER.txt', float)
    for i in range(6):
        dur_i = dict()
        for j in range(48):
            dur_i.update({str(j):data[i*48+j]})
        actdurDict.update({'dur_'+str(i):dur_i})

    # and return the final taskdurDict
    os.chdir(cdir)
    return actdurDict

def get_bathDict(cluster, **kwargs):
    '''
    Create the dictionary with startprobabilities of presence in bathroom.
    '''
    #first go the the correct location
    cdir = os.getcwd()
    DATA_PATH = '../Data/IndividualActivities'
    os.chdir(DATA_PATH)
    
    # create an empty dictionary
    bathDict = dict()

    # Second we load the start probability from the corresponding file
    # based on the day of the week and the employment type
    FILNAM = 'BathCluster'+str(cluster) +'.txt'
    data = np.loadtxt(FILNAM, float)
    bathDict.update({'bath':data.T[0]})
        
    # and return the final bathDict
    bathDict.update({'period':600, 'steps':144})
    os.chdir(cdir)
    return bathDict
