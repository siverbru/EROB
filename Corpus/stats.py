# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 11:57:02 2013

@author: Ruben Baetens
"""

import numpy as np
import random

import data


def get_probability(rnd, prob, p_type='cum'):
    '''
    Find the x-value in a given comulative probability 'prob_cum' based on a
    given random y-value 'rnd'.
    '''
    if p_type != 'cum':
        prob = np.cumsum(prob)
        prob /= max(prob)
    idx = 1
    while rnd >= prob[idx-1]:
        idx += 1
    return idx

def sum_dict(dict_a, dict_b):
    '''
    Sum the values stored under the same keys in python dictionarys.
    '''
    # loop through the keys and sum both values given
    if len(dict_a.keys()) == 0:
        sum_dict = dict_b
    elif len(dict_b.keys()) == 0:
        sum_dict = dict_a
    else:
        sum_dict = dict()
        for key in dict_a.keys():
            if dict_a[key] is None:
                sum_dict.update({key:None})
            elif key != 'time':
                sum_dict.update({key:dict_a[key]+dict_b[key]})
            else:
                sum_dict.update({key:dict_a[key]})
    # and return
    return sum_dict

class MCSA(object):
    '''
    The MCSA class defines a Monte Carlo Survival Analysis.
    '''
    def __init__(self, cluster, **kwargs):
        # load the dataset of the cluster into ds
        ds = data.get_occDict(cluster)
        self.OSS = ds['ss']
        self.OPM = {1:ds['os_1'], 2:ds['os_2'], 3:ds['os_3']}
        self.ODM = {1:ds['ol_1'], 2:ds['ol_2'], 3:ds['ol_3']}

    def startstate(self):
        '''
        Get the startstate for first simulation day at 4:00 AM.
        '''
        # we define the startstate based on the given probability
        probs = [self.OSS['1'], self.OSS['2'], self.OSS['3']]
        state = get_probability(random.random(), probs)

        return int(state)

    def transition(self, state, timebin):
        '''
        Get next occupancy state from current state ending at time.
        '''
        # we define the new state based on the given probability
        probs = self.OPM[state][str(timebin)]
        newoc = get_probability(random.random(), probs)

        return int(newoc)

    def duration(self, state, timebin):
        '''
        Get the duration of current state started at time.
        '''
        # we define the new duration based on the given probability
        probs = self.ODM[state][str(timebin)]
        durat = get_probability(random.random(), probs)

        return durat

class MCSA_task(object):
    '''
    The MCSA class defines a Monte Carlo Survival Analysis.
    '''
    def __init__(self, HHsize, **kwargs):
        # load the dataset of the cluster into ds
        ds = data.get_taskdurDict(HHsize)
        self.ADM = {0:ds['dur_0'], 1:ds['dur_1'], 2:ds['dur_2'], 3:ds['dur_3'], 4:ds['dur_4'], 5:ds['dur_5']}

    def duration(self, act, timebin):
        '''
        Get the duration of current activity started at time.
        '''
        # we define the new duration based on the given probability
        probs = self.ADM[act][str(timebin)]
        durat = get_probability(random.random(), probs)

        return durat
    
class MCSA_act(object):
    '''
    The MCSA class defines a Monte Carlo Survival Analysis.
    '''
    def __init__(self, employment, **kwargs):
        # load the dataset of the cluster into ds
        # 0-3 as in act lists, 4 = un(dress), 5 = bath
        ds = data.get_actdurDict(employment)
        self.ADM = {0:ds['dur_0'], 1:ds['dur_1'], 2:ds['dur_2'], 3:ds['dur_3'], 4:ds['dur_4'], 5:ds['dur_5']}

    def duration(self, act, timebin):
        '''
        Get the duration of current activity started at time.
        '''
        # we define the new duration based on the given probability
        probs = self.ADM[act][str(timebin)]
        durat = get_probability(random.random(), probs)

        return durat
    
class DTMC_bath(object):
    '''
    The DTMC class defines a Discrete-Time Markov Chain
    '''
    # All object parameters are given in kwargs
    def __init__(self, clusterDict, **kwargs):
        # load the dataset of the cluster into ds
        self.ds = dict()
        self.ds.update({0:data.get_bathDict(clusterDict['sun'])})
        self.ds.update({1:data.get_bathDict(clusterDict['mon'])})
        self.ds.update({2:data.get_bathDict(clusterDict['tue'])})
        self.ds.update({3:data.get_bathDict(clusterDict['wed'])})
        self.ds.update({4:data.get_bathDict(clusterDict['thu'])})
        self.ds.update({5:data.get_bathDict(clusterDict['fri'])})
        self.ds.update({6:data.get_bathDict(clusterDict['sat'])})
    def get_var(self, dow, act, step):
        # get the probability of the given activity for daytype dow at step
        return self.ds[dow][act][step]

class DTMC_task(object):
    '''
    The DTMC class defines a Discrete-Time Markov Chain
    '''
    # All object parameters are given in kwargs
    def __init__(self, HHsize, **kwargs):
        # load the dataset of the cluster into ds
        if HHsize > 5:
            HHsize = 5
        self.ds = dict()
        for i in range(5):
            self.ds.update({i+1:data.get_taskDict(0, HHsize)})
        self.ds.update({6:data.get_taskDict(1, HHsize)})
        self.ds.update({0:data.get_taskDict(1, HHsize)})
    def get_var(self, dow, act, step):
        # get the probability of the given activity for daytype dow at step
        return self.ds[dow][act][step]
    
class DTMC_act(object):
    '''
    The DTMC class defines a Discrete-Time Markov Chain
    '''
    # All object parameters are given in kwargs
    def __init__(self, employment, **kwargs):
        # load the dataset of the cluster into ds
        # 0 = weekday, 1= weekend
        self.ds = dict()
        for i in range(5):
            self.ds.update({i+1:data.get_actDict(0, employment)})
        self.ds.update({0:data.get_actDict(1, employment)})
        self.ds.update({6:data.get_actDict(1, employment)})
    def get_var(self, dow, act, step):
        # get the probability of the given activity for daytype dow at step
        return self.ds[dow][act][step]

