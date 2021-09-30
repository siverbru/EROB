# -*- coding: utf-8 -*-
"""
Created on Mon October 07 16:16:10 2013

@author of original model: Ruben Baetens
@author of changes: Silke Verbruggen
"""

import sys
from sys import exit
import random
import numpy as np
import time
import datetime
import calendar
import os
import csv
import _pickle as cPickle
import itertools
from scipy.stats import weibull_min, genextreme, loglaplace

import stats
import data
import windowhabits

sys.path.append("..")
from Data.Households import households
from Data.Appliances import set_appliances

class Household(object):
    '''
    The Household class allows simulation of the household 
    for building energy simulations.

    Main functions are:
        - __init__()
        - self.parameterize(), which determines the properties of the household
        - self.simulate(), main function which includes following functions:
            - Occupancy, determine the occupancy of each household member
            - Activity
            - ....
    '''

    def __init__(self, name, OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, i, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt, **kwargs):
        '''
        Initiation of Household object.
        '''
        # check and define inputs
        try:
            if not isinstance(name, str):
                raise ValueError('Given name %d is not a string' % str(name))
        except:
            raise TypeError('give another name')
        self.creation = time.asctime()
        self.name = name
        self.OccONFILE = OccONFILE
        self.ActONFILE = ActONFILE
        self.shtype = shtype
        self.shrooms = shrooms
        self.members = members
        self.apps = apps
        self.clusters = clusters
        self.nbedr = nbedr
        self.BEDR = BEDR
        self.nday = nday
        self.year = year
        self.nBui = i
        self.npers = npers
        self.habits = habits
        self.HHhabit = HHhabit
        self.SeCo = SeCo
        self.VentS = VentS
        self.DW = DW
        self.YearBuilt = YearBuilt
        self.parameterize(**kwargs)

    def parameterize(self, **kwargs):
        '''
        Determine the household (members, clusters, appliances, tappoints) 
        either based on inputs or based on statistics from the Belgian TUS 2013.
        '''
        def check(self, **kwargs):
            '''
            Check if the input parameters are valid and consequent
            '''
            if self.members:
                if self.npers != 0:
                    if len(self.members) != self.npers:
                        raise Error('Number of persons does not coincide with defined members')
                if self.clusters:
                    if len(self.members) != len(self.clusters):
                        raise Error('Number of members does not coincide with number of occupancy profiles')
                if self.BEDR:
                    if len(self.members) != len(self.BEDR):
                        raise Error('Number of members does not coincide with number of members assigned to bedroom')
            elif self.npers != 0:
                if self.clusters:
                    if self.npers != len(self.clusters):
                        raise Error('Number of persons does not coincide with number of occupancy profiles')
                if self.BEDR:
                    if self.npers != len(self.BEDR):
                        raise Error('Number of persons does not coincide with number of members assigned to bedroom')
            elif self.clusters:
                if self.BEDR:
                    if len(self.clusters) != len(self.BEDR):
                        raise Error('Number of occupancy profiles does not coincide with number of members assigned to bedroom')
            
            if self.nbedr != 0:
                if self.BEDR:
                    if self.nbedr < max(self.BEDR)+1:
                        raise Error('Number of bedrooms assigned to members exceeds the available number of bedrooms')
                        exit()
            return None
      
        def members(self, **kwargs):
            '''
            Define the employment type of all household members based on time
            use survey statistics (TUS 2005) or the given kwargs.
            '''
            members = []
            # First we check if the members are given as input. 
            if not self.members:
                # If no members are given, 
                # check if # persons, # bedrooms or bedroom assignment is given 
                # and assign based on TUS-data.
                if self.npers == 0 and self.nbedr == 0 and len(self.BEDR) == 0 and self.OccONFILE == 0:
                    key = random.randint(1, len(households))
                    members = households[key]
                    npers = len(members)
                elif self.npers != 0 or len(self.BEDR) > 0 or self.OccONFILE == 1:
                    if self.npers != 0:
                        npers = self.npers
                    elif len(self.BEDR) > 0:
                        npers = len(self.BEDR)
                    else:
                        cdir = os.getcwd()
                        os.chdir('../Data/InputFiles')
                        d = os.getcwd()
                        occ_ye = np.loadtxt(d + '/occ'+str(self.nBui)+'.txt', np.int32)
                        occ_year = occ_ye.tolist()
                        if len(occ_year)>20:
                            npers = 1
                        else:
                            npers = len(occ_year)
                        os.chdir(cdir)    
                    while len(members) != npers:
                        key = random.randint(1, len(households))
                        members = households[key]
                elif self.nbedr != 0:
                    data = np.loadtxt('assignPers.txt', float).T[self.nbedr-1]
                    rnd = np.random.random()
                    npers = stats.get_probability(rnd, data[0:], p_type='prob')
                    while len(members) != npers:
                        key = random.randint(1, len(households))
                        members = households[key]
            else:
                members = self.members
                npers = len(members)
                
            self.npers = npers
            print (members)
            return members
        
        def clusters(self, **kwargs):
            '''
            Allocate each household member to the correct cluster based on the
            members occupation in time use survey data.
            '''
            clusters = []
            
            # First we check if the clusters are given as input
            # If no clusters are given, assign based on employment-type
            if not self.clusters:
                members = self.members
                for ind in members:
                    clu_i = data.get_clusters(ind)
                    clusters.append(clu_i)
            else:
                clusters = self.clusters
            print (clusters)
            return clusters
        
        def bedrooms(self, **kwargs):
            '''
            Determine number of bedrooms and assign bedrooms to household members.
            if assignment of bedrooms not given: 
            determine in which bedroom each occupant sleeps, 
            couples sleep together
            and others sleep together when no more rooms are available
            '''
            BEDR = []
            
            # First we check if the Bedroom Assignment is given as input
            if not self.BEDR:
                if self.nbedr == 0:
                    if self.npers > 5:
                        PERS = 4
                    else:
                        PERS = self.npers-1
                    data = np.loadtxt('assignBedr.txt', float)[PERS]
                    rnd = np.random.random()
                    nbedr = stats.get_probability(rnd, data[0:], p_type='prob')
                else:
                    nbedr = self.nbedr
                # determine in which bedroom each occupant sleeps, couples sleep together and others sleep together when no more rooms are available (max 3 pers per room)
                # determine how many couples present
                couple = 0
                if len(self.members)>1:
                    if self.members[0] in ['FTE','PTE','Unemployed', 'Retired'] and self.members[1] in ['FTE','PTE','Unemployed', 'Retired']:
                        couple += 1
                if couple == 0:
                    for i in range(self.npers):
                        if i < nbedr:
                            BEDR.append(i)
                        elif i < nbedr*2-1:
                            BEDR.append(i-nbedr+1)
                        else:
                            BEDR.append(i-(nbedr*2)+1)
                else:
                    for i in range(2):
                        BEDR.append(0)
                    for i in range(self.npers-2):
                        if i < nbedr-1:
                            BEDR.append(i + 1)
                        elif i < (nbedr-1)*2:
                            BEDR.append(i - nbedr + 2)
                        else:
                            BEDR.append(i - 2*nbedr + 3)
            else:
                BEDR = self.BEDR
                if self.nbedr == 0:
                    nbedr = max(self.BEDR)+1
                else:
                    nbedr = self.nbedr
                    
            self.nbedr = nbedr
            return BEDR

        def appliances(self, **kwargs):
            '''
            Define the pressent household appliances based on average national
            statistics independent of household member composition.
            '''
            app_n = []
            
            # First we check if the appliances are given as input. 
            # If no appliances are given, assign based on the rate of ownership.
            if not self.apps:
                for app in set_appliances:
                    if set_appliances[app]['type'] == 'appliance':
                        obj = Equipment(**set_appliances[app])
                        owner = obj.owner >= random.random()
                        app_n.append(app) if owner else None
     
                    if not ('FridgeFreezer' in app_n) and not ('Refrigerator' in app_n):
                        prob=set_appliances['FridgeFreezer']['owner']/(set_appliances['FridgeFreezer']['owner']+set_appliances['Refrigerator']['owner'])
                        app_n.append('FridgeFreezer') if prob >= random.random()  else app_n.append('Refrigerator') 
                    
                    if not ('FridgeFreezer' in app_n) and not ('Freezer' in app_n):
                        app_n.append('Freezer')    
                        
                HHsize = self.npers-1
                if HHsize > 4:
                    HHsize = 4
                data = np.loadtxt('IndividualActivities/ownTV.txt', float)[HHsize]
                rnd = np.random.random()
                Ntv = stats.get_probability(rnd, data[0:], p_type='cum')-1
                for i in range(Ntv):
                    app_n.append('tv')
                data2 = np.loadtxt('IndividualActivities/ownPC.txt', float)[HHsize]
                rnd = np.random.random()
                Npc = stats.get_probability(rnd, data2[0:], p_type='cum')-1
                for i in range(Npc):
                    app_n.append('PC')
                data3 = np.loadtxt('IndividualActivities/ownSmartphone.txt', float)[HHsize]
                rnd = np.random.random()
                Nsp = stats.get_probability(rnd, data3[0:], p_type='cum')-1
                for i in range(Nsp):
                    app_n.append('Smartphone')
                    
            else:
                app_n = self.apps
            #print (app_n)
            return app_n

        def tappings():
            '''
            Define the present household tapping types.
            '''
            tap_n = ['showerFlow', 'bathFlow', 'OtherFlow', 'dishFlow']
            return tap_n
        
        # run all
        self.check = check(self, **kwargs)
        self.members = members(self, **kwargs)
        self.clusters = clusters(self, **kwargs)
        self.BEDR = bedrooms(self, **kwargs)
        self.apps = appliances(self, **kwargs)
        self.taps = tappings()
        self.check = check(self, **kwargs)

        return None

    def simulate(self, **kwargs):
        '''
        The simulate function includes the simulation of the household
        occupancies, internal heat gains, plug loads, lighting loads, 
        DHW-tappings, heating preferences, solar shading use and window use.
        '''
        ndays = self.nday
        year = self.year
        self.__chronology__(year, ndays)
        self.__occupancy__(**kwargs)
        self.__activity__()
        self.__metabolic__()
        self.__plugload__()
        self.__dhwload__()
        self.__shsetting__()
        self.__window__()

    def __chronology__(self, year, ndays):
        '''
        A basic internal calendar is made, storing the days and months of the
        simulated year. (Monday = 0, ..., Sunday = 6)
        '''
        # first we determine the first week of the year
        fdoy = datetime.datetime(year,1,1).weekday()
        fweek = list(range(7))[fdoy:]
        
        # whereafter we fill the complete year
        if ndays:
            nday = ndays
        else:
            nday = 366 if calendar.isleap(year) else 365
        day_of_week = (fweek+53*list(range(7)))[:nday]
        
        # and return the day_of_week for the entire year
        self.dow = day_of_week
        
        return None

    def __occupancy__(self, min_form = True, min_time = False, **kwargs):
        '''
        Simulation of a number of days based on cluster 'BxDict'.
        - Including weekend days,
        - starting from a regular monday at 4:00 AM.
        '''
        def check(occday, min_form = True, min_time = False):
            '''
            We set a check which becomes True if the simulated day behaves
            according to the cluster, as a safety measure for impossible
            solutions.
            '''
            # First we check if the simulated occ-chain has the same shape
            shape = True
            if min_form:
                location = np.zeros(1, dtype=int)
                reduction = occday[0]*np.ones(1, dtype=int)
                for i in range(len(occday)-1):
                    if occday[i+1] != occday[i]:
                        location = np.append(location, i+1)
                        reduction = np.append(reduction,occday[i+1])

            # And second we see if the chain has nu sub-30 min differences
            length = True
            if min_time:
                minlength = 99
                for i in location:
                    j = 0
                    while occday[i+j] == occday[i] and i+j < len(occday)-1:
                        j = j+1
                    if j < minlength:
                        minlength = j
                # and we neglect the very short presences of 20 min or less
                length = not minlength < 3

            # both have to be true to allow continuation, and we return boolean
            return shape and length

        def dayrun(start, cluster):
            '''
            Simulation of a single day according to start state and transition 
            probabilities in correspondence with the cluster-type.
            '''

            # we set the default daycheck at False for the first run and loop
            # creating days while False, meaning while the simulated day does
            # not correspond to the agreed-on rules in check().
            daycheck = False
            end = datetime.datetime.utcnow() + datetime.timedelta(seconds = 10)
            
            # create Monte Carlo Survival Analysis (MCSA) object in stats.py.
            SA = stats.MCSA(cluster)
            
            # and then keep simulating a day until daycheck is True
            while daycheck == False:
                # set start state conditions
                tbin = 0
                occs = np.zeros(144, dtype=int)
                occs[0] = start
                t48 = np.array(sorted(list(range(1, 49)) * 3))
                dt = SA.duration(start, t48[0])
                # and loop sequentially transition and duration functions
                while tbin < 143:
                    tbin += 1
                    if dt == 0:
                        occs[tbin] = SA.transition(occs[tbin-1], t48[tbin])
                        dt = SA.duration(occs[tbin], t48[tbin]) - 1
                        # -1 is necessary, as the occupancy state already started
                    else:
                        occs[tbin] = occs[tbin - 1]
                        dt += -1
                        
                # whereafer we control if this day is ok
                daycheck = check(occs)
                # and we include a break if the while-loop takes too long until
                # check()-conditions are fulfilled.
                if datetime.datetime.utcnow() > end:
                    break

            # return occupants array if daycheck is ok
            return occs

        def merge(occ):
            '''
            Merge the occupancy profiles of all household members to a single
            profile denoting the most active state of all members.
            '''
            # We start defining an array of correct length filled with the
            # least active state and loop to see if more-active people are
            # present at the depicted moment.
            occs = int(3)*np.ones(len(occ[0]))
            for member in occ:
                for to in range(len(member)):
                    if member[to] < occs[to]:
                        occs[to] = member[to]

            # return the merge occupancy states
            return occs

        # determine if the occupancy-chains are already defined on an external file
        # If ONFILE is set at 1 the external file will be used
        ONFILE = self.OccONFILE
        cdir = os.getcwd()
        os.chdir('../Data/InputFiles')
        d = os.getcwd()
        if ONFILE == 1:
            occ_ye = np.loadtxt(d + '/occ'+str(self.nBui)+'.txt', np.int32)
            if len(occ_ye) > 20:
                occ_year = np.zeros((1,len(occ_ye)))
                occ_year[0] = occ_ye
            else:
                occ_year = occ_ye.tolist()
            occ_merged = []
            occ_me = np.loadtxt(d + '/occ_m'+str(self.nBui)+'.txt', np.float64)
            occ_merged.append(occ_me)
            os.chdir(cdir)
            
        # If ONFILE is set at 0 the occupancy-chains are newly created
        else:
            os.chdir(cdir)
            occ_week = []
            # run the model for each daytype successively by which we can create a typical week.
            for member in self.clusters:
                SA = stats.MCSA(member['mon'])
                start = SA.startstate()
                mon = dayrun(start, member['mon'])
                tue = dayrun(mon[-1], member['tue'])
                wed = dayrun(tue[-1], member['wed'])
                thu = dayrun(wed[-1], member['thu'])
                fri = dayrun(thu[-1], member['fri'])
                sat = dayrun(fri[-1], member['sat'])
                sun = dayrun(sat[-1], member['sun'])
                week = np.concatenate((mon, tue, wed, thu, fri, sat, sun))
                occ_week.append(week)
                
            # A merge occupancy is created depicted the most active state of all household members.
            occ_merg = merge(occ_week)
            
            # combine the weekly occupancy states for the entire year by repeating them every week 
            # and correcting for the first day of year, and repeat for merged occ
            bins = 144
            tstart = bins*self.dow[0]
            tstop = tstart + bins*self.nday
            occ_year = []
            for line in range(len(occ_week)):
                occ_year.append(np.tile(occ_week,54)[line][tstart:tstop])
            occ_merged = []
            occ_merged.append(np.tile(occ_merg,54)[tstart:tstop])
            
            # print to file for re-use
            
            FILNAM = d + '/occ'+str(self.nBui)+'.txt'
            hea ='#1'+ '\n' + 'double data('+str(len(occ_year[0]))+','+str(len(occ_year))+')'
            np.savetxt(fname=FILNAM,header=hea, X = occ_year,fmt = '%d', newline="\n")
            FILNAM = d + '/occ_m'+str(self.nBui)+'.txt'
            hea ='#1'+ '\n' + 'double data('+str(len(occ_merged))+','+str(1)+')'
            np.savetxt(fname=FILNAM,header=hea, X = occ_merged,fmt = '%d')
            

        # return occupancies      
        self.occ = occ_year
        self.occ_m = occ_merged
        p = 0
        for t in range(len(occ_merged[0])):
            if occ_merged[0][t] < 3:
                p += 1
        presence = p/len(occ_merged[0])
        # print (presence)
        
    def __activity__(self):
        '''
        Activity model based on time use statistics from Belgium (2013). Modelling approach based on the work of Aerts (ref).
        Determines when different energy consuming activities and location specific activities are performed and by which householdmember.
        This results in activity profiles for each individual.
        
        There are two types of activities:
            - Household tasks: activities such as doing laundry, cooking,... that are usually performed by one person, but probability to perform the activity is based on household type
            - Personal activities: activities such as taking a shower, using a computer,... which can be independently performed from other household members
        
        In contrast to the occupancy model, the activity model is not continuous, this means that not all actions are followed directly by another action.
        This is because we only model the energy consuming and location specific activities, neglecting many other actions.
        '''
        
        nbin = 48
        nday = self.nday
        nstep = nday*nbin*3
        dow = self.dow
        npers = self.npers
        occ = self.occ
        occ_m = self.occ_m
        TASKS = ['cook','dishes','vacuum','iron','wash','dry']
        INDACT= ['pc','tv','audio','adm']
        order = ['FTE','PTE','Unemployed','Retired','School','Student']
        ACTtasks= np.zeros((6,npers,nstep+1))# order list TASKS
        ACTind= np.zeros((6,npers,nstep+1))# 0=pc, 1=tv, 2= audio, 3=adm, 4=(un)dress, 5=bath
        SA = stats.MCSA_task(npers)
        cdir = os.getcwd()
        os.chdir('../Data/InputFiles')
        d = os.getcwd()
        os.chdir(cdir)
        # If ONFILE is set at 1 the external file will be used
        if self.ActONFILE == 1:
            acttaskdata = np.loadtxt(d+'\ActTasks'+str(self.nBui)+'.txt', float)
            for i in range(6):
                acttask =  []
                for j in range(nstep):
                    acttask.append(acttaskdata[i*nstep+j])
                acttask.append(acttaskdata[-1])
                ACTtasks[i]=np.transpose(acttask)
            indactdata = np.loadtxt(d+'\ACTind'+str(self.nBui)+'.txt', float)
            for i in range(6):
                indact =  []
                for j in range(nstep):
                    indact.append(indactdata[i*nstep+j])
                indact.append(indactdata[-1])
                ACTind[i]=np.transpose(indact)
        else:
            # HOUSEHOLDTASKS
            # Determine when the household tasks are performed based on 
            # the housheold occupancy, size of the household and day of the week.
            calibrationfactor = np.loadtxt('../Data/HouseholdTasks/calibrationHHtasks.txt', float)
            avdurationHHtask = []
            for i in range(6):
                HHtask = TASKS[i]
                actdata = stats.DTMC_task(npers)
                cal = calibrationfactor[i]
                to = -1
                left = -1
                act_hh = np.zeros(nday*nbin*3+1)
                for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                    dow_i = dow[doy]
                    to += 1
                    if left <= 0:
                        occs = 1 if occ_m[0][to] == 1 else 0
                        prob = occs * cal * actdata.get_var(dow_i, HHtask, step)
                        if random.random() < prob:
                            left = SA.duration(i, step) - 1
                            act_hh[to] += 1
                    else:
                        left += -1
                        act_hh[to] += 1  
                avdurationHHtask.append(sum(act_hh)/(len(act_hh)))
                
                # assign household task to individual householdmembers based on 
                # each householdmembers availability and member probability 
                # (=prob that specific members of the household perform the task) 
                # (e.g. prob to cook higher for FTE than for school)
                Pmemb = np.loadtxt('Pmember.txt', float).T[i]
                dur = 0
                to = -1
                for doy, step in itertools.product(range(nday), range(nbin*3)):
                    to += 1
                    if act_hh[to-1] == 0 and act_hh[to] == 1:
                        dur += 1
                        start = to
                    elif act_hh[to] == 1:
                        dur += 1
                    elif act_hh[to-1] == 1 and act_hh[to] == 0:
                        end = to
                        k = -1
                        P_task = []
                        for member in self.members:
                            k += 1
                            emp_i = order.index(member)
                            P_memb = Pmemb[emp_i]
                            pres = 0
                            otheract = 0
                            for tz in range(start,end):
                                if occ[k][tz] == 1:
                                    pres += 1
                                    for task in range(i-1):
                                        if ACTtasks[task][k][tz] == 1:
                                            otheract += 1
                            if otheract > 0:
                                P_avi = 0
                            else:
                                P_av = pres/(end-start)
                                if P_av < 0.8:
                                    P_avi = 0
                                else:
                                    P_avi = P_av
                            P_assign = P_avi*P_memb
                            P_task.append(P_assign)
                        if sum(P_task) != 0:
                            Pmax = max(P_task)
                            ind = P_task.index(Pmax)
                            for tz in range(start,end):
                                ACTtasks[i][ind][tz] += 1
                        dur = 0
            # and print to calibrate
            FILE = open('../Data/HouseholdTasks/duration.txt',"a+")
            csv.writer(FILE, delimiter=' ',lineterminator = '\n').writerow(avdurationHHtask)
            FILE.close()
            # and print to file for re-use
            
            FILE = open(d + '/ACTtasks'+str(self.nBui)+'.txt',"a+")
            FILE.truncate(0)
            k = -1
            for task in TASKS:
                k += 1
                hea = '# Activity ' + task +'\n'
                FILE.write(hea)
                csv.writer(FILE, delimiter=' ',lineterminator = '\n').writerows(np.transpose(ACTtasks[k]))
            FILE.close()
            
            
            # Determine when individual activities are performed based on 
            # occupancy cluster, employment type and day of the week.
            calibrationfactor2 = np.loadtxt('../Data/IndividualActivities/calibrationIndAct.txt', float)
            p = -1
            for member in self.members:
                p += 1
                actdata = stats.DTMC_act(member)
                bathdata = stats.DTMC_bath(self.clusters[p])
                SIA = stats.MCSA_act(member)
                # LOCATION ACTIVITIES
                # Determine when individual activities are performed linked to 
                # the specific rooms.
                
                # BEDROOM (act = (un)dress = 4)
                act_bed = np.zeros(nday*nbin*3+1)
                to = -1
                left = -1
                for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                    dow_i = dow[doy]
                    to += 1
                    if occ[p][to] == 1 and occ[p][to-1]== 2:
                        left = SIA.duration(4, step) - 1
                        act_bed[to] += 1
                    elif left > 0:
                        for i in range(6):
                            otheract = 0
                            if ACTtasks[i][p][to] == 1:
                                otheract += 1
                        if otheract > 0 or occ[p][to] != 1:
                            left = 0
                        else:
                            left += -1
                            act_bed[to] += 1
                ACTind[4][p] = act_bed
                
                # BATHROOM (act = presence bathroom = 5)
                act_bath = np.zeros(nday*nbin*3+1)
                to = -1
                left = -1
                cal = calibrationfactor2[4]
                for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                    dow_i = dow[doy]
                    to += 1
                    if left <= 0:
                        if occ[p][to] == 1 and ACTind[4][p][to] == 0:
                            for i in range(6):
                                otheract = 0
                                if ACTtasks[i][p][to] == 1:
                                    otheract += 1
                            if otheract == 0:
                                occs = 1
                            else:
                                occs = 0
                        else:
                            occs = 0
                        prob = occs * cal * bathdata.get_var(dow_i, 'bath', step)
                        if random.random() < prob:
                            left = SIA.duration(5, step) - 1
                            act_bath[to] += 1
                    else:
                        if occ[p][to] == 1 and ACTind[4][p][to] == 0:
                            for i in range(6):
                                otheract = 0
                                if ACTtasks[i][p][to] == 1:
                                    otheract += 1
                            if otheract == 0:
                                left += -1
                                act_bath[to] += 1
                            else:
                                left = 0 
                ACTind[5][p] = act_bath
                
                # Determine occupancy in dayzone.
                occ_day = np.zeros(nday*nbin*3+1)
                occ_day_hhtask = np.zeros(nday*nbin*3+1)
                to = -1
                for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                    to += 1
                    otheract = 0
                    for i in range(6):                           
                        if ACTtasks[i][p][to] == 1:
                            otheract += 1
                    if occ[p][to] == 1 and ACTind[4][p][to] == 0 and ACTind[5][p][to] == 0 and otheract == 0:
                        occ_day[to] += 1
                    if occ[p][to] == 1 and ACTind[4][p][to] == 0 and ACTind[5][p][to] == 0:
                        occ_day_hhtask[to] += 1    
                
                # Determine when other individual activities are performed 
                # based on dayzone occupancy and activity compatibility
                for act in ['pc','adm']:
                    indexact = INDACT.index(act)
                    cal = calibrationfactor2[indexact]
                    act_ind = np.zeros(nday*nbin*3+1)
                    to = -1
                    left = -1
                    for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                        dow_i = dow[doy]
                        to += 1
                        if left <= 0:
                            occs = 1 if occ_day[to] == 1 else 0
                            prob = occs * cal * actdata.get_var(dow_i, act, step)
                            if random.random() < prob:
                                left = SIA.duration(indexact, step) - 1
                                act_ind[to] += 1
                        else:
                            if occ_day[to] == 1:
                                left += -1
                                act_ind[to] += 1
                            else:
                                left = 0
                    ACTind[indexact][p] = act_ind
                
                for act in ['tv']:
                    indexact = INDACT.index(act)
                    cal = calibrationfactor2[indexact]
                    act_ind = np.zeros(nday*nbin*3+1)
                    to = -1
                    left = -1
                    for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                        dow_i = dow[doy]
                        to += 1
                        if left <= 0:
                            if occ_day_hhtask[to] == 1 and ACTind[0][p][to] == 0 and ACTind[3][p][to] == 0 and ACTtasks[2][p][to]==0 and ACTtasks[4][p][to] == 0 and ACTtasks[5][p][to] == 0:
                                occs = 1
                            else:
                                occs = 0
                            prob = occs * cal * actdata.get_var(dow_i, act, step)
                            if random.random() < prob:
                                left = SIA.duration(indexact, step) - 1
                                act_ind[to] += 1
                        else:
                            if occ_day_hhtask[to] == 1 and ACTind[0][p][to] == 0 and ACTind[3][p][to] == 0 and ACTtasks[2][p][to]==0 and ACTtasks[4][p][to] == 0 and ACTtasks[5][p][to] == 0:
                                left += -1
                                act_ind[to] += 1
                            else:
                                left = 0
                    ACTind[indexact][p] = act_ind
                        
                for act in ['audio']:
                    indexact = INDACT.index(act)
                    cal = calibrationfactor2[indexact]
                    act_ind = np.zeros(nday*nbin*3+1)
                    to = -1
                    left = -1
                    for doy, step, step2 in itertools.product(range(nday), range(nbin), range(3)):
                        dow_i = dow[doy]
                        to += 1
                        if left <= 0:
                            if occ[p][to] == 1 and ACTind[1][p][to] == 0 and ACTtasks[2][p][to] == 0:
                                occs = 1
                            else:
                                occs = 0
                            prob = occs * cal * actdata.get_var(dow_i, act, step)
                            if random.random() < prob:
                                left = SIA.duration(indexact, step) - 1
                                act_ind[to] += 1
                        else:
                            if occ[p][to] == 1 and ACTind[1][p][to] == 0 and ACTtasks[2][p][to] == 0:
                                left += -1
                                act_ind[to] += 1
                            else:
                                left = 0
                    ACTind[indexact][p] = act_ind 
            
                # and print to calibrate
                durIndAct = np.zeros(5)
                k = -1
                for i in [0,1,2,3,5]:
                    k += 1
                    durIndAct[k] = sum(ACTind[i][p])/len(ACTind[i][p])
                FILE = open('../Data/IndividualActivities/duration.txt',"a+")
                csv.writer(FILE, delimiter=' ',lineterminator = '\n').writerow(durIndAct)
                FILE.close()
                
            # and print to file for re-use
            FILE = open(d + '/ACTind'+str(self.nBui)+'.txt',"a+")
            FILE.truncate(0)
            k = -1
            for act in ['pc','tv','audio','adm','dress','bath']:
                k += 1
                hea = '# Activity ' + act +'\n'
                FILE.write(hea)
                csv.writer(FILE, delimiter=' ',lineterminator = '\n').writerows(np.transpose(ACTind[k]))
            FILE.close()
            
            
        # return activityprofiles per activity per person   
        self.ACTind = ACTind
        self.ACTtasks = ACTtasks
        
        # return activityprofiles per activity 
        TASK = np.zeros((len(TASKS), nbin*3*nday+1))
        ACT = np.zeros((len(INDACT), nbin*3*nday+1)) 
        for i in range(len(TASKS)):
            for m in range(npers):
                TASK[i] += self.ACTtasks[i][m]
        for i in range(len(INDACT)):
            for m in range(npers):
                ACT[i] += self.ACTind[i][m]
        self.TASK = TASK
        self.ACT = ACT
                 

    def __metabolic__(self):
        '''
        Calculate metabolic rate and corresponding heat gains and CO2-production.
        Calculate room occupancy and assign load to corresponding room.
        '''        
        nbin = 144
        nday = self.nday
        npers = self.npers
        occ = self.occ
        ACTind = self.ACTind
        ACTtasks = self.ACTtasks
        nroom = self.nbedr + 3
        
        # lists per room with nroom: Dayzone(0), Kitchen(1), Bathroom(2), Bedroom(s)(3-nroom)
        CO2 = np.zeros((nroom,nday*nbin+1))
        MetHeat = np.zeros((nroom,nday*nbin+1))
        occ_room = np.zeros((nroom,nday*nbin+1))
        # properties for determination of metabolic rate and corresponding CO2- and heat-production.
        xp = [1,1.6]
        fp = [11.875,19]
        x = [0.8, 1.25]
        f = [70,110]
        
        for i in range(npers):
            memb = self.members[i]
            # properties for determination of metabolic rate and corresponding CO2- and heat-production.
            met_sl = random.gauss(0.95,0.1) # sleeping
            met_bath = random.gauss(1.5,0.1) # bathroom activities
            met_iron = random.gauss(3.8,0.1) # ironing, vacuuming
            met_tv = random.gauss(1.2,0.1) # watching tv, sitting
            met_cook = random.gauss(3.3,0.1) # cooking, working in kitchen
            met_other = random.gauss(1.5,0.1) # light effort, other activities
            met = [met_sl, met_bath, met_iron, met_tv, met_cook, met_other]
            CO2_i = np.interp(met, xp, fp)/6
            Heat_i = np.interp(met, x, f)
            B = self.BEDR[i]
            to = -1
            for doy, step in itertools.product(range(nday), range(nbin)):
                to += 1
                if occ[i][to] == 2:
                    occ_room[3+B][to] += 1
                    CO2[3+B][to] += CO2_i[0]
                    MetHeat[3+B][to] += Heat_i[0]
                elif ACTind[4][i][to] == 1:
                    occ_room[3+B][to] += 1
                    CO2[3+B][to] += CO2_i[5]
                    MetHeat[3+B][to] += Heat_i[5]
                elif ACTtasks[0][i][to] == 1 or ACTtasks[1][i][to] == 1:
                    occ_room[1][to] += 1
                    CO2[1][to] += CO2_i[4]
                    MetHeat[1][to] += Heat_i[4]
                elif ACTind[5][i][to] == 1:
                    occ_room[2][to] += 1
                    CO2[2][to] += CO2_i[1]
                    MetHeat[2][to] += Heat_i[1]
                elif ACTind[0][i][to] == 1 or ACTind[3][i][to] == 1:
                    if memb in ['Student','School']:
                        occ_room[3+B][to] += 1
                        CO2[3+B][to] += CO2_i[3]
                        MetHeat[3+B][to] += Heat_i[3]
                    else:
                        occ_room[0][to] += 1
                        CO2[0][to] += CO2_i[3]
                        MetHeat[0][to] += Heat_i[3]
                elif ACTtasks[2][i][to] == 1 or ACTtasks[3][i][to] == 1:
                    occ_room[0][to] += 1
                    CO2[0][to] += CO2_i[2]
                    MetHeat[0][to] += Heat_i[2]
                elif ACTind[1][i][to] == 1:
                    occ_room[0][to] += 1
                    CO2[0][to] += CO2_i[3]
                    MetHeat[0][to] += Heat_i[3]
                elif occ[i][to] == 1:
                    occ_room[0][to] += 1
                    CO2[0][to] += CO2_i[5]
                    MetHeat[0][to] += Heat_i[5]
        
        # return metabolic heat gains, co2 production and occupancies per room    
        self.occROOM = occ_room
        self.CO2 = CO2
        self.MetHeat = MetHeat
        
        return None

    def __plugload__(self):
        '''
        Simulation of the electric load of the appliances and lighting
        and the corresponding internal heat gains.
        '''

        def receptacles(self):
            '''
            Simulation of the use of appliances.
            Both for appliances linked to specific household tasks and individual 
            activities and for appliances not linked to specific activities.
            '''

            # define parameters
            nbin = 144
            nday = self.nday
            nroom = self.nbedr + 3
            npers = self.npers
            TASKS = ['cook','dishes','vacuum','iron','wash','dry']
            INDACT= ['pc','tv','audio','adm']
            TASK = self.TASK
            ACT = self.ACT
            appliances = [['Oven','Microwave','Hob','Kettle'],['DishWasher'],['Vacuum'],['Iron'],['WashingMachine'],['TumbleDryer']]
            
            power = np.zeros(nbin*nday+1)
            radi = np.zeros((nroom, nbin*nday+1))
            conv = np.zeros((nroom, nbin*nday+1))

            # for vacuum and iron assign load when performing activity
            for act in ['vacuum','iron']:
                ind = TASKS.index(act)
                APP = appliances[ind]
                for app in APP:
                    power_i = np.zeros(nbin*nday+1)
                    radi_i = np.zeros(nbin*nday+1)
                    conv_i = np.zeros(nbin*nday+1)
                    if app in self.apps:
                        load = set_appliances[app]['cycle_power']
                        standby_load = set_appliances[app]['standby_power']
                        frad = set_appliances[app]['frad']
                        fcon = set_appliances[app]['fconv']
                        to = -1
                        for doy, step in itertools.product(range(nday), range(nbin)):
                            to += 1
                            if TASK[ind][to] == 1:
                                power_i[to] += load
                                radi_i[to] += load*frad
                                conv_i[to] += load*fcon
                            else:
                                power_i[to] += standby_load
                                radi_i[to] += standby_load*frad
                                conv_i[to] += standby_load*fcon
                        power += power_i
                        radi[0] += radi_i
                        conv[0] += conv_i
                            
            # washmachine, dryer and dishwasher (when present in the household)
            # activated when the corresponding activity is finished.
            # e.g doing laundry mainly entails sorting the laundry and putting 
            # it in the wasmachine, so wasmachine starts at the end of doing laundry
            for act in ['dishes','wash','dry']:
                ind = TASKS.index(act)
                APP = appliances[ind]
                for app in APP:
                    power_i = np.zeros(nbin*nday+1)
                    radi_i = np.zeros(nbin*nday+1)
                    conv_i = np.zeros(nbin*nday+1)
                    if app in self.apps:
                        load = set_appliances[app]['cycle_power']
                        standby_load = set_appliances[app]['standby_power']
                        cycle_length = set_appliances[app]['cycle_length']/10
                        frad = set_appliances[app]['frad']
                        fcon = set_appliances[app]['fconv']
                        loc = set_appliances[app]['location']
                        cal = set_appliances[app]['cal']
                        to = -1
                        left = -1
                        for doy, step in itertools.product(range(nday), range(nbin)):
                            to += 1
                            if left < 0:
                                if TASK[ind][to-1] == 1 and TASK[ind][to] == 0 and random.random() < cal:
                                    power_i[to] += load
                                    radi_i[to] += load*frad
                                    conv_i[to] += load*fcon
                                    left = random.gauss(cycle_length, cycle_length/10)
                                else:
                                    power_i[to] += standby_load
                                    radi_i[to] += standby_load*frad
                                    conv_i[to] += standby_load*fcon
                            else:
                                if left < 1:
                                    power_i[to] += load*left
                                    radi_i[to] += load*left*frad
                                    conv_i[to] += load*left*fcon
                                    left += -1
                                else:
                                    power_i[to] += load
                                    radi_i[to] += load*frad
                                    conv_i[to] += load*fcon
                                    left += -1
                        power += power_i
                        radi[loc] += radi_i
                        conv[loc] += conv_i
                
            # for cooking not always appliances are used 
            # and different appliances can be used at the same time.
            for act in ['cook']:
                ind = TASKS.index(act)
                APP = appliances[ind]
                for app in APP:
                    power_i = np.zeros(nbin*nday+1)
                    radi_i = np.zeros(nbin*nday+1)
                    conv_i = np.zeros(nbin*nday+1)
                    standby_load = 100
                    if app in self.apps:
                        load_i = set_appliances[app]['cycle_power']
                        load = random.gauss(load_i, load_i/10)
                        standby_load = set_appliances[app]['standby_power']
                        cycle_length = set_appliances[app]['cycle_length']/10
                        frad = set_appliances[app]['frad']
                        fcon = set_appliances[app]['fconv']
                        cal = set_appliances[app]['cal']
                        to = -1
                        left = -1
                        for doy, step in itertools.product(range(nday), range(nbin)):
                            to += 1
                            if TASK[ind][to] == 1:
                                if left <= 0:
                                    if random.random() < cal:
                                        left = random.gauss(cycle_length, cycle_length/10) - 1
                                        power_i[to] += load
                                        radi_i[to] += load*frad
                                        conv_i[to] += load*fcon
                                    else:
                                        power_i[to] += standby_load
                                        radi_i[to] += standby_load*frad
                                        conv_i[to] += standby_load*fcon
                                else:
                                    if left < 1:
                                        power_i[to] += load*left
                                        radi_i[to] += load*left*frad
                                        conv_i[to] += load*left*fcon
                                        left += -1
                                    else:
                                        power_i[to] += load
                                        radi_i[to] += load*frad
                                        conv_i[to] += load*fcon
                                        left += -1
                            else:
                                power_i[to] += standby_load
                                radi_i[to] += standby_load*frad
                                conv_i[to] += standby_load*fcon
                                left = 0
                        
                    power += power_i
                    radi[1] += radi_i
                    conv[1] += conv_i

            # use pc when performing activity pc and 
            # when working/study/administration high chance of using pc as well
            pcmax = self.apps.count('PC')
            if pcmax >= npers:
                for act in ['pc', 'adm']:
                    ind = INDACT.index(act)
                    for m in range(npers):
                        power_i = np.zeros(nbin*nday+1)
                        radi_i = np.zeros(nbin*nday+1)
                        conv_i = np.zeros(nbin*nday+1)
                        load_i = set_appliances['PC']['cycle_power']
                        load = random.gauss(load_i, load_i/4)
                        frad = set_appliances['PC']['frad']
                        fcon = set_appliances['PC']['fconv']
                        standby_load = set_appliances['PC']['standby_power']
                        to = -1
                        usepc = 0
                        for doy, step in itertools.product(range(nday), range(nbin)):
                            to += 1
                            if act == 'pc':
                                if self.ACTind[ind][m][to] == 1:
                                    power_i[to] += load
                                    radi_i[to] += load*frad
                                    conv_i[to] += load*fcon
                                else:
                                    power_i[to] += standby_load
                                    radi_i[to] += standby_load*frad
                                    conv_i[to] += standby_load*fcon
                            else:
                                if self.ACTind[ind][m][to-1] == 0 and self.ACTind[ind][m][to-1] == 1 :
                                    prob = 0.5
                                    if random.random() < prob:
                                        usepc = 1
                                        power_i[to] += load
                                        radi_i[to] += load*frad
                                        conv_i[to] += load*fcon
                                    else:
                                        usepc = 0
                                elif self.ACTind[ind][m][to] == 1:
                                    if usepc == 1:
                                        power_i[to] += load
                                        radi_i[to] += load*frad
                                        conv_i[to] += load*fcon
                                    else:
                                        power_i[to] += standby_load
                                        radi_i[to] += standby_load*frad
                                        conv_i[to] += standby_load*fcon
                                else:
                                    power_i[to] += standby_load
                                    radi_i[to] += standby_load*frad
                                    conv_i[to] += standby_load*fcon
                        if self.members[m] in ['School', 'Student']:
                            B=self.BEDR[m]
                            power += power_i
                            radi[B+3] += radi_i
                            conv[B+3] += conv_i
                        else:
                            power += power_i
                            radi[0] += radi_i
                            conv[0] += conv_i
            else:
                pc = np.zeros(nbin*nday+1)
                for act in ['pc', 'adm']:
                    ind = INDACT.index(act)
                    for m in range(npers):
                        power_i = np.zeros(nbin*nday+1)
                        radi_i = np.zeros(nbin*nday+1)
                        conv_i = np.zeros(nbin*nday+1)
                        load_i = set_appliances['PC']['cycle_power']
                        load = random.gauss(load_i, load_i/4)
                        frad = set_appliances['PC']['frad']
                        fcon = set_appliances['PC']['fconv']
                        standby_load = set_appliances['PC']['standby_power']
                        to = -1
                        usepc = 0
                        for doy, step in itertools.product(range(nday), range(nbin)):
                            to += 1
                            if pc[to] < pcmax:
                                if act == 'pc':
                                    if self.ACTind[ind][m][to] == 1:
                                        power_i[to] += load
                                        radi_i[to] += load*frad
                                        conv_i[to] += load*fcon
                                        pc[to] += 1
                                else:
                                    if self.ACTind[ind][m][to-1] == 0 and self.ACTind[ind][m][to-1] == 1 :
                                        prob = 0.5
                                        if random.random() < prob:
                                            usepc = 1
                                            power_i[to] += load
                                            radi_i[to] += load*frad
                                            conv_i[to] += load*fcon
                                            pc[to] += 1
                                        else:
                                            usepc = 0
                                    elif self.ACTind[ind][m][to] == 1:
                                        if usepc == 1:
                                            power_i[to] += load
                                            radi_i[to] += load*frad
                                            conv_i[to] += load*fcon
                                            pc[to] += 1
                        if self.members[m] in ['School', 'Student']:
                            B=self.BEDR[m]
                            power += power_i
                            radi[B+3] += radi_i
                            conv[B+3] += conv_i
                        else:
                            power += power_i
                            radi[0] += radi_i
                            conv[0] += conv_i
                power_i = np.zeros(nbin*nday+1)
                radi_i = np.zeros(nbin*nday+1)
                conv_i = np.zeros(nbin*nday+1)
                to = -1
                for doy, step in itertools.product(range(nday), range(nbin)):
                    to += 1
                    if pc[to] < pcmax:
                        power_i[to] += standby_load*(pcmax-pc[to])
                        radi_i[to] += standby_load*(pcmax-pc[to])*frad
                        conv_i[to] += standby_load*(pcmax-pc[to])*fcon
                power += power_i
                radi[0] += radi_i
                conv[0] += conv_i
            
            for act in ['pc']:
                if 'Printer' in self.apps:
                    ind = INDACT.index(act)
                    power_i = np.zeros(nbin*nday+1)
                    radi_i = np.zeros(nbin*nday+1)
                    conv_i = np.zeros(nbin*nday+1)
                    load_i = set_appliances['Printer']['cycle_power']
                    load = random.gauss(load_i, load_i/10)
                    standby_load = set_appliances['Printer']['standby_power']
                    frad = set_appliances['Printer']['frad']
                    fcon = set_appliances['Printer']['fconv']
                    cycle_length = set_appliances['Printer']['cycle_length']/10
                    cal = set_appliances['Printer']['cal']
                    for m in range(npers):
                        to = -1
                        for doy, step in itertools.product(range(nday), range(nbin)):
                            to += 1
                            if ACT[ind][to] > 0:
                                if left <= 0:
                                    if random.random() < cal:
                                        left = random.gauss(cycle_length, cycle_length/10) - 1
                                        if left < 0:
                                            power_i[to] = load*(left+1)
                                            radi_i[to] = load*frad*(left+1)
                                            conv_i[to] = load*fcon*(left+1)
                                        else:
                                            power_i[to] = load
                                            radi_i[to] = load*frad
                                            conv_i[to] = load*fcon
                                    else:
                                        power_i[to] = standby_load
                                        radi_i[to] = standby_load*frad
                                        conv_i[to] = standby_load*fcon
                                else:
                                    if left < 1:
                                        power_i[to] += load*left
                                        radi_i[to] += load*left*frad
                                        conv_i[to] += load*left*fcon
                                        left += -1
                                    else:
                                        left += -1
                                        power_i[to] = load
                                        radi_i[to] = load*frad
                                        conv_i[to] = load*fcon
                            else:
                                power_i[to] = standby_load
                                radi_i[to] = standby_load*frad
                                conv_i[to] = standby_load*fcon
                                    
                    power += power_i
                    radi[0] += radi_i
                    conv[0] += conv_i
            
            # use hifi when performing activity listening to audio
            power_i = np.zeros(nbin*nday+1)
            radi_i = np.zeros(nbin*nday+1)
            conv_i = np.zeros(nbin*nday+1)
            for act in ['audio']:
                ind = INDACT.index(act)
                load_i = set_appliances['HiFi']['cycle_power']
                load = random.gauss(load_i, load_i/10)
                standby_load = set_appliances['HiFi']['standby_power']
                frad = set_appliances['HiFi']['frad']
                fcon = set_appliances['HiFi']['fconv']
                to = -1
                for doy, step in itertools.product(range(nday), range(nbin)):
                    to += 1
                    if ACT[ind][to] > 0:
                        power_i[to] += load
                        radi_i[to] += load*frad
                        conv_i[to] += load*fcon
                    else:
                        power_i[to] += standby_load
                        radi_i[to] += standby_load*frad
                        conv_i[to] += standby_load*fcon
                        
            power += power_i
            radi[0] += radi_i
            conv[0] += conv_i
            
            # use tv when performing activity tv, sharing occurs often
            # each tv is assumed to be equiped with digital tv requiring some 
            # kind of digibox/digicorder (which is included in the load of the tv)
            power_i = np.zeros(nbin*nday+1)
            radi_i = np.zeros(nbin*nday+1)
            conv_i = np.zeros(nbin*nday+1)
            for act in ['tv']:
                ind = INDACT.index(act)
                tvmax = self.apps.count('tv')
                load_i = set_appliances['tv']['cycle_power']
                load = random.gauss(load_i, load_i/4)
                standby_load = set_appliances['tv']['standby_power']
                frad = set_appliances['tv']['frad']
                fcon = set_appliances['tv']['fconv']
                to = -1
                for doy, step in itertools.product(range(nday), range(nbin)):
                    to += 1
                    if ACT[ind][to] <= tvmax:
                        power_i[to] += load*ACT[ind][to] + standby_load*(tvmax-ACT[ind][to])
                        radi_i[to] += load*frad*ACT[ind][to] + standby_load*(tvmax-ACT[ind][to])*frad
                        conv_i[to] += load*fcon*ACT[ind][to] + standby_load*(tvmax-ACT[ind][to])*fcon
                    else:
                        power_i[to] += load*tvmax
                        radi_i[to] += load*frad*tvmax
                        conv_i[to] += load*fcon*tvmax
                        
            power += power_i
            radi[0] += radi_i
            conv[0] += conv_i         
            
            # APPLIANCES NOT LINKED TO ACTIVITY
            # appliances such as fridge and freezer are continuously on
            # and work according to a fixed cycle.
            cyclingApp = ['FridgeFreezer', 'Freezer','Refrigerator']
            for app in cyclingApp:
                if app in self.apps:
                    cycledelay = set_appliances[app]['delay']/10
                    load_i = set_appliances[app]['cycle_power']
                    load = random.gauss(load_i, load_i/10)
                    standby_load = set_appliances[app]['standby_power']
                    frad = set_appliances[app]['frad']
                    fcon = set_appliances[app]['fconv']
                    cycle_length = set_appliances[app]['cycle_length']/10
                    power_i = np.zeros(nbin*nday+1)
                    radi_i = np.zeros(nbin*nday+1)
                    conv_i = np.zeros(nbin*nday+1)
                    to = -1
                    left = -1
                    delay = random.gauss(cycledelay, cycledelay/4)
                    for doy, step in itertools.product(range(nday), range(nbin)):
                        to += 1
                        if left <= 0:
                            power_i[to] += standby_load
                            radi_i[to] += standby_load*frad
                            conv_i[to] += standby_load*fcon
                            if delay <= 1:
                                left = random.gauss(cycle_length, cycle_length/4)
                                delay = random.gauss(cycledelay, cycledelay/4)
                            else:
                                delay += -1
                        else:
                            left += -1
                            power_i[to] += load
                            radi_i[to] += load*frad
                            conv_i[to] += load*fcon
                    power += power_i
                    radi[0] += radi_i
                    conv[0] += conv_i
        
            electr = ['Smartphone']
            for app in electr:
                smartphones = self.apps.count('Smartphone')
                for s in range(min(smartphones,npers)):
                    load_i = set_appliances[app]['cycle_power']
                    load = random.gauss(load_i, load_i/10)
                    standby_load = set_appliances[app]['standby_power']
                    frad = set_appliances[app]['frad']
                    fcon = set_appliances[app]['fconv']
                    cal = set_appliances[app]['cal']
                    cycle_length = set_appliances[app]['cycle_length']
                    power_i = np.zeros(nbin*nday+1)
                    radi_i = np.zeros(nbin*nday+1)
                    conv_i = np.zeros(nbin*nday+1)
                    to = -1
                    left = -1
                    n_s = 0
                    for doy, step in itertools.product(range(nday), range(nbin)):
                        to += 1
                        if self.occ[s][to] == 1:
                            if left <= 0:
                                if random.random() < cal:
                                    left = random.gauss(cycle_length, cycle_length/10) - 1
                                    power_i[to] += load
                                    radi_i[to] += load*frad
                                    conv_i[to] += load*fcon
                                    n_s += 1
                                else:
                                    power_i[to] += standby_load
                                    radi_i[to] += standby_load*frad
                                    conv_i[to] += standby_load*fcon
                            else:
                                left += -1
                                power_i[to] += load
                                radi_i[to] += load*frad
                                conv_i[to] += load*fcon
                    power += power_i
                    radi[0] += radi_i
                    conv[0] += conv_i 
        
            # output 
            load = int(np.sum(power)/6/1000)
            print (' - Appliance load is %s kWh' % str(load))
            
            # return total electricity use and internal heat gains per room    
            self.power = power
            self.radi = radi
            self.conv = conv

            return None

        def lightingload(self):
            '''
            Simulate use of lighting based on occupancy in the room and 
            solar radiation. Adaptive model in dayzone based on Widen (2009).
            '''

            # information on occupancy 
            occROOM = self.occROOM
            nroom = len(occROOM)
            
            # define load in the different rooms
            power_max_dayzone = 100
            power_max_room = 20
            prob_adj = 0.8 # hourly probability to adjust
            pow_adj = 25 # power by which is adjusted
            frad = 0.55
            fcon = 0.45
            
            nday = self.nday
            nbin = 144
            irr_max = 200
            nbedr = self.nbedr
            
            # load solar radiation levels which determine the need for lighting.
            # The loaded data represent the global horizontal radiation
            # at a time-step of 1-hour for Uccle, Belgium.
            cdir = os.getcwd()
            os.chdir('..')
            d = os.getcwd()
            SI = np.loadtxt(d + '/Data/Climate/SI.txt',str)
            os.chdir(cdir)
            
            irr = np.zeros(nday*nbin+1)
            to = -1
            tl = -1
            for doy, step in itertools.product(range(nday), range(24)):
                to += 1
                for step2 in range(6):
                    tl += 1
                    irr[tl] = SI[to]
            
            lightingload = np.zeros((nroom,nbin*nday+1))
            lightcon = np.zeros((nroom,nbin*nday+1))
            lightrad = np.zeros((nroom,nbin*nday+1))
            pow_id = np.zeros(nday*nbin+1)
            to = -1
            for doy, step in itertools.product(range(nday), range(nbin)):
                to += 1
                for room in range(nroom-nbedr):
                    if occROOM[room][to] > 0 and irr[to] <= irr_max:
                        if room == 0:
                            if random.random() <= prob_adj:
                                pow_id[to] = power_max_dayzone*(1 - irr[to]/irr_max)
                                delta = lightingload[room][to-1] - pow_id[to]
                                delta_min = np.abs(delta - pow_adj)
                                delta_plus = np.abs(delta + pow_adj)
                                if delta > 0 and delta_min < np.abs(delta) :
                                    lightday = lightingload[room][to-1]-pow_adj
                                elif delta < 0 and delta_plus < np.abs(delta):
                                    lightD = lightingload[room][to-1]+pow_adj
                                    if lightD > power_max_dayzone:
                                        lightday = power_max_dayzone
                                    else:
                                        lightday = lightD
                                else:
                                    lightday = lightingload[room][to-1]
                            else:
                                lightday = lightingload[room][to-1]
                            lightingload[room][to] += lightday
                            lightcon[room][to] = lightday*fcon
                            lightrad[room][to] = lightday*frad
                        else:
                            lightingload[room][to] += power_max_room
                            lightcon[room][to] = power_max_room*fcon
                            lightrad[room][to] = power_max_room*frad

                for m in range(self.npers):
                    room = self.BEDR[m]
                    if occROOM[3+room][to] > 0 and self.occ[m][to] != 2 and irr[to] <= irr_max:
                        lightingload[3+room][to] = power_max_room
                        lightcon[3+room][to] = power_max_room*fcon
                        lightrad[3+room][to] = power_max_room*frad                               
                        
            # output 
            load = int(np.sum(lightingload)/6/1000)
            print (' - Lighting load is %s kWh' % str(load))

            # return total electricity use and internal heat gains per room    
            self.lightload = lightingload
            self.lightRad = lightrad
            self.lightCon = lightcon
            
            return None

        receptacles(self)
        lightingload(self)

        return None

    def __dhwload__(self):
        '''
        Simulate use of domestic hot water based on data from Instal2020 (ref)
        DHW is determined per minute
        '''
        nmin = self.nday*24*60
        flow = np.zeros(nmin+1)
        bathroompresence = self.ACTind[5]
        tapcycles = []
        for tap in self.taps:
            eq = Equipment(**set_appliances[tap])
            r_tap, n_tap = eq.stochastic_flow(tap, self.nday, self.dow, self.occ_m, bathroompresence, self.npers, self.ACTtasks, self.apps)
            flow += r_tap['mDHW']
            tapcycles.append(n_tap)
        self.DHW = flow
        
        load = np.sum(flow)
        self.loadpppd = int(load/self.nday/self.npers)
        print (' - Draw-off is %s l/pp.day' % str(self.loadpppd))
        
        # and print to calibrate
        # FILE = open('../Data/DHW/DHW.txt',"a+")
        # csv.writer(FILE, delimiter=' ',lineterminator = '\n').writerow(tapcycles)
        # FILE.close()
        
        return None
    
    def __shsetting__(self):
        '''
        Simulation of the space heating setting points.
        - Including weekend days,
        - starting from a regular monday at 4:00 AM.
        '''

        # we define setting types based on their setpoint temperatures 
        # when being active (1), sleeping (2) or absent (3).
        types = dict()
        types.update({'1' : {1:17.0, 2:14.0, 3:15.0}})
        types.update({'2' : {1:18.5, 2:15.0, 3:18.5}})
        types.update({'3' : {1:20.0, 2:15.0, 3:19.5}})
        types.update({'4' : {1:20.0, 2:11.0, 3:19.5}})
        types.update({'5' : {1:20.0, 2:14.5, 3:15.0}})
        types.update({'6' : {1:21.0, 2:20.5, 3:21.0}})
        types.update({'7' : {1:21.5, 2:15.5, 3:21.5}})
        
        # and the probabilities these types occur based on Dutch research,
        # i.e. Leidelmeijer and van Grieken (2005).
        types.update({'prob' : [0.04, 0.16, 0.35, 0.08, 0.11, 0.05, 0.20]})
        
        # and given a type, denote which rooms are heated
        given = dict()
        given.update({'1' : [[],['dayzone'],['dayzone','nightzone'],['dayzone','bathroom','nightzone']]})
        given.update({'2' : [['dayzone','bathroom']]})
        given.update({'3' : [['dayzone'],['dayzone','bathroom'],['dayzone','nightzone']]})
        given.update({'4' : [['dayzone'],['dayzone','nightzone']]})
        given.update({'5' : [['dayzone']]})
        given.update({'6' : [['dayzone','bathroom','nightzone']]})
        given.update({'7' : [['dayzone','bathroom']]})

        # select a type from the given types and probabilities
        if int(self.shtype) <= 0:
            rnd = np.random.random()
            shtype = str(stats.get_probability(rnd, types['prob'], 'prob'))
        else:
            shtype = self.shtype
            
        if not self.shrooms:
            if len(given[shtype]) != 1:
                nr = np.random.randint(np.shape(given[shtype])[0])
                shrooms = given[shtype][nr]
            else:
                shrooms = given[shtype][0]
        else:
            shrooms = self.shrooms

        # create a profile for the heated rooms
        shnon = 12*np.ones(len(self.occ_m[0])+1)
        shset = np.hstack((self.occ_m[0],self.occ_m[0][-1]))
        for key in types[shtype].keys():
            for i in range(len(shset)):
                if int(shset[i]) == key:
                    shset[i] = types[shtype][key]

        # and couple to the heated rooms
        sh_settings = dict()
        for room in ['dayzone', 'nightzone', 'bathroom']:
            if room in shrooms:
                sh_settings.update({room:shset})
            else:
                sh_settings.update({room:shnon})
                
        # and store
        self.sh_settings = sh_settings
        print (' - Average comfort setting of dayzone is %s °C' % str(round(np.average(sh_settings['dayzone']),2)))
        
        return None
    
    def __window__(self, **kwargs):
    
        '''
        Simulation of window opening behaviour based on window use habits.
        For more information regarding the implemented model (REF)
        - 0: window is closed
        - 1: window is opened
        '''
        
        ######################################################################
        # Determine the window use habits if not already defined in inputs

        if not self.habits:
            habits = windowhabits.get_habits(self.VentS, self.DW, self.YearBuilt, self.members, self.HHhabit, self.SeCo)
        else:
            habits = self.habits

        ####################################################################################################################################
        # assign habit name to habit number
        # possible habits the household perform when going to bed in winter
        habits_sleep_w = dict()
        habits_sleep_w.update({0:'no change'})
        habits_sleep_w.update({1:'Close'})
        habits_sleep_w.update({2:'CloseAcc'})
        habits_sleep_w.update({3:'CloseExclBed'})
        habits_sleep_w.update({5:'CloseBed'})
        habits_sleep_w.update({6:'Open'})
        
        # possible habits the household perform when going to bed in summer
        habits_sleep_s = dict()
        habits_sleep_s.update({0:'no change'})
        habits_sleep_s.update({1:'Close'})
        habits_sleep_s.update({2:'CloseAcc'})
        habits_sleep_s.update({3:'CloseExclBed'})
        habits_sleep_s.update({4:'OpenBed'})
        habits_sleep_s.update({5:'Open'})
        habits_sleep_s.update({6:'CloseBed'})

        # possible habits the household perform when leaving the home in winter
        habits_away_w = dict()
        habits_away_w.update({0:'no change'})
        habits_away_w.update({1:'Close'})
        habits_away_w.update({2:'CloseAcc'})
        habits_away_w.update({3:'CloseExclBed'})
        habits_away_w.update({4:'OpenBed'})
        habits_away_w.update({5:'Open'})
        
        # possible habits the household perform when leaving the home in summer
        habits_away_s = dict()
        habits_away_s.update({0:'no change'})
        habits_away_s.update({1:'Close'})
        habits_away_s.update({2:'CloseAcc'})
        habits_away_s.update({3:'CloseExclBed'})
        habits_away_s.update({4:'OpenBed'})
        habits_away_s.update({5:'Open'})
        
        # determine which habits the household performs in the bedroom in winter
        habits_bedr_w = dict()
        habits_bedr_w.update({-1:'Closed'}) #assumption if no windows --> window is always closed
        habits_bedr_w.update({0:'Undefined'})
        habits_bedr_w.update({1:'Closed'})
        habits_bedr_w.update({2:'OpenMorningShort'})
        habits_bedr_w.update({3:'OpenMorningLong'})
        habits_bedr_w.update({4:'Open'})
        habits_bedr_w.update({5:'OpenNight'})
        habits_bedr_w.update({6:'OpenBeforeSleep'})
    
        # determine which habits the household performs in the bedroom in summer
        habits_bedr_s = dict()
        habits_bedr_s.update({-1:'Closed'})
        habits_bedr_s.update({0:'Undefined'})
        habits_bedr_s.update({1:'Closed'})
        habits_bedr_s.update({2:'OpenMorningShort'})
        habits_bedr_s.update({3:'OpenMorningLong'})
        habits_bedr_s.update({4:'Open'})
        habits_bedr_s.update({5:'OpenNight'})
        habits_bedr_s.update({6:'OpenBeforeSleep'})
        habits_bedr_s.update({7:'OpenDay'})
        
        # determine which habits the household performs in the bathroom in winter
        habits_bath_w = dict()
        habits_bath_w.update({-1:'Closed'})
        habits_bath_w.update({0:'Undefined'}) 
        habits_bath_w.update({1:'Closed'})
        habits_bath_w.update({2:'OpenPresent'})
        habits_bath_w.update({3:'OpenAfterPresent'})
        habits_bath_w.update({4:'OpenMorningShort'})
        habits_bath_w.update({5:'Open'})
        habits_bath_w.update({6:'ClosedPresent'})
        habits_bath_w.update({7:'OpenMorningLong'})
        
        # determine which habits the household performs in the bathroom in summer
        habits_bath_s = dict()
        habits_bath_s.update({-1:'Closed'})
        habits_bath_s.update({0:'Undefined'}) 
        habits_bath_s.update({1:'Closed'})
        habits_bath_s.update({2:'OpenPresent'})
        habits_bath_s.update({3:'OpenAfterPresent'})
        habits_bath_s.update({4:'OpenMorningShort'})
        habits_bath_s.update({5:'Open'})
        habits_bath_s.update({6:'ClosedPresent'})
        habits_bath_s.update({7:'OpenMorningLong'})
        habits_bath_s.update({8:'OpenNight'})
        
        # determine which habits the household performs in the livingroom in winter
        habits_living_w = dict()
        habits_living_w.update({-1:'Closed'})
        habits_living_w.update({0:'Undefined'}) 
        habits_living_w.update({1:'Closed'})
        habits_living_w.update({2:'OpenMorningShort'})
        habits_living_w.update({3:'OpenMorningLong'})
        habits_living_w.update({4:'OpenPresent'})
        habits_living_w.update({5:'OpenEvening'})
        habits_living_w.update({6:'ClosedPresent'})
        habits_living_w.update({7:'Open'})
        
        # determine which habits the household performs in the livingroom in summer
        habits_living_s = dict()
        habits_living_s.update({-1:'Closed'})
        habits_living_s.update({0:'Undefined'}) 
        habits_living_s.update({1:'Closed'})
        habits_living_s.update({2:'OpenMorningShort'})
        habits_living_s.update({3:'OpenMorningLong'})
        habits_living_s.update({4:'OpenPresent'})
        habits_living_s.update({5:'OpenNight'})
        habits_living_s.update({6:'Open'})
        habits_living_s.update({7:'OpenEvening'})
        habits_living_s.update({8:'ClosedPresent'})
        
        # determine which habits the household performs in the kitchen in winter
        habits_kitchen_w = dict()
        habits_kitchen_w.update({-1:'Closed'})
        habits_kitchen_w.update({0:'Undefined'}) 
        habits_kitchen_w.update({1:'Closed'})
        habits_kitchen_w.update({2:'OpenCooking'})
        habits_kitchen_w.update({3:'OpenAfterCooking'})
        habits_kitchen_w.update({4:'OpenMorningShort'})
        habits_kitchen_w.update({5:'OpenEvening'})
        habits_kitchen_w.update({6:'Open'})
        habits_kitchen_w.update({7:'OpenPresent'})
        habits_kitchen_w.update({8:'OpenMorningLong'})
        
        # determine which habits the household performs in the kitchen in summer
        habits_kitchen_s = dict()
        habits_kitchen_s.update({-1:'Closed'})
        habits_kitchen_s.update({0:'Undefined'}) 
        habits_kitchen_s.update({1:'Closed'})
        habits_kitchen_s.update({2:'OpenCooking'})
        habits_kitchen_s.update({3:'OpenAfterCooking'})
        habits_kitchen_s.update({4:'OpenMorningShort'})
        habits_kitchen_s.update({5:'OpenEvening'})
        habits_kitchen_s.update({6:'Open'})
        habits_kitchen_s.update({7:'OpenPresent'})
        habits_kitchen_s.update({8:'OpenMorningLong'})
        habits_kitchen_s.update({9:'OpenNight'})

        ######################################################################
        # Simulate window use
        # First determine the probabilities of opening a window for each room 
        # in case no habits are present. 
        # Model of Maeyens and Janssend based on weather.
        nday = self.nday
        nhour = nday*24
        nbin = 144
        nbedr = self.nbedr
        nroom = nbedr+3
        P = np.zeros((nroom, nhour))
        Pi = np.zeros((nroom, nhour))
        W = np.zeros((nroom, nhour))
        W_P = np.zeros((nroom, nbin*nday+1))
        
        # load weather variables
        cdir = os.getcwd()
        os.chdir('..')
        d = os.getcwd()
        SI_file = np.loadtxt(d + '/Data/Climate/SI.txt', np.float64)
        SI = SI_file.tolist() # Solar irradiation W/m²
                
        Te_file = np.loadtxt(d + '/Data/Climate/Te.txt', np.float64)
        Te = Te_file.tolist() # Outdoor temperature

        v_file = np.loadtxt(d + '/Data/Climate/v.txt', np.float64)
        v = v_file.tolist() # Wind velocity m/s
            
        F_file = np.loadtxt(d + '/Data/Climate/F.txt', np.float64)
        F = F_file.tolist() # Francostoro presence correction factor
        
        season = np.loadtxt(d + '/Data/Climate/season.txt',str)
        os.chdir(cdir)
        
        for i in range(nroom):
            to = -1
            for k in range(nhour):
                if SI[k] > 300: 
                    Pi[i][k] = (2.1*Te[k]+3)/100
                else:
                    Pi[i][k] = (1.16*Te[k]+3)/100
                if v[k] > 3: 
                    P[i][k] = (1.55-0.183*v[k])*Pi[i][k]*F[k]
                else:
                    P[i][k] = Pi[i][k]*F[k]
                r = random.random()
                if r < P[i][k]:
                    W[i][k] = 1 
                for run in range(6):
                    to += 1
                    if W[i][k] == 1:
                        W_P[i][to] = 1
                    else:
                        W_P[i][to] = 0
        
        ######################################################################
        # Simulate window use per 10 min.
        # Simulation of window use requires information on occupancy and activity 
        occROOM = self.occROOM
        occ_m = self.occ_m
        occ = self.occ
        cooking = self.TASK[0]
        window = -1 * np.ones((nroom, nbin*nday+1))
        to = -1
        k = 0 # time counters for morning habit
        l = 0
        m = 0
        n = 0
        for doy, step in itertools.product(range(nday), range(nbin)):
            to += 1
            if season[to] == 'Winter':
                habit_living = habits_living_w[habits[0]]
                habit_kitchen = habits_kitchen_w[habits[2]]
                habit_bedr = habits_bedr_w[habits[4]]
                habit_bath = habits_bath_w[habits[6]]
                habit_away = habits_away_w[habits[8]]
                habit_sleep = habits_sleep_w[habits[10]]

            else:
                habit_living = habits_living_s[habits[1]]
                habit_kitchen = habits_kitchen_s[habits[3]]
                habit_bedr = habits_bedr_s[habits[5]]
                habit_bath = habits_bath_s[habits[7]]
                habit_away = habits_away_s[habits[9]]
                habit_sleep = habits_sleep_s[habits[11]]

            # Simulate window use in the livingroom when habits are present
            if habit_living == 'Open':
                window[0][to] = 1
            elif habit_living == 'Closed':
                window[0][to] = 0
            elif habit_living == 'OpenNight':
                if occ_m[0][to] == 2:
                    window[0][to] = 1
                else:
                    window[0][to] = 0 
            elif habit_living == 'OpenMorningShort':# 30 minutes open after sleep
                if occ_m[0][to-1] == 2 and occ_m[0][to] == 1:
                    window[0][to] = 1
                    m = 3
                elif m > 0 and occ_m[0][to] == 1:
                    window[0][to] = 1
                    m += -1
                else:
                    window[0][to] = 0
                    m = 0
            elif habit_living == 'OpenMorningLong':# 3 hours open after sleep or until everyone leaves
                if occ_m[0][to-1] == 2 and occ_m[0][to] == 1:
                    window[0][to] = 1
                    m = 18
                elif m > 0 and occ_m[0][to] == 1:
                    window[0][to] = 1
                    m += -1
                else:
                    window[0][to] = 0
                    m = 0
            elif habit_living == 'OpenPresent':
                if occROOM[0][to] == 0:
                    window[0][to] = 0
                else:
                    window[0][to] = 1
            elif habit_living == 'ClosedPresent':
                if occROOM[0][to] == 0:
                    window[0][to] = 1
                else:
                    window[0][to] = 0
            elif habit_living == 'OpenEvening':# 3 hours open before sleep
                if to < nday*nbin-17:
                    if occ_m[0][to+17] == 2 and occ_m[0][to] == 1:
                        window[0][to] = 1
                    else:
                        window[0][to] = 0
                else:
                    if occ_m[0][nday*nbin] == 2 and occ_m[0][to] == 1:
                        window[0][to] = 1
                    else:
                        window[0][to] = 0
                            
            # Simulate window use in the kitchen when habits are present
            if habit_kitchen == 'Open':
                window[1][to] = 1
            elif habit_living == 'Closed':
                window[1][to] = 0
            elif habit_kitchen == 'OpenNight':
                if occ_m[0][to] == 2:
                    window[1][to] = 1
                else:
                    window[1][to] = 0 
            elif habit_kitchen == 'OpenMorningShort':# 30 minutes open after sleep
                if occ_m[0][to-1] == 2 and occ_m[0][to] == 1:
                    window[1][to] = 1
                    n = 3
                elif n > 0 and occ_m[0][to] == 1:
                    window[1][to] = 1
                    n += -1
                else:
                    window[1][to] = 0
                    n = 0
            elif habit_kitchen == 'OpenMorningLong':# 3 hours open after sleep
                if occ_m[0][to-1] == 2 and occ_m[0][to] == 1:
                    window[1][to] = 1
                    n = 18
                elif n > 0 and occ_m[0][to] == 1:
                    window[1][to] = 1
                    n += -1
                else:
                    window[1][to] = 0
                    n = 0
            elif habit_kitchen == 'OpenPresent':
                if occROOM[1][to] == 0:
                    window[1][to] = 0
                else:
                    window[1][to] = 1
            elif habit_kitchen == 'OpenEvening':# 3 hours open before sleep
                if to < nday*nbin-17:
                    if occ_m[0][to+17] == 2 and occ_m[0][to] == 1:
                        window[1][to] = 1
                    else:
                        window[1][to] = 0
                else:
                    if occ_m[0][nday*nbin] == 2 and occ_m[0][to] == 1:
                        window[1][to] = 1
                    else:
                        window[1][to] = 0
            elif habit_kitchen == 'OpenCooking':
                if cooking[to] == 1:
                    window[1][to] = 1
                else:
                    window[1][to] = 0
            elif habit_kitchen == 'OpenAfterCooking':
                if cooking[to] == 0:
                    if cooking[to-1] == 1:
                        window[1][to] = 1
                        n = 6
                    elif n > 0 and occ_m[0][to] == 1:
                        window[1][to] = 1
                        n += -1
                    else:
                        window[1][to] = 0
                        n = 0
                else:
                    window[1][to] = 0
                    n = 0
                    
            # Simulate window use in the bedroom when habits are present
            if habit_bedr == 'Open':
                for i in range(nbedr):
                    window[i+3][to] = 1
            elif habit_bedr == 'Closed':
                for i in range(nbedr):
                    window[i+3][to] = 0
            elif habit_bedr == 'OpenDay':
                for i in range(nbedr):
                    if occROOM[i+3][to] == 0:
                        window[i+3][to] = 1
                    else:
                        window[i+3][to] = 0
            elif habit_bedr == 'OpenNight':
                for i in range(self.npers):
                    if occ[i][to] == 2:
                        window[self.BEDR[i]+3][to] = 1
                    else:
                        if window[self.BEDR[i]+3][to] < 1:
                            window[self.BEDR[i]+3][to] = 0
                        else:
                            window[self.BEDR[i]+3][to] = 1
            elif habit_bedr == 'OpenMorningShort':# 3o minutes open after sleep
                for i in range(nbedr):
                    if occROOM[i+3][to-1] > 0 and occROOM[i+3][to] == 0:
                        window[i+3][to] = 1
                        k = 3
                    elif k > 0 and occ_m[0][to] == 1:
                        window[i+3][to] = 1
                        k += -1
                    else:
                        window[i+3][to] = 0
                        k = 0
            elif habit_bedr == 'OpenMorningLong':# 3 hours open after sleep
                for i in range(nbedr):
                    if occROOM[i+3][to-1] > 0 and occROOM[i+3][to] == 0:
                        window[i+3][to] = 1
                        k = 18
                    elif k > 0 and occ_m[0][to] == 1:
                        window[i+3][to] = 1
                        k += -1
                    else:
                        window[i+3][to] = 0
                        k = 0
            elif habit_bedr == 'OpenBeforeSleep':# 1 hour before sleep open
                for i in range(nbedr):
                    if to < nday*nbin-5:
                        if occROOM[i+3][to+5] > 0 and occROOM[i+3][to] == 0 and occ_m[0][to] == 1:
                            window[i+3][to] = 1
                        else:
                            window[i+3][to] = 0
                    else:
                        if occROOM[i+3][nday*nbin] > 0 and occROOM[i+3][to] == 0 and occ_m[0][to] == 1:
                            window[i+3][to] = 1
                        else:
                            window[i+3][to] = 0

            # Simulate window use in the bathroom when habits are present
            if habit_bath == 'Open':
                window[2][to] = 1
            elif habit_bath == 'Closed':
                window[2][to] = 0
            elif habit_bath == 'OpenPresent':
                if occROOM[2][to] == 0:
                    window[2][to] = 0
                else:
                    window[2][to] = 1
            elif habit_bath == 'ClosedPresent':
                if occROOM[2][to] == 0:
                    window[2][to] = 1
                else:
                    window[2][to] = 0
            elif habit_bath == 'OpenAfterPresent':
                if occROOM[2][to] == 0:
                    if occROOM[2][to-1] > 0:
                        window[2][to] = 1
                        l = 6
                    elif l > 0 and occ_m[0][to] == 1:
                        window[2][to] = 1
                        l += -1
                    else:
                        window[2][to] = 0
                        l = 0
                else:
                    window[2][to] = 0
                    l = 0
            elif habit_bath == 'OpenMorningShort':# 30 minutes open after sleep
                if occ_m[0][to-1] == 2 and occ_m[0][to] == 1:
                    window[2][to] = 1
                    l = 3
                elif l > 0 and occ_m[0][to] == 1:
                    window[2][to] = 1
                    l += -1
                else:
                    window[2][to] = 0
                    l = 0
            elif habit_bath == 'OpenMorningLong':# 3 hours open after sleep
                if occ_m[0][to-1] == 2 and occ_m[0][to] == 1:
                    window[2][to] = 1
                    l = 18
                elif l > 0 and occ_m[0][to] == 1:
                    window[2][to] = 1
                    l += -1
                else:
                    window[2][to] = 0
                    l = 0
            elif habit_bath == 'OpenNight':
                if occ_m[0][to] == 2:
                    window[2][to] = 1
                else:
                    window[2][to] = 0
                
            # Secondly, the behaviour when away and asleep are determined: OCCUPANCY-HABIT    
            # when all occupants are AWAY
            if occ_m[0][to] == 3:   
                if habit_away == 'Close':
                    for i in range(nroom):
                        window[i][to] = 0
                elif habit_away == 'CloseAcc':
                    window[0][to] = 0
                    window[1][to] = 0
                elif habit_away == 'CloseExclBed':
                    for i in range(3):
                        window[i][to] = 0
                elif habit_away == 'Open':
                    for i in range(nroom):
                        window[i][to] = 1
                elif habit_away == 'OpenBed':
                    for i in range(3,nroom):
                        window[i][to] = 1
                            
            # when all occupants are ASLEEP           
            if occ_m[0][to] == 2:
                if habit_sleep == 'Close' or habit_sleep == 'CloseExclBed':
                    for i in range(3):
                        window[i][to] = 0
                elif habit_sleep == 'Open':
                    for i in range(3):
                        window[i][to] = 1
                elif habit_sleep == 'CloseAcc':
                    window[0][to] = 0
                    window[1][to] = 0
            for i in range(3,nroom):
                if occROOM[i][to] > 0:
                    if habit_sleep == 'Close' or habit_sleep == 'CloseBed':
                        window[i][to] = 0
                    if habit_sleep == 'Open' or habit_sleep == 'OpenBed':
                        window[i][to] = 1
                
            # For the undefined periods (e.g. when no window habit is present)
            # a weather based window use model is applied (Maeyens & Janssens).
            for i in range(nroom):
                if window[i][to] == int(-1):
                    if occ_m[0][to] == 1:
                        if W_P[i][to] == 1:
                            window[i][to] = 1
                        else:
                            window[i][to] = 0
                    else:
                        if to == 0:
                            window[i][to] = 0
                        else:
                            window[i][to] = window [i][to-1]          

        # and store
        self.window = window

        return None                         

    def roundUp(self):
        '''
        Round the simulation by wrapping all data and reduce storage size.
        '''

        #######################################################################
        # first we move and sumarize data to the most upper level.
        nday = self.nday
        nbin = 144
        nbedr = self.nbedr
        nroom = nbedr + 3

        
        mDHW = np.zeros(nbin*nday + 1)
        tl = -1
        for to in range(nbin*nday):
            for i in range(6):
                tl += 1
                mDHW[to] += self.DHW[tl]

        sh_day = self.sh_settings['dayzone']
        sh_night = self.sh_settings['nightzone']
        sh_bath = self.sh_settings['bathroom']
        P = self.power
        for i in range(nroom):  
            P += self.lightload[i]
        
        QRad = np.zeros((11, nbin*nday + 1))
        QCon = np.zeros((11, nbin*nday + 1))

        Window = np.zeros((11, nbin*nday + 1))

        occ = np.zeros((11, nbin*nday + 1))
        CO2_prod = np.zeros((11, nbin*nday + 1))
        
        missing = -1*np.ones(nbin*nday + 1)
        for i in range(11):
            if i < nroom:
                QRad[i] += self.radi[i] + self.lightRad[i] + 0.5*self.MetHeat[i]/6
                QCon[i] += self.conv[i] + self.lightCon[i] + 0.5*self.MetHeat[i]/6
                occ[i] += self.occROOM[i]
                CO2_prod[i] += self.CO2[i]
                Window[i] = self.window[i]
            else:
                QRad[i] = missing
                QCon[i] = missing
                occ[i] = missing
                CO2_prod[i] = missing
                Window[i] = missing
                
        # we roll the outputs so the start-time is 0:00 instead of 4:00

        self.DHW = np.roll(mDHW,24)
        self.P = np.roll(P,24)
        self.sh_day = np.roll(sh_day,24)
        self.sh_night = np.roll(sh_night,24)
        self.sh_bath = np.roll(sh_bath,24)
        self.QRad = np.roll(QRad,24,axis=1)
        self.QCon = np.roll(QCon,24,axis=1)

        self.Window = np.roll(Window,24,axis=1)

        self.occ = np.roll(occ,24,axis=1)
        self.CO2 = np.roll(CO2_prod,24,axis=1)
        
        
        return None

    def pickle(self):
        '''
        Pickle the generated profile and its results for storing later.
        '''
        cPickle.dump(self, open(self.name+'.p','wb'))
        return None

class Equipment(object):
    '''
    Data records for appliance simulation based on generated activity and
    occupancy profiles
    '''
    # All object parameters are given in kwargs
    def __init__(self, **kwargs):
        # copy kwargs to object parameters
        for (key, value) in kwargs.items():
            setattr(self, key, value)
            
    def stochastic_flow(self, tap, nday, dow, occ, bathpres, pers, ACTtasks, apps):
        '''
        Simulate domestic hot water use based on data from Instal2020
        '''
        # determine duration distribution parameters
        cdir = os.getcwd()
        os.chdir('../Data/DHW')
        d = os.getcwd()
        if tap == "showerFlow":
            probShower = np.loadtxt(d + '/probShower.txt', np.float64)
            probFlow = probShower
        elif tap == "bathFlow":
            probBath = np.loadtxt(d + '/probBath.txt', np.float64)
            probFlow = probBath
        elif tap == "OtherFlow":
            probOther = np.loadtxt(d + '/probOther.txt', np.float64)
            probFlow = []
            probFlow.append(probOther)
        if pers < 5:
            persons = pers
        else:
            persons = 4
        
        # determine calibration factor
        cali = np.loadtxt(d + '/calibration.txt', float)
        calibration = []
        calibration.append(cali[persons-1])
        
        if tap == 'OtherFlow':
            cal = calibration[0][2]
        elif tap == 'showerFlow':
            cal = calibration[0][0]
        elif tap == 'dishFlow':
            cal = set_appliances['dishFlow']['cal']
        else:
            cal = calibration[0][1]

        nhour = 24
        nmin = 24*60*nday

        n_fl = 0
        flow = np.zeros(nmin+1)
        for p in range(pers):
            to = -1 # time counter for occupancy
            tl = -1
            left = -1 # time counter for appliance duration
            for doy, step, step2 in itertools.product(range(nday), range(nhour), range(6)):
                to += 1
                for i in range(10):
                    tl += 1
                    if left <= 0:
                        # person should be present in the bathroom
                        bpres = 1 if bathpres[p][to] == 1 else 0
                        occs = 1 if occ[0][to] == 1 else 0
                        dishwash = 1 if ACTtasks[1][p][to] == 1 and 'DishWasher' not in apps else 0
                        # determine probability of tapping in week and weekend
                        if tap == "OtherFlow":
                            prob = occs * probFlow[0][step]
                        elif tap == "dishFlow":
                            prob = dishwash
                        elif dow == 0 or dow == 6: 
                            prob = bpres * probFlow[1][step]
                        else:
                            prob = bpres * probFlow[0][step]
                        # check if there is a statechange in the appliance
                        if random.random() < prob * cal:
                            n_fl += 1
                            if tap == "showerFlow":
                                left = weibull_min.rvs(1.71,0.96,8.04)
                            elif tap == "OtherFlow": # correction needed as the time is in seconds,often less than 1 minute
                                left = genextreme.rvs(-1.12,2.39,2.16)/60
                            elif tap == "bathFlow":
                                left = random.gauss(self.cycle_length, self.cycle_length/2)
                            elif tap =="dishFlow":
                                left = loglaplace.rvs(1.32,-0.0042,3.00)
                                while left > 10:
                                    left = loglaplace.rvs(1.32,-0.0042,3.00)
                            flow[tl] += random.gauss(self.cycle_flow, self.cycle_flow/10)
                            left += -1
                    elif tap == "OtherFlow":
                        if occ[0][to] == 1:
                            flow[tl] += random.gauss(self.cycle_flow, self.cycle_flow/10)
                            left += -1
                        else:
                            left = -1
                    elif tap == "dishFlow":
                        if ACTtasks[1][p][to] == 1:
                            flow[tl] += random.gauss(self.cycle_flow, self.cycle_flow/10)
                            left += -1
                        else:
                            left = -1
                    else:
                        if bathpres[p][to] == 1:
                            flow[tl] += random.gauss(self.cycle_flow, self.cycle_flow/10)
                            left += -1
                        else:
                            left = -1

        r_fl = {'time':time, 'mDHW':flow}
        os.chdir(cdir)
        return r_fl, n_fl 
    
class Error(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg