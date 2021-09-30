# -*- coding: utf-8 -*-
"""
Created on April 6 2021

@author: Silke Verbruggen

Code to calibrate the activity probabilities.
"""

import os
import feeder
import numpy as np
import csv
from Data.Appliances import set_appliances

os.chdir('..')
directory = os.getcwd()
dir1 = (directory + '/Data/Inputfiles')
dir2 = (directory + '/Data/HouseholdTasks')
dir3 = (directory + '/Data/IndividualActivities')
dir4 = (directory + '/Data/DHW')
os.chdir(directory + '/Data')
d = os.getcwd()

rep=5 # amount of repeated calibration rounds to check convergence

for j in range(rep):
    nday = 365
    year = 2021
    nBui = 100
    members = []
    clusters = []
    npers = 1
    nbedr = 0
    BEDR = []
    apps = []
    OccONFILE = 0
    ActONFILE = 0
    shtype = '-1'
    shrooms = []
    nstep = nday*144
    habits = []
    HHhabit = -1
    SeCo = -1
    VentS = -1
    DW = -1
    YearBuilt = -1

    FILE = open(directory +'/Data/DHW/DHW.txt',"a+")
    FILE.truncate(0)
    FILE.close

    test = feeder.IDEAS_Feeder('Test', nBui, d, OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt)

    TASKS = ['cook','dishes','vacuum','iron','wash','dry']
    INDACT= ['pc','tv','audio','adm','bath']

    # determine how long householdtask performed per household
    '''
    DurTasks = []
    averagedurationHHtask = np.loadtxt(dir2 +'\duration.txt', float)
    for i in range(6):
        taskduration =  []
        for n in range(nBui):
            taskduration.append(averagedurationHHtask[n][i])
        avduration = sum(taskduration)/len(taskduration)
        DurTasks.append(avduration)
    

    # determine how long individual activities performed per person
    DurAct = []
    averagedurationAct = np.loadtxt(dir3 +'\duration.txt', float)
    for i in range(5):
        actduration =  []
        for n in range(len(averagedurationAct)):
            actduration.append(averagedurationAct[n][i])
        aveduration = sum(actduration)/len(actduration)
        DurAct.append(aveduration)

    # determine how long each appliance is used
    DurAct = []
    averagedurationAct = np.loadtxt(dir3 +'\cookduration.csv', float)
    for i in range(4):
        actduration =  []
        for n in range(nBui):
            actduration.append(averagedurationAct[n][i])
        aveduration = sum(actduration)/len(actduration)
        DurAct.append(aveduration)
    '''
    # determine average number of cycles
    cycles = []
    averagecycles = np.loadtxt(dir4 +'\DHW.txt', float)
    for i in range(4):
        cycl =  []
        for n in range(nBui):
            if averagecycles[n][i] > 0:
                cycl.append(averagecycles[n][i])
        cycles.append(sum(cycl)/len(cycl))

    '''
    # compare number of cycles with reference: "convergence factors"

    convHHtask = np.loadtxt(dir2 +'\convergence.txt', float)
    conv = []
    for i in range(6):
        conv.append(DurTasks[i]/convHHtask[i])
    print (conv)
    

    convIndAct = np.loadtxt(dir3 +'\convergence.txt', float)
    conv2 = []
    for i in range(5):
        conv2.append(DurAct[i]/convIndAct[i])
    print (conv2)


    convCook = np.loadtxt(dir3 +'\cookconvergence.txt', float)
    conv3 = []
    for i in range(4):
        conv3.append(DurAct[i]/convCook[i])
    print (conv3)
    '''
    convDHW = np.loadtxt(dir4 +"\cycles.txt", float)
    conDHW = convDHW[0]
    conv3 = []
    for i in range(3):
        conv3.append(cycles[i]/conDHW[i])
    print (conv3)
    '''
    conv4 = cycles[3]/266
    print (conv4)

    #calculate new calibration factors

    cal = []
    calibrationfactor = np.loadtxt(dir2+'\calibrationHHtasks.txt', float)
    for i in range(6):
        cal.append(((calibrationfactor[i]/conv[i])+calibrationfactor[i])/2)
    print (cal)
    file= open(dir2+'\calibrationHHtasks.txt', "w")
    csv.writer(file, delimiter=' ',lineterminator = '\n').writerow(cal)
    file.close()    
    

    cal2 = []
    calibrationfactor2 = np.loadtxt(dir3+'\calibrationIndAct.txt', float)
    for i in range(5):
        cal2.append(((calibrationfactor2[i]/conv2[i])+calibrationfactor2[i])/2)
    print (cal2)
    file2= open(dir3+'\calibrationIndAct.txt', "w")
    csv.writer(file2, delimiter=' ',lineterminator = '\n').writerow(cal2)
    file2.close()

    k=-1
    for app in ['Oven','Microwave','Hob','Kettle']:
        k += 1
        cal = set_appliances[app]['cal']
        set_appliances[app]['cal'] = ((cal/conv3[k])+cal)/2
        print (set_appliances[app]['cal'])
            
    with open('Appliances.py', 'w') as file:
        file.write('set_appliances = \\\n' + str(set_appliances))
    
    '''
    cal3 = []
    calibrationfactor3 = np.loadtxt(dir4 +'\calibration.txt', float)
    for i in range(3):
        cal3.append(((calibrationfactor3[0][i]/conv3[i])+calibrationfactor3[0][i])/2)
    print(cal3)
    file3= open(dir4+'\calibration.txt', "w")
    csv.writer(file3, delimiter=' ',lineterminator = '\n').writerow(cal3)
    csv.writer(file3, delimiter=' ',lineterminator = '\n').writerow(cal3)
    csv.writer(file3, delimiter=' ',lineterminator = '\n').writerow(cal3)
    csv.writer(file3, delimiter=' ',lineterminator = '\n').writerow(cal3)
    file3.close()

    '''
    cal = set_appliances['dishFlow']['cal']
    set_appliances['dishFlow']['cal'] = ((cal/conv4)+cal)/2
            
    with open('Appliances.py', 'w') as file:
        file.write('set_appliances = \\\n' + str(set_appliances))
    '''