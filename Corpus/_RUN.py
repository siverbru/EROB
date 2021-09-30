# -*- coding: utf-8 -*-
"""
Last update 20 April 2021

@author: Silke Verbruggen
"""
import os
import feeder

# give a name for the resultfile
name = 'Test'

# define simulation length (nday) and year (only used to define with which day of the year the simulation starts).
nday = 365
year = 2021

# define number of buildings to simulate
nBui = 1
# define HOUEHOLDMEMBERS.
# choice to not define anything, or to define just the number of persons, 
# or define the employment types, or even more detailed define occupancy types

# Number of persons:
# not predefined = 0
npers = 0
 
# Members choose from:
# FTE(fulltime empl.),PTE(parttime empl.),Unemployed,Retired,School(U18),Student
# e.g. ['FTE', 'FTE', 'Retired', 'School']
# if not predefined, leave empty []
members = []

# define CLUSTERS (occupancy profiles) for each weekday. 
# Choice from nightshift (1), early shift (2), short absence (+-3h) (3), dayshift (4), half day absence (5), long absence (6), always present (7) 
# e.g. [{'mon':4,'tue':4, 'wed':4, 'thu':4, 'fri':6, 'sat':3, 'sun':7}]
# if not predefined, leave empty []
clusters = []

# define BEDROOMS
# Choice to not define anything, define the number of bedrooms or 
# define the bedroom assignment

# define number of bedrooms (max = 8)
# not predefined = 0
nbedr = 0

# if both members and number of bedrooms are defined, option to define who sleeps together in one room
# start from roomnumber 0 with maximum of nbedr-1
# e.g. for household ['FTE','FTE','School'], [0,0,1] when the two FTE persons sleep together and the child sleeps in a seperate bedroom.
# if not predefined, leave empty []
BEDR = []

# define APPLIANCES. Choice from list defined in appliances.py
# if not predefined, leave empty []
apps = []


# define some building properties. if unknown = -1.
VentS = -1          # 0 = no system, 1 = exhaust ventilation, 2 = balanced ventilation
DW = -1             # 1 = apartment, 2 = house
YearBuilt = -1     # year of built/renovation 

# Indicators if the occupancy and actvity profiles are already defined on files (= 1), or not (= 0).
# Used for repeated simulations with exactly the same occupancy and appliance use.
OccONFILE = 0
ActONFILE = 0

# define heating profile, not predefined = -1
# define number that denotes the temperature to which is heated respectively when active, asleep and absent.
# 1={17,14,15}, 2={18.5,15,18.5}, 3={20,15,19.5}, 4={20,11,19.5}, 5={20,14.5,15}, 6={21,20.5,21}, 7={21.5,15.5,21.5}
shtype = -1
# define which zones are heated. e.g. ['dayzone','bathroom','nightzone'], unknown = []
shrooms = []

# define window use habits if known
# Either define the habits for each room and leaving and sleeping. unknown = []
# Order:
# living room winter: 0=no habit, 1=always closed, 2=open morning short, 3=open morning long, 4=open when present, 5=open evening, 6=closed when present, 7=always open
# living room summer: 0=no habit, 1=always closed, 2=open morning short, 3=open morning long, 4=open when present, 5=open night, 6=always open, 7=open evening, 8=closed when present
# kitchen winter: 0=no habit, 1=always closed, 2=open when cooking, 3=open after cooking, 4=open morning short, 5=open evening, 6=always open, 7=open when present, 8= open morning long
# kitchen summer: 0=no habit, 1=always closed, 2=open when cooking, 3=open after cooking, 4=open morning short, 5=open evening, 6=always open, 7=open when present, 8= open morning long, 9=open night
# bedroom winter: 0=no habit, 1=always closed, 2=open morning short, 3=open morning long, 4=always open, 5=open night, 6=open before sleep
# bedroom summer: 0=no habit, 1=always closed, 2=open morning short, 3=open morning long, 4=always open, 5=open night, 6=open before sleep, 7= open during the day
# bathroom winter: 0=no habit, 1=always closed, 2=open when present, 3=open after present, 4=open morning short, 5=always open, 6=closed when present, 7=open morning long
# bathroom summer: 0=no habit, 1=always closed, 2=open when present, 3=open after present, 4=open morning short, 5=always open, 6=closed when present, 7=open morning long, 8=open night
# leaving winter: 0=no change, 1=close all, 2=close accesible, 3=close all except bedroom, 4=open bedroom, 5= open all
# leaving summer: 0=no change, 1=close all, 2=close accesible, 3=close all except bedroom, 4=open bedroom, 5= open all
# sleeping winter: 0=no change, 1=close all, 2=close accesible, 3=close all except bedroom, 5= close bedroom, 6= open all
# sleeping summer: 0=no change, 1=close all, 2=close accesible, 3=close all except bedroom, 4= open bedroom, 5=open all, 6=close bedroom
habits = []
# or define household habit in winter and seasonality coherence. unknown = -1
HHhabit = -1        # 1=all closed, 2=all shortly opened, 3=day/night, 4=Bedr opened most, 5=short/long presence
SeCo = -1           # 1=same as in winter, 2=all opened more, 3=dayzone opened more, 4=nightzone opened more

##################### Run script
os.chdir('..')
directory = os.getcwd()
os.chdir(directory + '/Data')
d = os.getcwd()
test = feeder.IDEAS_Feeder(name, nBui, d, OccONFILE, ActONFILE, shtype, shrooms, members, apps, clusters, nbedr, BEDR, nday, year, npers, habits, HHhabit, SeCo, VentS, DW, YearBuilt)
