# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:48:51 2021

@author: siverbru
"""
import math
import numpy as np

'''
This model calculates which window use habits are present for each household in the simulation.
Based on the work presented in (REF)
'''
def get_habits(VentS, DW, Year, members, HHhabit, SeCo, **kwargs):
    # load inputdata if available
    if VentS == -1:
        probVS = [0.661,0.167,0.171]
        cumprobVS = np.cumsum(probVS)
        rnd = np.random.random()
        idz = 0
        while rnd >= cumprobVS[idz]:
            idz += 1
        VentS = idz
    Vent3_no = 1 if VentS == 0 else 0
    Vent3_C = 1 if VentS == 1 else 0
    if DW == -1:
        if np.random.random() < 0.204:
            DW = 1
        else:
            DW = 0
    DW2 = 1 if DW == 1 else 0
    YEAR_unk = 1 if Year == -1 else 0
    YEAR_1 = 1 if Year <= 1950 else 0
    YEAR_2 = 1 if Year > 1950 and Year <= 1975 else 0
    YEAR_3 = 1 if Year > 1975 and Year <= 2005 else 0
    YEAR_4 = 1 if Year > 2005 and Year <= 2015 else 0
    RET = 1 if 'Retired' in members else 0
    STU = 1 if 'Student' in members else 0
    UNE = 1 if 'Unemployed' in members else 0
    Fam_s = 1 if len(members) == 1 else 0
    Fam_c = 1 if len(members) == 2 and 'School' not in members else 0
    CHILD = 1 if 'School' in members else 0

    # step 1: determine which household winter habit each household has.
    if HHhabit != -1:
        HHtot = HHhabit
    else:
        logHH_O = -2.589 + 3.009*Vent3_no + 2.278*Vent3_C + 0.733*YEAR_unk + 0.809*YEAR_1 + 0.625*YEAR_2 + 0.770*YEAR_3 + 0.882*YEAR_4 - 1.108*Fam_s - 0.440*Fam_c + 2.062*STU + 1.370*UNE - 1.040*CHILD
        logHH_DN = -0.636 + 2.191*Vent3_no + 0.958*Vent3_C - 1.087*YEAR_unk - 0.302*YEAR_1 + 0.060*YEAR_2 + 0.291*YEAR_3 + 0.194*YEAR_4 - 0.495*Fam_s - 0.610*Fam_c + 0.717*STU + 1.253*UNE - 1.486*CHILD
        logHH_BED = 0.484 + 1.355*Vent3_no + 0.229*Vent3_C + 0.043*YEAR_unk - 0.533*YEAR_1 - 0.566*YEAR_2 + 0.960*YEAR_3 - 0.164*YEAR_4 - 1.848*Fam_s - 1.153*Fam_c + 1.409*STU + 0.598*UNE - 1.478*CHILD
        logHH_SL = -1.759 + 0.462*Vent3_no + 0.666*Vent3_C - 0.971*YEAR_unk - 0.763*YEAR_1 - 0.137*YEAR_2 - 0.748*YEAR_3 - 0.951*YEAR_4 + 1.422*Fam_s + 0.937*Fam_c + 2.828*STU + 3.017*UNE + 0.252*CHILD
            
        total = 1 + math.exp(logHH_O) + math.exp(logHH_DN) + math.exp(logHH_BED) + math.exp(logHH_SL)
        P_HH_C = 1/total
        P_HH_O = math.exp(logHH_O)/total
        P_HH_DN = math.exp(logHH_DN)/total
        P_HH_BED = math.exp(logHH_BED)/total
        P_HH_SL = math.exp(logHH_SL)/total
        probHH = [P_HH_C, P_HH_O, P_HH_DN, P_HH_BED, P_HH_SL]
            
        cumprobHH = np.cumsum(probHH)
        rnd = np.random.random()
        idx = 1
        while rnd >= cumprobHH[idx-1]:
            idx += 1
        HHtot = idx
        
    HH_C = HH_O = HH_DN = HH_BED = HH_SL = 0
    if HHtot == 1:
        HH_C = 1
    elif HHtot == 2:
        HH_O = 1
    elif HHtot == 3:
        HH_DN = 1
    elif HHtot == 4:
        HH_BED = 1
    else:
        HH_SL = 1

    # step 2: determine based on these householdhabits the window use habits in winter.
    if HH_C == 1:
        LIWtot = KIWtot = BEWtot = 1
        probBAW = [0.922,0.013,0.065]
        cumprobBAW = np.cumsum(probBAW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobBAW[idz-1]:
            idz += 1
        BAWtot = idz
            
    elif HH_O == 1:
        probLIW = [0.2820, 0.5490, 0.0985, 0.0705]
        cumprobLIW = np.cumsum(probLIW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobLIW[idz-1]:
            idz += 1
        if idz == 1:
            LIWtot = 0
        else:
            LIWtot = idz
                
        probKIW = [0.321, 0.143, 0.179, 0.054, 0.304]
        cumprobKIW = np.cumsum(probKIW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobKIW[idz-1]:
            idz += 1
        KIWtot = idz-1
                
        probBAW = [0.200, 0.020, 0.080, 0.560, 0.140]
        cumprobBAW = np.cumsum(probBAW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobBAW[idz-1]:
            idz += 1
        BAWtot = idz-1
                
        probBEW = [0.074, 0.642, 0.198, 0.025, 0.062]
        cumprobBEW = np.cumsum(probBEW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobBEW[idz-1]:
            idz += 1
        if idz == 1:
            BEWtot = 0
        elif idz == 5:
            BEWtot = 6
        else:
            BEWtot = idz
                
    elif HH_DN ==1:
        # function to determine day vs night habits
        # determine dayzone habits (liv = Kit) and nightzone habits (BED = BATH)
        probLIW = [0.043, 0.857, 0.100]
        cumprobLIW = np.cumsum(probLIW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobLIW[idz-1]:
            idz += 1
        if idz <= 2:
            LIWtot = idz
            KIWtot = idz
            
            probBEW = [0.017, 0.085, 0.119, 0.779]
            cumprobBEW = np.cumsum(probBEW)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBEW[idx-1]:
                idx += 1
            if idx == 1:
                BEWtot = BAWtot = idx - 1
            elif idx <= 3:
                BEWtot = idx + 2
                BAWtot = 5
            else:
                probBEW_sh = [0.544, 0.282, 0.174]
                cumprobBEW_sh = np.cumsum(probBEW_sh)
                rnd = np.random.random()
                idy = 1
                while rnd >= cumprobBEW_sh[idy-1]:
                    idy += 1
                if idy == 3:
                    BEWtot = 6
                else: 
                    BEWtot = idy + 1
                probBAW_sh = [0.151, 0.818, 0.031]
                cumprobBAW_sh = np.cumsum(probBAW_sh)
                rnd = np.random.random()
                idy = 1
                while rnd >= cumprobBAW_sh[idy-1]:
                    idy += 1
                BAWtot = idy + 1
        else:
            probLIW_sh = [0.14, 0.86]
            cumprobLIW_sh = np.cumsum(probLIW_sh)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIW_sh[idx-1]:
                idx += 1
            if idx == 1:
                LIWtot = 4
            else:
                LIWtot = 2
            probKIW_sh = [0.754, 0.246]
            cumprobKIW_sh = np.cumsum(probKIW_sh)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIW_sh[idx-1]:
                idx += 1
            if idx == 1:
                KIWtot = 2
            else:
                KIWtot = 4
            probBEW = [0.428, 0.572]
            cumprobBEW = np.cumsum(probBEW)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBEW[idx-1]:
                idx += 1
            BEWtot = idx + 3
            BAWtot = 5 
                
    elif HH_SL ==1:
        # function to determine long (LI/BE) vs short (KI/BA) presence habits
        # determine living and bedroom habits
        probBEW = [0.086, 0.343, 0.343, 0.114, 0.029, 0.086]
        cumprobBEW = np.cumsum(probBEW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobBEW[idz-1]:
            idz += 1
        if idz <= 2:
            LIWtot = idz - 1
            BEWtot = idz - 1
        elif idz == 5:
            BEWtot = 4
            LIWtot = 7
        elif idz == 6:
            BEWtot = 5
            LIWtot = 2
        elif idz == 4:
            BEWtot = 3
            probLIW = [0.33, 0.67]
            cumprobLIW = np.cumsum(probLIW)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIW[idx-1]:
                idx += 1
            LIWtot = idx + 1
        else:
            BEWtot = 2
            probLIW = [0.75, 0.125, 0.125]
            cumprobLIW = np.cumsum(probLIW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobLIW[idz-1]:
                idz += 1
            LIWtot = idz + 1
        # based on the living and bedroom habits determine habits in kitchen and bathroom
        if BEWtot == 0:
            KIWtot = 1
            BAWtot = 1
        elif BEWtot == 1:
            probKIWBAW = [0.25, 0.25, 0.5]
            cumprobKIWBAW = np.cumsum(probKIWBAW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobKIWBAW[idz-1]:
                idz += 1
            if idz == 1:
                KIWtot = 0
                BAWtot = 0
            else:
                KIWtot = 2
                BAWtot = idz 
        else:
            probKIWBAW = [0.125, 0.875]
            cumprobKIWBAW = np.cumsum(probKIWBAW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobKIWBAW[idz-1]:
                idz += 1
            KIWtot = idz - 1
            BAWtot = idz - 1

    else:
        # function to determine 'bed most'-habits
        # determine habit in bedroom
        logBEW_M = 3.689 - 2.303*DW2
        logBEW_O = 3.219 - 1.427*DW2
        logBEW_ON = 2.890 - 2.890*DW2
        logBEW_BS = 0.693 + 0.405*DW2
        
        total = 1 + math.exp(logBEW_M) + math.exp(logBEW_O) + math.exp(logBEW_ON) + math.exp(logBEW_BS)
        P_BEW_NH = 1/total
        P_BEW_M = math.exp(logBEW_M)/total
        P_BEW_O = math.exp(logBEW_O)/total
        P_BEW_ON = math.exp(logBEW_ON)/total
        P_BEW_BS = math.exp(logBEW_BS)/total
        probBEW = [P_BEW_NH, P_BEW_M, P_BEW_O, P_BEW_ON, P_BEW_BS]
        cumprobBEW = np.cumsum(probBEW)
        rnd = np.random.random()
        idx = 1
        while rnd >= cumprobBEW[idx-1]:
            idx += 1
        if idx == 1:
            BEWtot = 0
        elif idx == 2:
            probBEW_m = [0.795,0.205]
            cumprobBEW_m = np.cumsum(probBEW_m)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobBEW_m[idz-1]:
                idz += 1
            BEWtot = 1+idz
        else:
            BEWtot = idx+1
        # determine based on bedroom habit the habit in the living room
        if BEWtot in [0,2,3,6]:
            LIWtot = 1
        elif BEWtot == 4:
            probLIW = [0.862, 0.138]
            cumprobLIW = np.cumsum(probLIW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobLIW[idz-1]:
                idz += 1
            LIWtot = idz
        elif BEWtot == 5:
            probLIW = [0.158, 0.737, 0.0525, 0.0525]
            cumprobLIW = np.cumsum(probLIW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobLIW[idz-1]:
                idz += 1
            if idz == 4:
                LIWtot = 4
            else:
                LIWtot = idz - 1
        # habits in Kitchen and Bathroom based on probabilities in this group.
        probKIW = [0.084, 0.614, 0.169, 0.096, 0.037]
        cumprobKIW = np.cumsum(probKIW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobKIW[idz-1]:
            idz += 1
        KIWtot = idz - 1
        
        probBAW = [0.043, 0.574, 0.053, 0.298, 0.021, 0.011]
        cumprobBAW = np.cumsum(probBAW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobBAW[idz-1]:
            idz += 1
        BAWtot = idz - 1

            
    # determine the sleeping habits based on the household and room habits
    if HH_C == 1:
        probSLW = [0.3, 0.7]
        cumprobSLW = np.cumsum(probSLW)
        rnd = np.random.random()
        idz = 1
        while rnd >= cumprobSLW[idz-1]:
            idz += 1
        SLWtot = idz - 1
        
    elif HH_O == 1:
        if BEWtot == 4:
            SLWtot = 0
        else:
            probSLW = [0.169, 0.597, 0.091, 0.143]
            cumprobSLW = np.cumsum(probSLW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobSLW[idz-1]:
                idz += 1
            SLWtot = idz - 1
    else:
        if BEWtot in [4,5]:
            SLWtot = 3
        else:
            probSLW = [0.115, 0.755, 0.036, 0.094]
            cumprobSLW = np.cumsum(probSLW)
            rnd = np.random.random()
            idz = 1
            while rnd >= cumprobSLW[idz-1]:
                idz += 1
            SLWtot = idz - 1
        
    SLW_C = SLW_nc = SLW_ca = SLW_ceb = 0
    if SLWtot == 0:
        SLW_nc = 1
    elif SLWtot == 1:
        SLW_C = 1
    elif SLWtot == 2:
        SLW_ca = 1
    elif SLWtot == 3:
        SLW_ceb = 1
    
    # determine the leaving habits based on the household and sleeping habits 
    logLEW_nc = -0.799 + 1.601*HH_C - 0.655*HH_O + 0.449*HH_DN + 0.867*HH_BED + 1.994*SLW_nc + 0.024*SLW_C - 2.612*SLW_ca
    logLEW_C = 0.497 + 1.594*HH_C - 0.778*HH_O + 0.171*HH_DN + 0.327*HH_BED - 0.571*SLW_nc + 0.947*SLW_C - 2.973*SLW_ca
    total = 1 + math.exp(logLEW_nc) + math.exp(logLEW_C)
    P_LEW_ca = 1/total
    P_LEW_nc = math.exp(logLEW_nc)/total
    P_LEW_C = math.exp(logLEW_C)/total
    probLEW = [P_LEW_nc, P_LEW_C, P_LEW_ca]
    cumprobLEW = np.cumsum(probLEW)
    rnd = np.random.random()
    idx = 1
    while rnd >= cumprobLEW[idx-1]:
        idx += 1
    LEWtot = idx - 1
    
    LEW_nc = LEW_C = 0
    if LEWtot == 0:
        LEW_nc = 1
    elif LEWtot == 1:
        LEW_C = 1
        
    # Step 3: determine seasonality coherence based on household habits and employment
    if SeCo != -1:
        SCtot = SeCo
    else:
        logSC_O = 0.821 + 1.049*HH_C - 0.798*HH_O + 0.354*HH_DN + 0.207*HH_BED - 1.032*RET
        logSC_D = -1.288 + 0.415*HH_C + 0.432*HH_O + 2.392*HH_DN + 2.698*HH_BED - 0.046*RET
        logSC_N = -2.506 + 1.787*HH_C + 2.201*HH_O + 1.443*HH_DN + 2.186*HH_BED + 0.377*RET
            
        total = 1 + math.exp(logSC_O) + math.exp(logSC_D) + math.exp(logSC_N)
        P_SC_S = 1/total
        P_SC_O = math.exp(logSC_O)/total
        P_SC_D = math.exp(logSC_D)/total
        P_SC_N = math.exp(logSC_N)/total
        probSC = [P_SC_S, P_SC_O, P_SC_D, P_SC_N]
        cumprobSC = np.cumsum(probSC)
        rnd = np.random.random()
        idx = 1
        while rnd >= cumprobSC[idx-1]:
            idx += 1
        SCtot = idx
            
    # step 4: determine summer habits based on winter habits and seasonality coherence
    if SCtot == 1:
        LIStot = LIWtot
        KIStot = KIWtot
        BEStot = BEWtot
        BAStot = BAWtot
    elif SCtot == 2:
        # opened more --> each winter habited increased (prob)
        if BAWtot == 0:
            probBAS = [0.454, 0.091, 0, 0.091, 0, 0.364]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx - 1
        elif BAWtot == 1:
            probBAS = [0.123, 0.383, 0.110, 0.137, 0.020, 0.192, 0, 0.035]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx - 1
        elif BAWtot == 2:
            probBAS = [0.667, 0, 0, 0.333]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx + 1
        elif BAWtot == 3:
            probBAS = [0.139, 0.194, 0.020, 0.611, 0, 0.036]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx + 1
        else:
            probBAS = [0.089, 0.750, 0, 0.161]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx + 3
                
        if BEWtot == 0:
            probBES = [0.2, 0.4, 0.4]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx + 2
        elif BEWtot == 1:
            probBES = [0.100, 0, 0.100, 0.133, 0.233, 0.389, 0.045]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx - 1
        elif BEWtot == 2:
            probBES = [0.024, 0.048, 0.738, 0.190]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx + 1
        elif BEWtot == 3:
            probBES = [0.055, 0, 0.833, 0.112]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx + 1
        elif BEWtot == 4:
            probBES = [0.143, 0.857]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx + 3
        elif BEWtot == 5:
            BEStot = 4
        else:
            probBES = [0.667, 0.333]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx + 3
                
        if LIWtot == 0:
            probLIS = [0.091, 0.364, 0.545]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx + 2 
        elif LIWtot == 1:
            probLIS = [0.113, 0, 0.143, 0.038, 0.421, 0.113, 0.173]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx - 1
        elif LIWtot == 2:
            probLIS = [0.154, 0, 0, 0, 0.077, 0.231, 0.538]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx - 1
        elif LIWtot == 3:
            probLIS = [0.5, 0, 0, 0, 0, 0.5]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx - 1 
        else:
            LIStot = 6
            
        if KIWtot == 0:
            probKIS = [0.454, 0, 0, 0, 0.182, 0, 0.364]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        elif KIWtot == 1:
            probKIS = [0.147, 0.137, 0.137, 0.000, 0.118, 0.000, 0.402, 0, 0, 0.059]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        elif KIWtot == 2:
            probKIS = [0.077, 0.077, 0.077, 0, 0, 0, 0.769]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        elif KIWtot == 3:
            probKIS = [0.25, 0, 0, 0, 0.75]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx + 1
        else:
            probKIS = [0.25, 0 ,0, 0, 0.25, 0, 0.5]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
            
    elif SCtot == 3:
        # dayzone opened more
        BEStot = BEWtot
        BAStot = BAWtot
        if LIWtot == 0:
            probLIS = [0.8, 0.2]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx + 3 
        elif LIWtot == 1:
            probLIS = [0.243, 0, 0.077, 0.026, 0.397, 0.090, 0.167]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx - 1
        elif LIWtot == 2:
            probLIS = [0.167, 0, 0, 0, 0, 0, 0.833]
            cumprobLIS = np.cumsum(probLIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobLIS[idx-1]:
                idx += 1
            LIStot = idx - 1
        elif LIWtot == 3:
            LIStot = 6
        else:
            LIStot = 5
            
        if KIWtot == 0:
            probKIS = [0.167, 0, 0, 0, 0, 0, 0.667, 0, 0, 0.167]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        elif KIWtot == 1:
            probKIS = [0.278, 0.111, 0.130, 0.000, 0.222, 0.000, 0.241, 0, 0, 0.018]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        elif KIWtot == 2:
            probKIS = [0.167, 0, 0.5, 0, 0, 0, 0.333]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        elif KIWtot == 3:
            probKIS = [0.167, 0, 0.667, 0, 0, 0, 0, 0.167]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx - 1
        else:
            probKIS = [0.143, 0 , 0.857]
            cumprobKIS = np.cumsum(probKIS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobKIS[idx-1]:
                idx += 1
            KIStot = idx + 3
    else:
        # nightzone opened more
        LIStot = LIWtot
        KIStot = KIWtot
        if BEWtot == 0:
            BEStot = 3
        elif BEWtot == 1:
            probBES = [0.167, 0, 0.167, 0, 0.500, 0.167]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx - 1
        elif BEWtot == 2:
            probBES = [0.05, 0, 0, 0, 0.600, 0.350]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx - 1
        elif BEWtot in [3,5]:
            BEStot = 4
        elif BEWtot == 4:
            probBES = [0.5, 0, 0, 0, 0.5]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx
        else:
            probBES = [0.333, 0.667]
            cumprobBES = np.cumsum(probBES)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBES[idx-1]:
                idx += 1
            BEStot = idx + 3
            
        if BAWtot in [0,5]:
            BAStot = 5
        elif BAWtot == 1:
            probBAS = [0.555, 0, 0.333,0,0.112]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx
        elif BAWtot == 2:
            probBAS = [0.75, 0.25]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            if idx == 1:
                BAStot = 2
            else:
                BAStot = 5
        elif BAWtot == 3:
            probBAS = [0.308, 0.231, 0.028, 0.385, 0, 0.050]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx + 1
        else:
            probBAS = [0.357, 0, 0, 0.643]
            cumprobBAS = np.cumsum(probBAS)
            rnd = np.random.random()
            idx = 1
            while rnd >= cumprobBAS[idx-1]:
                idx += 1
            BAStot = idx + 3
            
    BES_C = BES_O = BES_ON = BES_ml = 0
    if BEStot == 1:
        BES_C = 1
    elif BEStot == 4:
        BES_O = 1
    elif BEStot == 5:
        BES_ON = 1
    elif BEStot == 3:
        BES_ml = 1
        
    # determine sleeping habits based on seasonality coherence and bedroom habits
    logSLS_C = 1.618 - 1.912*SLW_nc - 2.296*SLW_ceb - 0.249*BES_C - 2.118*BES_O - 2.003*BES_ON
    logSLS_ca = 1.429 - 1.954*SLW_nc - 0.733*SLW_ceb - 2.509*BES_C - 0.540*BES_O - 0.774*BES_ON
    logSLS_ceb = 1.398 - 1.833*SLW_nc + 0.794*SLW_ceb - 2.710*BES_C - 0.162*BES_O - 0.109*BES_ON
    logSLS_o = -1.613 - 0.001*SLW_nc + 0.070*SLW_ceb - 0.002*BES_C - 0.150*BES_O + 1.743*BES_ON
    total = 1 + math.exp(logSLS_C) + math.exp(logSLS_ca)+ math.exp(logSLS_ceb)+ math.exp(logSLS_o)
    P_SLS_nc = 1/total
    P_SLS_C = math.exp(logSLS_C)/total
    P_SLS_ca = math.exp(logSLS_ca)/total
    P_SLS_ceb = math.exp(logSLS_ceb)/total
    P_SLS_o = math.exp(logSLS_o)/total
    probSLS = [P_SLS_nc, P_SLS_C, P_SLS_ca, P_SLS_ceb, 0, P_SLS_o]
    cumprobSLS = np.cumsum(probSLS)
    rnd = np.random.random()
    idx = 1
    while rnd >= cumprobSLS[idx-1]:
        idx += 1
    SLStot = idx - 1
    
    SLS_nc = SLS_C = SLS_ca = SLS_o = 0
    if SLStot == 0:
        SLS_nc = 1
    elif SLStot == 1:
        SLS_C = 1
    elif SLStot == 2:
        SLS_ca = 1
    elif SLStot == 5:
        SLS_o = 1
        
    # dtermine leaving habits in summer based on seasonality coherence, 
    # bedroom habits and leaving habits in winter
    logLES_C = 1.540 - 2.095*SLS_nc + 1.271*SLS_C + 0.973*SLS_ca - 0.819*SLS_o - 1.644*LEW_nc + 1.114*LEW_C - 1.268*BES_ml - 1.644*BES_O
    logLES_ca = 2.756 - 1.390*SLS_nc + 0.143*SLS_C + 2.207*SLS_ca - 2.273*SLS_o - 2.830*LEW_nc - 1.381*LEW_C - 0.198*BES_ml - 0.387*BES_O
    total = 1 + math.exp(logLES_C) + math.exp(logLES_ca)
    P_LES_nc = 1/total
    P_LES_C = math.exp(logLES_C)/total
    P_LES_ca = math.exp(logLES_ca)/total
    probLES = [P_LES_nc, P_LES_C, P_LES_ca]
    cumprobLES = np.cumsum(probLES)
    rnd = np.random.random()
    idx = 1
    while rnd >= cumprobLES[idx-1]:
        idx += 1
    LEStot = idx - 1 
                
    habits = [LIWtot, LIStot, KIWtot, KIStot, BEWtot, BEStot, BAWtot, BAStot, LEWtot, LEStot, SLWtot, SLStot]   

    return habits
    

