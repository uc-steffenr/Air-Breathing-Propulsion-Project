import numpy as np
import matplotlib.pyplot as plt
from sys import exit

# NOTE might only need stage 1 and 2 values

class EngineDesign:
    def __init__(self,**kwargs):
        #########################################
        #            GIVEN VALUES               #
        #########################################
        self.BPR = 7
        self.m = 520 # kg/s
        self.p0 = 1.01 # bar
        self.T0 = 288 # K

        self.mc = self.m/(1+self.BPR)
        self.mf = self.m - self.mc
        print('mc = {0:.2f} kg/s'.format(self.mc)) # core airflow
        # print('mf = {0:.2f} kg/s'.format(self.mf))
        # self.mc = 20

        self.eta_inf_c = 0.90
        self.eta_inf_f = 0.91
        self.pi_f = 1.6
        self.pi_c = 18

        self.deHaller = 0.65

        #########################################
        #           NECESSARY CONSTANTS         #
        #########################################
        self.gamma_c = 1.4
        self.gamma_h = 1.333
        self.R = 287 # J/(kg-K)
        self.cpa = 1005 # J/(kg-K)
        self.cpg = 1148 # J/(kg-K)     

        ##########################
        #   ROOT RADIUS VALUES   #
        ##########################
        self.rr = []
        self.Cw1_r = []
        self.C1_r = []
        self.V1_r = []
        self.Vw1_r = []
        self.Cw2_r = []
        self.C2_r = []
        self.V2_r = []
        self.Vw2_r = []
        self.Cw3_r = []
        self.C3_r = []

        self.alpha1_r = []
        self.alpha2_r = []
        self.alpha3_r = []
        self.beta1_r = []
        self.beta2_r = []

        ##########################
        #   MEAN RADIUS VALUES   #
        ##########################
        self.Cw1_m = []
        self.C1_m = []
        self.V1_m = []
        self.Vw1_m = []
        self.Cw2_m = []
        self.C2_m = []
        self.V2_m = []
        self.Vw2_m = []
        self.Cw3_m = []
        self.C3_m = []

        self.alpha1_m = []
        self.beta1_m = []
        self.alpha2_m = []
        self.beta2_m = []
        self.alpha3_m = []

        ##########################
        #    TIP RADIUS VALUES   #
        ##########################
        self.Cw1_t = []
        self.C1_t = []
        self.V1_t = []
        self.Vw1_t = []
        self.Cw2_t = []
        self.C2_t = []
        self.V2_t = []
        self.Vw2_t = []
        self.Cw3_t = []
        self.C3_t = []

        self.alpha1_t = []
        self.alpha2_t = []
        self.alpha3_t = []
        self.beta1_t = []
        self.beta2_t = []    

        #########################################
        #            EXTRA ARGUMENTS            #
        #########################################
        self.showValues = False
        self.machTipExit = False
        self.deHallerExit = False
        self.constant = 'mean' # else tip

        for key,val in kwargs.items():
            if key == 'showValues':
                self.showValues = val
            elif key == 'deHallerExit':
                self.deHallerExit = val
            elif key == 'machTipExit':
                self.machTipExit = val
            elif key == 'constant':
                self.constant = val
        
        # decides how to store non-constant value
        if self.constant == 'mean':
            self.rt = []
        else:
            self.rm = []




    #####################################################################
    #                       COMPRESSOR METHODS                          #
    #####################################################################
    from compressorSpeedandDimensions import compressorSpeedandDimensions    
    from compressorStageEstimation import compressorStageEstimation    
    from compressorStageDesign import compressorStageDesign
    from compressorAirAngles import compressorAirAngles
   