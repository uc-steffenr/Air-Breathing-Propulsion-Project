import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from plotData import *
from VelocityTriangle import VelocityTriangle


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
        self.To4 = 2150 # K

        self.mc = self.m/(1+self.BPR)
        self.mf = self.m - self.mc
        print('mc = {0:.2f} kg/s'.format(self.mc)) # core airflow
        # print('mf = {0:.2f} kg/s'.format(self.mf))
        # self.mc = 20

        self.eta_inf_c = 0.90
        self.eta_inf_f = 0.91
        self.pi_f = 1.6
        self.pi_c = 18
        self.eta_b = 0.98
        self.eta_inf_t = 0.93

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

    def plotData(self):
        ar = [self.alpha1_r,self.alpha2_r]
        am = [self.alpha1_m,self.alpha2_m]
        at = [self.alpha1_t,self.alpha2_t]
        br = [self.beta1_r,self.beta2_r]
        bm = [self.beta1_m,self.beta2_m]
        bt = [self.beta1_t,self.beta2_t]

        plotAirAngles(ar,am,at,br,bm,bt)
        plotPoRatio(self.poRatio)
        plotPo(self.po)
        plotTo(self.To)
        return

    #####################################################################
    #                         TURBINE METHODS                           #
    #####################################################################
    # from turbineSpeedandDimensions import turbineSpeedandDimensions
   
   
   
   # Solve Velocity Triangles
    def CompleteTriangles(self):        
        
        for i in range(len(self.alpha1_r)):
            C_a    = 227.5
            # Root Angles
            alp1r = self.alpha1_r[i]
            bet1r = self.beta1_r[i]
            alp2r = self.alpha2_r[i]
            bet2r = self.beta2_r[i]
            # Mean Angles
            alp1m = self.alpha1_m[i]
            bet1m = self.beta1_m[i]
            alp2m = self.alpha2_m[i]
            bet2m = self.beta2_m[i]
            # Tip Angles
            alp1t = self.alpha1_t[i]
            bet1t = self.beta1_t[i]
            alp2t = self.alpha2_t[i]
            bet2t = self.beta2_t[i]

            # Initialize different Triangles
            VelTri1_r = VelocityTriangle(alpha=alp1r, beta=bet1r, Ca=C_a)
            VelTri2_r = VelocityTriangle(alpha=alp2r, beta=bet2r, Ca=C_a)
            VelTri1_m = VelocityTriangle(alpha=alp1m, beta=bet1m, Ca=C_a)
            VelTri2_m = VelocityTriangle(alpha=alp2m, beta=bet2m, Ca=C_a)
            VelTri1_t = VelocityTriangle(alpha=alp1t, beta=bet1t, Ca=C_a)
            VelTri2_t = VelocityTriangle(alpha=alp2t, beta=bet2t, Ca=C_a)
            # Solve Everything
            '''
            VelTri1_r.Solve(alpha=alp1r, beta=bet1r, Ca=C_a)
            VelTri2_r.Solve(alpha=alp2r, beta=bet2r, Ca=C_a)
            VelTri1_m.Solve(alpha=alp1m, beta=bet1m, Ca=C_a)
            VelTri2_m.Solve(alpha=alp2m, beta=bet2m, Ca=C_a)
            VelTri1_t.Solve(alpha=alp1t, beta=bet1t, Ca=C_a)
            VelTri2_t.Solve(alpha=alp2t, beta=bet2t, Ca=C_a)
            '''
            # Append everything
            self.C1_r.append(VelTri1_r.C)
            self.V1_r.append(VelTri1_r.V)
            self.Vw1_r.append(VelTri1_r.Vw)
            self.C2_r.append(VelTri2_r.C)
            self.V2_r.append(VelTri2_r.V)
            self.Vw2_r.append(VelTri2_r.Vw)

            self.C1_m.append(VelTri1_m.C)
            self.V1_m.append(VelTri1_m.V)
            self.Vw1_m.append(VelTri1_m.Vw)
            self.C2_m.append(VelTri2_m.C)
            self.V2_m.append(VelTri2_m.V)
            self.Vw2_m.append(VelTri2_m.Vw)

            self.C1_t.append(VelTri1_t.C)
            self.V1_t.append(VelTri1_t.V)
            self.Vw1_t.append(VelTri1_t.Vw)
            self.C2_t.append(VelTri2_t.C)
            self.V2_t.append(VelTri2_t.V)
            self.Vw2_t.append(VelTri2_t.Vw)
        # End Loop

        # Putting in missing Stuff
        self.alpha3_r = [self.alpha1_r[1], self.alpha1_r[2], self.alpha1_r[3], self.alpha1_r[4], self.alpha1_r[5], self.alpha1_r[6], self.alpha1_r[7], 0]
        self.alpha3_m = [self.alpha1_m[1], self.alpha1_m[2], self.alpha1_m[3], self.alpha1_m[4], self.alpha1_m[5], self.alpha1_m[6], self.alpha1_m[7], 0]
        self.alpha3_t = [self.alpha1_t[1], self.alpha1_t[2], self.alpha1_t[3], self.alpha1_t[4], self.alpha1_t[5], self.alpha1_t[6], self.alpha1_t[7], 0]
            
        self.Cw3_r = [self.Cw1_r[1], self.Cw1_r[2], self.Cw1_r[3], self.Cw1_r[4], self.Cw1_r[5], self.Cw1_r[6], self.Cw1_r[7], 0]
        self.Cw3_m = [self.Cw1_m[1], self.Cw1_m[2], self.Cw1_m[3], self.Cw1_m[4], self.Cw1_m[5], self.Cw1_m[6], self.Cw1_m[7], 0]
        self.Cw3_t = [self.Cw1_t[1], self.Cw1_t[2], self.Cw1_t[3], self.Cw1_t[4], self.Cw1_t[5], self.Cw1_t[6], self.Cw1_t[7], 0]

        self.C3_r = [self.C1_r[1], self.C1_r[2], self.C1_r[3], self.C1_r[4], self.C1_r[5], self.C1_r[6], self.C1_r[7], C_a]
        self.C3_m = [self.C1_m[1], self.C1_m[2], self.C1_m[3], self.C1_m[4], self.C1_m[5], self.C1_m[6], self.C1_m[7], C_a]
        self.C3_t = [self.C1_t[1], self.C1_t[2], self.C1_t[3], self.C1_t[4], self.C1_t[5], self.C1_t[6], self.C1_t[7], C_a]


            
    
    def TabData(self):


        


        ##########################
        #       ROOT VALUES      #
        ##########################
        #----------Angles----------
        a1r = [i*180/np.pi for i in self.alpha1_r]
        a2r = [i*180/np.pi for i in self.alpha2_r]
        b1r = [i*180/np.pi for i in self.beta1_r]
        b2r = [i*180/np.pi for i in self.beta2_r]

        AngData_r = {"a1r": a1r, 'b1r': b1r, 'a2r': a2r, 'b2r': b2r}
        AngFrame_r = pd.DataFrame(data=AngData_r)

        #----------Velocities----------
        Vel1Data_r = {"Cw1_r": self.Cw1_r, "C1_r": self.C1_r, "Vw1_r":self.Vw1_r, "V1_r":self.V1_r}
        Vel2Data_r = {"Cw2_r": self.Cw2_r, "C2_r": self.C2_r, "Vw2_r":self.Vw2_r, "V2_r":self.V2_r}
        Vel3Data_r = {"Cw3_r": self.Cw3_r, "C3_r": self.C3_r}

        Vel1Frame_r = pd.DataFrame(data=Vel1Data_r)
        Vel2Frame_r = pd.DataFrame(data=Vel2Data_r)
        Vel3Frame_r = pd.DataFrame(data=Vel3Data_r)

        # Store as CSV:
        AngFrame_r.to_csv('AngFrame_r.csv')
        Vel1Frame_r.to_csv('Vel1Frame_r.csv')
        Vel2Frame_r.to_csv('Vel2Frame_r.csv')
        Vel3Frame_r.to_csv('Vel3Frame_r.csv')


        ##########################
        #       MEAN VALUES      #
        ##########################
        #----------Angles----------
        a1m = [i*180/np.pi for i in self.alpha1_m]
        a2m = [i*180/np.pi for i in self.alpha2_m]
        b1m = [i*180/np.pi for i in self.beta1_m]
        b2m = [i*180/np.pi for i in self.beta2_m]

        AngData_m = {"a1m": a1m, 'b1m': b1m, 'a2m': a2m, 'b2m': b2m}
        AngFrame_m = pd.DataFrame(data=AngData_m)
       
        
        #----------Velocities----------
        Vel1Data_m = {"Cw1_m": self.Cw1_m, "C1_m": self.C1_m, "Vw1_m":self.Vw1_m, "V1_m":self.V1_m}
        Vel2Data_m = {"Cw2_m": self.Cw2_m, "C2_m": self.C2_m, "Vw2_m":self.Vw2_m, "V2_m":self.V2_m}
        Vel3Data_m = {"Cw3_m": self.Cw3_m, "C3_m": self.C3_m}

        Vel1Frame_m = pd.DataFrame(data=Vel1Data_m)
        Vel2Frame_m = pd.DataFrame(data=Vel2Data_m)
        Vel3Frame_m = pd.DataFrame(data=Vel3Data_m)

        # Store as CSV:
        AngFrame_m.to_csv('AngFrame_m.csv')
        Vel1Frame_m.to_csv('Vel1Frame_m.csv')
        Vel2Frame_m.to_csv('Vel2Frame_m.csv')
        Vel3Frame_m.to_csv('Vel3Frame_m.csv')

        ##########################
        #       TIP VALUES       #
        ##########################
        #----------Angles----------
        a1t = [i*180/np.pi for i in self.alpha1_t]
        a2t = [i*180/np.pi for i in self.alpha2_t]
        b1t = [i*180/np.pi for i in self.beta1_t]
        b2t = [i*180/np.pi for i in self.beta2_t]

        AngData_t = {"a1t": a1t, 'b1t': b1t, 'a2t': a2t, 'b2t': b2t}
        AngFrame_t = pd.DataFrame(data=AngData_t)

        #----------Velocities----------
        Vel1Data_t = {"Cw1_t": self.Cw1_t, "C1_t": self.C1_t, "Vw1_t":self.Vw1_t, "V1_t":self.V1_t}
        Vel2Data_t = {"Cw2_t": self.Cw2_t, "C2_t": self.C2_t, "Vw1_t":self.Vw2_t, "V1_t":self.V2_t}
        Vel3Data_t = {"Cw3_t": self.Cw3_t, "C3_t": self.C3_t}

        Vel1Frame_t = pd.DataFrame(data=Vel1Data_t)
        Vel2Frame_t = pd.DataFrame(data=Vel2Data_t)
        Vel3Frame_t = pd.DataFrame(data=Vel3Data_t)

        # Store as CSV:
        AngFrame_t.to_csv('AngFrame_t.csv')
        Vel1Frame_t.to_csv('Vel1Frame_t.csv')
        Vel2Frame_t.to_csv('Vel2Frame_t.csv')
        Vel3Frame_t.to_csv('Vel3Frame_t.csv')
        