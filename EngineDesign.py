import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from plotData import *

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
        self.eta_m = 0.99

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
    
    def compressorStageDesign(self,numStages):
        self.alpha1 = []
        self.beta1 = []
        self.Cw1 = []
        self.Vw1 = []
        self.V1 = []
        self.C1 = []

        self.alpha2 = []
        self.beta2 = []
        self.Cw2 = []
        self.Vw2 = []
        self.V2 = []
        self.C2 = []

        self.alpha3 = []
        self.beta3 = [] 
        self.Cw3 = []
        self.Vw3 = []
        self.C3 = []
        self.V3 = []

        self.dToS = (self.To3-self.To13)/numStages
        # need to estimate change in temperature for every stage based on dToS
        dT_first_last = self.dToS*0.85 # kinda guestimate for first and last delta T
        dT = self.dToS*1.1

        print('dToS = {0:.2f} K'.format(self.dToS))
        print('First and last dToS = {0:.2f} K'.format(dT_first_last))
        print('Other stage dToS - {0:.2f} K'.format(dT))

        # work done factors
        lam_1 = 0.98
        lam_2 = 0.93
        lam_3 = 0.88
        lam = 0.84

        # for stages 1 and 2
        dCw = Cw2 = (self.cpa*dT_first_last)/(lam_1*self.Um)
        a1 = 0
        b1 = np.arctan(self.Um/self.Ca)
        b2 = np.arctan((self.Um-Cw2)/self.Ca)
        a2 = np.arctan(Cw2/self.Ca)

        print('beta1 = {0:.3f}'.format(b1))
        print('beta2 = {0:.3f}'.format(b2))
        print('Cw2 = {0:.2f} m/s'.format(Cw2))

        deHaller_test = np.cos(b1)/np.cos(b2)
        print('V2/V1 = {0:.3f}'.format(deHaller_test))

        if deHaller_test < self.deHaller:
            print('FIRST STAGE FAILS DEHALLER TEST')
            print('V2/V1 = {0:.3f}'.format(deHaller_test))

        self.alpha1.append(a1)
        self.beta1.append(b1)
        self.alpha2.append(a2)
        self.beta2.append(b2)

        po_ratio = (1 + ((self.eta_inf_c*dT_first_last)/self.To13))**(self.gamma_c/(self.gamma_c-1))
        po3_1 = self.po13*po_ratio
        To3_1 = self.To13 + dT_first_last
        return
        


    def HPturbineab(self,Ut_t,phi,psi,rRatio):
        #########################################
        #              DESIGN INPUTS            #
        #########################################
        self.Ut_t = Ut_t 
        self.phi = phi
        self.psi = psi
        self.rRatio = rRatio # keeping the same as compressor
        #########################################

        # Carryover from compressor
        self.N = self.N

        # Inlet conditions assuming combustor efficiency as 0.98
        self.po4 = self.po3*self.eta_b

        # Calculate our axial velocity with our givens
        self.Ca_t = self.phi*self.Ut_t
        print(self.Ca_t)
        self.dTos_t = (self.psi*(self.Ut_t**2))/(2*self.cpg)

        # Now we will get conditions at turbine inlet assuming there is no whirl component here (Cw1=0) so C1=Ca1
        self.T4 = self.To4 - ((self.Ca_t**2)/(2*self.cpg))
        self.p4 = self.po4*((self.T4/self.To4)**(self.gamma_h/(self.gamma_h-1)))
        self.rho4 = (self.p4*1e5)/(self.R*self.T4)

        # Radii at tip, root, and mean
        self.rt_4 = self.Ut_t/(2*np.pi*self.N)
        self.rr_4 = np.sqrt((self.rt_4**2)-(self.mc/(np.pi*self.rho4*self.Ca_t)))
        self.rm_4 = (self.rt_4+self.rr_4)/2 # constant through all turbine stages

        # To get the height of the blades
        self.h4 = 2*(self.rt_4-self.rm_4)

        # Area at inlet
        self.A4 = 2*np.pi*self.rm_4*self.h4

        # Finding work required by compressor for dTo45_t
        self.w_req = self.eta_inf_c*self.cpa*(self.To3-self.To2)
        self.dTo45 = self.w_req/self.cpg

        # Conditions at outlet
        self.To5 = self.To4 - self.dTo45
        self.po5 = self.po4*(self.To5/self.To4)**(self.gamma_h/(self.eta_inf_t*(self.gamma_h-1)))
        self.T5 = self.To5 - ((self.Ca_t**2)/(2*self.cpg))
        self.p5 = self.po5*(self.T5/self.To5)**(self.gamma_h/(self.gamma_h-1))
        self.rho5 = self.p5*1e5/(self.R*self.T5)

        # Area at outlet
        self.A5 = self.mc/(self.rho5*self.Ca_t)

        # Height at outlet
        self.h5 = self.A5/(self.rm_4*2*np.pi)

        # Radii at outlet tip and root
        self.rt_5 = self.rm_4 + (self.h5/2)
        self.rr_5 = self.rm_4 - (self.h5/2)


        # Check Mach number at tip since this will be maximum
        self.V4t = np.sqrt((self.Ut_t**2)+(self.Ca_t**2))
        self.a4 = np.sqrt(self.gamma_h*self.R*self.T4)
        self.M4t = self.V4t/self.a4

        # Check Mach number at tip since this will be maximum
        self.V5t = np.sqrt((self.Ut_t**2)+(self.Ca_t**2))
        self.a5 = np.sqrt(self.gamma_h*self.R*self.T5)
        self.M5t = self.V5t/self.a5

        # Now we will get number of stages from the temperature change across the stage and then across the entire turbine
        self.t_stages = self.dTo45/self.dTos_t 


        # now check tip mach number
        if self.M4t > 1.2:
            print('###########################')
            print('TIP MACH NUMBER EXCEEDS 1.2')
            print('###########################')
            print('M4t = {0:.2f}'.format(self.M4t))
            # potentially add system exit here??

        # now check tip mach number
        if self.M5t > 1.2:
            print('###########################')
            print('TIP MACH NUMBER EXCEEDS 1.2')
            print('###########################')
            print('M5t = {0:.2f}'.format(self.M5t))
            # potentially add system exit here??

        self.PR1 = (1-(self.dTos_t/(self.eta_inf_t*self.To4)))**(self.gamma_h/(self.gamma_h-1))
        self.PR2 = (1-(self.dTos_t/(self.eta_inf_t*(self.To4-self.dTos_t))))**(self.gamma_h/(self.gamma_h-1))
        self.PR3 = (1-(self.dTos_t/(self.eta_inf_t*(self.To4-(2*self.dTos_t)))))**(self.gamma_h/(self.gamma_h-1))
        

        if self.showValues:
            print('###############################################')
            print('TURBINE SPEED, DIMENSIONS, AND NUMBER OF STAGES')
            print('###############################################')
            print('To4 = {0:.2f} K'.format(self.To4))
            print('po4 = {0:.2f} bar'.format(self.po4))
            print('T4 = {0:.2f} K'.format(self.T4))
            print('p4 = {0:.2f} bar'.format(self.p4))
            print('rho4 = {0:.3f} kg/m^3'.format(self.rho4))
            print('To5 = {0:.2f} K'.format(self.To5))
            print('po5 = {0:.2f} bar'.format(self.po5))
            print('T5 = {0:.2f} K'.format(self.T5))
            print('p5 = {0:.2f} bar'.format(self.p5))
            print('rho5 = {0:.3f} kg/m^3'.format(self.rho5))
            print()
            print('A4 = {0:.4f} m^2'.format(self.A4))
            print('h4 = {0:.4f} m'.format(self.h4))
            print('Turbine Inlet rt = {0:.4f} m'.format(self.rt_4))
            print('Turbine Inlet rr = {0:.4f} m'.format(self.rr_4))
            print()
            print('A5 = {0:.4f} m^2'.format(self.A5))
            print('h5 = {0:.4f} m'.format(self.h5))
            print('Turbine Outlet rt = {0:.4f} m'.format(self.rt_5))
            print('Turbine Outlet rr = {0:.4f} m'.format(self.rr_5))
            print('Turbine rm = {0:.4f} m'.format(self.rm_4))
            print('N = {0:.2f} rev/s'.format(self.N))
            print('V4t = {0:.2f} m/s'.format(self.V4t))
            print('a4 = {0:.2f} m/s'.format(self.a4))
            print('M4t = {0:.2f}'.format(self.M4t))
            print()
            print('# of stages: {0:.2f}',format(self.t_stages))
            print('Temperature change across stage: {0:.2f}',format(self.dTos_t))
            print('Pressure Ratio across stage 1: {0:.2f}',format(self.PR1))
            print('Pressure Ratio across stage 2: {0:.2f}',format(self.PR2))
            print('Pressure Ratio across stage 3: {0:.2f}',format(self.PR3))

       

        return 

    def TurbineAirAngles(self,phi,psi,DegReact1,DegReact2,DegReact3):
        self.phi = phi
        self.psi = psi
        self.DegReact1 = DegReact1
        self.DegReact2 = DegReact2
        self.DegReact3 = DegReact3
        
        
        # Gas stage 1 angles to find velocity triangles later in code 
        self.beta2t1 = (180/np.pi)*np.arctan((self.psi/2)-(2*self.DegReact1)/(2*self.phi))
        self.alpha2t1 = (180/np.pi)*np.arctan(np.tan(self.beta2t1)+(1/self.phi))
        self.beta3t1 = (180/np.pi)*np.arctan((self.psi/2)+(2*self.DegReact1)/(2*self.phi))
        self.alpha3t1 = (180/np.pi)*np.arctan(np.tan(self.beta3t1)-(1/self.phi))

        # Gas stage 2 angles to find velocity triangles later in code 
        self.beta2t2 = (180/np.pi)*np.arctan((self.psi/2)-(2*self.DegReact2)/(2*self.phi))
        self.alpha2t2 = (180/np.pi)*np.arctan(np.tan(self.beta2t2)+(1/self.phi))
        self.beta3t2 = (180/np.pi)*np.arctan((self.psi/2)+(2*self.DegReact2)/(2*self.phi))
        self.alpha3t2 = (180/np.pi)*np.arctan(np.tan(self.beta3t2)-(1/self.phi))

        # Gas stage 2 angles to find velocity triangles later in code 
        self.beta2t3 = (180/np.pi)*np.arctan((self.psi/2)-(2*self.DegReact3)/(2*self.phi))
        self.alpha2t3 = (180/np.pi)*np.arctan(np.tan(self.beta2t3)+(1/self.phi))
        self.beta3t3 = (180/np.pi)*np.arctan((self.psi/2)+(2*self.DegReact3)/(2*self.phi))
        self.alpha3t3 = (180/np.pi)*np.arctan(np.tan(self.beta3t3)-(1/self.phi))

        # Calculate the rest of the velocity triangle components
        self.Vw1n    = self.Ca_t*np.tan(self.beta2t1)
        self.Cw1n    = self.Ca_t*np.tan(self.alpha2t1)
        self.V1n     = np.sqrt(self.Vw1n**2 + self.Ca_t**2)
        self.C1n     = np.sqrt(self.Cw1n**2 + self.Ca_t**2)

        self.Vw1r    = self.Ca_t*np.tan(self.beta3t1)
        self.Cw1r    = self.Ca_t*np.tan(self.alpha3t1)
        self.V1r     = np.sqrt(self.Vw1r**2 + self.Ca_t**2)
        self.C1r     = np.sqrt(self.Cw1r**2 + self.Ca_t**2)

        self.Vw2n    = self.Ca_t*np.tan(self.beta2t2)
        self.Cw2n    = self.Ca_t*np.tan(self.alpha2t2)
        self.V2n     = np.sqrt(self.Vw2n**2 + self.Ca_t**2)
        self.C2n     = np.sqrt(self.Cw2n**2 + self.Ca_t**2)

        self.Vw2r    = self.Ca_t*np.tan(self.beta3t2)
        self.Cw2r    = self.Ca_t*np.tan(self.alpha3t2)
        self.V2r     = np.sqrt(self.Vw2r**2 + self.Ca_t**2)
        self.C2r     = np.sqrt(self.Cw2r**2 + self.Ca_t**2)

        self.Vw3n    = self.Ca_t*np.tan(self.beta2t3)
        self.Cw3n    = self.Ca_t*np.tan(self.alpha2t3)
        self.V3n     = np.sqrt(self.Vw3n**2 + self.Ca_t**2)
        self.C3n     = np.sqrt(self.Cw3n**2 + self.Ca_t**2)

        self.Vw3r    = self.Ca_t*np.tan(self.beta3t3)
        self.Cw3r    = self.Ca_t*np.tan(self.alpha3t3)
        self.V3r     = np.sqrt(self.Vw3r**2 + self.Ca_t**2)
        self.C3r     = np.sqrt(self.Cw3r**2 + self.Ca_t**2)


        if self.showValues:
            print('######################################')
            print('TURBINE ANGLES AND VELOCITY COMPONENTS')
            print('######################################')
            print()
            print('Turbine Stage 1 Angles:')
            print('alpha 1 = 0 degrees') #No swirl entering turbine
            print('alpha 2 = {0:.4f} degrees'.format(self.alpha2t1))
            print('beta 2 = {0:.4f} degrees'.format(self.beta2t1))
            print('alpha 3 = {0:.4f} degrees'.format(self.alpha3t1))
            print('beta 3 = {0:.4f} degrees'.format(self.beta3t1))
            print()
            print('Turbine Stage 2 Angles:')
            print('alpha 1 = {0:.4f} degrees'.format(self.alpha3t1))
            print('alpha 2 = {0:.4f} degrees'.format(self.alpha2t2))
            print('beta 2 = {0:.4f} degrees'.format(self.beta2t2))
            print('alpha 3 = {0:.4f} degrees'.format(self.alpha3t2))
            print('beta 3 = {0:.4f} degrees'.format(self.beta3t2))
            print()
            print('Turbine Stage 3 Angles:')
            print('alpha 1 = {0:.4f} degrees'.format(self.alpha3t2))
            print('alpha 2 = {0:.4f} degrees'.format(self.alpha2t3))
            print('beta 2 = {0:.4f} degrees'.format(self.beta2t3))
            print('alpha 3 = {0:.4f} degrees'.format(self.alpha3t3))
            print('beta 3 = {0:.4f} degrees'.format(self.beta3t3))
            print()
            print('Turbine Stage 1 Velocity Triangle Components:')
            print('Vw1n = {0:.4f} m/s'.format(self.Vw1n))
            print('Cw1n = {0:.4f} m/s'.format(self.Cw1n))
            print('V1n = {0:.4f} m/s'.format(self.V1n))
            print('C1n = {0:.4f} m/s'.format(self.C1n))
            print('Vw1r = {0:.4f} m/s'.format(self.Vw1r))
            print('Cw1r = {0:.4f} m/s'.format(self.Cw1r))
            print('V1r = {0:.4f} m/s'.format(self.V1r))
            print('C1r = {0:.4f} m/s'.format(self.C1r))
            print()
            print('Turbine Stage 2 Velocity Triangle Components:')
            print('Vw2n = {0:.4f} m/s'.format(self.Vw2n))
            print('Cw2n = {0:.4f} m/s'.format(self.Cw2n))
            print('V2n = {0:.4f} m/s'.format(self.V2n))
            print('C2n = {0:.4f} m/s'.format(self.C2n))
            print('Vw2r = {0:.4f} m/s'.format(self.Vw2r))
            print('Cw2r = {0:.4f} m/s'.format(self.Cw2r))
            print('V2r = {0:.4f} m/s'.format(self.V2r))
            print('C2r = {0:.4f} m/s'.format(self.C2r))
            print()
            print('Turbine Stage 3 Velocity Triangle Components:')
            print('Vw3n = {0:.4f} m/s'.format(self.Vw3n))
            print('Cw3n = {0:.4f} m/s'.format(self.Cw3n))
            print('V3n = {0:.4f} m/s'.format(self.V3n))
            print('C3n = {0:.4f} m/s'.format(self.C3n))
            print('Vw3r = {0:.4f} m/s'.format(self.Vw3r))
            print('Cw3r = {0:.4f} m/s'.format(self.Cw3r))
            print('V3r = {0:.4f} m/s'.format(self.V3r))
            print('C3r = {0:.4f} m/s'.format(self.C3r))


        return
