import numpy as np
import matplotlib.pyplot as plt

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
        print('mc = {0:.2f} kg/s'.format(self.mc))
        print('mf = {0:.2f} kg/s'.format(self.mf))
        # self.mc = 20

        self.eta_inf_c = 0.90
        self.eta_inf_f = 0.91
        self.pi_f = 1.6
        self.pi_c = 18
        self.eta_b = 0.98
        self.eta_inf_t = 0.93

        self.deHaller = 0.65
        #########################################


        #########################################
        #            ASSUMED VALUES             #
        #########################################


        #########################################

        #########################################
        #           NECESSARY CONSTANTS         #
        #########################################
        self.gamma_c = 1.4
        self.gamma_h = 1.333
        self.R = 287 # J/(kg-K)
        self.cpa = 1005 # J/(kg-K)
        self.cpg = 1148 # J/(kg-K)
        #########################################

        self.showValues = False

        for key,val in kwargs.items():
            if key == 'showValues':
                self.showValues = True


    
    def compressorSpeedandDimensions(self,Ca,rRatio,Ut):
        #########################################
        #              DESIGN INPUTS            #
        #########################################
        self.Ca = Ca
        self.rRatio = rRatio
        self.Ut = Ut
        #########################################

        # assuming not moving (C0=0) and no loss in the intake
        self.To2 = self.T0
        self.po2 = self.p0

        self.To13 = self.To2*(self.pi_f**((self.gamma_c-1)/(self.eta_inf_f*self.gamma_c)))
        self.po13 = self.pi_f * self.po2

        self.T13 = self.To13 - ((self.Ca**2)/(2*self.cpa))
        self.p13 = self.po13*(self.T13/self.To13)**(self.gamma_c/(self.gamma_c-1))
        self.rho13 = (self.p13*1e5)/(self.R*self.T13)

        self.rt_inlet = np.sqrt(self.mc/(np.pi*(self.rho13*self.Ca*(1-self.rRatio**2))))
        self.N = self.Ut/(2*np.pi*self.rt_inlet)

        self.rr_inlet = self.rRatio*self.rt_inlet
        self.rm = (self.rt_inlet + self.rr_inlet) / 2

        self.V13t = np.sqrt(self.Ut**2 + self.Ca**2)
        self.a13 = np.sqrt(self.gamma_c*self.R*self.T13)
        self.M13t = self.V13t/self.a13

        # now check tip mach number
        if self.M13t > 1.2:
            print('###########################')
            print('TIP MACH NUMBER EXCEEDS 1.2')
            print('###########################')
            print('M13t = {0:.2f}'.format(self.M13t))
            # potentially add system exit here??
        
        # assuming mean radius is constant for all stages
        # assuming exit velocity of compressor is axial and equal to Ca
        self.po3 = self.pi_c*self.pi_f*self.p0
        self.To3 = self.To13*(self.po3/self.po13)**((self.gamma_c-1)/(self.eta_inf_c*self.gamma_c))
        self.T3 = self.To3 - (self.Ca**2)/(2*self.cpa)
        self.p3 = self.po3*(self.T3/self.To3)**(self.gamma_c/(self.gamma_c-1))
        self.rho3 = (self.p3*1e5)/(self.R*self.T3)
        self.A3 = self.mc/(self.rho3*self.Ca)
        self.h = self.A3/(2*np.pi*self.rm)
        self.rt_outlet = self.rm + (self.h/2)
        self.rr_outlet = self.rm - (self.h/2)


        if self.showValues:
            print('###############################')
            print('COMPRESSOR SPEED AND DIMENSIONS')
            print('###############################')
            print('To13 = {0:.2f} K'.format(self.To13))
            print('po13 = {0:.2f} bar'.format(self.po13))
            print('T13 = {0:.2f} K'.format(self.T13))
            print('p13 = {0:.2f} bar'.format(self.p13))
            print('rho13 = {0:.3f} kg/m^3'.format(self.rho13))
            print()
            print('po3 = {0:.2f} bar'.format(self.po3))
            print('To3 = {0:.2f} K'.format(self.To3))
            print('T3 = {0:.2f} K'.format(self.T3))
            print('p3 = {0:.2f} K'.format(self.p3))
            print('rho3 = {0:.3f} kg/m^3'.format(self.rho3))
            print()
            print('A3 = {0:.4f} m^2'.format(self.A3))
            print('h = {0:.4f} m'.format(self.h))
            print('Inlet rt = {0:.4f} m'.format(self.rt_inlet))
            print('Inlet rr = {0:.4f} m'.format(self.rr_inlet))
            print('Outlet rt = {0:.4f} m'.format(self.rt_outlet))
            print('Outlet rr = {0:.4f} m'.format(self.rr_outlet))
            print('rm = {0:.4f} m'.format(self.rm))
            print('N = {0:.2f} rev/s'.format(self.N))
            print('V13t = {0:.2f} m/s'.format(self.V13t))
            print('a13 = {0:.2f} m/s'.format(self.a13))
            print('M13t = {0:.2f}'.format(self.M13t))
    
    def compressorStageEstimation(self):
        self.Um = 2*np.pi*self.N*self.rm
        # Ca1 = Ca2 = Ca
        b1 = np.arctan(self.Um/self.Ca)
        V1 = self.Ca/np.cos(b1)
        V2 = V1 * self.deHaller
        b2 = np.arccos(self.Ca/V2)
        # neglecting work-done factor
        self.dToS = (self.Ca*self.Um*(np.tan(b1)-np.tan(b2)))/self.cpa
        self.est_stages = (self.To3-self.To13)/self.dToS

        if self.showValues:
            print('###########################')
            print('COMPRESSOR STAGE ESTIMATION')
            print('###########################')
            print('Um = {0:.2f} m/s'.format(self.Um))
            print('beta1 = {0:.2f} degrees'.format(np.rad2deg(b1)))
            print('V1 = {0:.2f} m/s'.format(V1))
            print('V2 = {0:.2f} m/s'.format(V2))
            print('beta2 = {0:.2f} degrees'.format(np.rad2deg(b2)))
            print('delta ToS = {0:.2f} K'.format(self.dToS))
            print('Estimated Stages = {0:.2f} stages'.format(self.est_stages))

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
    
    def airAngles(self):
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
        self.dTos_t = (self.psi*(self.Ut_t**2))/(2*self.cpg)
        print(self.Ca_t)

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
        self.rt_5 = self.rm_4 + (self.h4/2)
        self.rr_5 = self.rm_4 - (self.h4/2)


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
            print('# of stages: ',self.t_stages)


        
        # Gas angles to find velocity triangles later in code 
        # self.beta2t = np.arctan((self.psi/2)-(2*self.DegReact)/(2*self.phi))
        # self.alpha2t = np.arctan(np.tan(self.beta2t)+(1/self.phi))
        # self.beta3t = np.arctan((self.psi/2)+(2*self.DegReact)/(2*self.phi))
        # self.alpha3t = np.arctan(np.tan(self.beta3t)-(1/self.phi))
        





        return 

