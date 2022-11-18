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

        self.mc = self.m/(1+self.BPR)
        self.mf = self.m - self.mc
        print('mc = {0:.2f} kg/s'.format(self.mc))
        print('mf = {0:.2f} kg/s'.format(self.mf))
        # self.mc = 20

        self.eta_inf_c = 0.90
        self.eta_inf_f = 0.91
        self.pi_f = 1.6
        self.pi_c = 18

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
        self.alpha1_m = []
        self.beta1_m = []
        # self.Cw1_m = []
        # self.Vw1_m = []
        # self.V1_ = []
        # self.C1 = []

        self.alpha2_m = []
        self.beta2_m = []
        # self.Cw2 = []
        # self.Vw2 = []
        # self.V2 = []
        # self.C2 = []

        self.alpha3_m = []
        self.beta3_m = [] 
        # self.Cw3 = []
        # self.Vw3 = []
        # self.C3 = []
        # self.V3 = []

        self.dToS = (self.To3-self.To13)/numStages
        # need to estimate change in temperature for every stage based on dToS
        dT_first_last = self.dToS*0.85 # kinda guestimate for first and last delta T
        dT = self.dToS*1.1

        print('dToS = {0:.2f} K'.format(self.dToS))
        print('First and last dToS = {0:.2f} K'.format(dT_first_last))
        print('Other stage dToS - {0:.2f} K'.format(dT))

        # work done factors
        # lam_1 = 0.98
        # lam_2 = 0.93
        # lam_3 = 0.88
        # lam = 0.84

        lam = [0.98, 0.93, 0.88, 0.84, 0.83]

        To_dl = 0
        po_dl = 0
        Lam_dl = 0

        for i in range(numStages):
            if i == 0:
                Cw2 = (self.cpa*dT_first_last)/(lam[0]*self.Um)
                a1 = 0
                b1 = np.arctan(self.Um/self.Ca)
                b2 = np.arctan((self.Um-Cw2)/self.Ca)
                a2 = np.arctan(Cw2/np.Ca)

                # check deHaller on rotor
                deHaller_test = np.cos(b1)/np.cos(b2)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} ROTOR FAILS DEHALLER, V1/V2 = {1:.2f}'.format(i,deHaller_test))
                
                self.alpha1_m.append(a1)
                self.alpha2_m.append(a2)
                self.beta1_m.append(b1)
                self.beta2_m.append(b2)

                poRatio = 1 + ((self.eta_inf_c*dT_first_last)/self.To13)**(self.gamma_c/(self.gamma_c-1))
                po_dl = self.po13*poRatio
                To_dl = self.To13 + dT_first_last
                Lam_dl = 1 - (Cw2/(2*self.Um))

            elif i == 1:
                b1 = np.arctan((self.cpa*dT)/(2*Lam_dl*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                b2 = np.arctan(-(self.cpa*dT)/(2*Lam_dl*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                a1 = a3_dl = np.arctan((self.Um/self.Ca) - np.tan(b1))
                self.alpha3_m.append(a3_dl)
                a2 = np.arctan((self.Um/self.Ca) - np.tan(b2))
                deHaller_test = np.cos(a2)/np.cos(a3_dl)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} STATOR FAILS DEHALLER, V3/V3 = {1:.2f}'.format(i-1,deHaller_test))
                Cwm = (Cw1+Cw2)/2
                Lam = 1 - Cwm/self.Um
                poRatio = 1 + ((self.eta_inf_c*self.dToS)/(To_dl+self.dToS))
                po_dl = po_dl*poRatio
                To_dl = To_dl + self.dToS
                
            elif i == 3:
                self.dToS = self.dTos - 2
                

            elif i == numStages-1:
                break
            else:
                break

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
            print('FIRST STAGE ROTOR FAILS DEHALLER TEST')
            print('V2/V1 = {0:.3f}'.format(deHaller_test))

        self.alpha1.append(a1)
        self.beta1.append(b1)
        self.alpha2.append(a2)
        self.beta2.append(b2)

        po_ratio = (1 + ((self.eta_inf_c*dT_first_last)/self.To13))**(self.gamma_c/(self.gamma_c-1))
        po3_1 = self.po13*po_ratio
        To3_1 = self.To13 + dT_first_last

        Lam = 1 - (Cw2)/(2*self.Um)
        b1 = np.arctan((self.cpa*dT)/(2*Lam*self.Um*self.Ca) + (self.Um*Lam)/(self.Ca))
        b2 = np.arctan(-(self.cpa*dT)/(2*Lam*self.Um*self.Ca) + (self.Um*Lam)/(self.Ca))
        a1 = np.arctan((self.Um/self.Ca) - np.tan(b1))
        a2 = np.arctan((self.Um/self.Ca) - np.tan(b2))

        Cw1 = self.Ca*np.tan(a1)
        Cw2 = self.Ca*np.tan(a2)

        deHaller_test = np.cos(a2)/np.cos(a1)
        print('V2/V1 = {0:.3f}'.format(deHaller_test))
        if deHaller_test < self.deHaller:
            print('FIRST STAGE STATOR FAILS DEHALLER TEST')
            print('V2/V1 = {0:.3f}'.format(deHaller_test))


        return
    
    def airAngles(self):
        return