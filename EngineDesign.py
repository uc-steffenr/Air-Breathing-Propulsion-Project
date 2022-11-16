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

        self.eta_inf_c = 0.90
        self.eta_inf_f = 0.91
        self.pi_f = 1.6
        self.pi_c = 18
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
        V2 = V1 * 0.72
        b2 = np.arccos(self.Ca/V2)
        # neglecting work-done factor
        dToS = (self.Ca*self.Um*(np.tan(b1)-np.tan(b2)))/self.cpa
        self.est_stages = (self.To3-self.To2)/dToS

        if self.showValues:
            print('###########################')
            print('COMPRESSOR STAGE ESTIMATION')
            print('###########################')
            print('Um = {0:.2f} m/s'.format(self.Um))
            print('beta1 = {0:.2f} degrees'.format(np.rad2deg(b1)))
            print('V1 = {0:.2f} m/s'.format(V1))
            print('V2 = {0:.2f} m/s'.format(V2))
            print('beta2 = {0:.2f} degrees'.format(np.rad2deg(b2)))
            print('delta ToS = {0:.2f} K'.format(dToS))
            print('Estimated Stages = {0:.2f} stages'.format(self.est_stages))

        return
    
    def compressorStageDesign(self):
        return
    
    def airAngles(self):
        return