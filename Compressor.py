import numpy as np
import matplotlib.pyplot as plt

class Compressor:
    def __init__(self,BPR,ma,p0,T0):
        self.BPR = BPR
        self.ma = ma
        self.p0 = p0
        self.T0 = T0

        self.mc = self.m/(1+self.BPR)
        self.mf = self.m - self.mc

        self.eta_inf_c = 0.90
        self.eta_inf_f = 0.91
        self.pi_f = 1.6

        self.gamma_c = 1.4
        self.gamma_h = 1.333
        self.R = 287 # J/(kg-K)
        self.cpa = 1005 # J/(kg-K)
        self.cpg = 1148 # J/(kg-K)


    def speedAndDimension(self,Ca,rRatio,Ut):
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

        self.T13 = self.Tt13 - ((self.Ca**2)/(2*self.cpg))
        self.p13 = self.pt13*(self.T13/self.Tt13)**(self.gamma_c/(self.gamma_c-1))
        self.rho13 = self.p13/(self.R*self.T13)

        self.rt = np.sqrt(self.mc/(np.pi(self.rho13*self.Ca*(1-self.rRatio**2))))
        self.N = self.Ut/(2*np.pi*self.rt)

        self.Vo13 = np.sqrt(self.Ut**2 + self.Ca**2)
    
    def numStageEstimation(self):
        return
    
    def stageDesign(self):
        return
    
    def airAngles(self):
        return