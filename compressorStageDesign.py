import numpy as np
import matplotlib.pyplot as plt
from sys import exit

def compressorStageDesign(self,numStages):
        self.alpha1_m = []
        self.beta1_m = []
        self.alpha2_m = []
        self.beta2_m = []
        self.alpha3_m = []

        self.dToS = (self.To3-self.To13)/numStages
        # need to estimate change in temperature for every stage based on dToS
        dT_first = self.dToS*0.85 # kinda guestimate for first delta T
        dT = self.dToS*1.1

        print('dToS = {0:.2f} K'.format(self.dToS))
        print('First dToS = {0:.2f} K'.format(dT_first))
        print('Other stage dToS - {0:.2f} K'.format(dT))

        lam = [0.98, 0.93, 0.88, 0.84, 0.83]

        To_dl = 0
        po_dl = 0
        Lam_dl = 0

        for i in range(numStages):
            if i == 0:
                Cw2 = (self.cpa*dT_first)/(lam[0]*self.Um)
                a1 = 0
                b1 = np.arctan(self.Um/self.Ca)
                b2 = np.arctan((self.Um-Cw2)/self.Ca)
                a2 = np.arctan(Cw2/self.Ca)

                # check deHaller on rotor
                deHaller_test = np.cos(b1)/np.cos(b2)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} ROTOR FAILS DEHALLER AT MEAN, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')

                poRatio = 1 + ((self.eta_inf_c*dT_first)/self.To13)**(self.gamma_c/(self.gamma_c-1))
                po_dl = self.po13*poRatio
                To_dl = self.To13 + dT_first
                Lam_dl = 1 - (Cw2/(2*self.Um))

                self.beta1_m.append(b1)
                self.beta2_m.append(b2)
                self.alpha1_m.append(a1)
                self.alpha2_m.append(a2)                

            elif i == 1:
                b1 = np.arctan((self.cpa*dT)/(2*lam[1]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                b2 = np.arctan(-(self.cpa*dT)/(2*lam[1]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                a1 = a3_dl = np.arctan((self.Um/self.Ca) - np.tan(b1))               
                a2 = np.arctan((self.Um/self.Ca) - np.tan(b2))

                # check deHaller for previous stage
                deHaller_test = np.cos(self.alpha2_m[i-1])/np.cos(a3_dl)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')

                deHaller_test = np.cos(b1)/np.cos(b2)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')
                
                Cw1 = self.Ca*np.tan(a1)
                Cw2 = self.Ca*np.tan(a2)
                Cwm = (Cw1+Cw2)/2
                Lam_dl = 1 - Cwm/self.Um
                poRatio = 1 + ((self.eta_inf_c*self.dToS)/(To_dl+self.dToS))
                po_dl = po_dl*poRatio
                To_dl = To_dl + self.dToS

                self.beta1_m.append(b1)
                self.beta2_m.append(b2)
                self.alpha1_m.append(a1)
                self.alpha2_m.append(a2)
                self.alpha3_m.append(a3_dl)
                
            elif i == 2:
                # self.dToS = self.dTos - 2 # lower the delta T_0s a bit so we can get Lambda = 0.5
                dT = dT - 2
                Lam_dl = 0.5
                b1 = np.arctan((self.cpa*dT)/(2*lam[2]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                b2 = np.arctan(-(self.cpa*dT)/(2*lam[2]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                a1 = a3_dl = b2
                a2 = b1

                deHaller_test = np.cos(self.alpha2_m[i-1])/np.cos(a3_dl)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')
                
                deHaller_test = np.cos(b1)/np.cos(b2)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')

                poRatio = 1 + ((self.eta_inf_c*self.dToS)/(To_dl + self.dToS))
                po_dl = po_dl * poRatio
                To_dl = To_dl + self.dToS

                self.beta1_m.append(b1)
                self.beta2_m.append(b2)
                self.alpha1_m.append(a1)
                self.alpha2_m.append(a2)
                self.alpha3_m.append(a3_dl)                

            elif i == numStages-1:
                poRatio = self.po3/po_dl
                dT = (poRatio**((self.gamma_c-1)/self.gamma_c) - 1)*(To_dl/self.eta_inf_c)

                b1 = np.arctan((self.cpa*dT)/(2*lam[4]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                b2 = np.arctan(-(self.cpa*dT)/(2*lam[4]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                a1 = a3_dl = b2
                a2 = b1
                a3 = 0 # because exit flow is axial

                deHaller_test = np.cos(self.alpha2_m[i-1])/np.cos(a3_dl)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')

                deHaller_test = np.cos(b1)/np.cos(b2)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')

                self.beta1_m.append(b1)
                self.beta2_m.append(b2)
                self.alpha1_m.append(a1)
                self.alpha2_m.append(a2)
                self.alpha3_m.append(a3_dl)  
                self.alpha3_m.append(a3)       

            else:
                b1 = np.arctan((self.cpa*dT)/(2*lam[3]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                b2 = np.arctan(-(self.cpa*dT)/(2*lam[3]*self.Um*self.Ca) + (self.Um*Lam_dl)/self.Ca)
                a1 = a3_dl = b2
                a2 = b1

                deHaller_test = np.cos(self.alpha2_m[i-1])/np.cos(a3_dl)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')
                
                deHaller_test = np.cos(b1)/np.cos(b2)
                if deHaller_test < self.deHaller:
                    print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                    if self.deHallerExit:
                        exit('Design is invalid')

                poRatio = 1 + ((self.eta_inf_c*self.dToS)/(To_dl + self.dToS))
                po_dl = po_dl * poRatio
                To_dl = To_dl + self.dToS

                self.beta1_m.append(b1)
                self.beta2_m.append(b2)
                self.alpha1_m.append(a1)
                self.alpha2_m.append(a2)
                self.alpha3_m.append(a3_dl)