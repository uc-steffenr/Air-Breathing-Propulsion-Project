import numpy as np
import matplotlib.pyplot as plt
from sys import exit


def compressorStageDesign(self,numStages,**kwargs):
        self.numStages = numStages
        self.To = [self.To13]
        self.po = [self.po13]
        self.poRatio = []

        self.dToS = (self.To3-self.To13)/numStages
        # need to estimate change in temperature for every stage based on dToS
        ### can tune these values
        dT_first = self.dToS*0.95 # kinda guestimate for first delta T
        # dT_first = 58
        dT = self.dToS*1.05

        print('dToS = {0:.2f} K'.format(self.dToS))
        print('First dToS = {0:.2f} K'.format(dT_first))
        print('Other stage dToS = {0:.2f} K'.format(dT))

        lam = [0.98, 0.93, 0.88, 0.84, 0.83]
        Lam2 = 0.7
        Lam_ = 0.55 # used for constant outer diameter

        To_dl = 0
        po_dl = 0

        for key,val in kwargs.items():
            if key == 'Lam2':
                Lam2 = val
            elif key == 'Lam':
                Lam_ = val

        if self.constant == 'mean':
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

                    poRatio = (1 + (self.eta_inf_c*dT_first)/self.To13)**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = self.po13*poRatio
                    To_dl = self.To13 + dT_first
                    Lam = 1 - (Cw2/(2*self.Um))

                    if self.showValues:
                        print('STAGE 1 DEGREE OF REACTION: ',Lam)

                    self.beta1_m.append(b1)
                    self.beta2_m.append(b2)
                    self.alpha1_m.append(a1)
                    self.alpha2_m.append(a2)
                    self.To.append(To_dl)
                    self.po.append(po_dl)                

                elif i == 1:
                    b1 = np.arctan((self.cpa*dT)/(2*lam[1]*self.Um*self.Ca) + (self.Um*Lam2)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[1]*self.Um*self.Ca) + (self.Um*Lam2)/self.Ca)
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

                    poRatio = (1 + ((self.eta_inf_c*dT)/(To_dl)))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = po_dl*poRatio
                    To_dl = To_dl + dT

                    self.beta1_m.append(b1)
                    self.beta2_m.append(b2)
                    self.alpha1_m.append(a1)
                    self.alpha2_m.append(a2)
                    self.alpha3_m.append(a3_dl)
                    self.To.append(To_dl)
                    self.po.append(po_dl) 
                    
                elif i == 2:
                    # self.dToS = self.dTos - 2 # lower the delta T_0s a bit so we can get Lambda = 0.5
                    dT = dT - 3 # can tune this too
                    Lam = 0.5
                    b1 = np.arctan((self.cpa*dT)/(2*lam[2]*self.Um*self.Ca) + (self.Um*Lam)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[2]*self.Um*self.Ca) + (self.Um*Lam)/self.Ca)
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

                    poRatio = (1 + (self.eta_inf_c*dT)/(To_dl))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = po_dl * poRatio
                    To_dl = To_dl + dT

                    self.beta1_m.append(b1)
                    self.beta2_m.append(b2)
                    self.alpha1_m.append(a1)
                    self.alpha2_m.append(a2)
                    self.alpha3_m.append(a3_dl)   
                    self.To.append(To_dl)
                    self.po.append(po_dl)              

                elif i == numStages-1:
                    # print(To_dl)
                    # print(po_dl)
                    poRatio = self.po3/po_dl
                    print('LAST PORATIO = ',poRatio)
                    dT = (poRatio**((self.gamma_c-1)/self.gamma_c) - 1)*(To_dl/self.eta_inf_c)
                    print('LAST DT = ',dT)
                    print('LAM = ',Lam)

                    b1 = np.arctan((self.cpa*dT)/(2*lam[4]*self.Um*self.Ca) + (self.Um*Lam)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[4]*self.Um*self.Ca) + (self.Um*Lam)/self.Ca)
                    a1 = a3_dl = b2
                    a2 = b1
                    a3 = 0 # because exit flow is axial

                    poRatio = (1 + (self.eta_inf_c*dT)/(To_dl))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    self.po3_real = po_dl*poRatio
                    self.To3_real = To_dl + dT

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
                    b1 = np.arctan((self.cpa*dT)/(2*lam[3]*self.Um*self.Ca) + (self.Um*Lam)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[3]*self.Um*self.Ca) + (self.Um*Lam)/self.Ca)
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

                    poRatio = (1 + (self.eta_inf_c*dT)/(To_dl))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = po_dl * poRatio
                    To_dl = To_dl + dT

                    self.beta1_m.append(b1)
                    self.beta2_m.append(b2)
                    self.alpha1_m.append(a1)
                    self.alpha2_m.append(a2)
                    self.alpha3_m.append(a3_dl)
                    self.To.append(To_dl)
                    self.po.append(po_dl) 

# NOTE pretty much just gonna change the Um's to Ut's and see how that goes...
        else:
            for i in range(numStages):
                if i == 0:
                    Cw2 = (self.cpa*dT_first)/(lam[0]*self.Ut)
                    a1 = 0
                    b1 = np.arctan(self.Ut/self.Ca)
                    b2 = np.arctan((self.Ut-Cw2)/self.Ca)
                    a2 = np.arctan(Cw2/self.Ca)

                    # check deHaller on rotor
                    deHaller_test = np.cos(b1)/np.cos(b2)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} ROTOR FAILS DEHALLER AT MEAN, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    poRatio = (1 + (self.eta_inf_c*dT_first)/self.To13)**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = self.po13*poRatio
                    To_dl = self.To13 + dT_first
                    Lam = 1 - (Cw2/(2*self.Ut))

                    if self.showValues:
                        print('STAGE 1 DEGREE OF REACTION: ',Lam)

                    self.beta1_t.append(b1)
                    self.beta2_t.append(b2)
                    self.alpha1_t.append(a1)
                    self.alpha2_t.append(a2)
                    self.To.append(To_dl)
                    self.po.append(po_dl)       

                elif i == 1:
                    b1 = np.arctan((self.cpa*dT)/(2*lam[1]*self.Ut*self.Ca) + (self.Ut*Lam2)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[1]*self.Ut*self.Ca) + (self.Ut*Lam2)/self.Ca)
                    a1 = a3_dl = np.arctan((self.Ut/self.Ca) - np.tan(b1))               
                    a2 = np.arctan((self.Ut/self.Ca) - np.tan(b2))

                    # check deHaller for previous stage
                    deHaller_test = np.cos(self.alpha2_t[i-1])/np.cos(a3_dl)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    deHaller_test = np.cos(b1)/np.cos(b2)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    poRatio = (1 + ((self.eta_inf_c*dT)/(To_dl)))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = po_dl*poRatio
                    To_dl = To_dl + dT

                    self.beta1_t.append(b1)
                    self.beta2_t.append(b2)
                    self.alpha1_t.append(a1)
                    self.alpha2_t.append(a2)
                    self.alpha3_t.append(a3_dl)
                    self.To.append(To_dl)
                    self.po.append(po_dl)     

                elif i == 2:
                    # self.dToS = self.dTos - 2 # lower the delta T_0s a bit so we can get Lambda = 0.5
                    dT = dT - 3 # can tune this too
                    Lam = Lam_
                    b1 = np.arctan((self.cpa*dT)/(2*lam[2]*self.Ut*self.Ca) + (self.Ut*Lam)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[2]*self.Ut*self.Ca) + (self.Ut*Lam)/self.Ca)
                    a1 = a3_dl = np.arctan((self.Ut/self.Ca) - np.tan(b1))               
                    a2 = np.arctan((self.Ut/self.Ca) - np.tan(b2))

                    deHaller_test = np.cos(self.alpha2_t[i-1])/np.cos(a3_dl)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')
                    
                    deHaller_test = np.cos(b1)/np.cos(b2)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    poRatio = (1 + (self.eta_inf_c*dT)/(To_dl))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = po_dl * poRatio
                    To_dl = To_dl + dT

                    self.beta1_t.append(b1)
                    self.beta2_t.append(b2)
                    self.alpha1_t.append(a1)
                    self.alpha2_t.append(a2)
                    self.alpha3_t.append(a3_dl)   
                    self.To.append(To_dl)
                    self.po.append(po_dl)     
                
                elif i == numStages-1:
                    # print(To_dl)
                    # print(po_dl)
                    poRatio = self.po3/po_dl
                    print('LAST PORATIO = ',poRatio)
                    dT = (poRatio**((self.gamma_c-1)/self.gamma_c) - 1)*(To_dl/self.eta_inf_c)
                    print('LAST DT = ',dT)
                    print('LAM = ',Lam)

                    b1 = np.arctan((self.cpa*dT)/(2*lam[4]*self.Ut*self.Ca) + (self.Ut*Lam)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[4]*self.Ut*self.Ca) + (self.Ut*Lam)/self.Ca)
                    a1 = a3_dl = np.arctan((self.Ut/self.Ca) - np.tan(b1))               
                    a2 = np.arctan((self.Ut/self.Ca) - np.tan(b2))
                    a3 = 0 # because exit flow is axial

                    poRatio = (1 + (self.eta_inf_c*dT)/(To_dl))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    self.po3_real = po_dl*poRatio
                    self.To3_real = To_dl + dT

                    deHaller_test = np.cos(self.alpha2_t[i-1])/np.cos(a3_dl)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    deHaller_test = np.cos(b1)/np.cos(b2)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    self.beta1_t.append(b1)
                    self.beta2_t.append(b2)
                    self.alpha1_t.append(a1)
                    self.alpha2_t.append(a2)
                    self.alpha3_t.append(a3_dl)  
                    self.alpha3_t.append(a3)       

                else:
                    b1 = np.arctan((self.cpa*dT)/(2*lam[3]*self.Ut*self.Ca) + (self.Ut*Lam)/self.Ca)
                    b2 = np.arctan(-(self.cpa*dT)/(2*lam[3]*self.Ut*self.Ca) + (self.Ut*Lam)/self.Ca)
                    a1 = a3_dl = np.arctan((self.Ut/self.Ca) - np.tan(b1))               
                    a2 = np.arctan((self.Ut/self.Ca) - np.tan(b2))

                    deHaller_test = np.cos(self.alpha2_t[i-1])/np.cos(a3_dl)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} STATOR FAILS DEHALLER, C3/C2 = {1:.2f}'.format(i,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')
                    
                    deHaller_test = np.cos(b1)/np.cos(b2)
                    if deHaller_test < self.deHaller:
                        print('STAGE {0} ROTOR FAILS DEHALLER, V2/V1 = {1:.2f}'.format(i+1,deHaller_test))
                        if self.deHallerExit:
                            exit('Design is invalid')

                    poRatio = (1 + (self.eta_inf_c*dT)/(To_dl))**(self.gamma_c/(self.gamma_c-1))
                    self.poRatio.append(poRatio)
                    po_dl = po_dl * poRatio
                    To_dl = To_dl + dT

                    self.beta1_t.append(b1)
                    self.beta2_t.append(b2)
                    self.alpha1_t.append(a1)
                    self.alpha2_t.append(a2)
                    self.alpha3_t.append(a3_dl)
                    self.To.append(To_dl)
                    self.po.append(po_dl) 
