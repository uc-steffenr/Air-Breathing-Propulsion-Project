import numpy as np
from sys import exit

def compressorAirAngles(self):        
    if self.constant == 'mean':
        # this part is assuming constant mean diameter
        # for stage 1
        # self.rr.append(self.rr_inlet)
        # self.rt.append(self.rt_inlet)
        Ur1 = self.rr_inlet * self.N * 2*np.pi
        Ut1 = self.rt_inlet * self.N * 2*np.pi

        self.beta1_r.append(np.arctan(Ur1/self.Ca))
        self.beta1_t.append(np.arctan(Ut1/self.Ca))

        self.Cw1_m.append(self.Ca * np.tan(self.alpha1_m[0]))
        self.Cw2_m.append(self.Ca * np.tan(self.alpha2_m[0]))
        self.C3_m.append(self.Ca/np.cos(self.alpha3_m[0]))
        T3 = self.To[0] - (self.C3_m[0]**2)/(2*self.cpa)
        p3 = self.po[0]*(T3/self.To[0])**(self.gamma_c/(self.gamma_c-1))
        rho3 = (p3*1e5)/(self.R*T3)
        A3 = self.mc/(rho3*self.Ca)
        h = A3/(2*np.pi*self.rm)
        print('h = ',h)
        # conditions at STATOR exit... assumed that radii at exit of rotor blades are
        # the mean of those at rotor inlet and stator exit (pg. 221; 240 digital)
        rt3 = self.rm + (h/2)
        rr3 = self.rm - (h/2)
        rr2 = (self.rr_inlet + rr3)/2
        rt2 = (self.rt_inlet + rt3)/2

        self.rr.append(rr2)
        self.rr.append(rr3)
        self.rt.append(rt2)
        self.rt.append(rt3)

        # NOTE maybe store radii values too
        
        self.Cw1_r.append(self.Cw1_m[0] * (self.rm/self.rr_inlet))
        self.Cw1_t.append(self.Cw1_m[0] * (self.rm/self.rt_inlet))
        self.Cw2_r.append(self.Cw2_m[0] * (self.rm/rr2))
        self.Cw2_t.append(self.Cw2_m[0] * (self.rm/rt2))

        self.alpha1_r.append(np.arctan(self.Cw1_r[0]/self.Ca))
        self.alpha2_r.append(np.arctan(self.Cw2_r[0]/self.Ca))
        self.beta2_r.append(np.arctan((self.N*2*np.pi*rr2 - self.Cw2_r[0])/self.Ca))

        self.alpha1_t.append(np.arctan(self.Cw1_t[0]/self.Ca))
        self.alpha2_t.append(np.arctan(self.Cw2_t[0]/self.Ca))
        self.beta2_t.append(np.arctan((self.N*2*np.pi*rt2 - self.Cw2_t[0])/self.Ca))

        rr1 = rr3
        rt1 = rt3

        for i in range(1,self.numStages):
            self.Cw1_m.append(self.Ca * np.tan(self.alpha1_m[i]))
            self.Cw2_m.append(self.Ca * np.tan(self.alpha2_m[i]))

            self.C3_m.append(self.Ca/np.cos(self.alpha3_m[i]))
            T3 = self.To[i] - (self.C3_m[i]**2)/(2*self.cpa)
            p3 = self.po[i]*(T3/self.To[i])**(self.gamma_c/(self.gamma_c-1))
            # print(T3)
            # print(p3)
            rho3 = (p3*1e5)/(self.R*T3)
            A3 = self.mc/(rho3*self.Ca)
            h = A3/(2*np.pi*self.rm)

            rt3 = self.rm + (h/2)
            rr3 = self.rm - (h/2)
            rr2 = (rr1 + rr3)/2
            rt2 = (rt1 + rt3)/2

            self.rt.append(rt2)
            self.rt.append(rt3)
            self.rr.append(rr2)
            self.rr.append(rr3)

            self.Cw1_r.append(self.Cw1_m[i] * (self.rm/rr1))
            self.Cw1_t.append(self.Cw1_m[i] * (self.rm/rt1))
            self.Cw2_r.append(self.Cw2_m[i] * (self.rm/rr2))
            self.Cw2_t.append(self.Cw2_m[i] * (self.rm/rt2))

            # if it is station 2, calculate mean angles, otherwise we know Lam = 0.5
            # so we can use a1=b2 and b1=a2

            self.alpha1_r.append(np.arctan(self.Cw1_r[i]/self.Ca))
            self.beta1_r.append(np.arctan((self.N*2*np.pi*rr1 - self.Cw1_r[i])/self.Ca))
            self.alpha2_r.append(np.arctan(self.Cw2_r[i]/self.Ca))
            self.beta2_r.append(np.arctan((self.N*2*np.pi*rr2 - self.Cw2_r[i])/self.Ca))

            self.alpha1_t.append(np.arctan(self.Cw1_t[i]/self.Ca))
            self.beta1_t.append(np.arctan((self.N*2*np.pi*rt1 - self.Cw1_t[i])/self.Ca))
            self.alpha2_t.append(np.arctan(self.Cw2_t[i]/self.Ca))
            self.beta2_t.append(np.arctan((self.N*2*np.pi*rt2 - self.Cw2_t[i])/self.Ca))
    
    ##### Outer Diameter Stuff ####
    else:
        Ur1 = self.rr_inlet * self.N * 2*np.pi
        Um1 = self.rm_inlet * self.N * 2*np.pi
        self.rm = []

        self.beta1_r.append(np.arctan(Ur1/self.Ca))
        self.beta1_m.append(np.arctan(Um1/self.Ca))

        self.Cw1_t.append(self.Ca * np.tan(self.alpha1_t[0]))
        self.Cw2_t.append(self.Ca * np.tan(self.alpha2_t[0]))
        self.C3_t.append(self.Ca/np.cos(self.alpha3_t[0]))
        T3 = self.To[0] - (self.C3_t[0]**2)/(2*self.cpa)
        p3 = self.po[0]*(T3/self.To[0])**(self.gamma_c/(self.gamma_c-1))
        rho3 = (p3*1e5)/(self.R*T3)

        rr3 = np.sqrt(self.rt**2 - self.mc/(np.pi*rho3*self.Ca))
        rm3 = (self.rt + rr3)/2
        rr2 = (self.rr_inlet + rr3)/2
        rm2 = (self.rm_inlet + rm3)/2
        
        self.rr.append(self.rr_inlet)
        self.rr.append(rr2)
        self.rr.append(rr3)
        self.rm.append(self.rm_inlet)
        self.rm.append(rm2)
        self.rm.append(rm3)

        # NOTE maybe store radii values too
        
        self.Cw1_r.append(self.Cw1_t[0] * (self.rt/self.rr_inlet))
        self.Cw1_m.append(self.Cw1_t[0] * (self.rt/self.rm_inlet))
        self.Cw2_r.append(self.Cw2_t[0] * (self.rt/rr2))
        self.Cw2_m.append(self.Cw2_t[0] * (self.rt/rm2))

        self.alpha1_r.append(np.arctan(self.Cw1_r[0]/self.Ca))
        self.alpha2_r.append(np.arctan(self.Cw2_r[0]/self.Ca))
        self.beta2_r.append(np.arctan((self.N*2*np.pi*rr2 - self.Cw2_r[0])/self.Ca))

        self.alpha1_m.append(np.arctan(self.Cw1_m[0]/self.Ca))
        self.alpha2_m.append(np.arctan(self.Cw2_m[0]/self.Ca))
        self.beta2_m.append(np.arctan((self.N*2*np.pi*rm2 - self.Cw2_m[0])/self.Ca))

        rr1 = rr3
        rm1 = rm3

        for i in range(1,self.numStages):
            self.Cw1_t.append(self.Ca * np.tan(self.alpha1_t[i]))
            self.Cw2_t.append(self.Ca * np.tan(self.alpha2_t[i]))

            self.C3_t.append(self.Ca/np.cos(self.alpha3_t[i]))
            T3 = self.To[i] - (self.C3_t[i]**2)/(2*self.cpa)
            p3 = self.po[i]*(T3/self.To[i])**(self.gamma_c/(self.gamma_c-1))
            # print(T3)
            # print(p3)
            rho3 = (p3*1e5)/(self.R*T3)

            rr3 = np.sqrt(self.rt**2 - self.mc/(np.pi*rho3*self.Ca))
            rm3 = (self.rt + rr3)/2
            rr2 = (rr1 + rr3)/2
            rm2 = (rm1 + rm3)/2

            self.rr.append(rr2)
            self.rr.append(rr3)
            self.rm.append(rm2)
            self.rm.append(rm3)

            self.Cw1_r.append(self.Cw1_t[i] * (self.rt/rr1))
            self.Cw1_m.append(self.Cw1_t[i] * (self.rt/rm1))
            self.Cw2_r.append(self.Cw2_t[i] * (self.rt/rr2))
            self.Cw2_m.append(self.Cw2_t[i] * (self.rt/rm2))

            # if it is station 2, calculate mean angles, otherwise we know Lam = 0.5
            # so we can use a1=b2 and b1=a2

            self.alpha1_r.append(np.arctan(self.Cw1_r[i]/self.Ca))
            self.beta1_r.append(np.arctan((self.N*2*np.pi*rr1 - self.Cw1_r[i])/self.Ca))
            self.alpha2_r.append(np.arctan(self.Cw2_r[i]/self.Ca))
            self.beta2_r.append(np.arctan((self.N*2*np.pi*rr2 - self.Cw2_r[i])/self.Ca))

            self.alpha1_m.append(np.arctan(self.Cw1_m[i]/self.Ca))
            self.beta1_m.append(np.arctan((self.N*2*np.pi*rm1 - self.Cw1_m[i])/self.Ca))
            self.alpha2_m.append(np.arctan(self.Cw2_m[i]/self.Ca))
            self.beta2_m.append(np.arctan((self.N*2*np.pi*rm2 - self.Cw2_m[i])/self.Ca))