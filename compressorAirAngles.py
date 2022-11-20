import numpy as np

def compressorAirAngles(self):
    ##########################
    #   ROOT RADIUS VALUES   #
    ##########################
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
    
    # for stage 1
    Ur1 = self.rr_inlet * self.N * 2*np.pi
    Ut1 = self.rt_inlet * self.N * 2*np.pi

    self.beta1_r.append(np.arctan(Ur1/self.Ca))
    self.beta1_t.append(np.arctan(Ut1/self.Ca))

    self.Cw1_m.append(self.Ca * np.tan(self.alpha1_m[0]))
    self.Cw2_m.append(self.Ca * np.tan(self.alpha2_m[0]))
    self.C3_m.append(self.Ca/np.cos(self.alpha3_m[0]))
    T3 = self.To[0] - (self.C3_m[0]**2)/(2*self.cpa)
    p3 = self.po[0]*(T3/self.To[0])**(self.gamma_c/(self.gamma_c-1))
    rho3 = p3/(self.R*T3)
    A3 = self.mc/(rho3*self.Ca)
    h = A3/(2*np.pi*self.rm)
    # conditions at STATOR exit... assumed that radii at exit of rotor blades are
    # the mean of those at rotor inlet and stator exit (pg. 221; 240 digital)
    rt3 = self.rm + (h/2)
    rr3 = self.rm - (h/2)
    rr2 = (self.rr_inlet + rr3)/2
    rt2 = (self.rt_inlet + rt3)/2

    # NOTE maybe store radii values too

    # Ur2 = rr2 * self.N * 2*np.pi
    # Ut2 = rt2 * self.N * 2*np.pi
    # Ur3 = rr3 * self.N * 2*np.pi
    # Ut3 = rt3 * self.N * 2*np.pi
    
    self.Cw2_r.append(self.Cw2_m[0] * (self.rm/rr2))
    self.Cw2_t.append(self.Cw2_m[0] * (self.rm/rt2))

    rr1 = rt3
    rt1 = rt3

    for i in range(1,self.numStages):
        self.Cw1_m.append(self.Ca * np.tan(self.alpha1_m[i]))
        self.Cw2_m.append(self.Ca * np.tan(self.alpha2_m[i]))

        self.C3_m.append(self.Ca/np.cos(self.alpha3_m[i]))
        T3 = self.To[i] - (self.C3_m[i]**2)/(2*self.cpa)
        p3 = self.po[i]*(T3/self.To[i])**(self.gamma_c/(self.gamma_c-1))
        rho3 = p3/(self.R*T3)
        A3 = self.mc/(rho3*self.Ca)
        h = A3/(2*np.pi*self.rm)

        rt3 = self.rm + (h/2)
        rr3 = self.rm - (h/2)
        rr2 = (rr1 + rr3)/2
        rt2 = (rt1 + rt3)/2

        self.Cw2_r.append(self.Cw2_m[i] * (self.rm/rr2))
        self.Cw2_t.append(self.Cw2_m[i] * (self.rm/rt2))