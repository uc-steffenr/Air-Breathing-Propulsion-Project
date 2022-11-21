import numpy as np


def HPturbineSpeedandDimensions(self,phi,psi,DegReact,rRatio):
        #########################################
        #              DESIGN INPUTS            #
        #########################################
        self.phi = phi
        self.psi = psi
        self.DegReact = DegReact
        self.rRatio = rRatio # keeping the same as compressor
        #########################################

        # Carryover from compressor
        self.N = self.N

        # Inlet conditions assuming combustor efficiency as 0.98
        self.po4 = self.po3*self.eta_b

        # Gas angles to find velocity triangles later in code 
        self.beta2t = np.arctan((self.psi/2)-(2*self.DegReact)/(2*self.phi))
        self.alpha2t = np.arctan(np.tan(self.beta2t)+(1/self.phi))
        self.beta3t = np.arctan((self.psi/2)+(2*self.DegReact)/(2*self.phi))
        self.alpha3t = np.arctan(np.tan(self.beta3t)-(1/self.phi))

        # Finding work required by compressor for dToS_t
        self.w_req = self.eta_inf_c*self.cpa*(self.To3-self.To2)
        self.dToS_t = self.w_req/self.cpg

        # Using dToS_t we can get U at mean radius and axial velo
        self.Um_t = ((2*self.cpg*self.dToS_t)/self.psi)**(1/2)
        self.Ca_t = self.phi*self.Um_t # Constant across the stage

        # Now we will get conditions at turbine inlet assuming there is no whirl component here (Cw1=0) so C1=Ca1
        self.T4 = self.To4 - ((self.Ca_t**2)/(2*self.cpg))
        self.p4 = self.po4*((self.T4/self.To4)**(self.gamma_h/(self.gamma_h-1)))
        self.rho4 = (self.p4*1e5)/(self.R*self.T4)

        # So area at turbine inlet:
        self.A4 = self.mf/(self.rho4*self.Ca_t)

        # To get height we need to get radii at tip, hub, and middle
        self.rt_4 = np.sqrt(self.mf/(np.pi*(self.rho4*self.Ca_t*(1-self.rRatio**2))))
        self.rr_4 = self.rt_4*self.rRatio
        self.rm_4 = (self.rt_4 + self.rr_4) / 2
        self.h4 = self.A4/(2*np.pi*self.rm_4)
        self.rt_5 = self.rm_4 + (self.h4/2)
        self.rr_5 = self.rm_4 - (self.h4/2)

        # Check Mach number at tip since this will be maximum
        self.Ut_4 = 2*np.pi*self.rt_4*self.N
        self.V4t = np.sqrt((self.Ut_4**2)+(self.Ca_t**2))
        self.a4 = np.sqrt(self.gamma_h*self.R*self.T4)
        self.M4t = self.V4t/self.a4


        if self.showValues:
            print('############################')
            print('TURBINE SPEED AND DIMENSIONS')
            print('############################')
            print('To4 = {0:.2f} K'.format(self.To4))
            print('po4 = {0:.2f} bar'.format(self.po4))
            print('T4 = {0:.2f} K'.format(self.T4))
            print('p4 = {0:.2f} bar'.format(self.p4))
            print('rho4 = {0:.3f} kg/m^3'.format(self.rho4))
            print()
            print('A4 = {0:.4f} m^2'.format(self.A4))
            print('h = {0:.4f} m'.format(self.h))
            print('Turbine Inlet rt = {0:.4f} m'.format(self.rt_4))
            print('Turbine Inlet rr = {0:.4f} m'.format(self.rr_4))
            print('Turbine Outlet rt = {0:.4f} m'.format(self.rt_5))
            print('Turbine Outlet rr = {0:.4f} m'.format(self.rr_5))
            print('Turbine rm = {0:.4f} m'.format(self.rm_4))
            print('N = {0:.2f} rev/s'.format(self.N))
            print('V4t = {0:.2f} m/s'.format(self.V4t))
            print('a4 = {0:.2f} m/s'.format(self.a4))
            print('M4t = {0:.2f}'.format(self.M4t))


        # We will assume a constant mean radius
        self.rm_t = self.Um_t/(2*np.pi*self.N)
    
def turbineSpeedandDimensions(self):
    return