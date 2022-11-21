import numpy as np

def compressorStageEstimation(self):
        if self.constant == 'mean':
            self.Um = 2*np.pi*self.N*self.rm
            # Ca1 = Ca2 = Ca
            b1 = np.arctan(self.Um/self.Ca)
            V1 = self.Ca/np.cos(b1)
            V2 = V1 * self.deHaller
            b2 = np.arccos(self.Ca/V2)
            # neglecting work-done factor
            self.dToS = (self.Ca*self.Um*(np.tan(b1)-np.tan(b2)))/self.cpa
            self.est_stages = (self.To3-self.To13)/self.dToS
        else:
            b1 = np.arctan(self.Ut/self.Ca)
            V1 = self.Ca/np.cos(b1)
            V2 = V1 * self.deHaller
            b2 = np.arccos(self.Ca/V2)
            # NOTE come back and look at this... not confident
            self.dToS = (self.Ca*self.Ut*(np.tan(b1)-np.tan(b2)))/self.cpa
            self.est_stages = (self.To3-self.To13)/self.dToS

        if self.showValues:
            print('###########################')
            print('COMPRESSOR STAGE ESTIMATION')
            print('###########################')
            if self.constant == 'mean':
                print('Um = {0:.2f} m/s'.format(self.Um))
            print('beta1 = {0:.2f} degrees'.format(np.rad2deg(b1)))
            print('V1 = {0:.2f} m/s'.format(V1))
            print('V2 = {0:.2f} m/s'.format(V2))
            print('beta2 = {0:.2f} degrees'.format(np.rad2deg(b2)))
            print('delta ToS = {0:.2f} K'.format(self.dToS))
            print('Estimated Stages = {0:.2f} stages'.format(self.est_stages))