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
    
    def speedAndDimension(self):
        return
    
    def numStageEstimation(self):
        return
    
    def stageDesign(self):
        return
    
    def airAngles(self):
        return