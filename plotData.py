import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

def plotAirAngles(ar,am,at,br,bm,bt):
    # ar[0] -> alpha1_r
    # ar[1] -> alpha2_r
    # similar for others, won't include alpha3 in these

    n = len(ar[0])
    xticks = ['root','mean','tip']
    x = [0,1,2]

    for i in range(n):
        a1 = [ar[0][i],am[0][i],at[0][i]]
        a2 = [ar[1][i],am[1][i],at[1][i]]
        b1 = [br[0][i],bm[0][i],bt[0][i]]
        b2 = [br[1][i],bm[1][i],bt[1][i]]

        a1 = [j*180/np.pi for j in a1]
        a2 = [j*180/np.pi for j in a2]
        b1 = [j*180/np.pi for j in b1]
        b2 = [j*180/np.pi for j in b2]

        fig,ax = plt.subplots()
        ax.plot(x,a1,label='a1')
        ax.plot(x,a2,label='a2')
        ax.plot(x,b1,label='b1')
        ax.plot(x,b2,label='b2')
        ax.set_xticks(x,xticks)
        ax.legend()

    return

def plotPo(po):
    return

def plotTo(To):
    return

def plotRadii(rr,rm,rt,type='mean'):
    return