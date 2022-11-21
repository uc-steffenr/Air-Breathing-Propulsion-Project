import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

def plotAirAngles(ar,am,at,br,bm,bt):
    # ar[0] -> alpha1_r
    # ar[1] -> alpha2_r
    # similar for others, won't include alpha3 in these

    n = len(ar[0])
    xticks = ['Root','Mean','Tip']
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
        ax.set_title('Stage {0} Air Angles'.format(i+1))
        ax.set_ylabel('Angle (Degrees)')
        ax.legend()

    return

def plotPoRatio(poRatio):
    n = len(poRatio)
    x = [i+1 for i in range(n)]

    fig,ax = plt.subplots()
    ax.plot(x,poRatio)
    ax.set_xlabel('Stages')
    ax.set_ylabel('Pressure Ratio')

def plotPo(po):
    n = len(po)
    x = [i+1 for i in range(n)]

    fig,ax = plt.subplots()
    ax.plot(x,po)
    ax.set_xlabel('Stages')
    ax.set_ylabel('Stagnation Pressure')


def plotTo(To):
    n = len(To)
    x = [i+1 for i in range(n)]

    fig,ax = plt.subplots()
    ax.plot(x,To)
    ax.set_xlabel('Stages')
    ax.set_ylabel('Stagnation Temperature')
    return

def plotRadii(rr,rm,rt,type='mean'):
    n = len(rr)
    if type == 'mean':
        rm_ = np.ones(n) * rm
        x = [i for i in range(n)]
        



    return