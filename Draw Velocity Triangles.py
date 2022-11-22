import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d as inter
from matplotlib.patches import ConnectionPatch

MeanAngs = pd.read_csv('AngFrame_m.csv')#, skiprows=[0])
alpha1_mean = np.array(MeanAngs.iloc[:,1])
beta1_mean = np.array(MeanAngs.iloc[:,2])
alpha2_mean = np.array(MeanAngs.iloc[:,3])
beta2_mean = np.array(MeanAngs.iloc[:,4])
alpha3_mean = np.array(MeanAngs.iloc[:,5])


RootAngs = pd.read_csv('AngFrame_r.csv')#, skiprows=[0])
alpha1_root = np.array(RootAngs.iloc[:,1])
beta1_root = np.array(RootAngs.iloc[:,2])
alpha2_root = np.array(RootAngs.iloc[:,3])
beta2_root = np.array(RootAngs.iloc[:,4])
alpha3_root = np.array(RootAngs.iloc[:,5])

TipAngs = pd.read_csv('AngFrame_t.csv')#, skiprows=[0])
alpha1_tip = np.array(TipAngs.iloc[:,1])
beta1_tip = np.array(TipAngs.iloc[:,2])
alpha2_tip = np.array(TipAngs.iloc[:,3])
beta2_tip = np.array(TipAngs.iloc[:,4])
alpha3_tip = np.array(TipAngs.iloc[:,5])

MeanVel1s = pd.read_csv('Vel1Frame_m.csv', skiprows=[0])
cw1_mean = np.array(MeanVel1s.iloc[:,1])
c1_mean  = np.array(MeanVel1s.iloc[:,2])
vw1_mean = np.array(MeanVel1s.iloc[:,3])
v1_mean = np.array(MeanVel1s.iloc[:,4])

RVel1s = pd.read_csv('Vel1Frame_r.csv', skiprows=[0])
cw1_root = np.array(RVel1s.iloc[:,1])
c1_root  = np.array(RVel1s.iloc[:,2])
vw1_root = np.array(RVel1s.iloc[:,3])
v1_root = np.array(RVel1s.iloc[:,4])

TVel1s = pd.read_csv('Vel1Frame_t.csv', skiprows=[0])
cw1_tip = np.array(TVel1s.iloc[:,1])
c1_tip  = np.array(TVel1s.iloc[:,2])
vw1_tip = np.array(TVel1s.iloc[:,3])
v1_tip = np.array(TVel1s.iloc[:,4])

MVel2s = pd.read_csv('Vel2Frame_m.csv', skiprows=[0])
cw2_mean = np.array(MVel2s.iloc[:,1])
c2_mean  = np.array(MVel2s.iloc[:,2])
vw2_mean = np.array(MVel2s.iloc[:,3])
v2_mean = np.array(MVel2s.iloc[:,4])

RVel2s = pd.read_csv('Vel2Frame_r.csv', skiprows=[0])
cw2_root = np.array(RVel2s.iloc[:,1])
c2_root  = np.array(RVel2s.iloc[:,2])
vw2_root = np.array(RVel2s.iloc[:,3])
v2_root = np.array(RVel2s.iloc[:,4])

TVel2s = pd.read_csv('Vel2Frame_t.csv', skiprows=[0])
cw2_tip = np.array(TVel2s.iloc[:,1])
c2_tip  = np.array(TVel2s.iloc[:,2])
vw2_tip = np.array(TVel2s.iloc[:,3])
v2_tip  = np.array(TVel2s.iloc[:,4])

MVel3s = pd.read_csv('Vel3Frame_m.csv', skiprows=[0])
cw3_mean = np.array(MVel3s.iloc[:,1])
c3_mean  = np.array(MVel3s.iloc[:,2])

RVel3s = pd.read_csv('Vel3Frame_r.csv', skiprows=[0])
cw3_root = np.array(RVel3s.iloc[:,1])
c3_root  = np.array(RVel3s.iloc[:,2])

TVel3s = pd.read_csv('Vel3Frame_t.csv', skiprows=[0])
cw3_tip = np.array(TVel3s.iloc[:,1])
c3_tip  = np.array(TVel3s.iloc[:,2])


print(beta1_mean)







# C1y = np.cos(alpha1)*C1
# C1x = np.sin(alpha1)*C1
# Cw1 = cw1_mean[1]

Ca = 227.5

# Draw a simple arrow between two points in axes coordinates
# within a single axes.

# fig1, ax1 = plt.subplots(1)
# fig2, ax2 = plt.subplots(1)
# fig3, ax3 = plt.subplots(1)
# fig4, ax4 = plt.subplots(1)
# fig5, ax5 = plt.subplots(1)
# fig6, ax6 = plt.subplots(1)
# fig7, ax7 = plt.subplots(1)
# fig8, ax8 = plt.subplots(1)
# fig9, ax9 = plt.subplots(1)
# axList = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]


fig1, (ax11, ax12, ax13) = plt.subplots(3,1)
fig2, (ax21, ax22, ax23) = plt.subplots(3,1)
fig3, (ax31, ax32, ax33) = plt.subplots(3,1)
fig4, (ax41, ax42, ax43) = plt.subplots(3,1)
fig5, (ax51, ax52, ax53) = plt.subplots(3,1)
fig6, (ax61, ax62, ax63) = plt.subplots(3,1)
fig7, (ax71, ax72, ax73) = plt.subplots(3,1)
fig8, (ax81, ax82, ax83) = plt.subplots(3,1)


axList = [(ax11, ax12, ax13),
          (ax21, ax22, ax23),
          (ax31, ax32, ax33),
          (ax41, ax42, ax43),
          (ax51, ax52, ax53),
          (ax61, ax62, ax63),
          (ax71, ax72, ax73),
          (ax81, ax82, ax83)]

figList = [fig1,
           fig2,
           fig3,
           fig4,
           fig5,
           fig6,
           fig7,
           fig8]


def preRotor(alpha1, beta1, ca, c1, v1, cw1, vw1, axIn):


    caStart = (0., ca)
    caEnd = (0., 0)
    triangle(caStart, caEnd, "r", axIn, labelIn ="Ca")
    
    c1Start = (0., ca)
    c1End = ( cw1, 0 )
    triangle(c1Start, c1End, "r", axIn, labelIn = "C1")

    cw1Start = (0., 0)
    cw1End = (cw1, 0)
    triangle(cw1Start, cw1End, "r", axIn, labelIn = "Cw1")
    
    v1Start = (0., ca)
    v1End = ( -vw1, 0) 
    triangle(v1Start, v1End, "b", axIn, labelIn = "V1")

    vw1Start = (0., 0)
    vw1End = (-vw1, 0)
    triangle(vw1Start, vw1End, "b", axIn, labelIn = "Vw1")
    


def preStator(alpha2, beta2, ca, c2, v2, cw2, vw2, axIn):
    caStart = (0., ca)
    caEnd = (0., 0)
    triangle(caStart, caEnd, "r", axIn, labelIn = "Ca")
    
    c2Start = (0., ca)
    c2End = ( cw2, 0 )
    triangle(c2Start, c2End,"r", axIn, labelIn = "C2")

    cw2Start = (0., 0)
    cw2End = (cw2, 0)
    triangle(cw2Start, cw2End, "r", axIn, labelIn = "Cw2")
    
    v2Start = (0., ca)
    v2End = ( -vw2, 0 )
    triangle(v2Start, v2End, "b", axIn, labelIn = "V2")

    vw2Start = (0., 0)
    vw2End = (-vw2, 0)
    triangle(vw2Start, vw2End, "b", axIn, labelIn = "Vw2")
    
def postStator(alpha3, ca, c3, cw3, axIn):
    
    caStart = (0., ca)
    caEnd = (0., 0)
    triangle(caStart, caEnd, "r", axIn, labelIn ="Ca")
    
    c3Start = (0., ca)
    c3End = ( cw3, 0 )
    triangle(c3Start, c3End,  "r", axIn, labelIn ="C3")

    cw3Start = (0., 0)
    cw3End = (cw3, 0)
    triangle(cw3Start, cw3End,  "r", axIn, labelIn ="Cw3")
    

def triangle(start, stop, colorVal, axIn, labelIn):
    
    coordsA = "data"
    coordsB = "data"
    
    con = ConnectionPatch(start, stop, coordsA, coordsB, arrowstyle="-|>", color = str(colorVal), label = labelIn)
    axIn.plot([start[0], stop[0]], [start[1], stop[1]], "o", markersize = 0)
    axIn.add_artist(con)
    axIn.legend()
    
    


# for i in range(9):
#     axTemp = axList[i]
#     alpha1 = np.deg2rad(alpha1_mean)[i]
#     beta1 = np.deg2rad(beta1_mean)[i]
#     c1 = c1_mean[i]
#     v1 = v1_mean[i]
#     cw1 = cw1_mean[i]
#     vw1 = vw1_mean[i]

#     preRotor(alpha1, beta1, Ca, c1, v1, cw1, vw1, axTemp)

# for i in range(9):
#     axTemp = axList[i]
#     alpha2 = np.deg2rad(alpha2_mean)[i]
#     beta2 = np.deg2rad(beta2_mean)[i]
#     c2 = c2_mean[i]
#     v2 = v2_mean[i]
#     cw2 = cw2_mean[i]
#     vw2 = vw2_mean[i]

#     preStator(alpha2, beta2, Ca, c2, v2, cw2, vw2, axTemp)

'''
for i in range(len(alpha1_mean)):
    axTemp = axList[i][0]
    figTemp = figList[i]
    axTemp.set_title("Mean Station 1 of Stage "+str(i+1))
    alpha1 = np.deg2rad(alpha1_mean)[i]
    beta1 = np.deg2rad(beta1_mean)[i]
    c1 = c1_mean[i]
    v1 = v1_mean[i]
    cw1 = cw1_mean[i]
    vw1 = vw1_mean[i]

    preRotor(alpha1, beta1, Ca, c1, v1, cw1, vw1, axTemp)
    
    axTemp = axList[i][1]
    axTemp.set_title("Mean Station 2 of Stage "+str(i+1))
    alpha2 = np.deg2rad(alpha2_mean)[i]
    beta2 = np.deg2rad(beta2_mean)[i]
    c2 = c2_mean[i]
    v2 = v2_mean[i]
    cw2 = cw2_mean[i]
    vw2 = vw2_mean[i]

    preStator(alpha2, beta2, Ca, c2, v2, cw2, vw2, axTemp)
    
    
    axTemp = axList[i][2]
    axTemp.set_title("Mean Station 3 of Stage "+str(i+1))
    alpha3 = np.deg2rad(alpha3_mean)[i]
    c3 = c3_mean[i]
    cw3 = cw3_mean[i]
    
    postStator(alpha3, Ca, c3, cw3, axTemp)

    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    figTemp = figList[i] # <--- Could be weird
    figTemp.tight_layout()
    figTemp.savefig("Mean Stage "+str(i+1))
    plt.close()
    print('stage {0} done'.format(i+1))
print('mean done') 
'''
'''
for i in range(len(alpha1_mean)):
    axTemp = axList[i][0]
    figTemp = figList[i]
    axTemp.set_title("Root Station 1 of Stage "+str(i+1))
    alpha1 = np.deg2rad(alpha1_root)[i]
    beta1 = np.deg2rad(beta1_root)[i]
    c1 = c1_root[i]
    v1 = v1_root[i]
    cw1 = cw1_root[i]
    vw1 = vw1_root[i]

    preRotor(alpha1, beta1, Ca, c1, v1, cw1, vw1, axTemp)
    
    axTemp = axList[i][1]
    axTemp.set_title("Root Station 2 of Stage "+str(i+1))
    alpha2 = np.deg2rad(alpha2_root)[i]
    beta2 = np.deg2rad(beta2_root)[i]
    c2 = c2_root[i]
    v2 = v2_root[i]
    cw2 = cw2_root[i]
    vw2 = vw2_root[i]

    preStator(alpha2, beta2, Ca, c2, v2, cw2, vw2, axTemp)
    
    
    axTemp = axList[i][2]
    axTemp.set_title("Root Station 3 of Stage "+str(i+1))
    alpha3 = np.deg2rad(alpha3_root)[i]
    c3 = c3_root[i]
    cw3 = cw3_root[i]
    
    postStator(alpha3, Ca, c3, cw3, axTemp)

    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    figTemp = figList[i] # <--- Could be weird
    figTemp.tight_layout()
    figTemp.savefig("Root Stage "+str(i+1))
    plt.close()
    print('Root stage {0} done'.format(i+1))
'''

for i in range(len(alpha1_mean)):
    axTemp = axList[i][0]
    figTemp = figList[i]
    axTemp.set_title("Tip Station 1 of Stage "+str(i+1))
    alpha1 = np.deg2rad(alpha1_tip)[i]
    beta1 = np.deg2rad(beta1_tip)[i]
    c1 = c1_tip[i]
    v1 = v1_tip[i]
    cw1 = cw1_tip[i]
    vw1 = vw1_tip[i]

    preRotor(alpha1, beta1, Ca, c1, v1, cw1, vw1, axTemp)
    
    axTemp = axList[i][1]
    axTemp.set_title("Tip Station 2 of Stage "+str(i+1))
    alpha2 = np.deg2rad(alpha2_tip)[i]
    beta2 = np.deg2rad(beta2_tip)[i]
    c2 = c2_tip[i]
    v2 = v2_tip[i]
    cw2 = cw2_tip[i]
    vw2 = vw2_tip[i]

    preStator(alpha2, beta2, Ca, c2, v2, cw2, vw2, axTemp)
    
    
    axTemp = axList[i][2]
    axTemp.set_title("Tip Station 3 of Stage "+str(i+1))
    alpha3 = np.deg2rad(alpha3_tip)[i]
    c3 = c3_tip[i]
    cw3 = cw3_tip[i]
    
    postStator(alpha3, Ca, c3, cw3, axTemp)

    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    figTemp = figList[i] # <--- Could be weird
    figTemp.tight_layout()
    figTemp.savefig("Tip Stage "+str(i+1))
    plt.close()
    print('Tip stage {0} done'.format(i+1))

'''    
print('done')
for i in range(len(alpha1_mean)):
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    figTemp = figList[7-i]
    figTemp.tight_layout()
    figTemp.savefig("Stage "+str(9-i))
    plt.close()
'''

