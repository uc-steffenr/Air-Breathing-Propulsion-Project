import numpy as np
import matplotlib.pyplot as plt
from trianglesolver import solve

# Create Class to solve velocity triangles
class VelocityTriangle:
    def __init__(self, alpha=0, beta=0, Ca=227.5):
        # Initialize all Triangle Values to be NoneType
        self.beta  = beta
        self.alpha = alpha
        self.Ca    = Ca

        self.Vw    = Ca*np.tan(beta)
        self.Cw    = Ca*np.tan(alpha)
        self.V     = np.sqrt(self.Vw**2 + Ca**2)
        self.C     = np.sqrt(self.Cw**2 + Ca**2)

    ############################################################### 
    def DrawCompTri(self):
        print('don\'t use this yet')
        
       
# Test 1
'''alpha1 = 30*np.pi/180
Ca1 = 180
beta1 = 58.7*np.pi/180

triangle1 = VelocityTriangle()
triangle1.Solve(alpha=alpha1, Ca=Ca1, beta=beta1)
print('Cw = ',triangle1.Cw)
print('C = ',triangle1.C)
print('V = ',triangle1.V)
print('Vw = ',triangle1.Vw)
# results: Test passes
'''
'''
# Test 2
C2 = 614.21
beta2 = 34.23 * np.pi/180
alpha2 = 65 * np.pi/180

triangle2 = VelocityTriangle()
triangle2.Solve(C=C2, alpha=alpha2, beta=beta2, )
# Results:
# triangle solver doesn't like this because in my code, I don't account for the common value of Ca between the two
# I could do something to make this more modular and fix this, but test 1 works and that is good enough for this project 
'''

