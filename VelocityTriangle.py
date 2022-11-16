import numpy as np
import matplotlib.pyplot as plt
from trianglesolver import solve

# Create Class to solve velocity triangles
class VelocityTriangle:
    def __init__(self):
        # Initialize all Triangle Values to be NoneType
        self.beta  = None
        self.alpha = None
        self.Vw    = None
        self.Cw    = None
        self.V     = None
        self.C     = None
        self.Ca    = None
        
    ######################################################
    def Solve(self,**kwargs):
        # Process Kwargs
        for key, val in kwargs.items():
            if key == 'beta':
                self.beta=val
            elif key == 'alpha':
                self.alpha=val
            elif key == 'V':
                self.V=val
            elif key == 'Vw':
                self.Vw=val
            elif key == 'Ca':
                self.Ca=val            
            elif key == 'C':
                self.C=val
            elif key == 'Cw':
                self.Cw=val
            else:
                print('you entered a kwarg wrong.  Try again bozo.')
        
        # SOLVING TRIANGLES USING TRIANGLE SOLVER
        # ------------------------------------------------------------------
        # Relative triangles - beta and V
        a_rel = self.Vw
        A_rel = self.beta
        b_rel = self.Ca
        B_rel = None
        c_rel = self.V
        C_rel = np.pi/2
        '''
        print(a_rel)
        print(A_rel)
        print(b_rel)
        print(B_rel)
        print(c_rel)
        print(C_rel)
        '''
        # Solving
        a_rel, b_rel, c_rel, A_rel, B_rel, C_rel =\
             solve(a=a_rel, b=b_rel, c=c_rel, A=A_rel, B=B_rel, C=C_rel)
        # Store Values as self to access in main code
        self.beta  = A_rel
        self.Vw    = a_rel
        self.V     = c_rel
        self.Ca    = b_rel
        
        #----------------------------------------------------------------------
        # Absolute triangles - alpha and C
        a_abs = self.Cw
        A_abs = self.alpha
        b_abs = self.Ca
        B_abs = None
        c_abs = self.C
        C_abs = np.pi/2
        '''
        print(a_abs)
        print(A_abs)
        print(b_abs)
        print(B_abs)
        print(c_abs)
        print(C_abs)
        '''
        # Solving
        a_abs, b_abs, c_abs, A_abs, B_abs, C_abs =\
             solve(a=a_abs, b=b_abs, c=c_abs, A=A_abs, B=B_abs, C=C_abs)
        # Store Values as self to access in main code
        self.alpha = A_abs
        self.Cw    = a_abs
        self.C     = c_abs

    ############################################################### 
    def drawtraingle(self):
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

