import cmath
import math
import numpy as np
import matplotlib.pyplot as plt

def bisection(trans_mat, x_min, x_max, N=10000, tol=1e-10, max_iter=100000):
    split = (x_max-x_min)/N
    ans = []
    for i in range(N):
        # initialize the values of x and the iteration counter
        x_i = x_min + float(i)*split
        x_ii = x_min + float(i+1)*split
        if trans_mat(x_i).imag != 0 or trans_mat(x_ii).imag != 0:
            continue
        elif trans_mat(x_i).real* trans_mat(x_ii).real > 0:
            continue
        elif trans_mat(x_i).real* trans_mat(x_ii).real < 0:
            x = (x_i + x_ii) / 2
            n_iter = 0
        # loop until the root is found or the maximum number of iterations is reached
            while abs(trans_mat(x)) > tol and n_iter < max_iter:
                #print(x)
                if trans_mat(x_i).real * trans_mat(x).real < 0:
                    x_ii = x
                elif trans_mat(x_i).real * trans_mat(x).real > 0:
                    x_i = x
                x = (x_i + x_ii) / 2
                n_iter += 1
        # check if the maximum number of iterations was reached
            if n_iter == max_iter:
                print("Maximum number of iterations reached.")
            ans.append(x)
    return ans

def trans_mat(n_eff):
    global k
    global N
    global X
    
    M = np.array([[1,0],
                  [0,1]])
    for i in range(1,len(X)):
        g_i = cmath.sqrt(((n_eff*k)**2-(N[i]*k)**2))
        if g_i == 0:
            continue
        w = X[i] - X[i-1]
        m = np.array([[cmath.cosh(g_i*w),cmath.sinh(g_i*w)/(1j*g_i)]
                      ,[1j*g_i*cmath.sinh(g_i*w),cmath.cosh(g_i*w)]])
        M = np.dot(m,M)
    g_N_1 = cmath.sqrt((n_eff*k)**2-(N[len(N)-1]*k)**2)
    g_0 = cmath.sqrt((n_eff*k)**2-(N[0]*k)**2)
    trans_mat = (1j*(M[0][0]*g_N_1 + M[1][1]*g_0) + M[1][0] - M[0][1]*g_N_1*g_0).imag
    return trans_mat

def heaviside(x):
    return np.heaviside(x, 0)


l = 1.55 #um (wavelength)
k = 2*math.pi/l


N = [1.445, 3.47, 1.445]
N_max = 0
N_min = 100

for i in range(len(N)):
    if N[i] > N_max:
        N_max = N[i]
    
    if N[i] < N_min:
        N_min = N[i]
num = 50
for i in range(num):
    x = 0.05*i
    X = [0, x]
    root = bisection(trans_mat, N_min+1e-9, N_max, 5000)
    for j in range(len(root)):
        plt.plot(x,root[j], "-or")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
