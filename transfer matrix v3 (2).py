import cmath
import math
import numpy as np
import matplotlib.pyplot as plt

def bisection(trans_mat, x_min, x_max, N=10000, tol=1e-10, max_iter=1000):
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
        #print(g_i)
        if g_i == 0:
            continue
        w = X[i] - X[i-1]
        m = np.array([[cmath.cosh(g_i*w),cmath.sinh(g_i*w)/(1j*g_i)]
                      ,[1j*g_i*cmath.sinh(g_i*w),cmath.cosh(g_i*w)]])
        #print(m)
        M = np.dot(m,M)
    g_N_1 = cmath.sqrt((n_eff*k)**2-(N[len(N)-1]*k)**2)
    g_0 = cmath.sqrt((n_eff*k)**2-(N[0]*k)**2)
    trans_mat = (1j*(M[0][0]*g_N_1 + M[1][1]*g_0) + M[1][0] - M[0][1]*g_N_1*g_0).imag
    return trans_mat

l = 1.55 #um (wavelength)
k = 2*math.pi/l


N = [1.445, 1.45, 1.445, 1.45, 1.445]
#N = [1.445, 3, 4, 1.45]
#N = [1.44, 3.47, 1.44]
X = [0,5,10,15]
#X = [0, 0.5,1]
#X = [0,0.25]
N_max = 0
N_min = 100

for i in range(len(N)):
    if N[i] > N_max:
        N_max = N[i]
    
    if N[i] < N_min:
        N_min = N[i]
X_p = (X[len(X)-1] - X[0])/2
    
root = bisection(trans_mat, N_min+1e-9, N_max, 500)
print(root)




fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

for i in range(len(root)):
    g_o = cmath.sqrt((root[i]*k)**2-(N[0]*k)**2)
    g_N_1 = cmath.sqrt((root[i]*k)**2-(N[len(N)-1]*k)**2)
    EH = np.array([[1],[g_o]])
    x_o = np.linspace(-X_p,X[0],500)
    E_o = np.exp(g_o.real*(x_o-X[0]))
    N_o = N[0] *x_o/x_o
    ax1.plot(x_o, E_o,label=f"n={root[i]}", color=f"C{i}")
    ax2.plot(x_o, N_o, color='green', linestyle='--')

    M = np.array([[1,0],
                  [0,1]])
    for j in range(len(X)-1):
        x_i = np.linspace(X[j],X[j+1], 500)
        g_i = cmath.sqrt((N[j+1]*k)**2-(root[i]*k)**2)
        N_i = N[j+1] *x_i/x_i
        if g_i.imag == 0:
            E_i = EH[0][0].real * np.cos(g_i.real*(x_i-X[j])) + EH[1][0].real * np.sin(g_i.real*(x_i-X[j]))/g_i.real
        else:
            E_i = EH[0][0].real * np.cosh(g_i.imag*(x_i-X[j])) + EH[1][0].real * np.sinh(g_i.imag*(x_i-X[j]))/g_i.imag
        ax1.plot(x_i, E_i, color=f"C{i}")
        ax2.plot(x_i, N_i, color='green', linestyle='--')
        w = X[j+1] - X[j]
        if g_i == 0:
            continue
        
        m = np.array([[cmath.cos(g_i*w),cmath.sin(g_i*w)/g_i]
                    ,[-g_i*cmath.sin(g_i*w),cmath.cos(g_i*w)]])
        EH = np.dot(m,EH)
        
        if j == len(X)-2:
            x_N_1 = np.linspace(X[j+1],X[j+1]+X_p,500)
            y = EH[0][0].real * np.exp(-g_N_1.real*(x_N_1 -X[len(X)-1]))
            ax1.plot(x_N_1, y, color=f"C{i}")
            N_i = N[j+2] *x_i/x_i
            ax2.plot(x_N_1, N_i, color='green', linestyle='--')

ax1.set_xlabel("x")
ax1.set_ylabel("Ey(x)")
ax2.set_ylabel('n')
ax2.tick_params(axis='y', labelcolor='green')
plt.title("Mode distribution of Ey(x)")
plt.show()

 
