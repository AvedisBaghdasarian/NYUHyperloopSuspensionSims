import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#track constants
plateheight = .0127 #m
platelength = 3.81 #m
kappa = 2*np.pi/(2*platelength) #1/m
A = plateheight/2 #m
vmax = 134.112/50 #m/s



#pod constants
kmax = 500000 ##N/m
cs = 1000 #Ns/m
m = 150 #kg
l = 1.524/2 #m
L = l
a = 1.524/2 #m
I = m*(l/2)**2 #max moment com is just two point masses, kgm^2


# Ftheta = 2*a*A*np.sin(kappa*L)*np.sqrt(k**2+(c*w)**2) #N

#com
# Res_com  = Fcom / m / Zcom / w

# angle oscillations

fig, ax = plt.subplots()
zetapoints = 3
zeta = np.linspace(0.25, 1, zetapoints)
npoints = 500

#theta function of k
ycom = np.zeros([zetapoints, npoints])
acom = np.zeros([zetapoints, npoints])
theta = np.zeros([zetapoints, npoints])
alphatheta = np.zeros([zetapoints, npoints])
ki = np.linspace(0, 5000, npoints)
vi = np.linspace(.1, vmax, npoints)
w = kappa*vi #1/s
for j in range(zetapoints):
    for i in range(npoints):
        k = ki[i]
        w_0 = a*np.sqrt(2*ki[i]/m)
        c = zeta[j]*a*I*w_0
        # Ftheta = 2*a*A*np.sin(kappa*L)*np.sqrt(k**2+(c*w)**2)
        Ftheta = 2*a*A*1*np.sqrt(k**2+(c*w)**2)
        Ztheta = np.sqrt((2*w_0*zeta[j])**2 + 1/w**2 * (w_0**2-w**2)**2) #N
        thetaw = Ftheta/m/Ztheta/w #amplitude of theta as a function of w
        alphathetaw = thetaw * w * w #amplitude of angular accel as a function of w
        theta[j, i] = np.amax(thetaw) #max response for any w
        alphatheta[j, i] = np.amax(alphathetaw) #max alpha for any w

    for i in range(npoints):
        k = ki[i]
        w_0 = np.sqrt(2*ki[i]/m)
        c = zeta[j]*m*w_0
        # Fcom = 2*A*np.cos(kappa*L)*np.sqrt(k**2+(c*w)**2)
        Fcom = 2*A*1*np.sqrt(k**2+(c*w)**2)
        Zcom = np.sqrt((2*w_0*zeta[j])**2 + 1/w**2 * (w_0**2-w**2)**2) #N
        ycomw = Fcom/m/Zcom/w #amplitude of response as function of w
        acomw = ycomw * w * w #amplitude of accel as a function of w
        ycom[j, i] = np.amax(ycomw) #max response for any w
        acom[j, i] = np.amax(acomw) #max accel for any w
alpha = w*w*theta
acccom = w*w*ycom #only makes sense on constant w
ytheta = a*theta
acctheta = a*alpha
accY2 = acctheta+acccom
DY2 = ytheta+ycom



# for k in range(zetapoints):
#     ax.set(xlabel='omega/omega_0', ylabel='acceleration (m/s/s)',
#            title='accel vs angular number')
    # ax.plot(y[k, :])
    # ax.plot(ki, ycom[k, :])
plt.show()
