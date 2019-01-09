import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#track constants
plateheight = .0127 #m
platelength = 3.81 #m
kappa = 2*np.pi/(2*platelength) #1/m
A = plateheight/2 #m
vmax = 134.112 #m/s


#pod constants
kmax = 500000 ##N/m
cs = 1000 #Ns/m
m = 100 #kg
l = 1.524/2 #m
L = l
a = l #m
I = m*(l)**2 #max moment com is just two point masses, kgm^2
b = .15 #m
I2 = m * (.15)**2

print(I)
# Ftheta = 2*a*A*np.sin(kappa*L)*np.sqrt(k**2+(c*w)**2) #N


fig, ax = plt.subplots()
zetapoints = 3
zeta = np.linspace(0.25, 1, zetapoints)
npoints = 500

#theta function of k
ycom = np.zeros([zetapoints, npoints])
acom = np.zeros([zetapoints, npoints])
theta = np.zeros([zetapoints, npoints])
alphatheta = np.zeros([zetapoints, npoints])
theta2 = np.zeros([zetapoints, npoints])
alphatheta2 = np.zeros([zetapoints, npoints])
ki = np.linspace(0, 33000, npoints)
vi = np.linspace(.1, vmax, npoints)
w = kappa*vi #1/s
for j in range(zetapoints):
    for i in range(npoints):
        k = ki[i]
        w_0 = a*np.sqrt(2*k/I)
        c = zeta[j]*a*I*w_0
        # Ftheta = 2*a*A*np.sin(kappa*L)*np.sqrt(k**2+(c*w)**2)
        Ftheta = 2*a*A*1*np.sqrt(k**2+(c*w)**2)
        Ztheta = np.sqrt((2*w_0*zeta[j])**2 + 1/w**2 * (w_0**2-w**2)**2) #N
        thetaw = Ftheta/I/Ztheta/w #amplitude of theta as a function of w
        alphathetaw = thetaw * w * w #amplitude of angular accel as a function of w
        theta[j, i] = np.amax(thetaw) #max response for any w
        alphatheta[j, i] = np.amax(alphathetaw) #max alpha for any w

    for i in range(npoints):
        k = ki[i]
        w_0 = np.sqrt(2*k/m)
        c = zeta[j]*m*w_0
        # Fcom = 2*A*np.cos(kappa*L)*np.sqrt(k**2+(c*w)**2)
        Fcom = 2*A*1*np.sqrt(k**2+(c*w)**2)
        Zcom = np.sqrt((2*w_0*zeta[j])**2 + 1/w**2 * (w_0**2-w**2)**2) #N
        ycomw = Fcom/m/Zcom/w #amplitude of response as function of w
        acomw = ycomw * w * w #amplitude of accel as a function of w
        ycom[j, i] = np.amax(ycomw) #max response for any w
        acom[j, i] = np.amax(acomw) #max accel for any w

    for i in range(npoints):
        k = ki[i]
        w_0 = b*np.sqrt(2*k/I2)
        c = zeta[j]*b*I2*w_0
        # Ftheta = 2*a*A*np.sin(kappa*L)*np.sqrt(k**2+(c*w)**2)
        Ftheta2 =2*b*A*np.sqrt(k**2+(c*w)**2)
        Ztheta2 = np.sqrt((2*w_0*zeta[j])**2 + 1/w**2 * (w_0**2-w**2)**2) #N
        thetaw2= Ftheta2/I2/Ztheta2/w #amplitude of theta as a function of w
        alphathetaw2 = thetaw2 * w * w #amplitude of angular accel as a function of w
        theta2[j, i] = np.amax(thetaw2) #max response for any w
        alphatheta2[j, i] = np.amax(alphathetaw2) #max alpha for any w


ytheta = a*theta
ytheta2 = b*theta2
accthetaY2 = a*alphatheta
accY2 = acctheta+acccom
DY2 = ytheta+ycom
DY1 = ycom-ytheta
yeet = ytheta + ytheta2 + ycom

for k in range(zetapoints):
#     ax.set(xlabel='omega/omega_0', ylabel='acceleration (m/s/s)',
#            title='accel vs angular number')
    # ax.plot(y[k, :])
    ax.plot(ki, ycom[k, :])
plt.show()
