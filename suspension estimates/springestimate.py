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


fig, ax = plt.subplots()
zetapoints = 3
zetas = np.linspace(0.25, 1, zetapoints)
npoints = 500

springconstants = np.linspace(0, 33000, npoints)
vi = np.linspace(.1, vmax, npoints)
frequencies = kappa*vi #1/s

#estimates vibrations given a set of spring constants, set of zeta, set of oscillations, inertial term, distance, iterations)
def estimate_vib(ki, zeta, w, m, a, npoints):
    zetapoints = len(zeta)
    response = np.zeros([zetapoints, npoints])
    aresponse = np.zeros([zetapoints, npoints])

    for j in range(zetapoints):
        for i in range(npoints):
            k = ki[i]
            w_0 = a*np.sqrt(2*k/I)
            c = zeta[j]*a*m*w_0
            F = 2*a*A*1*np.sqrt(k**2+(c*w)**2)
            Z = np.sqrt((2*w_0*zeta[j])**2 + 1/w**2 * (w_0**2-w**2)**2) #N
            responsew = F/I/Z/w #amplitude of theta as a function of w
            aresponsew = responsew * w * w #amplitude of angular accel as a function of w
            response[j, i] = np.amax(responsew) #max response for any w
            aresponse[j, i] = np.amax(aresponsew) #max alpha for any w

    return response, aresponse

theta, alphatheta = estimate_vib(springconstants, zetas, frequencies, I, a, npoints)
com, acccom = estimate_vib(springconstants, zetas, frequencies, m, 1, npoints)
theta2, alphatheta2 = estimate_vib(springconstants, zetas, frequencies, I2, b, npoints)

# calculation of locations of interest
ytheta = a*theta
ytheta2 = b*theta2
accthetaY2 = a*alphatheta
accY2 = accthetaY2+acccom
DY2 = ytheta+ycom
DY1 = ycom-ytheta
yeet = ytheta + ytheta2 + ycom

for k in range(zetapoints):
    ax.set(xlabel='omega/omega_0', ylabel='acceleration (m/s/s)',
           title='accel vs angular number')
    ax.plot(springconstants, yeet[k, :])
plt.show()
