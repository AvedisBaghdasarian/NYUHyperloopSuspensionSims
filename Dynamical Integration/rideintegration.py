import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

## pod constants
wheelD = .0762 #m
k = 25000 #N/m
c = 2500 #Ns/m
m = 97.2 #kg
L = 1.524/2 #m
a = L #m
b = .2 #m
I1 = 18 #max moment com is just two point masses, kgm^2
I2 = 2
print(I2, I1)



track2L = np.load("track.npy")
trackmeta = np.load("trackmeta.npy")

racetime = 10.775 #s



##the following codeblock calculates a track for eqch wheel by perturbing the initial track fucntion
deltax, nplates, tracklength, profilesamples= trackmeta[:]
profilesamples = int(profilesamples)
tracksamples = len(track2L)
podsamples = int(np.round(2*a/deltax))

track2sampler = np.arange(tracksamples)-podsamples
trackRsampler = np.arange(tracksamples)-profilesamples

track1L = np.take(track2L,track2sampler , mode = 'wrap')
track2R  = np.take(track2L,trackRsampler , mode = 'wrap')
track1R  = np.take(track1L,trackRsampler , mode = 'wrap')

trackcom = (track1L+track1R+track2L+track2R)
tracktheta1 = (track2L+track2R - track1L-track1R)
tracktheta2 = (track2R+track1R - track2L-track1L)
trackcomdx =np.gradient(trackcom)
tracktheta1dx =np.gradient(tracktheta1)
tracktheta2dx = np.gradient(tracktheta2)



##initial conditions
y0 = np.zeros(8)
y0[2] = trackcom[0]/4
y0[4] = tracktheta1[0]/4
y0[6] = tracktheta2[0]/4



## the following is the RHS of the governing ODE
def dy(t, y): ## this order: x, dx, ycom, dycom, theta1, dtheta1, theta2, dtheta2
    [x, dx, ycom, dycom, theta1, dtheta1, theta2, dtheta2] = y[:]
    positionsample = int(np.round(x/deltax))

    xdot = dx
    #dxdot
    if t <= 4.887:
        dxdot = 27.44 #m/s
    elif t <= 5.887:
        dxdot = 0 #m/s
    else:
        dxdot = -27.44
    #end dxdot
    ycomdot = dycom
    dycomdot = (-k*(4*ycom-trackcom[positionsample]) - c*(4*dycom-dx*trackcomdx[positionsample]))/m
    theta1dot = dtheta1
    dtheta1dot = (a*(-k*(a*4*theta1-tracktheta1[positionsample]) - c*(a*4*dtheta1-dx*tracktheta1dx[positionsample])))/I1
    theta2dot = dtheta2
    dtheta2dot = (b*(-k*(b*4*theta2-tracktheta2[positionsample]) - c*(b*4*dtheta2-dx*tracktheta2dx[positionsample])))/I2

    return [xdot, dxdot, ycomdot, dycomdot, theta1dot, dtheta1dot, theta2dot, dtheta2dot]



##the following block is a quick and dirty optimizer for track oscillations over some variables. The
## code must be manually edited to change the variable that it is opeimized over. Properly, a
##versatile opeimization feature would be created with object oriented code, but there is no reason to
##commit such an effort for the same results
# optpoints = 50
# resulter = np.zeros(optpoints)
# g = np.linspace(401, 2500, optpoints)
# for i in range(optpoints):
#     c = g[i]
#     sol = solve_ivp(dy, [0,8], y0, max_step = 1)
#     x, dx, ycom, dycom, theta1, dtheta1, theta2, dtheta2 = sol.y[:]
#     positionsample = np.rint(x/deltax).astype(int)
#
#     d2ycom = (-k*(4*ycom-trackcom[positionsample]) - c*(4*dycom-dx*trackcomdx[positionsample]))/m
#     d2theta1 = (a*(-k*(a*4*theta1-tracktheta1[positionsample]) - c*(a*4*dtheta1-dx*tracktheta1dx[positionsample])))/I1
#     d2theta2 = (b*(-k*(b*4*theta2-tracktheta2[positionsample]) - c*(b*4*dtheta2-dx*tracktheta2dx[positionsample])))/I2
#
#     resulter[i] = np.amax(ycom+a*theta1+b*theta2-track2R[positionsample])



##the following integrates the pods run
sol = solve_ivp(dy, [0,8], y0, max_step = 1)
x, dx, ycom, dycom, theta1, dtheta1, theta2, dtheta2 = sol.y[:]

positionsample = np.rint(x/deltax).astype(int)



##the following calculates accelerations, as they are not recorded by the the solve_ivp method
d2ycom = (-k*(4*ycom-trackcom[positionsample]) - c*(4*dycom-dx*trackcomdx[positionsample]))/m
d2theta1 = (a*(-k*(a*4*theta1-tracktheta1[positionsample]) - c*(a*4*dtheta1-dx*tracktheta1dx[positionsample])))/I1
d2theta2 = (b*(-k*(b*4*theta2-tracktheta2[positionsample]) - c*(b*4*dtheta2-dx*tracktheta2dx[positionsample])))/I2



#plots
fig, ax = plt.subplots()
ax.set(xlabel='Dampening Constant, Ns/m', ylabel='Height, m',title='Abs Max Height of coupling point of 1 Suspension Assembly vs Suspension distance')

ax.plot(d2theta2)
#ycom+a*theta1+b*theta2-trackcom[positionsample]/4)
plt.show()
