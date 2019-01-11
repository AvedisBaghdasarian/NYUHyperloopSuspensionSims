import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


wheelD = .0762 #m
k = 16000 ##N/m
c = 1500 #Ns/m
m = 100 #kg
L = 1.524/2 #m
a = L #m
b = .15 #m
I1 = m*(L)**2 #max moment com is just two point masses, kgm^2
I2 = m * (.15)**2
print(I2, I1)

track2L = np.load("track.npy")
trackmeta = np.load("trackmeta.npy")


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
# wheelsamplespacing = 15
# wheelsamples = int(np.round(wheelD/wheelsamplespacing/dx))
# if wheelsamples%2 == 0:
#     halfwheelsamples = wheelsamplespacing*wheelsamples/2
#     wheelsampler = np.round(np.linspace(-halfwheelsamples, halfwheelsamples, wheelsamples))
# else:
#     halfwheelsamples = wheelsamplespacing*(wheelsamples-1)/2
#     wheelsampler = np.round(np.linspace(-halfwheelsamples, halfwheelsamples, wheelsamples))

racetime = 10.775 #s

y0 = np.zeros(8)
y0[2] = trackcom[0]
y0[4] = tracktheta1[0]
y0[6] = tracktheta2[0]

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

sol = solve_ivp(dy, [0,8], y0, max_step = 1)


x, dx, ycom, dycom, theta1, dtheta1, theta2, dtheta2 = sol.y[:]
positionsample = np.rint(x/deltax).astype(int)

d2ycom = (-k*(4*ycom-trackcom[positionsample]) - c*(4*dycom-dx*trackcomdx[positionsample]))/m
d2theta1 = (a*(-k*(a*4*theta1-tracktheta1[positionsample]) - c*(a*4*dtheta1-dx*tracktheta1dx[positionsample])))/I1
d2theta2 = (b*(-k*(b*4*theta2-tracktheta2[positionsample]) - c*(b*4*dtheta2-dx*tracktheta2dx[positionsample])))/I2

fig, ax = plt.subplots()

ax.plot(sol.t, d2theta2)

plt.show()
