##this code generated the trackfunction

import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

plateheight = .0127 #m
platelength = 3.81 #m
flatheight = .00102 #m
flatlength = .363 #m
gaplength  = 0.003175 #m
fulltrack = 1264.92 #m

trackratio = 50/166 #first 100 plates are considered, this is enough
tracklength = trackratio * fulltrack

samples = 21*16600 #samples of track,must be high enough that gaplength is accounted for. I prefer at least 3 samples to a gaplength
dx = tracklength / samples
print(dx)
track = np.zeros(samples)

##initial track plate model
platesamples = int(np.round(platelength/dx))

track[0:platesamples] = np.linspace(0, plateheight, platesamples)
track[platesamples:2*platesamples] = np.linspace(plateheight, 0, platesamples)

track[:] = track[(np.array(range(samples))%(2*platesamples))]

##track flatness profile
profilesamples = int(np.round(flatlength/dx)) #it would be better to round during evaluation rather than rounding the number of bins then linspace, but whatevs
profile = np.zeros(2*profilesamples)

profile[0:profilesamples] = np.linspace(0, flatheight, profilesamples)
profile[profilesamples:2*profilesamples] = np.linspace(flatheight, 0, profilesamples)

track[:] = track[:] + profile[(np.array(range(samples))%(2*profilesamples))]

gapsamples = int(np.round(gaplength/dx))
print(gapsamples)
if gapsamples%2 == 0:
    a = gapsamples/2
    for i in  range(int(trackratio*166)):
        for j in range(gapsamples):
            track[int(i*platesamples-a+j)] = 0

else:
    a = (gapsamples-1)/2
    for i in  range(int(trackratio*332)):
        for j in range(gapsamples):
            b = (i*platesamples-a+j)
            track[int(b)] = 0


fig, ax = plt.subplots()
# for k in range(zetapoints):
#     ax.set(xlabel='omega/omega_0', ylabel='acceleration (m/s/s)',
#            title='accel vs angular number')
#     ax.plot(x, acc[k,:])
ax.plot(track)

# plt.show()
