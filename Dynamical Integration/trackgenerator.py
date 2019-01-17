##this code generated the track

import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl


##track parameters
plateheight = .0127 #m
platelength = 3.81 #m
flatheight = .00102 #m
flatlength = .363 #m
gaplength  = 0.003175 #m
fulltrack = 1264.92 #m

trackratio = 200/332 #first 100 plates are considered, this is enough
tracklength = trackratio * fulltrack

samples = 21*16600 #samples of track,must be high enough that gaplength is accounted for. I prefer at least 3 samples to a gaplength
dx = tracklength / (samples-1)
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



##gaps
gapsamples = int(np.round(gaplength/dx))

if gapsamples%2 == 0:
    a = gapsamples/2
    for i in  range(int(trackratio*166)):
        for j in range(gapsamples):
            track[int(i*platesamples-a+j)] = track[int(i*platesamples-a+j)]-.00003 #it is so shallow because the wheel filters the amplitude
else:
    a = (gapsamples-1)/2
    for i in  range(int(trackratio*332)):
        for j in range(gapsamples):
            b = (i*platesamples-a+j)
            track[int(b)] = track[int(b)]-.00003


##save track to numpy array file, save track metadata
np.save("track", track)


average = np.average(track)
trackmeta = np.zeros(4)
trackmeta[:] = dx, trackratio*332, tracklength, profilesamples #dx, number of plates, lenght of track considered, flatness samples
np.save("trackmeta", trackmeta)



##plotting
fig, ax = plt.subplots()
ax.set(xlabel='Sample', ylabel='Height, m',title='Track')
ax.plot(track)

plt.show()
