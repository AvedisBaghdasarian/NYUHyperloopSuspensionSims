import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

plateheight = .0127 #m
platelength = 3.81 #m
flatheight = .00102 #m
flatlength = .363 #m

tracklength = 1264.92#m

samples = 16600 #samples of track
dx = tracklength / samples

track = np.zeros(samples)
platesamples = np.round(2*platelength/dx)
print(platesamples, 2*platelength/dx, (tracklength/2/platelength))
