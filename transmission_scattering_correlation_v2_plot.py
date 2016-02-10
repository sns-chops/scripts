#!/usr/bin/env python

asciiout = "corr.dat"


import numpy as np, pylab

angle, mon1_v, mon1_err, mon2_v, mon2_err, alldet_v, alldet_err, maskeddet_v, maskeddet_err \
    = np.loadtxt(asciiout).T

# pylab.errorbar(angle, mon1_v, yerr=mon1_err)
# pylab.errorbar(angle, mon2_v, yerr=mon2_err)
# pylab.errorbar(angle, alldet_v, yerr=alldet_err)
pylab.errorbar(angle, maskeddet_v, yerr=maskeddet_err)
pylab.show()
