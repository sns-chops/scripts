#!/usr/bin/env python

import numpy as np

from fitfit_model import Instrument, time_distrib_thru_fc
class ARCS(Instrument):
    # this
    IE_params = -0.00268304038288988,0.863669894898031,100.230838479278,-3.99754529853209,200.040128814641,1.2365793368859,350.02744912904,0.0346067554662553,0.953181679989923,0.025,2.19241457833097,4.44502602237474,49.4568791178012,0.67954001944297

    # and these parameters are from a1Gw2-11-f5_fit_fit.dat
    A = "1.51546780992587 0.431080986615399 2288.09554767074 20216.7410430439 3.26599271523613 -0.75268083061942 3598.57437439047 15613.751735194 23.5079799676861 -0.500561774397004"
    B = "0.204226546952013 0.422693595124972 -592.788347939597 19205.2500865301 1.70107957244692 -2.45893589986053 -1350.11529608879 16837.2873626292 4.13403206659535 -2.05426275276519"
    R = "0.946792303879114  0.0 -6.55236227284847  -49.8115885416586  0.129646355044244  2.27663310920202  2.62171627724354  -335.009222176655  0.0647089854557834  2.27663310920202"
    to = "107.196504837589  0.177528212550242  -4859.15231611581  39474.0031453964  0.00021045946207969  1.01864575152697  22314.3619803924  138.991460677331  0.00096235218431338  2.6059824501576"
    It_params = A,B,R,to

model = ARCS()

def test_computeICParams():
    print model.computeICParams(0.70795)
    print model.computeICParams(0.63096)
    print model.computeICParams(0.686)
    return

def test_IC():
    # model.plot_IC_atE(100*1e-3)
    # model.plot_IC_atE(700*1e-3)
    model.check_IC_sumrule(100*1e-3)
    model.check_IC_sumrule(700*1e-3)
    return

def test_I_E():
    E = np.arange(0, 1, 0.001)
    I = model.I_E(E)
    from matplotlib import pyplot as plt
    plt.figure()
    plt.plot(E,I)
    plt.show()
    plt.close()
    return

def test_tdtf():
    L = 11.61
    dt_fc = 1.0 # microsecond
    E = 100*1e-3
    t = np.arange(0, 1000., 1)
    t0 = 10.
    I = time_distrib_thru_fc(t, t0, E, L, dt_fc, model.I_E, model.IC_atE)
    from matplotlib import pyplot as plt
    plt.figure()
    plt.plot(t,I)
    plt.show()
    plt.close()
    return
                
def test_tdtf2():
    E = 100*1e-3
    t = np.arange(0, 1000., 1)
    t0 = 10.
    from fermichopper import FermiChopper as base
    class FC(base):
        frequency = 600
        slit_thickness = 1.5
        radius = 50
    I = model.time_distrib_thru_fc(t, t0, E, FC())
    from matplotlib import pyplot as plt
    plt.figure()
    plt.plot(t,I)
    plt.show()
    plt.close()
    return

def main():
    # test_computeICParams()
    # test_IC()
    # test_I_E()
    # test_tdtf()
    test_tdtf2()
    return

if __name__ == '__main__': main()
