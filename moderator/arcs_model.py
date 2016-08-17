#!/usr/bin/env python

import numpy as np

def IC_atE(E):
    A,B,R,to = computeICParams(E)
    def _(t):
        return IC(t, A,B,R,to)
    return _

def IC(t, A, B, R, to):
    """t: microsecond"""
    r = (1.-R)*A/2.*(A*(t-to*10))**2 * np.exp(-A*(t-to*10)) \
        +R*B*(A/(A-B))**3 *(
            np.exp(-B*(t-to*10)) 
            - np.exp(-A*(t-to*10))*(1+(A-B)*(t-to*10)+0.5*(A-B)**2*(t-to*10)**2)
        )
    r[t<to*10] = 0
    return r

def plot_IC_atE(E):
    params = computeICParams(E)
    t = np.arange(0, 1000., 1)
    I = IC(t, *params)
    from matplotlib import pyplot as plt
    plt.figure()
    plt.plot(t,I)
    plt.show()
    plt.close()
    return

def check_IC_sumrule(E):
    params = computeICParams(E)
    t0 = params[-1]
    # print t0
    t = np.arange(0, 100., 0.1)
    I = IC(t, *params)
    dt = t[1]-t[0]
    assert np.isclose(I.sum()*dt, 1, rtol=2e-2)
    return

def time_distrib_thru_fc(t, t0, E, L, dt_fc, I_E, f_t_E):
    """
    t: microsecond
    t0: microsecond
    E: eV
    L: meter
    dt_fc: microsecond
    t: microsecond
    """
    from mcni.utils import conversion as Conv
    v = Conv.e2v(E*1000)
    t_fc = L/v + t0/1e6
    t_travel = t_fc - t/1e6 # second
    v1 = L/t_travel # m/s
    E1 = Conv.VS2E * v1*v1 # in meV
    E1 /= 1000 # in eV
    res = I_E(E1) * 2*E1*dt_fc*1.e-6/t_travel*f_t_E(E1)(t)
    res[t_travel<0] = 0
    return res

bl17_IE_params = -0.00268304038288988,0.863669894898031,100.230838479278,-3.99754529853209,200.040128814641,1.2365793368859,350.02744912904,0.0346067554662553,0.953181679989923,0.025,2.19241457833097,4.44502602237474,49.4568791178012,0.67954001944297
def bl17_I_E(E):
    c, R1, T1, R2, T2, R3, T3, a, b, Ecut, s, delta, g, I = bl17_IE_params
    return I_E(E, c, R1, T1, R2, T2, R3, T3, a, b, Ecut, s, delta, g, I)

def I_E(E, c, R1, T1, R2, T2, R3, T3, a, b, Ecut, s, delta, g, I):
    from numpy import exp, sqrt, power
    k = 8.617e-5 # eV/K
    B = 7.36e-3 # eV
    D = lambda E: 1/(1+(Ecut/E)**s)
    x = g*(E-2*B)
    x[E<2*B] = 0
    rho = 1 + delta*exp(-x)*(1 + x +0.5*x**2)
    # impl copied from SNS_source_analytic component
    arg1 = I*1.0e12 * exp(-c/sqrt(E));
    arg2 = R1*E/power(k*T1,2) *exp(-E/(k*T1));
    arg3 = R2*E/power(k*T2,2) *exp(-E/(k*T2));
    arg4 = R3*E/power(k*T3,2) *exp(-power(E/(k*T3),b));
    arg5 = D(E)*rho/power(E,1-a);
    arg6 =(arg1 * ( arg2 + arg3 + arg4 + arg5 ));
    return arg6
            
def computeICParams(E):
    "compute parameters in Ikeda Carpenter function for given E (eV)"
    return [func(E) for func in funcs]

# this function
def f(x, a,b,c,d,f,g,h,i,j,k):
    return a*x**b*(1+c*x+d*x**2+(x/f)**g)/(1+h*x+i*x**2+(x/j)**k)

# and these parameters are from a1Gw2-11-f5_fit_fit.dat
A = "1.51546780992587 0.431080986615399 2288.09554767074 20216.7410430439 3.26599271523613 -0.75268083061942 3598.57437439047 15613.751735194 23.5079799676861 -0.500561774397004"
B = "0.204226546952013 0.422693595124972 -592.788347939597 19205.2500865301 1.70107957244692 -2.45893589986053 -1350.11529608879 16837.2873626292 4.13403206659535 -2.05426275276519"
R = "0.946792303879114  0.0 -6.55236227284847  -49.8115885416586  0.129646355044244  2.27663310920202  2.62171627724354  -335.009222176655  0.0647089854557834  2.27663310920202"
to = "107.196504837589  0.177528212550242  -4859.15231611581  39474.0031453964  0.00021045946207969  1.01864575152697  22314.3619803924  138.991460677331  0.00096235218431338  2.6059824501576"


def build_func(param_str):
    params = map(float, param_str.split())
    return lambda x: f(x, *params)

funcs = map(build_func, (A,B,R,to))

def test_computeICParams():
    print computeICParams(0.70795)
    print computeICParams(0.63096)
    print computeICParams(0.686)
    return

def test_IC():
    plot_IC_atE(100*1e-3)
    # plot_IC_atE(700*1e-3)
    # check_IC_sumrule(100*1e-3)
    # check_IC_sumrule(700*1e-3)
    return

def test_I_E():
    E = np.arange(0, 1, 0.001)
    I = bl17_I_E(E)
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
    I = time_distrib_thru_fc(t, t0, E, L, dt_fc, bl17_I_E, IC_atE)
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
    test_tdtf()
    return

if __name__ == '__main__': main()
