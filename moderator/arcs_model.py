#!/usr/bin/env python

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

def computeICParams(E):
    "compute parameters in Ikeda Carpenter function for given E (eV)"
    return [func(E) for func in funcs]

def test_computeICParams():
    print computeICParams(0.70795)
    print computeICParams(0.63096)
    print computeICParams(0.686)
    return

def main():
    test_computeICParams()
    return

if __name__ == '__main__': main()