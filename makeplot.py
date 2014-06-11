import pylab as plt
import numpy as np

def wake(opt):
    plt.figure()
    plt.title("Free Wake")
    r = opt.aso.aero2.induced.r
    z = opt.aso.aero2.induced.z
    plt.plot(r,z, "b-")
    plt.plot(r,z, "bo")
    plt.plot(-r,z, "b-")
    plt.plot(-r,z, "bo")

    plt.xlabel("r(m)")

def aerodef(opt):
    plt.figure()
    yE = opt.aso.discrete.yE
    cE = opt.aso.discrete.cE
    Cl = opt.aso.discrete.Cl
    cd = opt.aso.aero2.Cd
    re = opt.aso.aero2.Re
    R = opt.aso.aero2.R
    y100 = np.linspace(0,R,100)
    c100 = opt.aso.discrete.c100

    plt.plot(yE, cE, label='Chord (m)')
    plt.plot(yE, Cl, label='C_l')
    plt.plot(yE, 100*cd, label='C_d (10^{-2})')
    plt.plot(yE, Cl/cd/100., label='L/D (10^{2})')
    plt.plot(yE, re/1000000.,label='Re (10^6)')
    plt.plot(y100, c100)
    plt.xlabel("r(m)")
    plt.legend()



def plot(opt):
    wake(opt)
    aerodef(opt)