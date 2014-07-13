import pylab as plt
import numpy as np

pi = np.pi

def wake(opt):
    plt.title("Free Wake")
    r = opt.aso.aero2.induced.r
    z = opt.aso.aero2.induced.z
    plt.plot(r,z, "b-")
    plt.plot(r,z, "bo")
    plt.plot(-r,z, "b-")
    plt.plot(-r,z, "bo")

    plt.xlabel("r(m)")

def aero_def(opt):
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
    plt.plot(yE, 100*cd, label='C_d ($\mathregular{10^{-2}}$)')
    plt.plot(yE, Cl/cd/100., label='L/D ($\mathregular{10^{2}}$)')
    plt.plot(yE, re/1000000.,label='Re ($\mathregular{10^{6}}$)')
    plt.plot(y100, c100)
    plt.xlabel("r(m)")
    plt.title('Aerodynamic Definition')
    plt.legend()

def struct_def(opt):
    yE = opt.aso.discrete.yE
    d = opt.aso.aero2.d
    theta = opt.aso.struc.theta
    nTube = opt.aso.struc.nTube
    nCap = opt.aso.struc.nCap
    lBiscuit = opt.aso.struc.lBiscuit
    yWire = opt.aso.struc.yWire
    zWire = opt.aso.struc.zWire
    dQuad = opt.aso.struc.dQuad
    thetaQuad = opt.aso.struc.thetaQuad
    nTubeQuad = opt.aso.struc.nTubeQuad
    lBiscuitQuad = opt.aso.struc.lBiscuitQuad
    hQuad = opt.aso.struc.hQuad

    plt.title('Structural Definition')
    plt.plot(yE,d*100/2.54, label='Spar diameter (in)')
    plt.plot(yE,theta*180/pi/10, label='Wrap angle ($\mathregular{10^{1}}$ deg)')
    plt.plot(yE,nTube, label='Tube layers')
    plt.plot(yE,nCap, label='Cap layers')
    plt.plot(yE,lBiscuit*100/2.54/12, label='Biscuit spacing (ft)')
    plt.plot(yWire,zWire, "o", color='red', label='Wire location (m)', linewidth=3, markersize=9., fillstyle='none')
    plt.plot(0,dQuad*100/2.54, "s", color='black', clip_on=False, label='Quad definition', linewidth=3, fillstyle='none', markersize=9.)

    plt.xlabel("r(m)")
    plt.legend(numpoints=1)

def struc_prop(opt):
    yE = opt.aso.discrete.yE
    EA = opt.aso.struc.spar.EA
    EIz = opt.aso.struc.spar.EIx
    EIx = opt.aso.struc.spar.EIx
    GJ = opt.aso.struc.spar.GJ
    mSpar = opt.aso.struc.spar.mSpar
    mChord = opt.aso.struc.chord.mChord
    EIQuad = opt.aso.struc.quad.EIx
    GJQuad = opt.aso.struc.quad.GJ
    dy = opt.aso.struc.spar.dy

    plt.plot(yE,EA/(10**7), label='EA ($\mathregular{10^{7}}$)')
    plt.plot(yE,EIz/(10**4), label='EIz ($\mathregular{10^{4}}$)')
    plt.plot(yE,EIx/(10**4), label='EIx ($\mathregular{10^{4}}$)')
    plt.plot(yE,GJ/(10**4), label='GJ ($\mathregular{10^{4}}$)')
    plt.plot(yE,mSpar/dy, label='mSpar (kg/m)')
    plt.plot(yE,mChord/dy,'k--', label='mChord (kg/m)')
    plt.plot(0,EIQuad/(10**4),'s', clip_on=False, color='blue', label='EIquad ($\mathregular{10^{4}}$)',  linewidth=3, markersize=8., fillstyle='none')
    plt.plot(0,GJQuad/(10**4),'s', clip_on=False, color='green', label='GJquad ($\mathregular{10^{4}}$)', linewidth=3, markersize=8., fillstyle='none')
    plt.title('Structural Properties')
    plt.xlabel('r(m)')
    plt.legend(numpoints=1)

def aero_velocity_and_angles(opt):
    yE = opt.aso.discrete.yE
    vi = opt.aso.aero2.vi
    phi = opt.aso.aero2.phi
    alphaJig = opt.aso.results.alphaJig

    zeros = np.linspace(-1.0, 11.0, 10000)
    plt.plot(zeros, np.zeros(zeros.shape), color='black')
    plt.plot(yE, vi/(10**-1), label='$\mathregular{v_i}$ ($\mathregular{10^{-1}}$ m/s)')
    plt.plot(yE, np.rad2deg(phi), label='$\mathregular{\phi}$ (deg)')
    plt.plot(yE, np.rad2deg(alphaJig), label='$\mathregular{\\alpha_{jig}}$ (deg)')
    plt.title('Aerodynamic Velocities and Angles')
    plt.xlabel('r(m)')
    plt.legend(loc=4)

    plt.axis([0.0, 10.0, -20.0, 20.0])

def aerodynamic_forces_and_power(opt):
    yE = opt.aso.discrete.yE
    drag = opt.aso.aero2.Fblade.Fx
    thrust = opt.aso.aero2.thrust.dT
    torque = opt.aso.aero2.Fblade.Q
    moment = opt.aso.aero2.Fblade.My
    power  = opt.aso.aero2.Fblade.P
    Pi = opt.aso.aero2.Fblade.Pi
    Pp = opt.aso.aero2.Fblade.Pp

    plt.plot(yE, thrust, label='Thrust (N/M)')
    plt.plot(yE, drag/(10**-1), label='Drag ($\mathregular{10^{-1}}$ N/M)')
    plt.plot(yE, torque, label='Torque (Nm/m)')
    plt.plot(yE, moment, color='green', label='Moment (Nm/m)')
    plt.plot(yE, power, color='black', label='Power (W/m)')
    plt.plot(yE, Pi, 'k-.', label='Pi (W/m)')
    plt.plot(yE, Pp, 'k--', label='Pp (W/m)')
    plt.title('Aerodynamic Forces and Power')
    plt.xlabel('r(m)')
    plt.legend(loc=2)

def structural_loads(opt):
    yN = opt.aso.discrete.yN
    Z = opt.aso.struc.Finternal[2,:]
    Z_prime = opt.aso.struc.Finternal[3,:]
    X = opt.aso.struc.Finternal[0,:]
    X_prime = opt.aso.struc.Finternal[5,:]
    Y = opt.aso.struc.Finternal[1,:]
    torque = opt.aso.struc.Finternal[4,:]

    plt.plot(yN, Z, 'k--', color='blue', label='Z (N)')
    plt.plot(yN, Z_prime, color='blue', label="Z' (Nm)")
    plt.plot(yN, X/(10**-1), 'k--', color='red', label='X ($\mathregular{10^{-1}}$ N)')
    plt.plot(yN, -X_prime/(10**-1), color='red', label="X' ($\mathregular{10^{-1}}$ Nm)")
    plt.plot(yN, Y, color='black', label="Y' (Nm)")
    plt.plot(yN, torque, color='green', label='Torque (Nm)')
    plt.title('Structural Loads')
    plt.xlabel('r(m)')
    plt.legend(loc=4)

def structural_deformation(opt):
    yN = opt.aso.discrete.yN

    plt.title('Structural Deformation')
    plt.xlabel('r(m)')
    plt.legend(loc=1)

def out_of_plane_failure(opt):
    yN = opt.aso.discrete.yN
    fail_top_cap_11 = opt.aso.struc.fail.top.cap[0, :]
    fail_top_cap_22 = opt.aso.struc.fail.top.cap[1, :]
    fail_top_cap_12 = opt.aso.struc.fail.top.cap[2, :]

    fail_top_plus_11 = opt.aso.struc.fail.top.plus[0, :]
    fail_top_plus_22 = opt.aso.struc.fail.top.plus[1, :]
    fail_top_plus_12 = opt.aso.struc.fail.top.plus[2, :]

    fail_top_minus_11 = opt.aso.struc.fail.top.minus[0, :]
    fail_top_minus_22 = opt.aso.struc.fail.top.minus[1, :]
    fail_top_minus_12 = opt.aso.struc.fail.top.minus[2, :]

    plt.plot(yN, fail_top_cap_11, 'k-',  label='Top$\mathregular{_c}$ 11')
    plt.plot(yN, fail_top_cap_22, 'k--', label='Top$\mathregular{_c}$ 22')
    plt.plot(yN, fail_top_cap_12, 'k-.', label='Top$\mathregular{_c}$ 12')

    plt.plot(yN, fail_top_plus_11, 'k-',  color='blue', label='Top$\mathregular{_+}$ 11')
    plt.plot(yN, fail_top_plus_22, 'k--', color='blue', label='Top$\mathregular{_+}$ 22')
    plt.plot(yN, fail_top_plus_12, 'k-.', color='blue', label='Top$\mathregular{_+}$ 12')

    plt.plot(yN, fail_top_minus_11, 'k-',  color='red', label='Top$\mathregular{_-}$ 11')
    plt.plot(yN, fail_top_minus_22, 'k--', color='red', label='Top$\mathregular{_-}$ 22')
    plt.plot(yN, fail_top_minus_12, 'k-.', color='red', label='Top$\mathregular{_-}$ 12')

    plt.title('Out-of-Plane Failure')
    plt.xlabel('r(m)')
    plt.legend(loc=1)

def in_plane_failure(opt):
    yN = opt.aso.discrete.yN

    fail_back_plus_11 = opt.aso.struc.fail.back.plus[0, :]
    fail_back_plus_22 = opt.aso.struc.fail.back.plus[1, :]
    fail_back_plus_12 = opt.aso.struc.fail.back.plus[2, :]

    fail_back_minus_11 = opt.aso.struc.fail.back.minus[0, :]
    fail_back_minus_22 = opt.aso.struc.fail.back.minus[1, :]
    fail_back_minus_12 = opt.aso.struc.fail.back.minus[2, :]

    plt.plot(yN, fail_back_plus_11, 'k-',  color='blue', label='back$\mathregular{_+}$ 11')
    plt.plot(yN, fail_back_plus_22, 'k--', color='blue', label='back$\mathregular{_+}$ 22')
    plt.plot(yN, fail_back_plus_12, 'k-.', color='blue', label='back$\mathregular{_+}$ 12')

    plt.plot(yN, fail_back_minus_11, 'k-',  color='red', label='back$\mathregular{_-}$ 11')
    plt.plot(yN, fail_back_minus_22, 'k--', color='red', label='back$\mathregular{_-}$ 22')
    plt.plot(yN, fail_back_minus_12, 'k-.', color='red', label='back$\mathregular{_-}$ 12')

    plt.title('In-Plane Failure')
    plt.xlabel('r(m)')
    plt.legend(loc=1)

def buckling_failure(opt):
    yN = opt.aso.discrete.yN
    yWire = opt.aso.discrete.yWire

    buckle_z = opt.aso.struc.fail.buckling.z
    buckle_x = opt.aso.struc.fail.buckling.x
    buckle_torsion = opt.aso.struc.fail.buckling.torsion
    quad_buckling = opt.aso.struc.fail.quad_buckling
    wire = opt.aso.struc.fail.wire

    plt.plot(yN, buckle_z, color='blue', label='Buckle$\mathregular{_Z}$')
    plt.plot(yN, buckle_z, color='red', label='Buckle$\mathregular{_X}$')
    plt.plot(yN, buckle_torsion, color='green', label='Buckle$\mathregular{_{Tor}}$')
    plt.plot(0, quad_buckling, 's', clip_on=False, color='black', label='Quad$\mathregular{_{Buckle}}$', markersize=9., fillstyle='none', linewidth=3)
    plt.plot(yWire, wire, 'o', color='blue', label='Wire', markersize=9., fillstyle='none', linewidth=3)
    plt.title('Buckling Failure')
    plt.xlabel('r(m)')
    plt.legend(loc=1, numpoints=1)

    plt.axis([0.0, 10.0, 0.0, 1.4])

def plot(opt):
    plt.figure()
    aero_def(opt)

    plt.figure()
    struct_def(opt)

    plt.figure()
    struc_prop(opt)

    plt.figure()
    wake(opt)

    plt.figure()
    aero_velocity_and_angles(opt)

    plt.figure()
    aerodynamic_forces_and_power(opt)

    plt.figure()
    structural_loads(opt)

    #plt.figure()
    #structural_deformation(opt)

    plt.figure()
    out_of_plane_failure(opt)

    plt.figure()
    in_plane_failure(opt)

    plt.figure()
    buckling_failure(opt)

    plt.show()

def plot_single(opt):

    params = {'size':15,
          'font.size' : 10.5,
          'legend.fontsize': 9,
          'legend.linewidth': 1,
          #'text.usetex': True
          }

    plt.rcParams.update(params)
    plt.figure()
    plt.subplot(3,4,1)

    plt.subplot(3,4,2)
    aero_def(opt)

    plt.subplot(3,4,3)
    struct_def(opt)

    plt.subplot(3,4,4)
    struc_prop(opt)

    plt.subplot(3,4,5)
    aero_velocity_and_angles(opt)

    plt.subplot(3,4,6)
    aerodynamic_forces_and_power(opt)

    plt.subplot(3,4,7)
    structural_loads(opt)

    plt.subplot(3,4,8)
    structural_deformation(opt)

    plt.subplot(3,4,9)
    wake(opt)

    plt.subplot(3,4,10)
    out_of_plane_failure(opt)

    plt.subplot(3,4,11)
    in_plane_failure(opt)

    plt.subplot(3,4,12)
    buckling_failure(opt)

    plt.show()


if __name__ == "__main__":
    import pylab
    from heli_opt import HeliOpt
    opt = HeliOpt()
    opt.aso.Omega_opt = 1.0512
    opt.driver.start_iteration()
    opt.driver.run_iteration()

    print 'Parameter:  Omega =', opt.aso.config.Omega

    print 'Constraint: Weight-Lift =', (opt.aso.Mtot*9.8-opt.aso.Ttot)

    print 'Objective:  Ptot =', opt.aso.Ptot

    # enable_trace()


    print 'Parameter:  Omega =', opt.aso.config.Omega

    print 'Constraint: Weight-Lift =', (opt.aso.Mtot*9.8-opt.aso.Ttot)

    print 'Objective:  Ptot =', opt.aso.Ptot

    #opt.run()

    plot_single(opt)




