import pylab as plt
import numpy as np

pi = np.pi


def wake(case):
    plt.title("Free Wake")
    r = np.array(case['aso.aero2.induced.r'])
    z = np.array(case['aso.aero2.induced.z'])
    plt.plot(r, z, "b-")
    plt.plot(r, z, "bo")
    plt.plot(-r, z, "b-")
    plt.plot(-r, z, "bo")

    plt.xlabel("r(m)")


def aero_def(case):
    yE = np.array(case['aso.discrete.yE'])
    cE = np.array(case['aso.discrete.cE'])
    Cl = np.array(case['aso.discrete.Cl'])
    cd = np.array(case['aso.aero2.Cd'])
    re = np.array(case['aso.aero2.Re'])
    R  = np.array(case['aso.config.R'])
    y100 = np.linspace(0, R, 100)
    c100 = np.array(case['aso.discrete.c100'])

    plt.plot(yE, cE, label='Chord (m)')
    plt.plot(yE, Cl, label='C_l')
    plt.plot(yE, 100*cd, label='C_d ($\mathregular{10^{-2}}$)')
    plt.plot(yE, Cl/cd/100., label='L/D ($\mathregular{10^{2}}$)')
    plt.plot(yE, re/1000000., label='Re ($\mathregular{10^{6}}$)')
    plt.plot(y100, c100)
    plt.xlabel("r(m)")
    plt.title('Aerodynamic Definition')
    plt.legend()


def struct_def(case):
    yE           = np.array(case['aso.discrete.yE'])
    d            = np.array(case['aso.discrete.d'])
    theta        = np.array(case['aso.discrete.theta'])
    nTube        = np.array(case['aso.discrete.nTube'])
    nCap         = np.array(case['aso.discrete.nCap'])
    lBiscuit     = np.array(case['aso.discrete.lBiscuit'])
    yWire        = np.array(case['aso.config.yWire'])
    zWire        = np.array(case['aso.config.zWire'])
    dQuad        = np.array(case['aso.config.dQuad'])
    thetaQuad    = np.array(case['aso.config.thetaQuad'])
    nTubeQuad    = np.array(case['aso.config.nTubeQuad'])
    lBiscuitQuad = np.array(case['aso.config.lBiscuitQuad'])
    hQuad        = np.array(case['aso.config.hQuad'])

    plt.title('Structural Definition')
    plt.plot(yE, d*100/2.54, label='Spar diameter (in)')
    plt.plot(yE, theta*180/pi/10, label='Wrap angle ($\mathregular{10^{1}}$ deg)')
    plt.plot(yE, nTube, label='Tube layers')
    plt.plot(yE, nCap, label='Cap layers')
    plt.plot(yE, lBiscuit*100/2.54/12, label='Biscuit spacing (ft)')
    plt.plot(yWire, zWire, "o", color='red', label='Wire location (m)', linewidth=3, markersize=9.)
    plt.plot(0, dQuad*100/2.54, "s", color='black', clip_on=False, label='Quad definition', linewidth=3, markersize=9.)

    plt.xlabel("r(m)")
    plt.legend(numpoints=1)


def struc_prop(case):
    yE     = np.array(case['aso.discrete.yE'])
    EA     = np.array(case['aso.struc.spar.EA'])
    EIz    = np.array(case['aso.struc.spar.EIx'])
    EIx    = np.array(case['aso.struc.spar.EIx'])
    GJ     = np.array(case['aso.struc.spar.GJ'])
    mSpar  = np.array(case['aso.struc.spar.mSpar'])
    mChord = np.array(case['aso.struc.chord.mChord'])
    EIQuad = np.array(case['aso.struc.quad.EIx'])
    GJQuad = np.array(case['aso.struc.quad.GJ'])
    dy     = np.array(case['aso.struc.spar.dy'])

    plt.plot(yE, EA/(10**7),  label='EA ($\mathregular{10^{7}}$)')
    plt.plot(yE, EIz/(10**4), label='EIz ($\mathregular{10^{4}}$)')
    plt.plot(yE, EIx/(10**4), label='EIx ($\mathregular{10^{4}}$)')
    plt.plot(yE, GJ/(10**4),  label='GJ ($\mathregular{10^{4}}$)')
    plt.plot(yE, mSpar/dy,    label='mSpar (kg/m)')
    plt.plot(yE, mChord/dy, 'k--', label='mChord (kg/m)')
    plt.plot(0, EIQuad/(10**4), 's', clip_on=False, color='blue', label='EIquad ($\mathregular{10^{4}}$)',  linewidth=3, markersize=8.)
    plt.plot(0, GJQuad/(10**4), 's', clip_on=False, color='green', label='GJquad ($\mathregular{10^{4}}$)', linewidth=3, markersize=8.)
    plt.title('Structural Properties')
    plt.xlabel('r(m)')
    plt.legend(numpoints=1)


def aero_velocity_and_angles(case):
    yE       = np.array(case['aso.discrete.yE'])
    vi       = np.array(case['aso.aero2.vi'])
    phi      = np.array(case['aso.aero2.phi'])
    alphaJig = np.array(case['aso.results.alphaJig'])

    zeros = np.linspace(-1.0, 11.0, 10000)
    plt.plot(zeros, np.zeros(zeros.shape), color='black')
    plt.plot(yE, vi/(10**-1), label='$\mathregular{v_i}$ ($\mathregular{10^{-1}}$ m/s)')
    plt.plot(yE, np.rad2deg(phi), label='$\mathregular{\phi}$ (deg)')
    plt.plot(yE, np.rad2deg(alphaJig), label='$\mathregular{\\alpha_{jig}}$ (deg)')
    plt.title('Aerodynamic Velocities and Angles')
    plt.xlabel('r(m)')
    plt.legend(loc=4)

    plt.axis([0.0, 10.0, -20.0, 20.0])


def aerodynamic_forces_and_power(case):
    yE     = np.array(case['aso.discrete.yE'])
    thrust = np.array(case['aso.aero2.thrust.dT'])
    drag   = np.array(case['aso.aero2.Fblade']['Fx'])
    torque = np.array(case['aso.aero2.Fblade']['Q'])
    moment = np.array(case['aso.aero2.Fblade']['My'])
    power  = np.array(case['aso.aero2.Fblade']['P'])
    Pi     = np.array(case['aso.aero2.Fblade']['Pi'])
    Pp     = np.array(case['aso.aero2.Fblade']['Pp'])

    plt.plot(yE, thrust, label='Thrust (N/M)')
    plt.plot(yE, drag/(10**-1), label='Drag ($\mathregular{10^{-1}}$ N/M)')
    plt.plot(yE, torque, label='Torque (Nm/m)')
    plt.plot(yE, moment, color='green', label='Moment (Nm/m)')
    plt.plot(yE, power,  color='black', label='Power (W/m)')
    plt.plot(yE, Pi, 'k-.', label='Pi (W/m)')
    plt.plot(yE, Pp, 'k--', label='Pp (W/m)')
    plt.title('Aerodynamic Forces and Power')
    plt.xlabel('r(m)')
    plt.legend(loc=2)


def structural_loads(case):
    yN      = np.array(case['aso.discrete.yN'])
    Z       = np.array(case['aso.struc.Finternal'])[2, :]
    Z_prime = np.array(case['aso.struc.Finternal'])[3, :]
    X       = np.array(case['aso.struc.Finternal'])[0, :]
    X_prime = np.array(case['aso.struc.Finternal'])[5, :]
    Y       = np.array(case['aso.struc.Finternal'])[1, :]
    torque  = np.array(case['aso.struc.Finternal'])[4, :]

    plt.plot(yN, Z, 'k--', color='blue', label='Z (N)')
    plt.plot(yN, Z_prime, color='blue', label="Z' (Nm)")
    plt.plot(yN, X/(10**-1), 'k--', color='red', label='X ($\mathregular{10^{-1}}$ N)')
    plt.plot(yN, -X_prime/(10**-1), color='red', label="X' ($\mathregular{10^{-1}}$ Nm)")
    plt.plot(yN, Y, color='black', label="Y' (Nm)")
    plt.plot(yN, torque, color='green', label='Torque (Nm)')
    plt.title('Structural Loads')
    plt.xlabel('r(m)')
    plt.legend(loc=4)


def structural_deformation(case):
    yN = np.array(case['aso.discrete.yN'])
    qq = np.array(case["aso.struc.fem.q"]).reshape((6,11), order='F')
    a = np.array(case["aso.config.anhedral"])
    # qh = qq(3,:)-yN'*anhedral;
    qh = qq[2] - yN * a
    #plot(yN,qh,yN,qq(1,:),yN,qq(5,:)*180/pi);
    plt.plot(yN, qh, label="Z deflection")
    plt.plot(yN, qq[0], label="X deflection")
    plt.plot(yN, qq[-2]*180./np.pi, label="Twist")

    plt.title('Structural Deformation')
    plt.xlabel('r(m)')
    plt.legend(loc=1)


def out_of_plane_failure(case):
    yN = np.array(case['aso.discrete.yN'])
    fail_top_cap_11   = np.array(case['aso.struc.fail']['top']['cap'])[0, :]
    fail_top_cap_22   = np.array(case['aso.struc.fail']['top']['cap'])[1, :]
    fail_top_cap_12   = np.array(case['aso.struc.fail']['top']['cap'])[2, :]

    fail_top_plus_11  = np.array(case['aso.struc.fail']['top']['plus'])[0, :]
    fail_top_plus_22  = np.array(case['aso.struc.fail']['top']['plus'])[1, :]
    fail_top_plus_12  = np.array(case['aso.struc.fail']['top']['plus'])[2, :]

    fail_top_minus_11 = np.array(case['aso.struc.fail']['top']['minus'])[0, :]
    fail_top_minus_22 = np.array(case['aso.struc.fail']['top']['minus'])[1, :]
    fail_top_minus_12 = np.array(case['aso.struc.fail']['top']['minus'])[2, :]

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


def in_plane_failure(case):
    yN = np.array(case['aso.discrete.yN'])

    fail_back_plus_11  = np.array(case['aso.struc.fail']['back']['plus'])[0, :]
    fail_back_plus_22  = np.array(case['aso.struc.fail']['back']['plus'])[1, :]
    fail_back_plus_12  = np.array(case['aso.struc.fail']['back']['plus'])[2, :]

    fail_back_minus_11 = np.array(case['aso.struc.fail']['back']['minus'])[0, :]
    fail_back_minus_22 = np.array(case['aso.struc.fail']['back']['minus'])[1, :]
    fail_back_minus_12 = np.array(case['aso.struc.fail']['back']['minus'])[2, :]

    plt.plot(yN, fail_back_plus_11, 'k-',  color='blue', label='back$\mathregular{_+}$ 11')
    plt.plot(yN, fail_back_plus_22, 'k--', color='blue', label='back$\mathregular{_+}$ 22')
    plt.plot(yN, fail_back_plus_12, 'k-.', color='blue', label='back$\mathregular{_+}$ 12')

    plt.plot(yN, fail_back_minus_11, 'k-',  color='red', label='back$\mathregular{_-}$ 11')
    plt.plot(yN, fail_back_minus_22, 'k--', color='red', label='back$\mathregular{_-}$ 22')
    plt.plot(yN, fail_back_minus_12, 'k-.', color='red', label='back$\mathregular{_-}$ 12')

    plt.title('In-Plane Failure')
    plt.xlabel('r(m)')
    plt.legend(loc=1)


def buckling_failure(case):
    yN    = np.array(case['aso.discrete.yN'])
    yWire = np.array(case['aso.config.yWire'])

    buckle_z       = np.array(case['aso.struc.fail']['buckling']['z'])
    buckle_x       = np.array(case['aso.struc.fail']['buckling']['x'])
    buckle_torsion = np.array(case['aso.struc.fail']['buckling']['torsion'])
    quad_buckling  = np.array(case['aso.struc.fail']['quad_buckling'])
    wire           = np.array(case['aso.struc.fail']['wire'])

    plt.plot(yN, buckle_z, color='blue', label='Buckle$\mathregular{_Z}$')
    plt.plot(yN, buckle_z, color='red', label='Buckle$\mathregular{_X}$')
    plt.plot(yN, buckle_torsion, color='green', label='Buckle$\mathregular{_{Tor}}$')
    plt.plot(0, quad_buckling, 's', clip_on=False, color='black', label='Quad$\mathregular{_{Buckle}}$', markersize=9., linewidth=3)
    plt.plot(yWire, wire, 'o', color='blue', label='Wire', markersize=9., linewidth=3)
    plt.title('Buckling Failure')
    plt.xlabel('r(m)')
    plt.legend(loc=1, numpoints=1)

    plt.axis([0.0, 10.0, 0.0, 1.4])


def plot(case):
    plt.figure()
    aero_def(case)

    plt.figure()
    struct_def(case)

    plt.figure()
    struc_prop(case)

    plt.figure()
    wake(case)

    plt.figure()
    aero_velocity_and_angles(case)

    plt.figure()
    aerodynamic_forces_and_power(case)

    plt.figure()
    structural_loads(case)

    plt.figure()
    structural_deformation(case)

    plt.figure()
    out_of_plane_failure(case)

    plt.figure()
    in_plane_failure(case)

    plt.figure()
    buckling_failure(case)

    plt.show()


def plot_single(case):

    params = {
        'size': 15,
        'font.size': 10.5,
        'legend.fontsize': 9,
        'legend.linewidth': 1,
        #'text.usetex': True
    }

    plt.rcParams.update(params)
    plt.figure()
    plt.subplot(3, 4, 1)

    plt.subplot(3, 4, 2)
    aero_def(case)

    plt.subplot(3, 4, 3)
    struct_def(case)

    plt.subplot(3, 4, 4)
    struc_prop(case)

    plt.subplot(3, 4, 5)
    aero_velocity_and_angles(case)

    plt.subplot(3, 4, 6)
    aerodynamic_forces_and_power(case)

    plt.subplot(3, 4, 7)
    structural_loads(case)

    plt.subplot(3, 4, 8)
    structural_deformation(case)

    plt.subplot(3, 4, 9)
    wake(case)

    plt.subplot(3, 4, 10)
    out_of_plane_failure(case)

    plt.subplot(3, 4, 11)
    in_plane_failure(case)

    plt.subplot(3, 4, 12)
    buckling_failure(case)

    plt.show()


if __name__ == "__main__":
    from openmdao.lib.casehandlers.api import CaseDataset
    dataset = CaseDataset('heli_opt.json', 'json')
    data = dataset.data.by_case().fetch()
    case = data[-1]
    plot_single(case)
