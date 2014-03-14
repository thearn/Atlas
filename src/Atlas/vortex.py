import numpy as np
# from openmdao.main.api import Component
# from openmdao.lib.datatypes.api import Float, Array, Int
#

# class vortexWakeCover(Component):
#     """
#     Computes vortexWakeCover
#     """
#     yN = Array(np.zeros(2), iotype='in', desc='description')
#     rho = Float(0, iotype='in', desc='description')
#     dT = Array(np.zeros(2), iotype='in', desc='description')
#     vc = Float(0, iotype='in', desc='description')
#     Omega = Float(0, iotype='in', desc='description')
#     b = Float(0, iotype='in', desc='description')
#     h = Float(0, iotype='in', desc='description')
#     Nw = Int(0, iotype='in', desc='description')
#     Ntt = Int(0, iotype='in', desc='description')
#     Ntheta = Int(0, iotype='in', desc='description')
#     qh = Array(np.zeros(2), iotype='in', desc='description')
#     ycmax = Float(0, iotype='in', desc='description')
#     flag_cover = Int(0, iotype='in', desc='description')
#     flag_plot = Int(0, iotype='in', desc='description')
#     flag_dynamic_climb = Int(0, iotype='in', desc='description')
#     vi = Array(np.zeros(2), iotype='out', desc='description')
#     ring_gamma = Array(np.zeros(2), iotype='out', desc='description')
#     ring_z = Array(np.zeros(2), iotype='out', desc='description')
#     ring_r = Array(np.zeros(2), iotype='out', desc='description')
#     ring_vz = Array(np.zeros(2), iotype='out', desc='description')
#     ring_vr = Array(np.zeros(2), iotype='out', desc='description')
#     def execute(self):
#        pass


yN = np.arange(11)
rho = 1.18
dT = np.array([0.0649, 1.5887, 6.8049, 11.2556, 15.7854, 19.8941, 23.2617,
               25.7216, 27.2537, 18.6446])
vc = 0
Omega = 1.0367
b = 2
h = 1.5
Nw = 8
Ntt = 3
Ntheta = 20
qh = np.array([0, -0.0179, -0.0427, -0.0727, -0.1049, -0.1331, -0.1445,
               -0.1247, -0.0789, -0.0181, 0.0480])
ycmax = 1.4656

"""
-----------------
"""

Ns = len(yN) - 1
dy = np.diff(yN)
yE = np.zeros(Ns)
for i in xrange(Ns):
    yE[i] = 0.5 * (yN[i + 1] + yN[i])

R = yN[-1]
cr = 0.5*np.mean(dy)
dtheta = np.pi/Ntheta
thetaArray = np.linspace(dtheta/2., np.pi - dtheta/2., Ntheta)

ring = {}
ring["Gamma"] = np.zeros((Nw+1, Ns+1))
ring["z"] = np.zeros((Nw+1, Ns+1))
ring["r"] = np.zeros((Nw+1, Ns+1))
ring["vz"] = np.zeros((Nw+1, Ns+1))
ring["vr"] = np.zeros((Nw+1, Ns+1))
vi = np.zeros(Ns)

GammaBound = dT / (rho/(Omega*yE)*dy)
ring["Gamma"][0,0] = -GammaBound[0]
for i in xrange(1, Ns):
    ring["Gamma"][0,i] = GammaBound[i-1] - GammaBound[i]
ring["Gamma"][0, -1] = GammaBound[-1]

ring["r"][0,:] = yN
ring["z"][0,:] = qh

# free wake time step
for t in xrange(Nw):

    # proceed with substeps
    for tt in xrange(Ntt):
        # Compute induced velocity on all ix(Ns+1) rings from all iix(Ns+1) rings
        ring["vz"] = np.zeros((Nw+1, Ns+1))
        ring["vr"] = np.zeros((Nw+1, Ns+1))

        # for each disc
        for i in xrange(t):
            # and each ring
            for s in xrange(Ns+1):

                # add velocity induced from each disc
                for ii in xrange(t):
                    # and each ring on each disc (inner ring cancels itself out)
                    for ss in xrange(1, Ns+1):

                        zr = ring["z"][ii,ss]
                        r = ring["r"][ii,ss]
                        zp = ring["z"][i,s]
                        yp = ring["r"][i,s]

                        M = ring["Gamma"][ii, ss]*r*dtheta/(2*np.pi)
                        X2 = (-r*np.sin(thetaArray))**2
                        Y2 = (yp - r*np.cos(thetaArray))**2
                        Z2 = (zp - zr)**2
                        Normal = np.sqrt(X2 + Y2 + Z2)

                        Normal[Normal < cr] = cr

                        Norm3 = Normal**3

                        ring["vr"][i, s] += sum(-np.cos(thetaArray) * (zp - zr)/Norm3)*M
                        ring["vz"][i, s] += sum((np.cos(thetaArray)*yp-r)/Norm3)*M

                        # ground effect ring
                        zr = -2*h - ring["z"][ii, ss]
                        Z2 = (zp - zr)**2
                        Normal = np.sqrt(X2 + Y2 + Z2)
                        Normal[Normal < cr] = cr
                        Norm3 = Normal**3

                        ring["vr"][i, s] += -sum(-np.cos(thetaArray) * (zp - zr)/Norm3)*M
                        ring["vz"][i, s] += -sum((np.cos(thetaArray)*yp-r)/Norm3)*M

        # compute altitude and time power approximation
        if tt == 0:
            PiApprox = 8*sum(dT*vi)
        realtime = 2*np.pi/Omega/b*((t - 1)*Ntt + tt)/Ntt
        altitude = vc*realtime

        # convect rings downstream
        dt = 2*np.pi/Omega/b/Ntt























