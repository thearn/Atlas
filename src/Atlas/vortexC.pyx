from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int
import numpy as np
from numpy import pi, cos, sin, mean, linspace, sqrt


# class VortexRingC(Component):
#     """
#     Vortex ring calculations
#     Computes the induced velocity on the rotor blades given the
#     thrust distribution on the rotor
#     """
#     # inputs
#     b        = Int(iotype='in', desc='number of blades')
#     Ns       = Int(iotype='in', desc='number of elements')
#     yN       = Array(iotype='in', desc='node locations')
#     rho      = Float(iotype='in', desc='air density')
#     vc       = Float(iotype='in', desc='vertical velocity')
#     Omega    = Float(iotype='in', desc='rotor angular velocity')
#     h        = Float(iotype='in', desc='height of rotor')
#     dT       = Array(iotype='in', desc='thrust')
#     q        = Array(iotype='in', desc='deformation')
#     anhedral = Float(iotype='in')

#     # outputs
#     vi       = Array(iotype='out', desc='induced velocity')
#     Gamma    = Array(iotype='out', desc='vortex strength')
#     z        = Array(iotype='out', desc='')
#     r        = Array(iotype='out', desc='')
#     vz       = Array(iotype='out', desc='')
#     vr       = Array(iotype='out', desc='')

#     def execute(self):
#         dy = np.zeros((self.Ns, 1))
#         yE = np.zeros((self.Ns, 1))
#         for s in range(self.Ns):
#             dy[s] = self.yN[s+1] - self.yN[s]  # length of each element
#             yE[s] = 0.5 * (self.yN[s] + self.yN[s+1])  # radial location of each element

#         # set fidelity
#         if self.Ns == 15:
#             Nw = 15
#             Ntt = 5
#             Ntheta = 40
#         elif self.Ns == 10:
#             Nw = 8
#             Ntt = 3
#             Ntheta = 20
#         else:
#             Nw = 5
#             Ntt = 1
#             Ntheta = 15

#         # Break out deformations
#         qq = np.zeros((6, self.Ns+1))
#         for s in range(1, self.Ns+1):
#             qq[:, s] = self.q[s*6:s*6+6].T
#         qh = (qq[2, :] - self.yN.T * self.anhedral).flatten()

#         cr = 0.5 * mean(dy)
#         dtheta = pi / Ntheta
#         self.thetaArray = linspace(dtheta/2, pi - dtheta/2, Ntheta)

#         # pre-allocate
#         self.Gamma = np.zeros((Nw+1, self.Ns+1))
#         self.z     = np.zeros((Nw+1, self.Ns+1))
#         self.r     = np.zeros((Nw+1, self.Ns+1))
#         self.vz    = np.zeros((Nw+1, self.Ns+1))
#         self.vr    = np.zeros((Nw+1, self.Ns+1))
#         self.vi    = np.zeros((self.Ns, 1))

#         # create nacent vortex rings
#         GammaBound = self.dT / (self.rho*(self.Omega*yE)*dy)

#         self.Gamma[0, 0] = -GammaBound[0]
#         for s in range(1, self.Ns):
#             self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
#         self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]

#         self.r[0, :] = self.yN.T
#         self.z[0, :] = qh[:]

#         # free-wake time stepping
#         for t in range(Nw+1):
#             # proceed with substeps
#             for tt in range(Ntt):
#                 # Compute induced velocity on all ix(Ns+1) rings from all iix(Ns+1) rings
#                 self.vz = np.zeros((Nw+1, self.Ns+1))
#                 self.vr = np.zeros((Nw+1, self.Ns+1))
#                 for i in range(t):                       # for each disk
#                     for s in range(self.Ns+1):           # and for each ring on each disk
#                         for ii in range(t):              # add the velocity induced from each disk
#                             for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
#                                 zr = self.z[ii, ss]
#                                 r  = self.r[ii, ss]
#                                 zp = self.z[i, s]
#                                 yp = self.r[i, s]

#                                 M  = self.Gamma[ii, ss] * r * dtheta / (2*pi)
#                                 X2 = (-r*sin(self.thetaArray))**2
#                                 Y2 = (yp - r*cos(self.thetaArray))**2
#                                 Z2 = (zp - zr)**2
#                                 Normal = sqrt(X2 + Y2 + Z2)
#                                 for iii in range(max(self.thetaArray.shape)):
#                                     if Normal[iii] < cr:
#                                         Normal[iii] = cr
#                                 Norm3 = Normal**3

#                                 self.vr[i, s] = self.vr[i, s] + np.sum(-cos(self.thetaArray) * (zp - zr) / Norm3) * M
#                                 self.vz[i, s] = self.vz[i, s] + np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M

#                                 # ground effect ring
#                                 zr = -2*self.h - self.z[ii, ss]
#                                 Z2 = (zp - zr)**2
#                                 Normal = sqrt(X2 + Y2 + Z2)
#                                 for iii in range(max(self.thetaArray.shape)):
#                                     if Normal[iii] < cr:
#                                         Normal[iii] = cr
#                                 Norm3 = Normal**3

#                                 self.vr[i, s] = self.vr[i, s] - np.sum(-cos(self.thetaArray) * (zp - zr) / Norm3) * M
#                                 self.vz[i, s] = self.vz[i, s] - np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M

#                 # Compute altitude and time power approximation
#                 # if tt == 0:
#                 #     PiApprox = 8 * np.sum(self.dT.dot(vi))
#                 # realtime = 2*pi / self.Omega / self.b * ((t-1) * Ntt + tt) / Ntt
#                 # altitude = self.vc * realtime

#                 # Convect rings downstream
#                 dt = 2*pi / self.Omega / self.b / Ntt
#                 for i in range(t):              # for each disk
#                     for s in range(self.Ns+1):  # and for each ring on each disk
#                         self.z[i, s] = self.z[i, s] + self.vz[i, s]*dt
#                         self.r[i, s] = self.r[i, s] + self.vr[i, s]*dt
#                         self.Gamma[i, s] = self.Gamma[i, s]

#             # Shift elements in ring array
#             for i in range(t)[::-1]:
#                 self.Gamma[i+1, :] = self.Gamma[i, :]
#                 self.r[i+1, :] = self.r[i, :]
#                 self.z[i+1, :] = self.z[i, :]

#             # Create nacent vortex rings
#             GammaBound = self.dT / (self.rho*(self.Omega*yE)*dy)
#             self.Gamma[0, 0] = -GammaBound[0]
#             for s in range(1, self.Ns):
#                 self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
#             self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]
#             self.r[0, :] = self.yN.T
#             self.z[0, :] = qh[:]

#         # Compute induced velocity on rotor (rp = [0 r(s) 0])
#         for s in range(self.Ns):                # for each element
#             self.vi[s] = 0
#             for ii in range(Nw):                # add the velocity induced from each disk
#                 for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
#                     if (ii == 0):
#                         ringFrac = 0.675
#                     else:
#                         ringFrac = 1

#                     zr = self.z[ii, ss]
#                     r  = self.r[ii, ss]
#                     zp = (qh[s] + qh[s+1]) / 2
#                     yp = yE[s]

#                     M  = self.Gamma[ii, ss] * r * dtheta / (2*pi)
#                     X2 = (-r*sin(self.thetaArray))**2
#                     Y2 = (yp - r*cos(self.thetaArray))**2
#                     Z2 = (zp - zr)**2
#                     Normal = sqrt(X2 + Y2 + Z2)
#                     for iii in range(max(self.thetaArray.shape)):
#                         if Normal[iii] < cr:
#                             Normal[iii] = cr
#                     Norm3 = Normal**3

#                     self.vi[s] = self.vi[s] + ringFrac * np.sum((cos(self.thetaArray)*yp - r) / Norm3) * M

#                     # ground effect ring
#                     zr = -2*self.h - self.z[ii, ss]
#                     Z2 = (zp - zr)**2
#                     Normal = sqrt(X2 + Y2 + Z2)
#                     for iii in range(max(self.thetaArray.shape)):
#                         if Normal[iii] < cr:
#                             Normal[iii] = cr
#                     Norm3 = Normal**3

#                     self.vi[s] = self.vi[s] - ringFrac * np.sum((cos(self.thetaArray)*yp - r) / Norm3) * M

#         # vi is positive downwards
#         self.vi = -self.vi


cimport cython
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

from libc.math cimport sin, cos, sqrt

#from libc.math cimport sqrt # that is for scalars


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def main_loop(
    double h,
    double rho,
    np.ndarray[DTYPE_t, ndim=2] yE,
    np.ndarray[DTYPE_t, ndim=2] dy,
    double[:] qh,
    int Nw,
    int Ntt,
    int Ns,
    np.ndarray[DTYPE_t, ndim=2] z,
    np.ndarray[DTYPE_t, ndim=2] r,
    np.ndarray[DTYPE_t, ndim=2, mode='c']  Gamma,
    #double[:, ::1]  Gamma,
    double Omega,
    np.ndarray[DTYPE_t, ndim=2] dT,
    np.ndarray[DTYPE_t, ndim=2] yN,
    double b,
    double dtheta,
    unsigned int Ntheta,
    double[:] thetaArray,
    double cr):

    '''return vz, vr, z, r, Gamma'''

    cdef unsigned int t, tt, i, s, ii, ss, iii
    cdef double zr,r_scalar, zp, yp, M, Z2
    cdef np.ndarray[DTYPE_t, ndim=2] vz,vr
    cdef np.ndarray[np.double_t, ndim=1] Normal, X2, Y2
    cdef double dt

    cdef double pi_c = pi
    cdef int max_thetaArray_shape = Ntheta

    cdef double one_over_two_pi = 1.0 / 2 * pi

    cdef double acc = 0.0
    cdef double[:] tmp
    cdef int j
    cdef double v, n
    cdef np.ndarray[np.double_t, ndim=1] inv_Norm3

    cdef double[:] sin_thetaArray, cos_thetaArray

    X2 = np.empty(Ntheta)
    Y2 = np.empty(Ntheta)
    Normal = np.empty(Ntheta)
    inv_Norm3 = np.empty(Ntheta)
    sin_thetaArray = np.empty(Ntheta)
    cos_thetaArray = np.empty(Ntheta)

    #sin_thetaArray = sin(thetaArray)
    #cos_thetaArray = cos(thetaArray)
    for j in range(Ntheta):
        sin_thetaArray[j] = sin(thetaArray[j])
        cos_thetaArray[j] = cos(thetaArray[j])

    #sin_thetaArray = sin( thetaArray)
    #cos_thetaArray = cos( thetaArray)

    print "inside main_loop: NW,Ntt,Ns,thetaArray.shape", Nw, Ntt, Ns, Ntheta

     # free-wake time stepping
    for t in range(Nw+1):
        # proceed with substeps
        for tt in range(Ntt):
            # Compute induced velocity on all ix(Ns+1) rings from all iix(Ns+1) rings
            vz = np.zeros((Nw+1, Ns+1))
            vr = np.zeros((Nw+1, Ns+1))
            for i in range(t):                       # for each disk
                for s in range(Ns+1):           # and for each ring on each disk
                    for ii in range(t):              # add the velocity induced from each disk
                        for ss in range(1, Ns+1):  # and each ring on each disk (inner ring cancels itself out)
                            zr = z[ii, ss]
                            r_scalar  = r[ii, ss]
                            zp = z[i, s]
                            yp = r[i, s]
                            M  = Gamma[ii, ss] * r_scalar * dtheta * one_over_two_pi

                            #X2 = (-r_scalar*sin(thetaArray))**2
                            #X2 = (-r_scalar*sin_thetaArray)**2
                            for j in range(Ntheta):
                                v = (-r_scalar*sin_thetaArray[j])
                                X2[j] = v * v

                            #Y2 = (yp - r_scalar*cos(thetaArray))**2
                            #Y2 = (yp - r_scalar*cos_thetaArray)**2
                            for j in range(Ntheta):
                                v = (yp - r_scalar*cos_thetaArray[j])
                                Y2[j] = v * v

                            Z2 = (zp - zr)**2
                            #Normal = sqrt(X2 + Y2 + Z2)
                            for j in range(Ntheta):
                                Normal[j] = sqrt(X2[j] + Y2[j] + Z2)

                            for iii in range(Ntheta):
                                if Normal[iii] < cr:
                                    Normal[iii] = cr
                            #Normal[ Normal < cr ] = cr

                            #Norm3 = Normal**3
                            #Norm3 = Normal*Normal*Normal
                            #nv_Norm3 = 1.0/ (Normal*Normal*Normal)
                            for j in range(Ntheta):
                                n = Normal[j]
                                inv_Norm3[j] = 1.0 / (n*n*n)
                            #vr[i, s] = vr[i, s] + np.sum(-cos(thetaArray) * (zp - zr) / Norm3) * M
                            #tmp = -cos_thetaArray * (zp - zr) * inv_Norm3
                            acc = 0.0
                            for j in range(Ntheta):
                                acc += -cos_thetaArray[j] * (zp - zr) * inv_Norm3[j]
                            vr[i, s] +=  acc * M

                            #vz[i, s] = vz[i, s] + np.sum((cos(thetaArray) * yp - r_scalar) / Norm3) * M
                            #tmp = (cos_thetaArray * yp - r_scalar) * inv_Norm3
                            acc = 0.0
                            for j in range(Ntheta):
                                acc += (cos_thetaArray[j] * yp - r_scalar) * inv_Norm3[j]
                            vz[i, s] += acc * M

                            zr = -2*h - z[ii, ss]
                            Z2 = (zp - zr)**2
                            #Normal = sqrt(X2 + Y2 + Z2)
                            for j in range(Ntheta):
                                Normal[j] = sqrt(X2[j] + Y2[j] + Z2)

                            for iii in range(Ntheta):
                                if Normal[iii] < cr:
                                    Normal[iii] = cr
                            #Normal[ Normal < cr ] = cr

                            #Norm3 = Normal**3
                            #inv_Norm3 = 1.0/ (Normal*Normal*Normal)
                            for j in range(Ntheta):
                                n = Normal[j]
                                inv_Norm3[j] = 1.0 / (n*n*n)

                            #vr[i, s] = vr[i, s] - np.sum(-cos(thetaArray) * (zp - zr) / Norm3) * M
                            #tmp = -cos_thetaArray * (zp - zr) * inv_Norm3
                            acc = 0.0
                            for j in range(Ntheta):
                                acc += -cos_thetaArray[j] * (zp - zr) * inv_Norm3[j]
                            vr[i, s] -= - acc * M

                            #vz[i, s] = vz[i, s] - np.sum((cos(thetaArray) * yp - r_scalar) / Norm3) * M
                            #tmp = (cos_thetaArray * yp - r_scalar) * inv_Norm3
                            acc = 0.0
                            for j in range(Ntheta):
                                acc += (cos_thetaArray[j] * yp - r_scalar) * inv_Norm3[j]
                            vz[i, s] -= - acc * M

            # Compute altitude and time power approximation
            # if tt == 0:
            #     PiApprox = 8 * np.sum(dT.dot(vi))
            # realtime = 2*pi / Omega / self.b * ((t-1) * Ntt + tt) / Ntt
            # altitude = self.vc * realtime

            # Convect rings downstream
            dt = 2*pi_c / Omega / b / Ntt
            for i in range(t):              # for each disk
                for s in range(Ns+1):  # and for each ring on each disk
                    z[i, s] += vz[i, s]*dt
                    r[i, s] += vr[i, s]*dt
                    #Gamma[i, s] = Gamma[i, s]

        # Shift elements in ring array
        for i in range(t)[::-1]:
            Gamma[i+1, :] = Gamma[i, :]
            r[i+1, :] = r[i, :]
            z[i+1, :] = z[i, :]

        # Create nacent vortex rings
        GammaBound = dT / (rho*(Omega*yE)*dy)
        Gamma[0, 0] = -GammaBound[0]
        for s in range(1, Ns):
            Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
        Gamma[0, Ns] = GammaBound[Ns-1]
        r[0, :] = yN.T
        z[0, :] = qh[:]

    return vz, vr, z, r, Gamma


class VortexRingC(Component):
    """
    Vortex ring calculations
    Computes the induced velocity on the rotor blades given the
    thrust distribution on the rotor
    """

    # inputs
    b        = Int(iotype='in', desc='number of blades')
    Ns       = Int(iotype='in', desc='number of elements')
    yN       = Array(iotype='in', desc='node locations')
    rho      = Float(iotype='in', desc='air density')
    vc       = Float(iotype='in', desc='vertical velocity')
    Omega    = Float(iotype='in', desc='rotor angular velocity')
    h        = Float(iotype='in', desc='height of rotor')
    dT       = Array(iotype='in', desc='thrust')
    q        = Array(iotype='in', desc='deformation')
    anhedral = Float(iotype='in')

    # outputs
    vi       = Array(iotype='out', desc='induced velocity')
    Gamma    = Array(iotype='out', desc='vortex strength')
    z        = Array(iotype='out', desc='')
    r        = Array(iotype='out', desc='')
    vz       = Array(iotype='out', desc='')
    vr       = Array(iotype='out', desc='')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def execute(self):
        print "inside cython version"

        cdef unsigned int t, tt, i, s, ii, ss, iii
        cdef unsigned int Nw, Ntt, Ntheta, max_thetaArray_shape

        cdef np.ndarray[DTYPE_t, ndim=2] dy, yE
        cdef np.ndarray[DTYPE_t, ndim=2] qq
        #cdef np.ndarray[DTYPE_t, ndim=1] qh
        cdef double[:] qh  # this work also. Possibly faster too!

        # could try this form
        #cdef np.ndarray[np.double_t, ndim=2] some_array

        cdef double cr, dtheta

        cdef np.ndarray[DTYPE_t, ndim=2] GammaBound

        cdef double zr, r, zp, yp, M, Z2
        cdef np.ndarray[DTYPE_t, ndim=1] X2, Y2, Normal
        #cdef double[:] X2, Y2, Normal, Norm3 # does not seem to work with array operations like '+' and 'sqrt'. Get cython errors

        cdef double dt, ringFrac

        cdef double[:] thetaArray

        # how do we type attributes?????????????? like yN

        #cdef int k
        cdef np.intp_t k
        cdef np.ndarray[DTYPE_t, ndim=2] A

        k = 3
        A = np.zeros((k, 1))

        # Ns = self.Ns
        # dy = np.zeros((Ns, 1))
        dy = np.zeros((self.Ns, 1))
        yE = np.zeros((self.Ns, 1))
        for s in range(self.Ns):
            dy[s] = self.yN[s+1] - self.yN[s]  # length of each element
            yE[s] = 0.5 * (self.yN[s] + self.yN[s+1])  # radial location of each element

        # set fidelity
        if self.Ns == 15:
            Nw = 15
            Ntt = 5
            Ntheta = 40
        elif self.Ns == 10:
            Nw = 8
            Ntt = 3
            Ntheta = 20
        else:
            Nw = 5
            Ntt = 1
            Ntheta = 15

        max_thetaArray_shape = Ntheta

        # Break out deformations
        qq = np.zeros((6, self.Ns+1))
        for s in range(1, self.Ns+1):
            qq[:, s] = self.q[s*6:s*6+6].T
        qh = (qq[2, :] - self.yN.T * self.anhedral).flatten()

        #print qq.shape, qh.shape

        cr = 0.5 * mean(dy)
        dtheta = pi / Ntheta
        thetaArray = linspace(dtheta/2, pi - dtheta/2, Ntheta)

        # import pprint
        # pprint.pprint(thetaArray)

        # pre-allocate
        self.Gamma = np.zeros((Nw+1, self.Ns+1))
        self.z     = np.zeros((Nw+1, self.Ns+1))
        self.r     = np.zeros((Nw+1, self.Ns+1))
        self.vz    = np.zeros((Nw+1, self.Ns+1))
        self.vr    = np.zeros((Nw+1, self.Ns+1))
        self.vi    = np.zeros((self.Ns, 1))

        # create nacent vortex rings
        GammaBound = self.dT / (self.rho*(self.Omega*yE)*dy)

        self.Gamma[0, 0] = -GammaBound[0]
        for s in range(1, self.Ns):
            self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
        self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]

        self.r[0, :] = self.yN.T
        self.z[0, :] = qh[:]

        # t0 = time.time()
        self.vz, self.vr, self.z, self.r, self.Gamma = main_loop(
            self.h,
            self.rho,
            yE,
            dy,
            qh,
            Nw,
            Ntt,
            self.Ns,
            self.z,
            self.r,
            self.Gamma,
            self.Omega,
            self.dT,
            self.yN,
            self.b,
            dtheta,
            Ntheta,
            thetaArray,
            cr)

        # print "main_loop time", time.time() - t0

        # t0 = time.time()

        X2 = np.empty(Ntheta)
        Y2 = np.empty(Ntheta)
        Normal = np.empty(Ntheta)
        inv_Norm3 = np.empty(Ntheta)
        sin_thetaArray = np.empty(Ntheta)
        cos_thetaArray = np.empty(Ntheta)

        # Compute induced velocity on rotor (rp = [0 r(s) 0])
        for s in range(self.Ns):                # for each element
            self.vi[s] = 0
            for ii in range(Nw):                # add the velocity induced from each disk
                for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
                    if (ii == 0):
                        ringFrac = 0.675
                    else:
                        ringFrac = 1

                    zr = self.z[ii, ss]
                    r  = self.r[ii, ss]
                    zp = (qh[s] + qh[s+1]) / 2
                    yp = yE[s]

                    M  = self.Gamma[ii, ss] * r * dtheta / (2*pi)
                    #X2 = (-r*sin(thetaArray))**2
                    for j in range(Ntheta):
                        X2[j] = (-r*sin_thetaArray[j])**2
                    #Y2 = (yp - r*cos(thetaArray))**2
                    for j in range(Ntheta):
                        Y2[j] = (yp - r*cos_thetaArray[j])**2
                    Z2 = (zp - zr)**2
                    #Normal = sqrt(X2 + Y2 + Z2)
                    for j in range(Ntheta):
                        Normal[j] = sqrt(X2[j] + Y2[j] + Z2)
                    for iii in range(Ntheta):
                        if Normal[iii] < cr:
                            Normal[iii] = cr
                    for j in range(Ntheta):
                        n = Normal[j]
                        inv_Norm3[j] = 1.0 / (n*n*n)

                    #self.vi[s] = self.vi[s] + ringFrac * np.sum((cos(thetaArray)*yp - r) / Norm3) * M
                    acc = 0.0
                    for j in range(Ntheta):
                        acc += (cos_thetaArray[j]*yp - r) * inv_Norm3[j]
                    self.vi[s] += ringFrac * M * acc
                    # ground effect ring
                    zr = -2*self.h - self.z[ii, ss]
                    Z2 = (zp - zr)**2
                    #Normal = sqrt(X2 + Y2 + Z2)
                    for j in range(Ntheta):
                        Normal[j] = sqrt(X2[j] + Y2[j] + Z2)
                    for iii in range(Ntheta):
                        if Normal[iii] < cr:
                            Normal[iii] = cr
                    for j in range(Ntheta):
                        n = Normal[j]
                        inv_Norm3[j] = 1.0 / (n*n*n)
                    #Norm3 = Normal**3

                    #self.vi[s] = self.vi[s] - ringFrac * np.sum((cos(thetaArray)*yp - r) / Norm3) * M
                    acc = 0
                    for j in range(Ntheta):
                        acc += (cos_thetaArray[j]*yp - r) * inv_Norm3[j]
                    self.vi[s] -= ringFrac * acc * M

        # print "other loop time", time.time() - t0

        # vi is positive downwards
        self.vi = -self.vi
        print 'vi:\n', self.vi
