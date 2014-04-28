import timeit

setup = '''
import os
from scipy.io import loadmat

import Atlas
from Atlas import %s as VortexRing

# populate inputs
path = os.path.join(os.path.dirname(Atlas.__file__), 'test/vortex.mat')
data = loadmat(path, struct_as_record=True, mat_dtype=True)

comp = VortexRing()

comp.b        = int(data['b'][0][0])
comp.yN       = data['yN'].flatten()
comp.Ns       = max(comp.yN.shape) - 1

comp.rho      = data['rho'][0][0]
comp.vc       = data['vc'][0][0]
comp.Omega    = data['Omega'][0][0]
comp.h        = data['h'][0][0]
comp.dT       = data['dT']
comp.q        = data['q']
comp.anhedral = data['anhedral'][0][0]
'''

t = timeit.Timer("comp.run()", setup % 'VortexRing')
print "Python function", t.timeit(1), "sec"

t = timeit.Timer("comp.run()", setup % 'VortexRingC')
print "Cython function", t.timeit(1), "sec"
