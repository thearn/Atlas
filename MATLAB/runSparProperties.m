% specifications

yN = [ 0, 14.5057 ]
d = 0.1016
theta = 0.6109
nTube = 4
nCap = [ 0, 0 ]
lBiscuit = 0.3048
flags.CFRPType = 1

[EIx, EIz, EA, GJ, mSpar] = SparProperties(yN, d, theta, nTube, nCap, lBiscuit, flags)
