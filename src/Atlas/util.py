import math

def abs_complex(c):
    x, y = c.real, c.imag
    if x < 0:
        return -x - y*1j
    else:
        return x + y*1j

def arctan2(x,y):
    a=x.real
    b=x.imag
    c=y.real
    d=y.imag
    return complex(math.atan2(a,c),(c*b-a*d)/(a**2+c**2))