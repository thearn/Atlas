import math

def arctan2(x,y):
    a=x.real
    b=x.imag
    c=y.real
    d=y.imag
    return complex(math.atan2(a,c),(c*b-a*d)/(a**2+c**2))