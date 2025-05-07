import numpy as np
from numpy import rad2deg as r2d
from numpy import deg2rad as d2r
from numpy import tan, sqrt
from colorama import Fore


# BLOCH DEFAULT
ALPHA1 = 1.300943487
ALPHA2 = 7.428215271
ALPHA3 = 17.97942756
B = 21

def parameters(alpha1 = None, alpha2 = None, alpha3 = None, b = None, shup = False):
    """
    Arguments:
        alpha1, alpha2, alpha3 (deg), and b (mm)
    Returns:
        R, X0-L, H
    """
    global ALPHA1, ALPHA2, ALPHA3, B
    dargs = []
    if type(alpha1) is type(None):
        alpha1 = ALPHA1; dargs.append(["\u03B11",alpha1,'°'])
    if type(alpha2) is type(None):
        alpha2 = ALPHA2; dargs.append(["\u03B12",alpha2,'°'])
    if type(alpha3) is type(None):
        alpha3 = ALPHA3; dargs.append(["\u03B13",alpha3,'°'])
    if type(b) is type(None):
        b = B; dargs.append(["2b",2*b,' mm'])
    if len(dargs) > 0:
        if not shup: print(Fore.BLUE + "parameters():\n  Missing args set to default values: " + ", ".join([f"{a[0]} = {a[1]}{a[2]}" for a in dargs]) + Fore.BLACK)
    #
    D,R,X0_L,H = 1,1,1,1
    a1, a2, a3 = tan(d2r(alpha1)), tan(d2r(alpha2)), tan(d2r(alpha3))
    D = (a2-a1)*sqrt(1+a3**2)+(a3-a2)*sqrt(1+a1**2)-(a3-a1)*sqrt(1+a2**2)   # eq. 11
    R = (b/D)*(a3-a2)*(a3-a1)*(a2-a1)
    X0_L = (b/D)*( (a3**2-a1**2)*sqrt(1+a2**2) - (a3**2-a2**2)*sqrt(1+a1**2) - (a2**2-a1**2)*sqrt(1+a3**2))
    H = (b/D)*( (a3-a2)*(1+a2*a3)*sqrt(1+a1**2) - (a3-a1)*(1+a1*a3)*sqrt(1+a2**2) + (a2-a1)*(1+a1*a2)*sqrt(1+a3**2) )
    if not shup:
        print(Fore.BLUE + f"parameters():\n  R = {R}\n  X0-L = {X0_L}\n  H = {H}" + Fore.BLACK)
    return R, X0_L, H


def deviation(alpha = 0, R = None, X0_L = None, H = None, b = None, shup = False):
    """
    Arguments:
        alpha       An angle in degrees (number) or a list defining a range of angles [start, stop, num. points].
        R, X0_L, H  The mono parameters, calculated from alpha1, alpha2, and alpha3 (and b).
    """
    if not (type(alpha) is list or type(alpha) is np.ndarray):
        try: alpha = float(alpha)
        except:
            print(Fore.RED + "deviation(): Argument alpha must be a number or an array (len 3)." + Fore.BLACK)
            return None
    else:
        if not len(alpha) == 3:
            print(Fore.RED + "deviation(): If argument alpha is passed as an array it must have the form [start,stop,num. of steps]." + Fore.BLACK)
            return None
        alpha = np.linspace(alpha[0], alpha[1], alpha[2])
    #
    try:
        R, X0_L, H, b  = float(R), float(X0_L), float(H), float(b)
    except:
        r, x0_l, h = parameters(shup = True)
        if type(R) is type(None): R = r
        if type(X0_L) is type(None): X0_L = x0_l
        if type(H) is type(None): H = h
        if type(b) is type(None):
            global B
            b = B
        if not shup: print(Fore.BLUE + "deviation():\n  Some missing or invalid arguments are set to default values." + Fore.BLACK)
    #
    def devi(R, X0_L, H, b, alpha):
        """alpha in deg"""
        a = tan(d2r(alpha))
        return (b*(1-a**2)+R*sqrt(1+a**2)-H)/a - X0_L
    #
    if type(alpha) is float:
        return devi(R, X0_L, H, b, d2r(alpha))
    else:
        dev = np.zeros([len(alpha)])*np.NaN
        for i, alp in enumerate(alpha):
            dev[i] = devi(R, X0_L, H, b, alp)
        return alpha, dev

        

    

        





    
    

