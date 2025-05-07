"""
Methods for manipulating Stokes vectors and Mueller matrices.

This is a new, slimmed down version. Double and tripple checked for bugs and errors. Created 23.06.22.

Versions:

25.05.07    Updated stokes2polarization() / polarizationFromStokes() to also return eccentricity and ellipticity.
            Removed method stokesInfo().
25.05.06    Bugfix in inclinedStokes (from when introducing 'inclination' instead of 'angle' as attribute).
            Added some redundant method calls to have a more unified naming convention.
25.03.05    Had to specify argument 'angle' in Ellipse in method plotEllipse() due to up or downgrade of matplotlib.
24.04.12    Added stokesInfo(), which accepts a Stokes vector and returns a dict with polarizations, inclinations, E-fields, etc.
            Updated inclinedStokes() to be much smarter than the old, stupid version.
24.03.10    Minor bugfixes, added stokes2fields() and fields2stokes() (note difficulty to distinguish signs for S2 and S3 if both non-zero)
23.07.29    As a part of a scripts-wide thing: added gloal variable default_shup.
23.07.08    Minor bug fixes
23.07.05    Added custom labels (for the legend) in plotEllipse() (as argument ownlabel = ['label1', 'label2',...])
23.06.25    Cleaning up the code
23.06.22    Created.
"""


__version__ = "25.07.07"
__author__  = "Mats Leandersson"



import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from matplotlib import cm 
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
from colorama import Fore
from random import uniform


default_shup = False



# ====================================================================================
# ==================================================================================== S T O K E S
# ====================================================================================


def normalizeStokes(stokes = [1,0,0,1], norm = True, shup = None):
    """
    Returns a normalized Stokes vector.
    Normalizes s1, s2, and s3 to 1 if norm == True, else it normalizes them to s0.
    The argument stokes must be a list or array of size 4.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.normalizeStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return list(np.ones([4]) * np.NaN)
    if type(stokes) is list: stokes = np.array(stokes)
    if not len(stokes) == 4:
        print(Fore.RED + "myStokes.normalizeStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return list(np.ones([4]) * np.NaN)
    #
    if norm:
        s0 = np.sqrt(stokes[1]**2 + stokes[2]**2 + stokes[3]**2)
    else:
        s0 = np.sqrt(stokes[0])
    s1 = np.round(stokes[1]/s0,9)
    s2 = np.round(stokes[2]/s0,9)
    s3 = np.round(stokes[3]/s0,9)
    S = np.array([1., s1, s2, s3])
    if not shup:
        print("Stokes vector: {0}".format(stokes))
        print("Normalized:    {0}".format(S))
    return S




def getS0(stokes = [], shup = None):
    """
    Returns s0 calculated from s1, s2, and s3.
    The argument stokes must be a list or array of size 4.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.getS0(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.NaN
    if type(stokes) is list: stokes = np.array(stokes)
    if not len(stokes) == 4:
        print(Fore.RED + "myStokes.getS0(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.NaN
    #
    return np.sqrt(stokes[1]**2 + stokes[2]**2 + stokes[3]**2) 


def inclinedStokes(inclination = None, angle = None, shup = None):
    """
    Returns a Stokes vector for linear inclined polarization (deg.).

    Note: 
    Keeping argument angle for legacy reasons. Pass either inclination or angle.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    try: inclination = float(inclination)
    except:
        try: inclination = float(angle)
        except:
            print(Fore.RED + "myStokes.inclinedStokes(): Argument 'angle' must be a number in degrees." + Fore.RESET)
            return list(np.ones([4]) * np.NaN)
    #
    stokes = np.array([1,1,0,0])
    return rotationMatrix(angle = -inclination, shup = True).dot(stokes)


# Old, stupid version of inclinedStokes()
#def inclinedStokes(angle = 0., shup = None):
#    """
#    Returns a Stokes vector for linear inclined polarization at an angle given by angle (deg.).
#    """
#    if type(shup) is type(None):
#        global default_shup
#        shup = default_shup
#    #
#    try:
#        a = np.mod(angle,180.)
#    except:
#        print(Fore.RED + "myStokes.inclinedStokes(): Argument 'angle' must be a number in degrees." + Fore.RESET)
#        return list(np.ones([4]) * np.NaN)
#    a = np.deg2rad(a)
#    z1 = 1. / np.sqrt( (np.tan(2.*a))**2. + 1. );   S1 = z1
#    z2 = np.sqrt(1. - z1**2.);                      S2 = z2
#    if a >= 0. and a <= np.pi/4.:           S1 =  z1;   S2 =   z2
#    if a > np.pi/4.1 and a <= np.pi/2.:     S1 = -z1;   S2 =   z2
#    if a > np.pi/2. and a <= np.pi * 3/4:   S1 = -z1;   S2 =  -z2
#    if a >  np.pi * 3/4 and a <= np.pi:     S1 =  z1;   S2 =  -z2
#    if abs(S1)<1e-9: S1 = 0.0
#    if abs(S2)<1e-9: S2 = 0.0
#    S = np.array([1,S1,S2,0])
#    if not shup:
#        print("Inclination: {0}°".format(angle))
#        print("Stokes = {0}".format(S))
#    return S




def stokes2polarization(stokes = [1, 0, 0, 1], norm = True, shup = None):
    """Same as polarizationFromStokes()"""
    return polarizationFromStokes(stokes = stokes, norm = norm, shup = shup)

def polarizationFromStokes(stokes = [1, 0, 0, 1], norm = True, shup = None):
    """
    Returns degree of linear and circular polarization, inclination, and eccentricity and ellipticity.
    The argument stokes must be a list or array of size 4.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.polarizationFromStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.NaN, np.NaN, np.NaN
    if type(stokes) is list: stokes = np.array(stokes)
    if not len(stokes) == 4:
        print(Fore.RED + "myStokes.polarizationFromStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.NaN, np.NaN, np.NaN
    #
    if norm: Stokes = normalizeStokes(stokes, shup = True)
    else: Stokes = stokes
    #
    PL = round( np.sqrt(Stokes[1]**2 + Stokes[2]**2)/Stokes[0], 4)
    PC = abs(round( Stokes[3]/Stokes[0], 4))
    Psi = round( 1./2. * np.arctan2(Stokes[2], Stokes[1])*np.rad2deg(1), 1)
    #
    e = np.sqrt( 2 * np.sqrt(Stokes[1]**2 + Stokes[2]**2) / (1 + np.sqrt(Stokes[1]**2 + Stokes[2]**2)) )
    f = 1 - np.sqrt(1 - e**2)
    #
    if not shup:
        print(f"Stokes vector: {stokes}")
        print(f"PC  = {PC:.4f}")
        print(f"PL  = {PL:.4f}")
        print(f"psi = {Psi:.2f}° (inclination)")
        print(f"e   = {e:.4f} (eccentricity)")
        print(f"f   = {f:.4f} (ellipticity or 'flattening')")
    #
    return PL, PC, Psi, e, f



def stokes2optical(stokes = [1,0,0,1], k2 = 1., shup = None):
    """Same as opticalFromStokes()"""
    return opticalFromStokes(stokes = stokes, k2 = k2, shup = shup)

def opticalFromStokes(stokes = [1,0,0,1], k2 = 1., shup = None):
    """
    Returns the optical parameters rs, rp, and delta (deg.) for a mirror which reflected
    circular(+) radiation that resulted in a Stokes vector given by 'stokes'.
    The argument stokes must be a list or array of size 4.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.opticalFromStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.NaN, np.NaN, np.NaN
    if type(stokes) is list: stokes = np.array(stokes)
    if not len(stokes) == 4:
        print(Fore.RED + "myStokes.opticalFromStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.NaN, np.NaN, np.NaN
    #
    rs = np.sqrt(k2) * np.sqrt(stokes[0] + stokes[1])
    rp = np.sqrt(k2) * np.sqrt(stokes[0] - stokes[1])
    #delta =  np.rad2deg( np.arcsin(stokes[2] / (rs * rp)) )
    delta = np.rad2deg( np.arccos(stokes[3] / (rs * rp / k2)) )
    #
    if not shup:
        print("Stokes: {0}".format(stokes))
        print("rs = {0:.3f}, rp = {1:.3f}, delta = {2:.2f}°".format(rs, rp, delta))
    return rs, rp, delta





def randomStokes(shup = None):
    """
    Returns a random but normalized Stokes vector.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    S3 = uniform(-1, 1)
    PL = np.sqrt(1 - S3**2)
    S1 = np.sin(uniform(0, 2*np.pi)) * PL
    S2 = np.sqrt(PL**2 - S1**2)
    stokes = [1, S1, S2, S3]
    if not shup:
        print("Random Stokes vector: {0}".format(stokes))
        PL, PC, Psi, _, _ = polarizationFromStokes(stokes, shup = True)
        print("PL = {0:.3f}, PC = {1:.3f}, Psi = {2:.2f}°".format(PL, PC, Psi))
    return stokes



def stokes2string(stokes = [], decimals = 3, includes0 = True, shup = None):
    """
    The argument stokes must be a list or array of size 4.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.stokes2string(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return "(nan, nan, nan, nan)"
    if type(stokes) is list: stokes = np.array(stokes)
    if not len(stokes) == 4:
        print(Fore.RED + "myStokes.stokes2string(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return "(nan, nan, nan, nan)"
    #
    if not type(decimals) is int: decimals = 3
    #
    if includes0:
        retstr =  "({1:.{0}f}, {2:.{0}f}, {3:.{0}f}, {4:.{0}f})".format(decimals, stokes[0], stokes[1], stokes[2], stokes[3])
    else:
        retstr =  "({1:.{0}f}, {2:.{0}f}, {3:.{0}f})".format(decimals, stokes[1], stokes[2], stokes[3])
    return retstr




    



# ====================================================================================
# ==================================================================================== M U E L L E R
# ====================================================================================


def optical2matrix(rs = np.NaN, rp = np.NaN, delta = np.NaN, shup = None):
    """Same as matrixFromOptical()"""
    return matrixFromOptical(rs = rs, rp = rp, delta = delta, shup = shup)

def matrixFromOptical(rs = np.NaN, rp = np.NaN, delta = np.NaN, shup = None):
    """
    Returns a Mueller matrix made from the optical parameters
    rs, rp, and delta (deg.).
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    try:
        rs = float(rs); rp = float(rp); d = float(delta)
    except:
        print(Fore.RED + "myStokes.matrixFromOptical(): Arguments 'rs', 'rp', and 'delta' must all be numbers." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    d = np.deg2rad(delta)
    m11 = (rs**2 + rp**2)
    m12 = (rs**2 - rp**2)
    m33 = 2*rs*rp*np.cos(d)
    m34 = 2*rs*rp*np.sin(d)
    M = np.array([ [  m11,  m12,    0,    0 ],
                   [  m12,  m11,    0,    0 ],
                   [    0,    0,  m33,  m34 ],
                   [    0,    0, -m34,  m33 ] ]) * 0.5
    return M



def stokes2matrix(stokes = [1,0,0,1], shup = None):
    """Same as matrixFromStokes()"""
    return matrixFromStokes(stokes = stokes, shup = shup)

def matrixFromStokes(stokes = [1,0,0,1], shup = None):
    """
    Returns a mirror matrix. The argument 'stokes' is a Stokes vector (type list or array)
    of size 4, and it is the vector that is the result from circular polarized light (+)
    being reflected by the mirror.
    The argument stokes must be a list or array of size 4.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.matrixFromStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    if type(stokes) is list: stokes = np.array(stokes)
    if not len(stokes) == 4:
        print(Fore.RED + "myStokes.matrixFromStokes(): Argument 'stokes' must be a list or array of length 4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    rs, rp, delta = opticalFromStokes(stokes = stokes, shup = True)
    return np.around(matrixFromOptical(rs = rs, rp = rp, delta = delta), 12)




def matrix2optical(matrix = [], shup = None):
    """Same as opticalFromMatrix()"""
    return opticalFromMatrix(matrix = matrix, shup = shup)

def opticalFromMatrix(matrix = [], shup = None):
    """
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(matrix) is list or type(matrix) is np.ndarray):
        print(Fore.RED + "myStokes.opticalFromMatrix(): Argument 'matrix' must be a list or array of dimension 4x4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    if type(matrix) is list: matrix = np.array(matrix)
    if not len(matrix) == 4:
        print(Fore.RED + "myStokes.opticalFromMatrix(): Argument 'matrix' must be a list or array of dimension 4x4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    rs = np.sqrt(matrix[0][0] + matrix[0][1])
    rp = np.sqrt(matrix[0][0] - matrix[0][1])
    delta = np.rad2deg(np.arccos(matrix[2][2] / (rs*rp)))
    if not shup:
        print("The optical parameters making up the matrix are: rs = {0:.3f}, rp = {1:.3f}, delta = {2:.2f}°".format(rs, rp, delta))
    return rs, rp, delta




def rotateMatrix( matrix = [], angle = 0., shup = None):
    """
    Rotates a Mueller matrix (argument 'matrix' is a 4x4 list or array) by an angle
    given by argument 'eta' (in degrees).
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(matrix) is list or type(matrix) is np.ndarray):
        print(Fore.RED + "myStokes.rotateMatrix(): Argument 'matrix' must be a list or array of dimension 4x4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    if type(matrix) is list: matrix = np.array(matrix)
    if not len(matrix) == 4:
        print(Fore.RED + "myStokes.rotateMatrix(): Argument 'matrix' must be a list or array of dimension 4x4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    try:
        #angle = -angle
        angle = np.mod(angle, 360.)
    except:
        print(Fore.RED + "myStokes.rotateMatrix(): Argument 'angle' must be a number in degrees." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    return np.around(rotationMatrix(-angle).dot( matrix.dot(rotationMatrix(angle)) ), 12)


def rotationMatrix(angle = 0., shup = None):
    """
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    try:
        angle = np.mod(angle, 360.)
    except:
        print(Fore.RED + "myStokes.rotationMatrix(): Argument 'angle' must be a number in degrees." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    eta = np.deg2rad(angle)
    return np.array([[1., 0., 0., 0.], 
                     [0., np.cos(2.*eta), np.sin(2.*eta), 0.], 
                     [0., -np.sin(2.*eta), np.cos(2.*eta), 0.], 
                     [0., 0., 0., 1.]])
    
    


def inverseMatrix(matrix = [], shup = None):
    """
    Returns the inverse Mueller matrix.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    if not (type(matrix) is list or type(matrix) is np.ndarray):
        print(Fore.RED + "myStokes.inverseMatrix(): Argument 'matrix' must be a 2d list or array of size 4x4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    if type(matrix) is list: matrix = np.array(matrix)
    if not len(matrix) == 4:
        print(Fore.RED + "myStokes.inverseMatrix(): Argument 'matrix' must be a 2d list or array of size 4x4." + Fore.RESET)
        return np.ones([4,4]) * np.NaN
    #
    if np.linalg.det(matrix) != 0:
        return np.around(np.linalg.inv(matrix),12)
    else:
        print(Fore.RED + "myStokes.inverseMatrix(): The passed matrix has no inverse." + Fore.RESET)
        return np.ones([4,4]) * np.NaN

def invertMatrix(matrix = [], shup = False):
    """
    Returns the inverse Mueller matrix.
    """
    return inverseMatrix(matrix = matrix, shup = shup)





def randomMatrix(shup = None):
    """
    Returns a Mueller matrix with random (but mostly reasonable) optical parameters
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    r = uniform(0.2, 0.99)
    rs = 1
    rp = r * rs
    delta = (1 - r) * 60 + uniform(-5, 5)
    matrix = matrixFromOptical(rs, rp, delta, shup = shup)
    if not shup:
        print('Mirror with optical parameters rs = {0:.3f}, rp = {1:.3f}, delta = {2:.2f}°'.format(rs, rp, delta))
    return matrix






# ====================================================================================
# ==================================================================================== F I E L D S
# ====================================================================================





def stokes2fields(stokes = [1,1,0,0], periods = 3):
    """
    At the moment I can't figure out how to distinguish the signs for S2 and S3 if they are both non-zero.
    """
    if not (type(stokes) is list or type(stokes) is np.ndarray):
        print(Fore.RED + "myStokes.stokes2fields(): Argument stokes must be a Stokes vector (length 4)." + Fore.RESET); return
    S = np.array(stokes)
    if not len(S) == 4:
        print(Fore.RED + "myStokes.stokes2fields(): Argument stokes must be a Stokes vector (length 4)." + Fore.RESET); return
    S = normalizeStokes(S, shup = True)
    #
    if S[3] == 0:
        #Ex0 = np.sqrt(0.5 * (S[0] + S[1]))
        #Ey0 = np.sqrt(0.5 * (S[0] - S[1])) * np.sign(S[2])
        #d = 0
        Ex0 = np.sqrt(0.5 * (S[0] + S[1]))
        Ey0 = np.sqrt(0.5 * (S[0] - S[1]))
        if np.sign(S[2]) >= 0.: d = 0.
        else: d = np.pi
    else:
        sgn2, sgn3 = np.sign(S[2]), np.sign(S[3])
        if sgn2 == 0: sgn2 = 1
        if sgn3 == 0: sgn3 = 1
        if S[1] == 0:
            Ex0 = 1/np.sqrt(2)
            Ey0 = 1/np.sqrt(2) * sgn2*sgn3           
        else:
            Ex0 = np.sqrt(0.5 * (S[0] + S[1]))
            Ey0 = np.sqrt(0.5 * (S[0] - S[1])) * sgn2*sgn3
        #Ex0 = np.sqrt(0.5 * (S[0] + S[1]))
        #Ey0 = np.sqrt(0.5 * (S[0] - S[1])) * np.sign(S[2])
        d2 = np.arccos(S[2]/(round(2*Ex0*Ey0,6)))
        d3 = np.arcsin(S[3]/(round(2*Ex0*Ey0,6)))
        d = d2

    try: periods = float(abs(periods))
    except:
        print(Fore.RED + "myStokes.stokes2fields(): Argument must be a positive integer. Setting periods = 5." + Fore.RESET)
        periods = 5
    z = np.linspace(0, periods * 360, int(periods * 360 + 1))
    Ex = Ex0 * np.sin(np.deg2rad(z))
    Ey = Ey0 * np.sin(np.deg2rad(z) - d)
    return {"Ex0": Ex0, "Ey0": Ey0, "d": np.rad2deg(d), "z": z, "Ex": Ex, "Ey": Ey ,"S": S}


def fields2stokes(Ex0 = 1., Ey0 = 0., d = 0):
    try:
        Ex0, Ey0, d = float(Ex0), float(Ey0), float(d)
    except:
        print(Fore.RED + "myStokes.fields2stokes(): Arguments Ex0, Ey0, and d (degrees) must be floats." + Fore.RESET)
        return {}
    S0 = Ex0**2 + Ey0**2
    S1 = Ex0**2 - Ey0**2
    S2 = 2 * Ex0 * Ey0 * np.cos(np.deg2rad(d))
    S3 = 2 * Ex0 * Ey0 * np.sin(np.deg2rad(d))
    return np.array([S0, S1, S2, S3])



# ====================================================================================
# ==================================================================================== P L O T   S T U F F
# ====================================================================================





def plotEllipse(stokes = [1,0,0,1], ax = None, norm = True, color = True, alpha = True,
               lw = 1, legend = True, figsize = (3,3), title = '', stokesastitle = True, 
               stokesaslabel = False, ownlabel = [],
               title_fontsize = 12, legend_fontsize = 8, tick_fontsize = 10):
    """
    stokes:          A Stokes vector or a list of Stokes vectors. A Stokes vector has length 4. 
                     The first element is ignored if norm == True.
    ax:              AxesSubplot. Will be created if ax == None.
    norm:            Default True. Normalizes s1, s2, and s3 to 1.
    color:           Different colors on the ellipses if True.
    alpha:           Different alpha on the...
    lw:              Line width.
    legend:          Shows legend if True. Labels are 1, 2, 3,... 
    figsize:         Size of the figure if ax == None.
    title:           Title of the AxesSubplot.
    stokesastitle:   Sets the Stokes vector as title if True (and title == '')
    stokesaslabel:   Use the Stokes vectors as labels for the legend.
    ownlabel:        A list of strings to be used as labels. Ignored if ownlabel == []. Arg. stokesaslabel should be False.
    title_fontsize:  Fontsize for title
    legend_fontsize:
    tick_fontsize:   
    """
    #
    errMsg1 = Fore.RED + "myStokes.plotEllipse(): Argument 'stokes' must be a list of Stokes vectors, each with size 4." + Fore.RESET
    errMsg2 = Fore.RED + "myStokes.plotEllipse(): Every Stokes vector in the list 'stokes' must be a list with size 4." + Fore.RESET
    if not type(stokes) is list:
        print(errMsg1); return ax
    if not len(np.shape(stokes)) in [1,2]:
        print(errMsg1); return ax
    if len(stokes) == 4 and len(np.shape(stokes)) == 1:
        stokes = [stokes]
    if not len(np.shape(stokes)) >= 2:
        print(errMsg1); return ax
    for s in stokes:
        if not len(s) == 4:
            print(errMsg2); return ax
    #
    if type(ax) is type(None):
        fig, ax = plt.subplots(figsize = figsize)
    #
    ax.set_aspect('equal')
    ax.spines['left'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('center')
    ax.spines['top'].set_color('none')
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xticks([-1, -0.5, 0.5, 1])
    ax.set_yticks([-1, -0.5, 0.5, 1])
    #
    colors = [ cm.jet(x) for x in np.linspace(0.2, 0.8, len(stokes)) ]
    alphas = np.linspace(0.9, 1.0, len(stokes))
    #
    Stokes = []
    for i, s in enumerate(stokes):
        s0_ = np.sqrt(s[1]**2 + s[2]**2 + s[3]**2)
        if norm: s0 = s0_
        else:
            if s[0] < s0_:
                print(Fore.RED + 'myStokes.plotEllipse(): The first element in Stokes vector has an unphysical value. Correcting it.' + Fore.RESET)
                s0 = s0_
            else:
                s0 = np.sqrt(s[0])
        S = list(np.array(s)/s0)
        S[0] = 1
        Stokes.append(S)
    #
    if not type(ownlabel) is list:
        print(Fore.RED + "myStokes.plotEllipse(): The argument 'ownlabel' must be a list of strings. Ignoring this argument." + Fore.RESET)
        ownlabel = []
    if len(ownlabel) > 0:
        olok = True
        for ol in ownlabel:
            if not type(ol) is str: olok = False
        if not len(ownlabel) == len(stokes): olok = False
        if not olok:
            print(Fore.RED + "myStokes.plotEllipse(): The argument 'ownlabel' must be a list of strings (w. the same length as 'stokes'). Ignoring this argument." + Fore.RESET)
            ownlabel = []
    #
    for i, s in enumerate(Stokes):
        PL = np.sqrt(s[1]**2 + s[2]**2)
        A = np.sqrt(0.5 * (1. + PL))
        B = np.sqrt(0.5 * (1. - PL))
        Psi = 0.5 * np.arctan2(s[2],s[1]) * np.rad2deg(1.)
        if color: co = colors[i]
        else: co = 'tab:blue'
        if alpha: al = alphas[i]
        else: al = 1
        if stokesaslabel:
            label = stokes2string(s, 3, False)
        elif len(ownlabel) > 0:
            label = ownlabel[i]
        else:
            label = f"{i+1}"
        ell = Ellipse((0, 0), 2*A, 2*B, angle = Psi, facecolor = 'none', alpha = al, edgecolor = co, lw = lw, label = label)
        ell.set_clip_box(ax.bbox)
        ax.add_artist(ell)
    if legend and len(stokes) > 1 and color or stokesaslabel:
        ax.legend(fontsize = legend_fontsize)
    if stokesastitle and len(stokes) == 1:
        ax.set_title( "({0:6.3f}, {1:6.3f}, {2:6.3f}, {3:6.3f}, )".format(S[0], S[1], S[2], S[3]), fontsize = title_fontsize )
    if not title == '':
        ax.set_title(title, fontsize = title_fontsize)
    plt.tight_layout(pad=2.5)
    plt.yticks(fontsize = tick_fontsize)
    plt.xticks(fontsize = tick_fontsize)
    #
    return ax
    








