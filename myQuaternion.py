"""
This module was written to get pitch, roll, yaw, and translation values for
the mirror and grating when misaligning the whole mono tank, using quaterions
instead of euler angles. It can of course be used for anything.

class:

quaterion()     mimics a quaterion
                methods: conjugate, inverse, norm, normalize, scalar, vector, copy,
                            rotate, rotationMatrix, RxRyRz, axisAngle
                operations: +, -, *

methods:

rquaternion()                   returns a rotation quaternion.
Rx(), Ry(), Rz()                rotation matrices.
Rzyx()                          roation matrix for Rz*Ry*Rx
rotationMatrixFromVectors()     returns the rotation matrix for rotating a vector from v to u.
anglesFromRotationMatrix()      get the angles for Rz*Ry*Rx from a rotation matrix.
decimals()                      clean up leftover decimals for quaternions, matrices, lists, and floats.
"""

_version_ = "2023.09.30"
_author_  = "Mats Leandersson"


import numpy as np
from numpy import sqrt, pi, cos, sin, tan, arccos, arcsin, arctan, arctan2
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from copy import deepcopy


class quaternion():
    """
    Quaternion class.
    
    Supports the addition, subtraction, and multipplication operators (+,-,*).
    
    Can return:
        conjugate()         the conjugate quaternion, as quaternion
        norm()              the norm of the quaternion, as scalar
        normalize()         the normalized quaternion, as quaternion
        inverse()           the inverse quaternion, as quaternion, (the normalized conjugate)
        scalar()            the scalar component of the quaternion, as scalar
        vector()            the vector component of the quaternion, as array
        copy()              a copy of itself, as quaternion
        rotate()            the rotated quaternion, as quaternion, (rotated around axis (arg.) by angle (arg.))
        rotationMatrix()    the 4x4 rotation matrix, as array
        RxRyRz()            the pitch, roll, and yaw angles, as floats
    
    The quaternion values are stored in:
        q                   an array of length 4, with elements (w,x,y,z), where [x,y,z] is the vector part.

    Note:
        Create a rotation quaternion by using the stand-alone method rquaternion(axis, angle). The method returns
        quaternion(w, x, y, z) where:
            w = cos(angle/2)
            x, y, z = axis/|axis|*sin(angle/2)
    Note 2:
        Quaternions, lists, and arrays can be "cleaned" up a bit by using the stand alone method decimals().
    """

    def __init__(self, w = 0, x = 0, y = 0, z = 0, shup = True):
        """
        Arguments w, x, y, z are the 4 elements in the quaternion. Must be scalars.
        """
        self._is_q = True
        try: w = float(w)
        except: w = 0.
        try: x = float(x)
        except: x = 0.
        try: y = float(y)
        except: x = 0.
        try: z = float(z)
        except: z = 0.
        self.q = np.array([w, x, y, z])
        if not shup:
            print(f"q = ({w}, {x}, {y}, {z})")
    
    # multiplication with another "quaternion"
    def __mul__(self, other):
        """
        Returns the result from Q * q, where Q is this quaternion
        and q is the argument quaternion.
        return quaternion
        """
        try:
            if other._is_q:
                w1, x1, y1, z1 = self.q
                w2, x2, y2, z2 = other.q
                w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
                x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
                y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
                z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
            return quaternion(w, x, y, z)
        except:
            print("Is the second object really a quaternion??")
            return None
    
    # addition with another "quaternion"
    def __add__(self, other):
        """
        Returns the result from Q + q, where Q is this quaternion
        and q is the argument quaternion.
        return quaternion
        """
        try:
            if other._is_q:
                w1, x1, y1, z1 = self.q
                w2, x2, y2, z2 = other.q
                w = w1 + w2
                x = x1 + x2
                y = y1 + y2
                z = z1 + z2
            return quaternion(w, x, y, z)
        except:
            print("Is the second object really a quaternion??")
            return None
    
    # subtraction with another "quaternion"
    def __sub__(self, q):
        """
        Returns the result from Q - q, where Q is this quaternion
        and q is the argument quaternion.
        return quaternion
        """
        try:
            if q._is_q:
                w1, x1, y1, z1 = self.q
                w2, x2, y2, z2 = q.q
                w = w1 - w2
                x = x1 - x2
                y = y1 - y2
                z = z1 - z2
            return quaternion(w, x, y, z)
        except:
            print("Is the second object really a quaternion??")
            return None
    
    # returns the conjugate of the quaternion
    def conjugate(self):
        """
        Returns the conjugate of the quaternion.
        return quaternion
        """
        return quaternion(self.q[0], -self.q[1], -self.q[2], -self.q[3])

    # returns the quaternion's norm
    def norm(self):
        """
        Returns the norm of the quaternion.
        return float
        """
        return sqrt(self.q[0]**2 + self.q[1]**2 + self.q[2]**2 + self.q[3]**2)
    
    # returns the normalized quaternion
    def normalize(self):
        """
        Returns the normalized quaternion.
        return quaternion
        """
        n = self.norm()
        return quaternion(self.q[0]/n, self.q[1]/n, self.q[2]/n, self.q[3]/n)
    
    # returns the inverse of the quaternion
    def inverse(self):
        """
        Returns the inverse of the quaternion (i.e. the noramlized conjugate).
        return quaternion
        """
        q = self.conjugate()
        q = q.normalize()
        return q
        
    # returns the "squalar" part of the quaternion
    def scalar(self):
        """
        Returns the scalar part of the quaternion.
        return float
        """
        return self.q[0]
    
    # returns the "vector" part of the quaternion
    def vector(self):
        """
        Returns the vector part of the quaternion.
        return numpy.array
        """
        return self.q[1:]
    
    # returns a copy of the quaternion
    def copy(self):
        """
        Returns a copy of itself.
        return quaternion
        """
        return deepcopy(self)
    
    # return rotated quaternion
    """
    Returns the quaternion rotated angle radians around axis.
    return quaternion
    """
    def rotate(self, axis = np.array([0,0,1]), angle = 0):
        if not (type(axis) is list or type(axis) is np.ndarray):
            print("axis must be an array or list of length 3"); return None
        axis = np.array(axis)
        if not len(axis) == 3:
            print("axis must be an array or list of length 3"); return None
        try: angle = float(angle)
        except: print("angle must be a number (radians)"); return None
        #
        qrscalar = cos(angle/2)
        n = sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
        qrvector = axis/n * np.sin(angle/2)
        qr = quaternion(qrscalar, qrvector[0], qrvector[1], qrvector[2])
        return qr * self * qr.conjugate()
    
    # return rotation matrix
    def rotationMatrix(self, d = 4, r = None):
        """
        Returns the 4x4 (or 3x3 if d = 3) rotation matrix represented by the quaternion.
        Pass r (integer, number of decimal places) to clean up the result.
        return numpy.array
        """
        try: d = int(abs(d))
        except:
            print("must be 4 or 3 (to return the 4x4 or 3x3 rotation matrix). setting default d = 4.")
            d = int(4)
        if not d in [3,4]:
            print("must be 4 or 3 (to return the 4x4 or 3x3 rotation matrix). setting default d = 4.")
            d = int(4)
        if type(r) is type(None): pass
        elif type(r) is int: r = int(abs(r))
        else:
            print("r must be an integer or None. setting default r = None.")
            r = None
        #
        q = self.normalize()
        w, x, y, z = q.q[0], q.q[1], q.q[2], q.q[3]
        m = [[1-2*y**2-2*z**2, 2*x*y-2*w*z,     2*x*z+2*w*y,     0],
             [2*x*y+2*w*z,     1-2*x**2-2*z**2, 2*y*z-2*w*x,     0],
             [2*x*z-2*w*y,     2*y*z+2*w*x,     1-2*x**2-2*y**2, 0],
             [0,               0,               0,               1]]
        m = np.array(m)
        if not type(r) is type(None):
            m = m.flatten()
            for i, e in enumerate(m): m[i] = round(e,r)
            m = m.reshape([4,4]) 
        if d == 4: return np.array(m)
        else: return np.array(m)[:3,:3]
    
    # return the pitch, roll, and yaw
    def RxRyRz(self):
        """
        Returns the Rx, Ry, and Rz angles represented by the quaternion.
        return float, float, float
        """
        w, x, y, z = self.q[0], self.q[1], self.q[2], self.q[3]
        Rx   = np.arctan2(2.0*(y*z + w*x), w*w - x*x - y*y + z*z)
        Ry = np.arcsin(-2.0*(x*z - w*y))
        Rz  = np.arctan2(2.0*(x*y + w*z), w*w + x*x - y*y - z*z)
        return Rx, Ry, Rz
    
    # get rotation axis and angle
    def axisAngle(self):
        """
        Returns the rotation axis and the rotation angle of the quaternion.
        """
        q = self.normalize()
        if q.q[0] == 1:
            print("this quaternion is not representing a rotaion (w = 1)"); return
        angle = 2*np.arccos(q.q[0])
        axis = q.q[1:]/np.sin(angle)
        axis = axis/np.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
        return axis, angle

def rquaternion(axis = [0,0,1], angle = 0):
    """
    Make a "rotaion" quaternion. Returns a quaternion according to:
        quaternion(w, x, y, z)
    where:
        w = cos(angle/2)
        x, y, z = axis/|axis|*sin(angle/2)
    """
    try: angle = float(angle)
    except:
        print("argument angle must be a number (radians). setting default angle = 0.")
        angle = 0
    if not (type(axis) is list or type(axis) is np.ndarray):
        print("argument axis must be a list or array. setting default axis = [0,0,1].")
        axis = [0,0,1]
    if not len(axis) == 3:
        print("argument axis must have length 3. setting default axis = [0,0,1].")
        axis = [0,0,1]
    #
    axis = np.array(axis)
    norm = sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    axis = axis/norm
    #
    qrscalar = cos(angle/2)
    qrvector = axis * np.sin(angle/2)
    return quaternion(qrscalar, qrvector[0], qrvector[1], qrvector[2])

def decimals(o = None, d = 6):
    """
    Use to clean up values in quaternions and arrays after calculations where there
    are e.g. left 1E-16 and stuff in the elements.
    Argument o is an object, either a quaternion, array, or list.
    Returns a new object. The passed object is unaffected.
    """
    try: d = int(abs(d))
    except:
        print("d must be an integer. setting default d = 6.")
        d = int(6)
    if type(o) is quaternion:
        q = o.copy()
        #for i, e in enumerate(q.q): q.q[i] = round(e, d)
        return np.around(q.q, d)
    elif type(o) is float:
        return round(o, d)
    elif type(o) is list:
        return list(np.around(o, d))
    elif type(o) is np.ndarray:
        return np.around(o, d)
    else:
        print("this method works for quaternions, lists, arrays, and floats.")
        print(f"the object in the passed argument o is a {type(o)}.")
        return
    
def Rx(a = 0):
    """3x3 matrix for rotation around x (pitch)"""
    m = [[1, 0,          0        ],
         [0, np.cos(a), -np.sin(a)],
         [0, np.sin(a),  np.cos(a)]]
    return np.array(m)

def Ry(a = 0):
    """3x3 matrix for rotation around y (roll)"""
    m = [[ np.cos(a), 0,  np.sin(a)],
         [ 0,         1,  0        ],
         [-np.sin(a), 0,  np.cos(a)]]
    return np.array(m)

def Rz(a = 0):
    """3x3 matrix for rotation around z (yaw)"""
    m = [[np.cos(a), -np.sin(a), 0],
         [np.sin(a),  np.cos(a), 0],
         [0,          0,         1]]
    return np.array(m)

def Rzyx(ax = 0, ay = 0, az = 0):
    """
    Rotation matrix for RzRyRx
    """
    R11 = cos(ay) * cos(az)
    R12 = sin(ax) * sin(ay) * cos(az) - cos(ax) * sin(az)
    R13 = cos(ax) * sin(ay) * cos(az) + sin(ax) * sin(az)
    R21 = cos(ay) * sin(az)
    R22 = sin(ax) * sin(ay) * sin(az) + cos(ax) * cos(az)
    R23 = cos(ax) * sin(ay) * sin(az) - sin(ax) * cos(az)
    R31 = -sin(ay)
    R32 = sin(ax) * cos(ay)
    R33 = cos(ax) * cos(ay)
    m = [[R11, R12, R13],
         [R21, R22, R23],
         [R31, R32, R33]]
    return np.array(m)

def rotationMatrixFromVectors(vec1, vec2, shup = False):
    """ 
    vec1: "source" vector
    vec2: "destination" vector
    return mat: matrix which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    if not shup:
        print("Rotation matrix for v = R * u")
        print(f"u = {vec1}")
        print(f"v = {vec2}")
        print(f"R = \n{rotation_matrix}\n")
    return rotation_matrix

def anglesFromRotationMatrix(M, shup = False):
    """
    From a rotation matrix, get the angles for Rx, Ry, and Rz to produce the same rotation.
    rotated_vector = Rz * Ry * Rx * vector
    """
    ay = -arcsin(M[2][0])
    ax = arctan2(M[2][1]/cos(ay), M[2][2]/cos(ay))
    az = arctan2(M[1][0]/cos(ay), M[0][0]/cos(ay))
    if not shup:
        print("Rotation angles for Rx, Ry, Rz:")
        print(f"ax = {ax:9.6f} / {np.rad2deg(ax):11.8f}°")
        print(f"ay = {ay:9.6f} / {np.rad2deg(ay):11.8f}°")
        print(f"az = {az:9.6f} / {np.rad2deg(az):11.8f}°\n")
    return ax, ay, az