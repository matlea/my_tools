"""
Cleaning up my old MonoSafety stuff.
Removing unimportant stuff.
Mats Leandersson, 2021-10-19


Versions
    2023-09-02  Added the possibility to get params from drawMono() w.o. plotting.
    2023-08-31  Added zero order light to drawMono().
                Moved blaze() to here from myMono().
                Deleted the old drawMono_().
    2023-08-07  New version of drawMono(). The old version exists as drawMono_().
    2021-10-19  Copied and cleaned up stuff from old scripts.
"""


import numpy as np
from numpy import pi

from matplotlib import pyplot as plt

#from PyQt5.QtWidgets import QDesktopWidget
#from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHBoxLayout
#from PyQt5.QtWidgets import QLabel, QLineEdit, QCheckBox, QComboBox

from colorama import Fore



__version__ = "2023-09-02"
__author__  = "Mats Leandersson"

# ==============================================================================================
print(Fore.LIGHTWHITE_EX + f"--- myMonoExtra, version {__version__}" + Fore.BLACK)
# ==============================================================================================

# ===================================== Mono safety

class BlochMono():
    myType = 'BlochMono'

    # Constamnts
    Planck = 4.135667e-15
    LightSpeed = 299792458
    
    # Parameters defining the monochromator.
    # Notations accordig to R.Sakari (A-G).
    monoParameters = {
                    'mpA' :  64.819,
                    'mpB' :  -0.125,
                    'mpC' :  43.821,
                    'mpD' :  20.000,
                    'mpE' : 570.000,    # mirror length
                    'mpF' :  42.000,    # beam vertical translation
                    'mpG' : 140.000,
                    'minM' : 0.,
                    'maxM' : 20.,
                    'minG' : 0.,
                    'maxG' : 30.,
                    'ES'   : 10000.,
                    'LD'   : 800.
                    }
    # Collision swith parameters.
    # Notation according to me.
    collisionSwitchParameters = {
                    'CSL1'  : 175.82,
                    'CSL2'  : 41.92,
                    'CSL3'  : 50,
                    'CSL4'  : 64+80, 
                    'CSR'   : 0.,   
                    'CSlim' : 0.
                    }
    # Motor parameters
    # Notation according to me
    motorParameters = {
                        'Mo1_Y' : 944.82,       # y-distance from the mirror rotation center to motor 1's rotation center
                        'Mo1_X' : 425,          # x-distance from...
                        'Mo1_R' : 600,          # The mirror's rotation radius
                        'Mo1_min' : 742.18,     # Minimum length of the motor axis
                        'Mo1_max' : 960.86,     # Maximun...
                        'Mo1_offset' : 0,       # Value for angle = 0
                        'Mo2_Y' : 922,          # y-distance from the grating rotation center to tmotor 2's rotation center
                        'Mo2_X' : 425,          # x-distance from...
                        'Mo2_R' : 600,          # The grating's rotation radius
                        'Mo2_min' : 629.15,     # Minimum length of the motor axis
                        'Mo2_max' : 938.46,     # Maximun...
                        'Mo2_offset' : 0        # Value for angle = 0
                        }

    monoParameters_bkup = {}
    collisionSwitchParameters_bkup = {}
    motorParameters_bkup = {}


    def __init__(self, line_density = 800):

        self.monoParameters_bkup = dict(self.monoParameters)
        self.collisionSwitchParameters_bkup = dict(self.collisionSwitchParameters)
        self.motorParameters_bkup = dict(self.motorParameters)

        self.motorParameters['Mo1_offset'] = self.motorParameters['Mo1_min']
        self.motorParameters['Mo2_offset'] = self.motorParameters['Mo2_min']
        self._initStuff()

    
    def _initStuff(self):
        pass
    

    def setParams(self, shup = False, **kwargs):
        """
        Set the mono, collision switch, and motor parameters.
        E.g. .setParams(LD = 92)

        """
        if len(kwargs) == 0:
            print("BlochMono.setParams(): pass at least one argument...")
            return

        for item in kwargs.items():
            
            pItems = self.monoParameters.items()
            Items = []
            for pI in pItems: Items.append(pI[0])
            if item[0] in Items:
                self.monoParameters[item[0]] = item[1]
                if not shup:
                    print("BlochMono.set(): .{0} = {1}".format(item[0], item[1]))
            
            pItems = self.collisionSwitchParameters.items()
            Items = []
            for pI in pItems: Items.append(pI[0])
            if item[0] in Items:
                self.collisionSwitchParameters[item[0]] = item[1]
                if not shup:
                    print("BlochMono.set(): .{0} = {1}".format(item[0], item[1]))
            
            pItems = self.motorParameters.items()
            Items = []
            for pI in pItems: Items.append(pI[0])
            if item[0] in Items:
                self.motorParameters[item[0]] = item[1]
                if not shup:
                    print("BlochMono.set(): .{0} = {1}".format(item[0], item[1]))
            
            self._initStuff()

    def showParams(self):
        print('monoParameters:')
        items = self.monoParameters.items()
        for item in items: print("{0}\t{1}".format(item[0],item[1]))
        print('collisionSwitchParameters:')
        items = self.collisionSwitchParameters.items()
        for item in items: print("{0}\t{1}".format(item[0],item[1]))
        print('motorParameters:')
        items = self.motorParameters.items()
        for item in items: print("{0}\t{1}".format(item[0],item[1]))

    def getParams(self, *args):
        """
        Get a selection of the mono, collision switch, and motor parameters as a dict.
        E.g. .getParams('LD', 'CSlim') returns {'LD' : 800, 'CSlim' : 1}
        If no argument/arguments is/are passed the method prints a list of all parameters.
        
        """
        if len(args) == 0:
            self.showParams()
            return

        rets = {}
        for item in args:            
            pItems = self.monoParameters.items()
            Items = []
            for pI in pItems: Items.append(pI[0])
            if item in Items:
                rets[item] = self.monoParameters[item]
            pItems = self.collisionSwitchParameters.items()
            Items = []
            for pI in pItems: Items.append(pI[0])
            if item in Items:
                rets[item] = self.collisionSwitchParameters[item]            
            pItems = self.motorParameters.items()
            Items = []
            for pI in pItems: Items.append(pI[0])
            if item in Items:
                rets[item] = self.motorParameters[item]
            
        return rets
    
    def resetParams(self, shup = False):
        self.monoParameters = dict(self.monoParameters_bkup)
        self.collisionSwitchParameters = dict(self.collisionSwitchParameters_bkup)
        self.motorParameters = dict(self.motorParameters_bkup)
        self._initStuff()
        if not shup:
            print('BlochMono.reset(): All parameters set to default values (monoParameters, collisionSwitchParameters, motorParameters)')
        
        


    # ===================== Energy, cff, diffraction, ...

    def EC2A(self, energy = None, cff = None, ld = None, shup = False, silence = False):
        """
        Energy and cff -> Mirror and grating angles (returns [mirror, grating], True/False)

        (silence is a legacy name that is changed to shup for consistency.)

        """
        if silence: shup = True
                
        if type(energy) is type(None) or type(cff) is type(None):
            if not shup: print('Pass energy and cff.')
            return [None, None], False
        if not type(ld) is type(None): 
            self.monoParameters['LD'] = ld
        LD = self.monoParameters['LD'] *1e3

        try:
            Lambda = self.Planck * self.LightSpeed / energy
            Alpha = np.arccos( Lambda * LD / (cff**2 - 1) * ( np.sqrt( cff**2 + (cff**2 -1)**2 / Lambda**2 / LD**2 )-1 ))
            Beta = np.arcsin( cff * np.sin(Alpha) )
            Theta = 0.5  * (Alpha + Beta)
        except:
            if not shup: print('Probably invalid first order diffraction conditions.')
            return [None, None], False
        
        M = np.rad2deg(Theta)
        G = np.rad2deg(Beta)
        ok = True

        minM = self.monoParameters['minM']
        maxM = self.monoParameters['maxM']
        minG = self.monoParameters['minG']
        maxG = self.monoParameters['maxG']
        if M < minM or M > maxM or G < minG or G > maxG: 
            ok = False
        if not shup and not ok:
            print('Angle(s) out of range.')
        return [M, G], ok
    

    def A2EC(self, mirror = None, grating = None, ld = None, shup = False, silence = False):
        """
        Mirror and grating angles -> Energy and cff

        (silence is a legacy name that is changed to shup for consistency.)

        """
        if silence: shup = True

        minM = self.monoParameters['minM']
        maxM = self.monoParameters['maxM']
        minG = self.monoParameters['minG']
        maxG = self.monoParameters['maxG']
        if type(mirror) is type(None) or type(grating) is type(None) or mirror < minM or mirror > maxM: # fix this mess
            if mirror < minM or mirror > maxM:
                if not shup: print('Mirror angle out of range.')
            if grating < minG or grating > maxG:
                if not shup: print('Grating angle out of range.')
            return [np.NaN, np.NaN], False
        
        if not type(ld) is type(None):
            self.monoParameters['LD'] = ld
        else:
            ld = self.monoParameters['LD']
        LD = ld * 1e3

        try:
            a = np.deg2rad(2*mirror - grating)
            b = np.deg2rad(grating)
            #t = np.deg2rad(mirror) #not used
            if a == b:
                # zero order
                return [np.NaN, np.NaN], False
            if np.sin(a) == 0:
                return [np.NaN, np.NaN], False
            if (np.cos(a) - np.cos(b)) == 0:
                return [np.NaN, np.NaN], False
            Lambda = (1/LD) * (np.cos(a) - np.cos(b))
            Energy = self.Planck * self.LightSpeed / Lambda
            Cff = np.sin(b) / np.sin(a)
        except:
            if not shup: print('Probably invalid first order diffraction conditions.')
            return [np.NaN, np.NaN], False
        
        if Energy <= 0 or Cff <=1:
            if not shup: 
                print('Invalid first order diffraction conditions.')
                print('( energy = {0:}eV, cff = {1})'.format(round(Energy,2), round(Cff,2)))
            return [Energy, Cff], False

        return [Energy, Cff], True


    def M2A(self, M = 0, L = 0, shup = False):
        """
        Motor position (axis length) to angle.
        Arguments: M = motor number (1=mirror, 2=grating), L = axis length
        """
        if M == 1: # mirror
            a = self.motorParameters['Mo1_Y']
            b = self.motorParameters['Mo1_X']
            r = self.motorParameters['Mo1_R']
            offset = self.motorParameters['Mo1_offset']
            Lmin = self.motorParameters['Mo1_min'] - offset
            Lmax = self.motorParameters['Mo1_max'] - offset
        elif M == 2: # grating
            a = self.motorParameters['Mo2_Y']
            b = self.motorParameters['Mo2_X']
            r = self.motorParameters['Mo2_R']
            offset = self.motorParameters['Mo2_offset']
            Lmin = self.motorParameters['Mo2_min'] - offset
            Lmax = self.motorParameters['Mo2_max'] - offset
        else:
            if not shup:
                print("BlochMono.M2A(): Argument 'M' should be 1 or 2 (mirror or grating)")
            return np.NaN
        if not (Lmin <= L and L <= Lmax ):
            if not shup:
                print("BlochMono.M2A(): Argument 'L' is outside the parameter limits. Lmin = {0} mm, Lmax = {1} mm".format(Lmin, Lmax))
            return np.NaN
        c = 1/(2*r) * (b**2 + a**2 + r**2 - (L + offset)**2)
        angle =  -np.arccos(c/np.sqrt(a**2+b**2)) + np.arctan(a/b)
        return np.rad2deg(angle)
    
    def A2M(self, M = 0, A = 0, shup = False):
        """
        Angle to motor position (axis length).
        Arguments: M = motor number (1=mirror, 2=grating), A = angle

        """       
        if M == 1: # Mirror
            a = self.motorParameters['Mo1_Y']
            b = self.motorParameters['Mo1_X']
            r = self.motorParameters['Mo1_R']
            offset = self.motorParameters['Mo1_offset']
            Lmin = self.motorParameters['Mo1_min'] - offset
            Lmax = self.motorParameters['Mo1_max'] - offset
            Amin = self.monoParameters['minM']
            Amax = self.monoParameters['maxM']
        elif M == 2: # Grating
            a = self.motorParameters['Mo1_Y']
            b = self.motorParameters['Mo1_X']
            r = self.motorParameters['Mo1_R']
            offset = self.motorParameters['Mo1_offset']
            Lmin = self.motorParameters['Mo1_min'] - offset
            Lmax = self.motorParameters['Mo1_max'] - offset
            Amin = self.monoParameters['minG']
            Amax = self.monoParameters['maxG']
        else:
            if not shup:
                print("BlochMono.A2M(): Argument 'M' should be 1 or 2 (mirror or grating)")
            return np.NaN
        if not ( Amin <= A and A <= Amax):
            if not shup:
                print("BlochMono.A2M(): Argument 'A' is outside the parameter limits. Amin = {0} mm, Amax = {1} mm".format(Amin, Amax))
            return np.NaN

        Ar = np.radians(A)
        L = np.sqrt(b**2 + a**2 + r**2 - 2*r * (b * np.cos(Ar) + a * np.sin(Ar)))
        if not ( Lmin <= (L - offset) and (L - offset) <= Lmax ):
            if not shup:
                print("BlochMono.A2M(): the result {0:5.1f} is outside the parameter limits. Lmin = {1} mm, Lmax = {2} mm".format(L-offset, Lmin, Lmax))
            return np.NaN
        else:
            return L - offset

    

    # ===================== Geometry

    def _point1(self, angle = 0):    # mirror, point 1
        t = np.deg2rad(angle)
        mpC = self.monoParameters['mpC']
        x = mpC * np.cos(0.75 * 2*pi + t)
        y = mpC * np.sin(0.75 * 2*pi + t)
        return x, y

    def _point2(self, angle = 0):    # mirror, point 2
        t = np.deg2rad(angle)
        CSL1 = self.collisionSwitchParameters['CSL1']
        x1, y1 = self._point1(angle = angle)
        x2 = x1 + CSL1 * np.cos(0.5 * 2*pi + t) 
        y2 = y1 + CSL1 * np.sin(0.5 * 2*pi + t)
        return x2, y2

    def _point3(self, angle = 0):    # mirror, point 3
        t = np.deg2rad(angle)
        CSL2 = self.collisionSwitchParameters['CSL2']
        x2, y2 = self._point2(angle = angle)
        x3 = x2 + CSL2 * np.cos(0.25 * 2*pi + t)
        y3 = y2 + + CSL2 * np.sin(0.25 * 2*pi + t)
        return x3, y3
    
    def _point4(self):               # grating, point 4 (center, fixed)
        return self.monoParameters['mpB'], -self.monoParameters['mpA'] + self.monoParameters['mpF']

    def _point5(self, angle = 0):    # grating, point 5
        b = np.deg2rad(angle)
        x4, y4 = self._point4()
        mpG = self.monoParameters['mpG']
        x5 = x4 + 0.5 * mpG * np.cos(0.5 * 2*pi + b)
        y5 = y4 + 0.5 * mpG * np.sin(0.5 * 2*pi + b)
        return x5, y5
    
    def _point5b(self, angle = 0):    # grating, point 5b (not important, just for rawings)
        b = np.deg2rad(angle)
        x4, y4 = self._point4()
        mpG = self.monoParameters['mpG']
        x5b = x4 + 0.5 * mpG * np.cos(0.0 * 2*pi + b)
        y5b = y4 + 0.5 * mpG * np.sin(0.0 * 2*pi + b)
        return x5b, y5b

    def _point6(self, angle = 0):    # grating, point 6
        b = np.deg2rad(angle)
        CSL3 = self.collisionSwitchParameters['CSL3']
        x5, y5 = self._point5(angle = angle)
        x6 = x5 + CSL3 * np.cos(0.25 * 2*pi + b)
        y6 = y5 + CSL3 * np.sin(0.25 * 2*pi + b)
        return x6, y6
    
    def _point7(self, angle = 0):    # grating, point 7 (not in use)
        b = np.deg2rad(angle)
        CSL4 = self.collisionSwitchParameters['CSL4']
        x6, y6 = self._point6(angle = angle)
        x7 = x6 + CSL4 * np.cos(0.5 * 2*pi + b)
        y7 = y6 + CSL4 * np.sin(0.5 * 2*pi + b)
        return x7, y7
    
    def _point8(self, mirror = 0, grating = 0):  # grating, point 8 (the point on the girder closest to the lim. sw.)
        x3, y3 = self._point3(angle = mirror)
        x6, y6 = self._point6(angle = grating)
        k  = np.tan( np.deg2rad(grating) )
        kN = np.tan( np.deg2rad(grating + 90) )
        x8 = 1/(kN-k) * ( (y6-k*x6) - (y3-kN*x3))
        y8 = y6 + k*(x8-x6)
        return x8, y8
    
    def _point9(self, mirror = 0, grating  = 0): # mirror, point 9 (closest to grating)
        x1, y1 = self._point1(angle = mirror)
        x5, y5 = self._point5(angle = grating)
        # Line coefficient
        k = np.tan(np.deg2rad(mirror))
        kN = np.tan(np.deg2rad(mirror+90))
        # Point on the mirror
        x9 = 1/(kN-k) * ((y1-k*x1) - (y5-kN*x5))
        y9 = k * x9 + y1 - k * x1
        return x9, y9
    
    def _point10(self, mirror = 0, grating = 0): # mirror, point 10 (closest to girder)
        x1, y1 = self._point1(mirror)
        x7, y7 = self._point7(grating)
        # Line coefficient
        k = np.tan(np.deg2rad(mirror))  #
        kN = np.tan(np.deg2rad(mirror+90))  #
        # Point on the mirror
        x10 = 1/(kN-k) * ((y1-k*x1) - (y7-kN*x7))
        y10 = k * x10 + y1 - k * x1
        return x10, y10
    
    def _point11(self, mirror = 0): # beginning of the mirror
        t = np.deg2rad(mirror)
        mpE = self.monoParameters['mpE']
        x1, y1 = self._point1(angle = mirror)
        x2 = x1 + mpE * np.cos(0.5 * 2*pi + t) 
        y2 = y1 + mpE * np.sin(0.5 * 2*pi + t)
        return x2, y2


    def _csgap(self, mirror = 0, grating = 0):
        x3, y3 = self._point3(angle = mirror)
        x8, y8 = self._point8(mirror = mirror, grating = grating)
        dx = x3-x8
        dy = y8-y3
        CSR = self.collisionSwitchParameters['CSR']
        CSlim = self.collisionSwitchParameters['CSlim']
        gap = np.sqrt(dx**2 + dy**2)
        if (dx <= 0 or dy <= 0): gap = -gap
        OK = True
        if gap <= CSlim + CSR:        
            OK = False    
        return gap, OK
    
    def _gap(self, mirror = 0, grating = 0):
        """
        gap is a legacy method name, together with _gap. Use .csgap() and ._csgap()
        """
        return self._csgap(mirror = mirror, grating = grating)
    
    def csgap(self, energy = np.NaN, cff = np.NaN, linedensity = np.NaN, shup = False):
        """
        Returns value, [bool, bool] for a given energy, cff, and line density, where
        the value is the distance between the collision switch and grating girder, the
        first bool is the state of the angles values, and the second bool is the state
        of the collision switch.

        If no line density is passed the method uses the defult value.

        
        """
        if not type(shup) is type(True):
            shup = False
        reterr = np.NaN, [False,False]
        if np.isnan(energy):
            if not shup:
                print("BlochMono.csgap() - Pass (1st) attribute 'energy' as a value > 0 (eV).")
            return reterr
        if energy <= 0:
            if not shup:
                print("BlochMono.csgap() - Pass (1st) attribute 'energy' as a value > 0 (eV).")
            return reterr
        if np.isnan(cff):
            if not shup:
                print("BlochMono.csgap() - Pass (2nd) attribute 'cff' as a value > 1.")
            return reterr
        if cff <= 1:
            if not shup:
                print("BlochMono.csgap() - Pass (2nd) attribute 'cff' as a value > 1.")
            return reterr
        if np.isnan(linedensity):
            ld = self.monoParameters['LD']
        if ld <= 1:
            if not shup:
                print("BlochMono.csgap() - Pass (3rd) attribute 'linedensity' as a value > 1 (l/mm).")
            return reterr

        angles, ok1 = self.EC2A(energy = energy, cff = cff, ld = ld, silence = True)
        gap, ok2    = self._gap(mirror = angles[0], grating = angles[1])

        if not shup:
            if not ok1:
                print('BlochMono.csgap() - Angles out of range (mir = {0:6.3f}°, gr = {1:6.3f}°)'.format(angles[0], angles[1]))
            if not ok2:
                CSlim = self.collisionSwitchParameters['CSlim']
                print('BlochMono.csgap() - The col.sw. distance is too small ({0:4.1f} mm, with safety margin {1} mm).'.format(gap, CSlim))
        
        return gap, [ok1,ok2]
    

    def gap(self, energy = np.NaN, cff = np.NaN, linedensity = np.NaN, shup = False):  # LEGACY METHOD NAME
        """
        gap is a legacy method name, together with _gap. Use .csgap() and ._csgap()
        """
        return self.csgap(energy = energy, cff = cff, linedensity = linedensity, shup = shup)


    def quickPlotLimSw(self, energy = np.NaN, cff = np.NaN, fullmirror = False, shup = False):
        """
        This is a legacy method name
        """
        self.plotMonoLims(energy = energy, cff = cff, fullmirror = fullmirror, shup = shup)
    
    def plotMonoLims(self, energy = np.NaN, cff = np.NaN, mirror = np.NaN, grating = np.NaN, fullmirror = False, figscale = 1, shup = False):
        """
        To illustrate the mirror and grating positions relative to eachother, as well as the distances:
        1. closest mirror - grating
        2. collision switch
        3. closest mirror - grating girder

        Note: Draws a figure even if the angles are outside their limits but reports it (unless shup = True)

        Arguments:
            energy = number (>0)
            cff = number (>1)
            fullmirror = False, only draws the nesseccary part of the mirror if False
            figscale = 1 (2 makes the figure twice as big, etc.)
            shup = False, prints error messages if "shut up" is False.
        
            mirror, grating: use angles instead of energy and cff

        """
        if np.isnan(mirror) and np.isnan(grating):
            if np.isnan(energy) or np.isnan(cff):
                print("BlochMono.plotMonoLims() - Pass attributes 'energy' and 'cff' as values (>0,>1)."); return
            if energy <= 0 or cff <= 1:
                print("BlochMono.plotMonoLims() - Pass attributes 'energy' and 'cff' as values (>0,>1)."); return
            gap1, ok1 = self.gap(energy = energy, cff = cff, shup = shup)
            angles, oka = self.EC2A(energy = energy, cff = cff, silence = True)
            if not oka and not shup:
                print("BlochMono.plotMonoLims() - Angles outside limits ({0:5.2f}°, {1:5.2f}°".format(angles[0],angles[1]))
            gap0, ok0 = self._mggap(angles[0],angles[1])
            gap2, ok2 = self._mgirgap(angles[0], angles[1])
        elif np.isnan(mirror) ^ np.isnan(grating):
            if not not shup:
                print("BlochMono.plotMonoLims() - Pass BOTH 'mirror' and 'grating' (as numbers)."); return
        else:
            gap0, _ = self._mggap(mirror,grating)
            gap1, _ = self._csgap(mirror,grating)
            gap2, _ = self._mgirgap(mirror,grating)
            angles = [mirror, grating]
            ld = self.getParams('LD')['LD']
            ec, okec = self.A2EC(mirror = mirror, grating = grating, ld = ld, shup = True)
            if okec:
                energy = ec[0]
                cff = ec[1]
        
        fig = plt.figure(figsize=(12 * figscale,12 * figscale))
        ax = fig.add_subplot(aspect='equal')

        ax.scatter(0, 0, s = 20, color = 'k')
        x1, y1 = self._point1(angle = angles[0]) 
        ax.plot([0,x1], [0,y1], color = 'k', linewidth = 0.75)
        x2, y2 = self._point2(angle = angles[0])
        ax.plot([x1,x2], [y1,y2], color = 'tab:orange', linewidth = 3)
        x3, y3 = self._point3(angle = angles[0])
        ax.plot([x2,x3], [y2,y3], color = 'k', linewidth = 0.75)
        ax.scatter(x3, y3, s = 30, color = 'tab:red')

        x4, y4 = self._point4()
        ax.scatter(x4, y4, s = 20, color = 'k')
        x5, y5 = self._point5(angle = angles[1])
        ax.plot([x4,x5], [y4,y5], color = 'tab:orange', linewidth = 3)
        x5b, y5b = self._point5b(angle = angles[1])
        ax.plot([x4,x5b], [y4,y5b], color = 'tab:orange', linewidth = 3)
        x6, y6 = self._point6(angle = angles[1])
        ax.plot([x5,x6], [y5,y6], color = 'k', linewidth = 0.75)
        x7, y7 = self._point7(angle = angles[1])
        ax.plot([x6,x7], [y6,y7], color = 'k', linewidth = 0.75)
        x8, y8 = self._point8(mirror = angles[0], grating = angles[1])
        ax.scatter(x8, y8, s = 30, color = 'tab:red')
        ax.plot([x3,x8], [y3,y8], color = 'tab:red', linewidth = 0.75, linestyle='--')

        x9, y9 = self._point9(angles[0], angles[1])
        ax.scatter(x5, y5, s = 30, color = 'tab:red')
        ax.scatter(x9, y9, s = 30, color = 'tab:red')
        ax.plot([x5,x9], [y5,y9], color = 'tab:red', linewidth = 0.75, linestyle='--')

        x10, y10 = self._point10(angles[0], angles[1])
        ax.plot([x2,x10], [y2,y10], color = 'tab:orange', linewidth = 3)
        ax.scatter(x7, y7, s = 30, color = 'tab:red')
        ax.scatter(x10, y10, s = 30, color = 'tab:red')
        ax.plot([x7,x10], [y7,y10], color = 'tab:red', linewidth = 0.75, linestyle='--')

        if fullmirror:
            x11, y11 = self._point11(mirror = angles[0])
            ax.plot([x10,x11], [y10,y11], color = 'tab:orange', linewidth = 3)

        ttl = 'E = {0:5.2f} eV,  cff = {1:5.2f},  mir = {2:5.2f}°,  gra = {3:5.2f}°,  mir.-gra. = {4:5.1f} mm,  c.s. = {5:5.1f} mm,  mir.-gird. = {6:5.1f}'.format(energy, cff, angles[0], angles[1], gap0, gap1, gap2)
        ax.set_title(ttl)
    
    def _mggap(self, mirror = 0, grating = 0):
        """
        This method calculates the smallest distance between the mirror and the
        grating as a function of their angles.
        """
        x5, y5 = self._point5(grating)
        x9, y9 = self._point9(mirror, grating)
        # distance between the mirror and grating
        dY = y5 - y9
        dX = x9 - x5
        gap = np.sqrt(dY**2 + dX**2)
        OK = True
        if y5 < y9: 
            OK = False
            gap = - abs(gap)
        return gap, OK
    
    def _mgirgap(self, mirror = 0, grating = 0):
        """
        This method calculates the smallest distance between the mirror and the
        grating girder as a function of their angles.
        """
        #x1, y1 = self._point1(mirror)
        #x7, y7 = self._point7(grating)
        ## Line coefficient
        #k = np.tan(np.deg2rad(mirror))
        #kN = np.tan(np.deg2rad(mirror+90))
        # Point on the mirror
        #x10 = 1/(k-kN) * ((y7-kN*x7) - (y1-k*x1))  # Update this to use of the ._point10() method.
        #y10 = y1 + k*x7 - k*x1 

        x7, y7 = self._point7(grating)
        x10, y10 = self._point10(mirror = mirror, grating = grating)
        # distance between the mirror and grating
        dY = y7 - y10
        dX = x10 - x7
        gap = np.sqrt(dY**2 + dX**2)
        OK = True
        if y7<=y10: 
            OK = False
            gap = - abs(gap)
        return gap, OK

    # old stuff, leave it
    def Safe(self, energy = None, cff = None, ld = None, silence = True):
        a, ok = self.EC2A(energy = energy, cff = cff, ld = ld)
        if not ok:
            if not silence:
                print('Angles out of range')
            return False
        csok, csgap = self._gap(theta = a[0], beta = a[1])
        if csok and not silence: print('Gap at ls: {0}mm'.format(round(csgap,1)))
        return csok and ok

    
    # ===================== minimum cff and energy 

    def MinCff(self, energy = 20, ld = None):
        if not type(ld) is type(None):
            self.monoParameters['LD'] = ld
        cffs = np.concatenate([ np.arange(1.05, 4, 0.01), np.arange(4, 10, 0.05) ])
        cffs = np.concatenate([ cffs, np.arange(10, 20, 0.1) ])
        cffs = np.concatenate( [ cffs, np.arange(20, 50, 0.25) ])
        cffs = np.concatenate( [ cffs, np.arange(50, 200, 0.5) ])
        for i, cff in enumerate(cffs):
            a, aok = self.EC2A(energy = energy, cff = cff)
            csok = self.Safe(energy = energy, cff = cff)
            if aok and csok:
                return round(cffs[i], 2)
        return np.nan
    
    def MinEnergy(self, cff = 2.25, ld = None):
        if not type(ld) is type(None):
            self.monoParameters['LD'] = ld
        energies = np.concatenate([ np.arange(5, 8, 0.01), np.arange(8, 12, 0.02) ])
        energies = np.concatenate([ energies, np.arange(12, 20, 0.05) ])
        energies = np.concatenate( [ energies, np.arange(20, 50, 0.1) ])
        energies = np.concatenate( [ energies, np.arange(50, 80, 0.25) ])
        energies = np.concatenate( [ energies, np.arange(80, 200, 0.5) ])
        energies = np.concatenate( [ energies, np.arange(200, 1000, 1) ])
        for i, energy in enumerate(energies):
            a, aok = self.EC2A(energy = energy, cff = cff)
            csok = self.Safe(energy = energy, cff = cff)
            if aok and csok:
                return round(energies[i], 2)
        return np.nan


    

    
# =========================================================


def drawMono(energy = 20, cff = 2.25, line_density = 800, ax = None,
             offset = 0, size = 5,
             figsize = (10,4), xlim = (None,None), ylim = (None, None),
             title = '_none_', xlabel = '_none_', ylabel = '_none_', fontsize = 8,
             auto_title = 1, zol = None,
             plot = True):
    """
    New version of drawMono(). Will hopefully be better designed. Work in progress.
    """
    result = {}
    # -------------------------------- Check args
    note_msg = Fore.BLUE + "(drawMono() was updated. The error might be casued by that. Use the old version drawMono_() instead.)" + Fore.BLUE
    err_msg1 = Fore.RED + "Invalid arguments.\n'energy' as a number > 0, 'cff' as a number > 1, 'line_density' as a number > 0." + Fore.BLACK
    try:
        energy = float(energy); cff = float(cff); line_density = float(line_density)
    except:
        print(err_msg1); return ax, result
    if not (energy > 0 and cff > 1 and line_density > 0):
        print(err_msg1); return ax, result
    if not type(figsize) is tuple:
        print(Fore.RED + "Argument 'figsize' must be a 2-tuple. Setting to default." + Fore.BLACK)
        print(note_msg)
        figsize = (10,4)
    if not len(figsize) == 2:
        print(Fore.RED + "Argument 'figsize' must be a 2-tuple. Setting to default." + Fore.BLACK)
        print(note_msg)
        figsize = (10,4)
    if not (type(title) is str and type(xlabel) is str and type(ylabel) is str):
        print(Fore.RED + "Arguments 'title', 'xlabel', and 'ylabel' must be strings." + Fore.BLACK)
        print(note_msg)
        return ax, result
    try: fontsize = float(fontsize)
    except:
        print(Fore.RED + "Argument 'fontsize' must be a number. Setting to default 8." + Fore.BLACK)
        print(note_msg)
        fontsize = 8
    #
    mono = BlochMono(line_density = line_density)
    #
    mpA = mono.monoParameters['mpA']; mpB = mono.monoParameters['mpB']; mpC = mono.monoParameters['mpC']
    mpD = mono.monoParameters['mpD']; mpE = mono.monoParameters['mpE']; mpF = mono.monoParameters['mpF']
    mpG = mono.monoParameters['mpG']
    #
    a, anglesok = mono.EC2A(energy, cff, line_density)
    theta = np.deg2rad(a[0])
    beta = np.deg2rad(a[1])
    alpha = 2 * theta - beta
    #
    if not type(zol) is type(None):
        try:
            zol = abs(float(zol))
        except:
            print(Fore.RED + "Pass argument zol as a number (mirror/grating angle)." + Fore.BLACK); return ax, result
        alpha = beta = theta = np.deg2rad(zol)
        energy = cff = np.NaN
    #

    # --------------------------------
    def intersect(L1,L2):
        # y1(x) = a1 + b1*x
        b1 = (L1[3]-L1[1])/(L1[2]-L1[0])
        a1 = L1[1]-b1*L1[0]
        # y2(x) = a2 + b2*x
        b2 = (L2[3]-L2[1])/(L2[2]-L2[0])
        a2 = L2[1]-b2*L2[0]
        # intersection y1(xi) = y2(x1)
        xi = (a2-a1)/(b1-b2)
        yi = a1 + b1*(a2-a1)/(b1-b2)
        a = np.arctan2( (b1-b2), (1+b1*b2))
        return xi, yi, a
    # --------------------------------
    def getangle(Line):
        a = np.arctan2((Line[3]-Line[1]),(Line[2]-Line[0]))
        return a
    # ----------------------------------
    PointG0 = [0., mpF]
    PointG1 = [ PointG0[0] - 0.5 * mpG * np.cos(beta),
                PointG0[1] - 0.5 * mpG * np.sin(beta) ]
    PointG2 = [ PointG0[0] + 0.5 * mpG * np.cos(beta),
                PointG0[1] + 0.5 * mpG * np.sin(beta) ]
    PointMR = [ PointG0[0] + mpB, 
                mpA ]
    PointMH = [ PointMR[0] + mpC * np.sin(theta), 
               PointMR[1] - mpC * np.cos(theta) ]
    PointM2 = [ PointMH[0] - mpD * np.cos(theta),
                PointMH[1] - mpD * np.sin(theta) ]
    PointM1 = [ PointMH[0] - (mpD + mpE) * np.cos(theta),
                PointMH[1] - (mpD + mpE) * np.sin(theta) ]
    Grating = [ PointG1[0], PointG1[1], PointG2[0], PointG2[1] ]
    Mirror  = [ PointM1[0], PointM1[1], PointM2[0], PointM2[1] ]
    # ----------------------------------
    def tracer(IN,MIRROR):
        x,y,a = intersect(MIRROR, IN)
        IN = [ IN[0], IN[1], x, y ]
        a = a + getangle(MIRROR)
        OUT = [ IN[2], IN[3], IN[2] + 500*np.cos(a), IN[3] + 500*np.sin(a) ]
        return IN,OUT
    # ---------------------------------- define the rays
    RAY1 = [-650, offset, -650 + 1000, offset ]
    RAY1u = [ RAY1[0], RAY1[1]+size/2., RAY1[2], RAY1[3]+size/2]
    RAY1d = [ RAY1[0], RAY1[1]-size/2., RAY1[2], RAY1[3]-size/2]
    LS = [ ['-','-','-'], ['-','-','-'], ['-','-','-'] ] # center, upper, lower
    # ---------------------------------- Center ray hitting M2
    RAY1, RAY2 = tracer(RAY1, Mirror)
    if not( RAY1[2] >= PointM1[0] and RAY1[2] <= PointM2[0]):
        print(Fore.RED + "The center of the beam does not hit M2." + Fore.RED); return ax
    if not( RAY2[2] >= PointG1[0] and RAY1[2] <= PointG2[0]):
        print('drawMono() - Center ray outside grating!')
        return ax, result
    # ---------------------------------- Center ray hitting the grating
    RAY2, RAY3 = tracer(RAY2, Grating)
    RAY3[2] = PointG0[0] + (mpG/2)*1.5
    RAY3[3] = RAY3[1]
    # ---------------------------------- Upper beam hitting M2
    RAY1u, RAY2u = tracer(RAY1u, Mirror)
    if not( RAY1u[2] >= PointM1[0] and RAY1u[2] <= PointM2[0]):
        LS[0][1] = ':'
        RAY1u = [RAY1u[0], PointM2[1], RAY1u[3], PointM2[1]]
        RAY1u, RAY2u = tracer(RAY1u, Mirror)
    # ---------------------------------- Upper beam hitting the grating
    RAY2u, RAY3u = tracer(RAY2u, Grating)
    if not( RAY2u[2] >= PointG1[0] and RAY2u[2] <= PointG2[0]):
        LS[1][1] = ':'
        a = getangle(RAY2u)
        RAY2u = [PointG2[0], PointG2[1],  PointG2[0]-2000*np.cos(a), PointG2[1] -2000*np.sin(a)]
        RAY2u, tmp = tracer(RAY2u, Mirror)
        RAY2u = [ RAY2u[2], RAY2u[3], RAY2u[0], RAY2u[1]]
    RAY2u, RAY3u = tracer(RAY2u, Grating)
    RAY3u[2] = PointG0[0] + (mpG/2)*1.5
    RAY3u[3] = RAY3u[1]
    # ---------------------------------- Lower beam hitting the M2
    RAY1d, RAY2d = tracer(RAY1d, Mirror)
    if not( RAY1d[2] >= PointM1[0] and RAY1d[2] <= PointM2[0]):
        LS[0][2] = ':'
        RAY1d = [RAY1d[0], PointM1[1], RAY1d[3], PointM1[1]]
        RAY1d, RAY2d = tracer(RAY1d, Mirror)
    # ---------------------------------- Lower beam hitting the grating
    RAY2d, RAY3d = tracer(RAY2d, Grating)
    if not( RAY2d[2] >= PointG1[0] and RAY2d[2] <= PointG2[0]):
        LS[1][2] = ':'
        a = getangle(RAY2d)
        RAY2d = [PointG1[0], PointG1[1],  PointG1[0]-2000*np.cos(a), PointG1[1] -2000*np.sin(a)]
        RAY2d, tmp = tracer(RAY2d, Mirror)
        RAY2d = [ RAY2d[2], RAY2d[3], RAY2d[0], RAY2d[1]]
    RAY2d, RAY3d = tracer(RAY2d, Grating)
    RAY3d[2] = PointG0[0] + (mpG/2)*1.5
    RAY3d[3] = RAY3d[1]
    # ----------------------------------
    if plot:
        fig = None
        if not hasattr(ax, 'plot'):
            fig, ax = plt.subplots(figsize = figsize)
        #
        ax.plot([PointMR[0],PointMH[0]], [PointMR[1],PointMH[1]], marker = '', color='silver')
        ax.plot([PointMH[0],Mirror[2]], [PointMH[1],Mirror[3]], marker = '', color='silver')
        
        ax.plot([Grating[0],Grating[2]], [Grating[1],Grating[3]], marker = '', color='goldenrod', linewidth=4)
        ax.plot([Mirror[0],Mirror[2]], [Mirror[1],Mirror[3]], marker = '', color='goldenrod', linewidth=4)
        
        ax.plot([RAY1d[0], RAY1d[2]],  [RAY1d[1], RAY1d[3]],  marker = '', color='b', linestyle=LS[0][2])
        ax.plot([RAY2d[0], RAY2d[2]],  [RAY2d[1], RAY2d[3]],  marker = '', color='b', linestyle=LS[1][2])
        ax.plot([RAY3d[0], RAY3d[2]],  [RAY3d[1], RAY3d[3]],  marker = '', color='b', linestyle=LS[2][2])

        ax.plot([RAY1u[0], RAY1u[2]],  [RAY1u[1], RAY1u[3]],  marker = '', color='r', linestyle=LS[0][1])
        ax.plot([RAY2u[0], RAY2u[2]],  [RAY2u[1], RAY2u[3]],  marker = '', color='r', linestyle=LS[1][1])
        ax.plot([RAY3u[0], RAY3u[2]],  [RAY3u[1], RAY3u[3]],  marker = '', color='r', linestyle=LS[2][1])

        ax.plot([RAY1[0], RAY1[2]],  [RAY1[1], RAY1[3]],  marker = '', color='g', linestyle=LS[0][0])
        ax.plot([RAY2[0], RAY2[2]],  [RAY2[1], RAY2[3]],  marker = '', color='g', linestyle=LS[1][0])
        ax.plot([RAY3[0], RAY3[2]],  [RAY3[1], RAY3[3]],  marker = '', color='g', linestyle=LS[2][0])
        #
        plt.axis("scaled")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        #
        if title == '_none_':
            if auto_title == 0:
                ttl = f"E = {energy} eV, cff = {cff}, G = {line_density} l/mm, $\u03B1$ = {90-alpha*180/np.pi:.3f}°, $\u03B2$ = {beta*180/np.pi-90:.3f}°, $\Theta$ = {theta*180/np.pi:.3f}°"
            elif auto_title == 1:
                ttl = f"$\u03B1$ = {90-alpha*180/np.pi:.3f}°, $\u03B2$ = {beta*180/np.pi-90:.3f}°, $\Theta$ = {theta*180/np.pi:.3f}°"
            elif auto_title == 2:
                ttl = f"mirror = {theta*180/np.pi:.3f}°, grating = {beta*180/pi:.3f}°"
            else:
                ttl = f"E = {energy} eV, cff = {cff}, G = {line_density} l/mm"
            ax.set_title(ttl, fontsize = 9)
        else:
            ax.set_title(title, fontsize = fontsize)
        if not xlabel == '_none_': ax.set_xlabel(xlabel, fontsize = fontsize)
        if not ylabel == '_none_': ax.set_ylabel(ylabel, fontsize = fontsize)
        #
        if not type(fig) is type(None): fig.tight_layout()
    
    #
    result = {'energy': energy,
              'cff': cff,
              'line_density': line_density,
              'alpha': 90-alpha*180/np.pi,
              'beta': beta*180/pi-90,
              'theta': theta*180/np.pi,
              'mirror_angle': theta*180/np.pi,
              'grating_angle': beta*180/pi,
              'size': size,
              'size_diffracted': RAY3u[3] - RAY3d[3]}
    m2_center = -np.sqrt((Mirror[3]-RAY1[3])**2 + (Mirror[2]-RAY1[2])**2)
    result.update({'mirror_center': m2_center})
    m2_footprint = np.sqrt((RAY1u[3] - RAY1d[3])**2 + (RAY1u[2] - RAY1d[2])**2)
    result.update({'mirror_footprint': m2_footprint})
    g_footprint = np.sqrt((RAY2u[3] - RAY2d[3])**2 + (RAY2u[2] - RAY2d[2])**2)
    result.update({'grating_footprint': g_footprint})
    #
    if plot: return ax, result
    else: return result


    

# =============================================================================================    blaze 


def blaze(energy = 100, cff = 2.25, blaze_angle = 2.17, line_density = 800, shup = False, plot = True, title = 1):
    """
    Draw the diffraction and blaze condition. Returns an pyplot axis and a dict with values, 
    or just a dict if argument plot is False.

    Notes:
    reflected_beam: the angle which the on-blaze light comes out above (or below) the diffracted light
    """
    try:
        energy, cff, blaze_angle, line_density = abs(float(energy)), float(cff), abs(float(blaze_angle)), abs(float(line_density))
    except:
        print(Fore.RED + "Arguments alpha, beta, blaze_angle, and line_density must be numbers (deg. and l/mm)." + Fore.BLACK); return
    if cff <= 1:
        print(Fore.RED + "Argument cff must be larger than 1. Setting default 2.25." + Fore.BLACK)
    #
    if not(type(title) is int or type(title) is str):
        print(Fore.RED + "Argument title must be a string or an integer (1-3: auto-generated titles). Setting default." + Fore.BLACK)
        title = 1
    if type(title) is int:
        if title <= 0 or title >= 4:
            print(Fore.RED + "If argument title is passed as an integer it must be 1, 2, or 3. Setting default." + Fore.BLACK)
            title = 1
    #
    #mono_config = cpgm_e2a(energy = energy, cff = cff, order = 1, line_density = line_density, shup = True, return_dict = True)
    #alpha = mono_config['alpha']
    #beta = mono_config['beta']
    mono = BlochMono(line_density = line_density)
    a, anglesok = mono.EC2A(energy, cff, line_density)
    mirror = a[0]
    theta = mirror
    grating = a[1]
    beta = grating - 90
    alpha = 180 - 2 * theta + beta
    alpha_r, beta_r = np.deg2rad(alpha), np.deg2rad(beta)
    #
    alpha_r, beta_r, blaze_angle_r = np.deg2rad(alpha), np.deg2rad(beta), np.deg2rad(blaze_angle)
    d = 1/line_density # line spacing
    h = d * np.tan(blaze_angle_r) # groove height
    beta_reflected = 2 * blaze_angle - alpha
    beta_reflected_r = np.deg2rad(beta_reflected)
    #
    q = {'energy': energy, 'cff': cff, 'line_density': line_density, 'blaze_angle': blaze_angle, 'alpha': alpha, 'beta': beta}
    #q.update({'beta_on_blaze': beta_reflected, 'mirror': mono_config['mirror'], 'grating': mono_config['grating']})
    q.update({'beta_on_blaze': beta_reflected, 'mirror': mirror, 'grating': grating})
    q.update({'reflected_beam': -(beta_reflected - beta)})
    #
    if plot:
        X = np.linspace(-2.5*d, 1.5*d, 5)
        fig, ax = plt.subplots(figsize = (10,6))
        Nh = 8; Nr = 100
        for x in X:
            ax.plot([x, x+d],[-h/2,h/2], color = "k")
            ax.plot([x+d, x+d],[h/2,-h/2], color = "k")
        ax.plot([0,0], [0,(Nh-1)*h], linestyle = "--", linewidth = 0.75, color = "k")
        # incidence
        ax.plot([0, -np.sin(alpha_r) * Nr*h], [0, np.cos(alpha_r) * Nr*h], color = "tab:red", linestyle = "-", label = "incidence")
        # diffracted
        ax.plot([0, np.sin(-beta_r) * Nr*h], [0, np.cos(-beta_r) * Nr*h], color = "tab:orange", linestyle = "-", label = "diffracted" )
        # reflected
        ax.plot([0, np.sin(-beta_reflected_r) * Nr*h], [0, np.cos(-beta_reflected_r) * Nr*h], color = "tab:orange", linestyle = "--", label = "reflected" )
        plt.axis("scaled")
        ax.set_ylim(-h,Nh*h)
        ax.set_xlim(-1.75*d, 1.75*d)
        ax.legend(fontsize = 8)
        #
        if type(title) is str:
            ax.set_title(title, fontsize = 9)
        else:
            if title == 1:
                ax.set_title(f"\u03B1 = {alpha:6.3f}°, \u03B2 = {beta:7.3f}°, \u03B2(blaze) = {beta_reflected:7.3f}°", fontsize = 9)
            elif title == 2:
                ax.set_title(f"energy = {energy}eV, cff = {cff}, G = {line_density:.0f}l/mm, \u03B1 = {alpha:6.3f}°, \u03B2 = {beta:7.3f}°, \u03B2(blaze) = {beta_reflected:7.3f}°", fontsize = 9)
        #
        ax.invert_yaxis()
        ax.set_yticks([])
        fig.tight_layout()
    #
    if not shup:
        print(f"energy:         {energy} eV")
        print(f"cff:            {cff}")
        print(f"line density:   {line_density} l/mm")
        print(f"blaze angle:    {blaze_angle:9.5f}°")
        print(f"alpha:          {alpha:9.5f}°")
        print(f"beta:           {beta:9.5f}°")
        print(f"beta on-blaze:  {beta_reflected:9.5f}°")
        print(f"mirror:         {mirror:9.5f}°")
        print(f"grating:        {grating:9.5f}°")
        print(f"reflected beam: {q['reflected_beam']:9.5f}° (from diffracted)")
    #
    if plot: return ax, q
    else: return q

    
