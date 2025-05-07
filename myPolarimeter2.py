"""

Versions
    24.07.05    Added fontsizes to plotScan()
    24.03.10    Bugfixes
    24.03.03    Introduced logging (to learn about logging)
    24.02.28    minor update of plotScan() and add2ax()
    24.02.05    touch-ups
    24.02.01    rewritten, simplified
    23.10.11    pass path = "mats" to load default data on my own laptop
    23.08.23    added stuff (to use for a conference)
    23.07.09    the data can now be either from 4 mirrors or from a retarder and an analyzer.
    23.07.08    fixed minor bugs
    23.07.03    finished compiling and rewriting

"""

__version__ = "24.07.05"
__author__  = "Mats Leandersson"


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from colorama import Fore
from copy import deepcopy
from scipy.optimize import curve_fit

import os   # used for the ideal polarimeter

import logging
logging.basicConfig(level = logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers = [logging.StreamHandler(), logging.FileHandler(filename = "logfile.log", mode = "a")])


# --------
try:
    import myStokes as mys
except:
    logging.critical("myPolarimeter.py: Warning: The 'myStokes' module could not be imported. Most or all of the code in this module is depending on that module. Bummer.")

try:
    from scipy.optimize import curve_fit
except:
    logging.critical("myPolarimeter.py: Warning: Could not import 'curve_fit' from 'scipy.optimize'. Some methods will not work. Bummer.")
# --------


my_mirror_mirror_files = ["M1_optical.dat", "M2_optical.dat", "M1_optical.dat", "M4_optical.dat"]
my_path = "/Users/matlea/Modules/optical_parameters_polarimeter/"




class Polarimeter():
    """
    Methods:    
        analyzerScan()      1d, scanning 1 axis
        retarderScan()      1d, scanning 1 axis
        polarimeterScan()   2d, scanning 2 axes
        specialCaseScan()   1d, scanning 2 axes
        specialCaseScans()  1d, scanning 2 axes

        (The above methods returns dicts with all scan input and output. Can be plotted with plotScan(dict, ax, **kwargs).)
    
    Properties:
        energy              get, set    float
        energies            get         array
        stokes              get, set    array (get size 4, set size 3 or 4)
        retarder_matrix0    get         array (4x4), the non-rotated retarder matrix
        analyzer_matrix0    get         array (4x4), the non-rotated analyzer matrix
        retarder_matrix     get         array (4x4), the rotated retarder matrix
        analyzer_matrix     get         array (4x4), the rotated analyzer matrix
        eta1                get, set    float, retarder angle
        eta2                get, set    float, analyzer angle
        eta1_noise          get, set    float
        eta2_noise          get, set    float
        intensity           get         float, polarimeter intensity
        intensity_noise     get, set    float
        retarder_present    get, set    bool


    """
    def __init__(self, mirror_files = [], shup = False, my_setup = False):
        """"""
        self._ok = False
        self._shup = shup
        #
        self._arg_descr = {}
        self._arg_descr.update({"analyzerScan": {"args": ["eta_start", "eta_stop", "eta_steps", "eta1"], 
                                                "descr": ["start angle, deg.", "stop angle, deg.", "number of steps", "retarder angle, deg. (if not set value)"]}})
        self._arg_descr.update({"retarderScan": {"args": ["eta_start", "eta_stop", "eta_steps", "eta2"], 
                                                "descr": ["start angle, deg.", "stop angle, deg.", "number of steps", "analyzer angle, deg. (if not set value)"]}})
        self._arg_descr.update({"polarimeterScan": {"args": ["eta1_start", "eta1_stop", "eta1_steps", "eta2_start", "eta2_stop", "eta2_steps"],
                                                    "descr": ["retarder start angle, deg.", "retarder stop angle, deg.", "number of steps, retarder", 
                                                              "analyzer start angle, deg.", "analyzer stop angle, deg.", "analyzer of steps, retarder"]}})
        self._arg_descr.update({"specialCaseScan": {"args": ["special_case", "delta", "eta_start", "eta_stop", "eta_steps"], 
                                                    "descr": ["special case number, int., 1-4", "if other case, angle difference, deg.", 
                                                              "retarder start angle, deg.", "retarder stop angle, deg.", "number of steps, retarder"]}})
        self._arg_descr.update({"specialCaseScans": {"args": ["eta_start", "eta_stop", "eta_steps"], 
                                                    "descr": ["retarder start angle, deg.", "retarder stop angle, deg.", "number of steps"]}})
        
        #
        if my_setup:
            global my_mirror_mirror_files, my_path
            mirror_files = []
            for mf in my_mirror_mirror_files: mirror_files.append(my_path + mf)
        #
        if not type(mirror_files) is list:
            logging.error("Polarimeter(): Argument 'mirror_files' must be a list of size 4.")
            return
        if not len(mirror_files) == 4:
            if len(mirror_files) == 0:
                logging.warning("Polarimeter(): No optical data file list supplied. Generating fake files.")
                mirror_files = self._fakeData()
            else:
                logging.error("Polarimeter(): Argument 'mirror_files' must be a list of size 4.")
                return
        self._mirror_files = mirror_files
        self._loadMirrorData()
        if not self._ok:
            logging.critical("Polarimeter(): Missing data. The polarimeter will not work.")
            return
        self._energies = np.array(self._mirror_data[0][0])
        self._energy = None
        self._stokes = np.array([1,0,0,1])
        self._analyzer_matrix0 = None
        self._retarder_matrix0 = None
        self._analyzer_matrix = None
        self._retarder_matrix = None
        self._eta1 = 0.
        self._eta2 = 0.
        self._eta1_noise = 0.
        self._eta2_noise = 0.
        self._intensity_noise = 0.
        self._retarder_present = True
        #
        self.energy = self._energies[0]
        #
    
    # ------------------
    
    def _printArgs(self, method = ""):
        d = self._arg_descr.get(method, {})
        if d == {}: return
        args, descr = d.get("args"), d.get("descr")
        print(Fore.BLUE)
        print(f"{method}():")
        for i, arg in enumerate(args):
            print(f"  {arg:<28}{descr[i]}")
        print(Fore.RESET)

        
    # ------------------
    
    def _loadMirrorData(self):
        ok = True
        Data = []
        for i, file in enumerate(self._mirror_files):
            try:
                data = np.loadtxt(file).transpose()
                msg = f"Polarimeter(): Loaded data from file {i+1} ({file})."
                logging.info(msg)
            except:
                data = np.array([])
                logging.error(f"Polarimeter(): Could not load data from file {i+1} ({file}).")
                ok = False
            Data.append(data)
        self._mirror_data = np.copy(Data)
        self._ok = ok
    
    def _calculateMatrices(self):
        indx = abs(self._energy - self._energies).argmin()
        mm1 = mys.matrixFromOptical(rs = self._mirror_data[0][1][indx], rp = self._mirror_data[0][2][indx], delta = self._mirror_data[0][3][indx])
        mm2 = mys.matrixFromOptical(rs = self._mirror_data[1][1][indx], rp = self._mirror_data[1][2][indx], delta = self._mirror_data[1][3][indx])
        mm3 = mys.matrixFromOptical(rs = self._mirror_data[2][1][indx], rp = self._mirror_data[2][2][indx], delta = self._mirror_data[2][3][indx])
        mm4 = mys.matrixFromOptical(rs = self._mirror_data[3][1][indx], rp = self._mirror_data[3][2][indx], delta = self._mirror_data[3][3][indx])
        mm2 = mys.rotateMatrix(mm2, 180)
        self._retarder_matrix0 = mm3.dot(mm2.dot(mm1))
        self._analyzer_matrix0 = np.copy(mm4)
        self._rotateRetarderMatrix()
        self._rotateAnalyzerMatrix()
    
    def _rotateRetarderMatrix(self):
        eta1 = self.eta1
        if self._eta1_noise > 0.: eta1 = eta1 + (0.5 - np.random.uniform(0, 1)) * self._eta1_noise
        self._retarder_matrix = mys.rotateMatrix(self._retarder_matrix0, eta1, True)
    
    def _rotateAnalyzerMatrix(self):
        eta2 = self.eta2
        if self._eta2_noise > 0.: eta2 = eta2 + (0.5 - np.random.uniform(0, 1)) * self._eta2_noise
        self._analyzer_matrix = mys.rotateMatrix(self._analyzer_matrix0, eta2, True)
    
    # ---------

    @property
    def shup(self):
        return self._shup
    @shup.setter
    def shup(self, value):
        if not type(value) is bool:
            logging.warning("Polarimeter(): shup must be a bool.")
            return
        self._shup = value

    @property
    def energy(self):
        return self._energy
    @energy.setter 
    def energy(self, value):
        try: value = abs(float(value))
        except:
            logging.warning("Polarimeter(): energy must be a positive float."); return
        indx = abs(value - self._energies).argmin()
        e = self._energies[indx]
        if value < self._energies.min() or value > self._energies.max():
            logging.warning(f"Polarimeter(): the energy range of the polarimeter is {self._energies.min()} to {self._energies.max()} eV.")
        self._energy = e
        self._calculateMatrices()
        if not self._energy == value and not self.shup:
            logging.info(f"Polarimeter(): set energy to {self._energy} eV (the stored value closest to {value}).")
    
    @property
    def energies(self):
        return self._energies
    @energies.setter
    def energies(self, value):
        logging.warning("Polarimeter(): energies can not be set (loaded from data files).")

    @property
    def stokes(self):
        return self._stokes
    @stokes.setter
    def stokes(self, value):
        if not (type(value) is list or type(value) is np.ndarray):
            logging.error("Polarimeter(): stokes must be a list or array of length 3 (or 4)."); return
        if not len(value) in [3, 4]:
            logging.error("Polarimeter(): stokes must be a list or array of length 3 (or 4)."); return
        if len(value) == 3:
            norm = np.sqrt(value[0]**2 + value[1]**2 + value[2]**2)
            value = np.array([1., value[0]/norm, value[1]/norm, value[2]/norm])
        else:
            norm = np.sqrt(value[1]**2 + value[2]**2 + value[3]**2)
            value = np.array([1., value[1]/norm, value[2]/norm, value[3]/norm])
        self._stokes = value
        self._calculateMatrices()
    
    @property
    def analyzer_matrix0(self):
        return self._analyzer_matrix0
    @analyzer_matrix0.setter
    def analyzer_matrix0(self, *args):
        logging.warning("Polarimeter(): analyzer_matrix0 is calculated and can not be set manually.")
    
    @property
    def retarder_matrix0(self):
        return self._retarder_matrix0
    @retarder_matrix0.setter
    def retarder_matrix0(self, *args):
        logging.warning("Polarimeter(): retarder_matrix0 is calculated and can not be set manually.")
    
    @property
    def analyzer_matrix(self):
        return self._analyzer_matrix
    @analyzer_matrix.setter
    def analyzer_matrix(self, *args):
        logging.warning("Polarimeter(): analyzer_matrix is calculated and can not be set manually.")
    
    @property
    def retarder_matrix(self):
        return self._retarder_matrix
    @retarder_matrix.setter
    def retarder_matrix(self, *args):
        logging.warning("Polarimeter(): retarder_matrix is calculated and can not be set manually.")
    
    @property
    def eta1_noise(self):
        return self._eta1_noise
    @eta1_noise.setter
    def eta1_noise(self, value):
        try: value = abs(float(value))
        except: logging.error("Polarimeter(): eta1_noise must be a float."); return
        self._eta1_noise = value
        self._rotateRetarderMatrix()
    
    @property
    def eta2_noise(self):
        return self._eta2_noise
    @eta2_noise.setter
    def eta2_noise(self, value):
        try: value = abs(float(value))
        except: logging.error("Polarimeter(): eta2_noise must be a float."); return
        self._eta2_noise = value
        self._rotateAnalyzerMatrix()

    @property
    def eta1(self):
        return self._eta1
    @eta1.setter
    def eta1(self, value):
        try: value = float(value)
        except: logging.error("Polarimeter(): eta1 must be a float (degrees)."); return
        self._eta1 = value
        self._rotateRetarderMatrix()
    
    @property
    def eta2(self):
        return self._eta2
    @eta2.setter
    def eta2(self, value):
        try: value = float(value)
        except: logging.error("Polarimeter(): eta2 must be a float (degrees)."); return
        self._eta2 = value
        self._rotateAnalyzerMatrix()
    
    @property
    def intensity_noise(self):
        return self._intensity_noise
    @intensity_noise.setter
    def intensity_noise(self, value):
        try: value = abs(float(value))
        except: logging.error("Polarimeter(): intensity_noise must be a float."); return
        self._intensity_noise = value
    
    @property
    def intensity(self):
        if self.retarder_present:
            I = self._analyzer_matrix.dot(self._retarder_matrix).dot(self._stokes)[0]
        else:
            I = self._analyzer_matrix.dot(self._stokes)[0]
        if self._intensity_noise > 0:
            I = I * (1 + np.random.uniform(0, self._intensity_noise))
        return I
    @intensity.setter
    def intensity(self, *args):
        logging.warning("Polarimeter(): intensity is calculated and can not be set manually.")
    
    # ------------------------

    @property 
    def retarder_present(self):
        return self._retarder_present
    @retarder_present.setter 
    def retarder_present(self, value):
        if not type(value) is bool:
            logging.error("Polarimeter(): retarder_present must be a float."); return
        self._retarder_present = value

    # ------------------------
    
    @property
    def rs_m1(self):
        return self._mirror_data[0][1]
    @property
    def rp_m1(self):
        return self._mirror_data[0][2]
    @property
    def d_m1(self):
        return self._mirror_data[0][2]
    @property
    def rs_m2(self):
        return self._mirror_data[1][1]
    @property
    def rp_m2(self):
        return self._mirror_data[1][2]
    @property
    def d_m2(self):
        return self._mirror_data[1][2]
    @property
    def rs_m3(self):
        return self._mirror_data[2][1]
    @property
    def rp_m3(self):
        return self._mirror_data[2][2]
    @property
    def d_m3(self):
        return self._mirror_data[2][2]
    @property
    def rs_m4(self):
        return self._mirror_data[3][1]
    @property
    def rp_m4(self):
        return self._mirror_data[3][2]
    @property
    def d_m4(self):
        return self._mirror_data[3][2]

    # ------------------------

    def analyzerScan(self, eta_start = 0., eta_stop = 360., eta_steps = 361, eta1 = None):
        """
        Returns an intensity curve for an analyzer scan.
        Arguments eta_start, eta_stop, and eta_steps defines the analyzer axis. Argument eta1 sets the retarder angle (if it should be changed from current).
        """
        result = {"eta": np.array([]), "intensity": np.array([]), "energy": self.energy, "stokes": self.stokes, "eta1": self.eta1}
        try:
            eta_start, eta_stop, eta_steps = float(eta_start), float(eta_stop), int(eta_steps)
        except:
            logging.error("Polarimeter.analyzerScan(): args eta_start and eta_stop must be floats or ints and arg eta_steps must be an int.")
            self._printArgs("analyzerScan")
            return result
        if not type(eta1) is type(None):
            try:
                eta1 = float(eta1)
                self.eta1 = eta1
                result.update({"eta1": eta1})
            except:
                logging.error(f"Polarimeter.analyzerScan(): arg eta1 must be a number (degrees). using already set {self.eta1}.")
                self._printArgs("analyzerScan")
        ETA = self.eta2
        eta = np.linspace(eta_start, eta_stop, eta_steps)
        intensity = np.zeros(len(eta)) * np.NaN
        for i, e in enumerate(eta):
            self.eta2 = e
            intensity[i] = self.intensity
        result.update({"eta": eta, "intensity": intensity, "scan_type": "analyzer_scan"})
        self.eta2 = ETA
        return result
    
    def retarderScan(self, eta_start = 0., eta_stop = 360., eta_steps = 361, eta2 = None):
        """
        Returns an intensity curve for a retarder scan.
        Arguments eta1_start, eta1_stop, and eta1_steps defines the retarder axis. Argument eta2 sets the analyzer angle (if it should be changed from current).
        """
        result = {"eta": np.array([]), "intensity": np.array([]), "energy": self.energy, "stokes": self.stokes, "eta2": self.eta2}
        try:
            eta_start, eta_stop, eta_steps = float(eta_start), float(eta_stop), int(eta_steps)
        except:
            logging.error("Polarimeter.retarderScan(): args eta_start and eta_stop must be floats or ints and arg eta_steps must be an int.")
            self._printArgs("retarderScan")
            return result
        if not type(eta2) is type(None):
            try:
                eta2 = float(eta2)
                self.eta2 = eta2
                result.update({"eta2": eta2})
            except:
                logging.error(f"Polarimeter.retarderScan(): arg eta2 must be a number (degrees). using already set {self.eta2}.")
                self._printArgs("retarderScan")
        ETA = self.eta1
        eta = np.linspace(eta_start, eta_stop, eta_steps)
        intensity = np.zeros(len(eta)) * np.NaN
        for i, e in enumerate(eta):
            self.eta1 = e
            intensity[i] = self.intensity
        result.update({"eta": eta, "intensity": intensity, "scan_type": "retarder_scan"})
        self.eta1 = ETA
        return result
    
    def polarimeterScan(self, eta1_start = 0, eta1_stop = 360, eta1_steps = 181, eta2_start = 0, eta2_stop = 360, eta2_steps = 181):
        """
        Returns the intensity map for a polarimeter scan (retarder, analyzer axes).
        Arguments:  eta1_start, eta1_stop, eta1_steps (retarder axis)
                    eta2_start, eta2_stop, eta2_steps (analyzer axis)
        """
        result = {"eta1": np.array([]), "eta2": np.array([]), "intensity": np.array([]), "energy": self.energy, "stokes": self.stokes}
        try:
            eta1_start, eta1_stop, eta1_steps =  float(eta1_start), float(eta1_stop), int(eta1_steps)
            eta2_start, eta2_stop, eta2_steps =  float(eta2_start), float(eta2_stop), int(eta2_steps)
        except:
            logging.error("Polarimeter.polarimeterScan(): args eta1_start, eta1_stop, eta2_start, and eta2_stop")
            print("must be floats and eta1_steps and eta2_steps must be ints.")
            self._printArgs("polarimeterScan")
            return result
        ETA1, ETA2 = self.eta1, self.eta2
        eta1 = np.linspace(eta1_start, eta1_stop, eta1_steps)
        eta2 = np.linspace(eta2_start, eta2_stop, eta2_steps)
        result.update({"eta1": eta1, "eta2": eta2})
        imap = np.zeros([len(eta1), len(eta2)])
        for ie1, e1 in enumerate(eta1):
            self.eta1 = e1
            for ie2, e2 in enumerate(eta2):
                self.eta2 = e2
                imap[ie1][ie2] = self.intensity
        result.update({"intensity": imap, "scan_type": "polarimeter_scan"})
        self.eta1, self.eta2 = ETA1, ETA2
        return result
    
    def specialCaseScan(self, special_case = None, delta = None, eta_start = 0, eta_stop = 360, eta_steps = 181):
        """
        Returns the intensity for a special case if arg special_case is passed (int., 1-4) or an arbitrary case if arg delta is passed (number, deg.).
        Other args: eta_start (number), eta_stop (number), eta_steps (int.), defines the retarder angle axis.
        """
        result = {"eta": np.array([]), "intensity": np.array([]), "special_case": None, "delta": None, "energy": self.energy, "stokes": self.stokes}
        if not type(special_case) is type(None) and not type(delta) is type(None):
            logging.error(f"Polarimeter.specialCaseScan(): pass EITHER arg special_case as an integer (1-4) or arg delta as an angle (deg.)"); return result
        #
        if not type(special_case) is type(None):
            if not type(special_case) is int:
                logging.error(f"Polarimeter.specialCaseScan(): arg special_case must be an int (1-4).")
                self._printArgs("specialCaseScan"); return result
            if not special_case in [1,2,3,4]:
                logging.error(f"Polarimeter.specialCaseScan(): arg special_case must be an int (1-4).")
                self._printArgs("specialCaseScan"); return result
            if special_case == 1: delta = 0.
            elif special_case == 2: delta = 90.
            elif special_case == 3: delta = 45.
            else: special_case = -45.
            result.update({"special_case": special_case, "delta": delta})
        else:
            try: delta = float(delta)
            except:
                logging.error(f"Polarimeter.specialCaseScan(): arg delta must be a number (deg.).")
                self._printArgs("specialCaseScan"); return result
            result.update({"delta": delta})
        #
        try:
            eta_start, eta_stop, eta_steps = float(eta_start), float(eta_stop), int(eta_steps)
        except:
            logging.error("Polarimeter.specialCaseScan(): args eta_start and eta_stop must be floats or ints and arg eta_steps must be an int.")
            self._printArgs("specialCaseScan"); return result
        #
        eta = np.linspace(eta_start, eta_stop, eta_steps)
        intensity = np.zeros(len(eta)) * np.NaN
        ETA1, ETA2 = self.eta1, self.eta2
        for ie, e in enumerate(eta):
            self.eta1, self.eta2 = e, e + delta
            intensity[ie] = self.intensity
        result.update({"eta": eta, "intensity": intensity})
        result.update({"scan_type": "special_case_scan"})
        return result
            
    def specialCaseScans(self, eta_start = 0, eta_stop = 360, eta_steps = 181):
        """
        Returns intensity curves for the 4 special cases. Arguments: eta_start, eta_stop, eta_steps (retarder).
        """
        result = {"eta": np.array([]), "intensity1": np.array([]), "intensity2": np.array([]), "intensity3": np.array([]), "intensity4": np.array([]), 
                  "deltas": np.array([0, 90, 45, -45]), "energy": self.energy, "stokes": self.stokes}
        try:
            eta_start, eta_stop, eta_steps = float(eta_start), float(eta_stop), int(eta_steps)
        except:
            logging.error("Polarimeter.specialCaseScans(): args eta_start and eta_stop must be floats or ints and arg eta_steps must be an int.")
            self._printArgs("specialCaseScans"); return result
        #
        eta = np.linspace(eta_start, eta_stop, eta_steps)
        intensity1, intensity2, intensity3, intensity4 = np.zeros(len(eta)) * np.NaN, np.zeros(len(eta)) * np.NaN, np.zeros(len(eta)) * np.NaN, np.zeros(len(eta)) * np.NaN
        ETA1, ETA2 = self.eta1, self.eta2
        for ie, e in enumerate(eta):
            self.eta1 = e
            self.eta2 = e
            intensity1[ie] = self.intensity
            self.eta2 = e + 90
            intensity2[ie] = self.intensity
            self.eta2 = e + 45
            intensity3[ie] = self.intensity
            self.eta2 = e - 45
            intensity4[ie] = self.intensity
        result.update({"eta": eta, "intensity1": intensity1, "intensity2": intensity2, "intensity3": intensity3, "intensity4": intensity4})
        result.update({"scan_type": "special_case_scans"})
        return result
    
    # ----

    def _fakeData(self):
        """Generates fake (but reasonable) data files in case real ones are missing."""
        data1, data2, data3, data4 = [], [], [], []
        for i, e in enumerate(np.linspace(10, 100, 91)):
            rs = 1 - i/400
            rp = np.sqrt(1 - rs**2)
            data1.append([e, rs, rp, -50 + i/4])
            data2.append([e, rs, rp, -55 + i/4])
            data3.append([e, rs, rp, -50 + i/4])
            rs = 1 - i/200
            rp = np.sqrt(1 - rs**2)
            data4.append([e, rs, rp, -20 + i/10])
        fns, files, data = [], [], [data1, data2, data3, data4]
        for i in [1,2,3,4]:
            fn = f"fake_optical_data_mirror{i}.dat"
            fns.append(fn)
            files.append(open(fn, "w"))
        for i in range(len(data1)):
            for j, file in enumerate(files):
                file.write(f"{data[j][i][0]:5.1f}\t{data[j][i][1]:5.3f}\t{data[j][i][2]:5.3f}\t{data[j][i][3]:5.3f}\n")
        for file in files: file.close()
        logging.info("Polarimeter._fakeData(): Generated fake (but reasonable) optical data files for the four mirrors.")
        return fns


        




def plotScan(scan = {}, ax = None, **kwargs):
    """
    Plot scan results from the polarimeter class.

    scan    a dict with scan results
    ax      pyplot axis, optional

    kwargs: aspect (str), ticks (number), ticks_minor (number), title (str), cmap (str), vmin (number), vmax (number), cbar (bool),
            xlabel (str), ylabel (str), color (str), linestyle (str), linewidth (number), legend (bool)
    
    """
    if not type(scan) is dict:
        logging.error("plotScan(): args scan must be a dict from a polarimeter scan.")
        return ax
    #
    eta1, eta2, intensity = scan.get("eta1", None), scan.get("eta2", None), scan.get("intensity", None)
    energy, stokes = scan.get("energy", None), scan.get("stokes", None)
    scan_type = scan.get("scan_type", "")

    fig = None
    aspect = kwargs.get("aspect", "equal")
    if not aspect in ["auto", "equal"]: aspect = "equal"
    ticks = kwargs.get("ticks", 90)
    ticks_minor = kwargs.get("ticks_minor", 0)
    title = kwargs.get("title", "")
    cmap = kwargs.get("cmap", "viridis")
    vmin = kwargs.get("vmin", None)
    vmax = kwargs.get("vmax", None)
    cbar = kwargs.get("cbar", False)
    xlabel = kwargs.get("xlabel", "")
    ylabel = kwargs.get("ylabel", "")
    color = kwargs.get("color", "k")
    linestyle = kwargs.get("linestyle", "-")
    linewidth = kwargs.get("linewidth", 0.8)
    label = kwargs.get("label", "")
    legend = kwargs.get("legend", True)
    title_fontsize = kwargs.get("title_fontsize", 12)
    label_fontsize = kwargs.get("label_fontsize", 10)
    legend_fontsize = kwargs.get("legend_fontsize", 9)

    if scan_type == "polarimeter_scan":
        if type(ax) is type(None):
            fig, ax = plt.subplots(figsize = kwargs.get("figsize", (4,4)))
        ims = ax.imshow(intensity, extent = [eta2[0], eta2[-1], eta1[-1], eta1[0]], aspect = aspect, cmap = cmap, vmin = vmin, vmax = vmax)
        if cbar: _ = plt.colorbar(ims, ax = ax)
        ax.invert_yaxis()
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticks))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ticks))
        if ticks_minor > 0:
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
        if xlabel == "": xlabel = "$\eta_{2}$, °"
        if ylabel == "": ylabel = "$\eta_{1}$, °"
        ax.set_xlabel(xlabel, fontsize = label_fontsize)
        ax.set_ylabel(ylabel, fontsize = label_fontsize)
        if title == "":
            s1, s2, s3 = round(stokes[1],3), round(stokes[2],3), round(stokes[3],3)
            title = f"E = {energy} eV, S = [{s1}, {s2}, {s3}]"
            ax.set_title(title, fontsize = title_fontsize)
        else:
            ax.set_title(title, fontsize = title_fontsize)
    
    if scan_type in ["analyzer_scan", "retarder_scan", "special_case_scan"]:
        if type(ax) is type(None):
            fig, ax = plt.subplots(figsize = kwargs.get("figsize", (6,4)))
        x = scan.get("eta")
        ax.plot(x, intensity, color = color, linestyle = linestyle, linewidth = linewidth, label = label)
        if xlabel == "":
            if "analyzer" in scan_type: xlabel = "$\eta_{2}$, °"
            elif "retarder" in scan_type: xlabel = "$\eta_{1}$, °"
            elif "special" in scan_type: xlabel = "$\eta_{1}$, °" 
        ax.set_xlabel(xlabel, fontsize = label_fontsize)
        if ylabel == "": ylabel = "Intensity, a.u."
        ax.set_ylabel(ylabel, fontsize = label_fontsize)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticks))
        if ticks_minor > 0:
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
        if title == "":
            s1, s2, s3 = round(stokes[1],3), round(stokes[2],3), round(stokes[3],3)
            title = f"E = {energy} eV, S = [{s1}, {s2}, {s3}]"
            if "analyzer" in scan_type: title = f"Analyzer scan, {title}, $\eta_{1}$ = {scan.get('eta1')}°"
            elif "retarder" in scan_type: title = f"Retarder scan, {title}, $\eta_{2}$ = {scan.get('eta2')}°"
            elif "special" in scan_type:
                special_case = scan.get("special_case", None)
                if not type(special_case) is type(None):
                    title = f"Special case scan, case = {special_case} ($\Delta$ = {scan.get('delta')}°), {title}"
                else:
                    title = f"Diagonal scan, $\Delta$ = {scan.get('delta')}°, {title}"
            ax.set_title(title, fontsize = title_fontsize)
        if legend and not label == "": ax.legend(fontsize = legend_fontsize)
    
    elif scan_type == "special_case_scans":
        if type(ax) is type(None):
            fig, ax = plt.subplots(figsize = kwargs.get("figsize", (6,4)))
        x = scan.get("eta")
        ax.plot(x, scan.get("intensity1"), linestyle = linestyle, linewidth = linewidth, label = "case 1 ($\Delta$ = 0°)")
        ax.plot(x, scan.get("intensity2"), linestyle = linestyle, linewidth = linewidth, label = "case 2 ($\Delta$ = 90°)")
        ax.plot(x, scan.get("intensity3"), linestyle = linestyle, linewidth = linewidth, label = "case 3 ($\Delta$ = 45°)")
        ax.plot(x, scan.get("intensity4"), linestyle = linestyle, linewidth = linewidth, label = "case 4 ($\Delta$ = -45°)")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticks))
        if ticks_minor > 0:
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
        if xlabel == "": xlabel = "$\eta_{1}$, °"
        ax.set_xlabel(xlabel, fontsize = label_fontsize)
        if ylabel == "": ylabel = "Intensity, a.u."
        ax.set_ylabel(ylabel, fontsize = label_fontsize)
        if legend: ax.legend(fontsize = legend_fontsize)
        if title == "":
            s1, s2, s3 = round(stokes[1],3), round(stokes[2],3), round(stokes[3],3)
            title = f"Special case scans, E = {energy} eV, S = [{s1}, {s2}, {s3}]"
        ax.set_title(title, fontsize = title_fontsize)
    
    if not type(fig) is type(None): fig.tight_layout()
    return ax



def add2ax(scan = {}, ax = None, **kwargs):
    """
    """
    if not type(scan) is dict:
        logging.error("plotScan(): args scan must be a dict from a polarimeter scan.")
        return ax
    #
    eta, eta1, eta2, intensity = scan.get("eta", None), scan.get("eta1", None), scan.get("eta2", None), scan.get("intensity", None)
    energy, stokes = scan.get("energy", None), scan.get("stokes", None)
    scan_type = scan.get("scan_type", "")
    #
    if scan_type == "polarimeter_scan":
        logging.warning("plotScan(): can not add 2d data. Yet. Sorry.")
        return ax
    #
    ticks = kwargs.get("ticks", None)
    ticks_minor = kwargs.get("ticks_minor", None)
    title = kwargs.get("title", None)
    xlabel = kwargs.get("xlabel", None)
    ylabel = kwargs.get("ylabel", None)
    color = kwargs.get("color", "k")
    linestyle = kwargs.get("linestyle", "-")
    linewidth = kwargs.get("linewidth", 0.8)
    label = kwargs.get("label", "")
    legend = kwargs.get("legend", False)

    #
    if not type(title) is type(None): ax.set_title(title)
    if not type(ticks) is type(None): 
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticks))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ticks))
    if not type(ticks_minor) is type(None):
        if ticks_minor > 0:
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
    if not type(xlabel) is type(None): ax.set_xlabel(xlabel)
    if not type(ylabel) is type(None): ax.set_ylabel(ylabel)
    #
    dothekwargstuff = False
    if scan_type in ["analyzer_scan", "retarder_scan", "special_case_scan"]:
        ax.plot(eta, scan["intensity"], color = color, linestyle = linestyle, linewidth = linewidth, label = label)
        dothekwargstuff = True
        if legend: ax.legend()
    elif scan_type == "special_case_scans":
        dothekwargstuff = True
        ax.plot(eta, scan.get("intensity1"), linestyle = linestyle, linewidth = linewidth, label = "case 1 ($\Delta$ = 0°)")
        ax.plot(eta, scan.get("intensity2"), linestyle = linestyle, linewidth = linewidth, label = "case 2 ($\Delta$ = 90°)")
        ax.plot(eta, scan.get("intensity3"), linestyle = linestyle, linewidth = linewidth, label = "case 3 ($\Delta$ = 45°)")
        ax.plot(eta, scan.get("intensity4"), linestyle = linestyle, linewidth = linewidth, label = "case 4 ($\Delta$ = -45°)")
        if legend: ax.legend()
    #
    if dothekwargstuff:
        if not type(title) is type(None): ax.set_title(title)
        if not type(ticks) is type(None): 
            ax.xaxis.set_major_locator(ticker.MultipleLocator(ticks))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(ticks))
        if not type(ticks_minor) is type(None):
            if ticks_minor > 0:
                ax.xaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
                ax.yaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
        if not type(xlabel) is type(None): ax.set_xlabel(xlabel)
        if not type(ylabel) is type(None): ax.set_ylabel(ylabel)
    #
    return ax





def fitSpecialCase(**kwargs):
    """
    """
    x = kwargs.get("x", None)
    y = kwargs.get("y", None)
    if type(x) is type(None) or type(y) is type(None):
        scan = kwargs.get("scan", None)
        if type(scan) is dict:
            if len(scan.get("eta", np.array([]))) > 0: x = scan.get("eta", np.array([]))
            if len(scan.get("intensity", np.array([]))) > 0: y = scan.get("intensity", np.array([]))
        else:
            logging.error("fitSpecialCase(): Pass x and y as arrays (or a scan dict from Polarimeter).")
            return {}
    if type(x) is type(None) or type(y) is type(None):
        logging.error("fitSpecialCase(): Pass x and y as arrays (or a scan dict from Polarimeter).")
        return {}
    #
    I, V, Psi = kwargs.get("I", 5), kwargs.get("V", 0.5), kwargs.get("Psi", 0)
    par, cov = _fitSpecialCase(x = x, y = y, p = (I, V, Psi))
    #
    result = {"I": par[0],
            "V": par[1],
            "Psi": par[2],
            "par": tuple(par),
            "cov": tuple(cov),
            "method": _specialCase,
            "x": x,
            "y": y}
    
    return result
    



    #result = {"eta": np.array([]), "intensity": np.array([]), "special_case": None, "delta": None, "energy": self.energy, "stokes": self.stokes}
        
def _specialCase(x, *p):
    I, V, Psi = p
    return I * (1 + V * np.sin(2 * np.deg2rad(x) + np.deg2rad(Psi)))

def _fitSpecialCase(x, y, p):
    try:
        par, cov = curve_fit(_specialCase, x, y, p0 = p)
    except Exception as e:
        logging.warning("_fitSpecialCase(): Curve fitting failed.", exc_info = True)
        par = [np.NaN, np.NaN, np.NaN]
        cov = par.copy()
    return par, cov




# ===================================================================================================================
# ===================================================================================================================
# ===================================================================================================================























# ===================================================================================================================
# ===================================================================================================================
# ===================================================================================================================

class polarimeter_():
    """
    NOTE THAT: The data files have 4 columns from now on (energy, rs, rp, delta) instead of 5 as earlier.

    NOTE THAT: The method specialCases() is still under development. Note to self: check 'compare_walan_and_i.ipynb'. It works.

    Args for __init__():
        mirror_files:       list of 2 or 4 mirror files.
        path                path to the mirror files
        etc...


    """

    def __init__(self, mirror_files = ['mirror1.txt', 'mirror2.txt', 'mirror3.txt', 'mirror4.txt'], path = '', make_ideal = False, shup = False, **kwargs):
        #
        if kwargs.get('help', False):
            help(polarimeter)
            return
        #

        if not type(mirror_files) is list and not type(path) is str:
            print(Fore.RED + "myPolarimeter.polarimeter(): Argument 'mirror_files' must be a list of size 4  or 2 and argument 'path' must be a string." + Fore.BLACK)
            return
        if not(len(mirror_files) == 4 or len(mirror_files) == 2):
            print(Fore.RED + "myPolarimeter.polarimeter(): Argument 'mirror_files' must be a list containing 4 or 2 file names." + Fore.BLACK)
            return
        #
        # When running on my own computer...
        if path.lower() == "mats":
            global my_mirror_mirror_files, my_path
            mirror_files, path = my_mirror_mirror_files, my_path
            print(mirror_files, path)
        #
        self.mirror_files = []
        self.mirror_data = []
        for file in mirror_files:
            if path == "": self.mirror_files.append(path + file)
            else:
                if not path.endswith('/'): pafi = f"{path}/{file}"
                else: pafi = f"{path}{file}"
                self.mirror_files.append(pafi)
        #
        if not self._loadMirrorData():
            if not shup:
                print("\nPass file names for either four files (one for each mirror) in a list of strings (length 4) or two files")
                print("(one for the retarder and one for the analyzer) in a list of stings (length 2). The argument is: mirror_files.")
                print("Pass argument path as a sting.")
                print(Fore.BLUE)
                print("Example: p = polarimeter(mirror_files = ['retarder.dat', 'analyzer.dat'], path = 'mydata/polarimeter')")
                print(Fore.BLACK)
                print("For now: creating an ideal polarimeter (no attenuation in the retarder, no retardation in the analyzer)")
                print("if argument 'make_ideal' was True, or a typical polarimeter if False. Default is False.\n")
            _ = self._idealPolarimeter(ideal = make_ideal)
        #
        print("Energy range: {0} to {1} ({2} values)".format(min(self.mirror_data[0][0]), max(self.mirror_data[0][0]), len(self.mirror_data[0][0])))
        self.energies = self.mirror_data[0][0]
        self.setEnergy(self.energies[0], shup = shup)
        self.setStokes([1,0,0,1], shup = shup)
        if not shup:
            print(Fore.GREEN + "\nReady to go!" + Fore.BLACK)

        self.polarimeter_scan = {'energy': np.NaN, 'stokes': np.array([]), 'eta1': np.array([]), 'eta2': np.array([]), 'intensity': np.array([])}
        self.special_cases = {}

        if not shup:
            print(Fore.BLUE)
            print("Main methods:")
            print("  .setEnergy(energy = -1, shup = False)")
            print("  .setStokes(stokes = [1,0,0,1], norm = True, shup = False)")
            print("  .polarimeterIntensity(eta1 = 0, eta2 = 0, noise = 0, shup = False)")
            print("  .polarimeterScan(eta1 = [0, 360, 5], eta2 = [0, 360, 5], noise = 0, plot = True, shup = False, **kwargs)")
            print(Fore.BLUE)


    # -----------------
    def _loadMirrorData(self):
        """
        Load mirror data from 2 or 4 files.
        """
        ok = True
        Data = []
        for i, file in enumerate(self.mirror_files):
            try:
                data = np.loadtxt(file).transpose()
                print(Fore.GREEN + f"Loaded data from file {i+1} ({file})." + Fore.BLACK)
            except:
                data = np.array([])
                print(Fore.RED + f"Could not load data from file {i+1} ({file})." + Fore.BLACK)
                ok = False
            Data.append(data)
        self.mirror_data = np.copy(Data)
        return ok
        
    # -----------------
    def _idealPolarimeter(self, ideal = True):
        """
        ideal == True:
        create an ideal polarimeter with no attenuation in the retarder and no phase shift in the analyzer.
        ideal == False:
        create a kind of 'typical' polarimeter.
        """
        data1, data2, data3, data4 = [], [], [], []
        if ideal:
            for i, e in enumerate(np.linspace(10, 100, 91)):
                data1.append([e, 1, 1, -50 + i/4])
                data2.append([e, 1, 1, -55 + i/4])
                data3.append([e, 1, 1, -60 + i/4])
                rs = 1 - i/200
                rp = np.sqrt(1 - rs**2)
                data4.append([e, rs, rp, 0])
        else:
            for i, e in enumerate(np.linspace(10, 100, 91)):
                rs = 1 - i/400
                rp = np.sqrt(1 - rs**2)
                data1.append([e, rs, rp, -50 + i/4])
                data2.append([e, rs, rp, -55 + i/4])
                data3.append([e, rs, rp, -60 + i/4])
                rs = 1 - i/200
                rp = np.sqrt(1 - rs**2)
                data4.append([e, rs, rp, -20 + i/10])
        self.mirror_data = np.array([np.transpose(data1), np.transpose(data2), np.transpose(data3), np.transpose(data4)])
        self.mirror_files = ['(no file)', '(no file)', '(no file)', '(no file)']
    

    # -----------------
    def setEnergy(self, energy = -1, shup = False):
        """
        Sets an energy and calculates the mirror matrices, and/or the retarder and analyzer matrices.
        """
        try:
            energy = float(energy)
        except:
            print(Fore.RED + "myPolarimeter.polarimeter.setEnergy(): Argument 'energy' must be a positive number in the range {0} to {1}.".format(min(self.energies), max(self.energies)) + Fore.BLACK)
            return np.NaN
        if energy <= 0:
            print(Fore.RED + "myPolarimeter.polarimeter.setEnergy(): Argument 'energy' must be a positive number." + Fore.BLACK)
            return np.NaN
        indx = abs(energy - self.energies).argmin()
        energy_set = self.energies[indx]
        if not shup:
            if not energy == energy_set:
                print(f"The closest available energy to {energy} is {energy_set}.")
        self.energy = energy_set
        if not shup:
            print(f"Setting energy = {self.energy}.")
        #
        if len(self.mirror_files) == 2:
            self.retarder_matrix = mys.matrixFromOptical(rs = self.mirror_data[0][1][indx], rp = self.mirror_data[0][2][indx], delta = self.mirror_data[0][2][indx])
            self.analyzer_matrix = mys.matrixFromOptical(rs = self.mirror_data[1][1][indx], rp = self.mirror_data[1][2][indx], delta = self.mirror_data[1][2][indx])
        elif len(self.mirror_files) == 4:
            mm1 = mys.matrixFromOptical(rs = self.mirror_data[0][1][indx], rp = self.mirror_data[0][2][indx], delta = self.mirror_data[0][2][indx])
            mm2 = mys.matrixFromOptical(rs = self.mirror_data[1][1][indx], rp = self.mirror_data[1][2][indx], delta = self.mirror_data[1][2][indx])
            mm3 = mys.matrixFromOptical(rs = self.mirror_data[2][1][indx], rp = self.mirror_data[2][2][indx], delta = self.mirror_data[2][2][indx])
            mm4 = mys.matrixFromOptical(rs = self.mirror_data[3][1][indx], rp = self.mirror_data[3][2][indx], delta = self.mirror_data[3][2][indx])
            mm2 = mys.rotateMatrix(mm2, 180)
            self.retarder_matrix = mm3.dot(mm2.dot(mm1))
            self.analyzer_matrix = np.copy(mm4)

        return self.energy

    # -----------------
    def setStokes(self, stokes = [1,0,0,1], norm = True, shup = False):
        """
        Sets the Stokes vector.
        """
        if not (type(stokes) is list or type(stokes) is np.ndarray):
            print(Fore.RED + "myPolarimeter.polarimeter.setStokes(): Argument 'stokes' must be a list or array of length 4. (1)" + Fore.BLACK)
            return
        if not len(stokes) == 4:
            print(Fore.RED + "myPolarimeter.polarimeter.setStokes(): Argument 'stokes' must be a list or array of length 4. (2)" + Fore.BLACK)
            return
        stokes = np.array(stokes)
        if norm:
            stokes = mys.normalizeStokes(stokes, shup = shup)
        self.stokes = stokes
        if not shup:
            print(f"Setting Stokes = {stokes}") 
        return stokes
    
    # -------------------
    def retarderMatrix(self, angle = 0, shup = False):
        """
        Calculates the retarder matrix for a given angle.
        """
        return mys.rotateMatrix(self.retarder_matrix, angle, shup = shup)
    
    def analyzerMatrix(self, angle = 0, shup = False):
        """
        Calculates the analyzer matrix for a given angle.
        """
        return mys.rotateMatrix(self.analyzer_matrix, angle, shup = shup)
    
    def polarimeterMatrix(self, eta1 = 0, eta2 = 0, shup = False):
        """
        Calculated the full polarimeter matrix for two given angles.
        """
        R = self.retarderMatrix(angle = eta1)
        A = self.analyzerMatrix(angle = eta2)
        try:
            PM = A.dot(R)
        except Exception as e:
            print(Fore.RED + "myPolarimeter.polarimeter.polarimeterMatrix(): The polarimeter matrix could not be calculated." + Fore.BLACK)
            PM = np.ones([4,4])*np.NaN
        return PM
    
    def polarimeterIntensity(self, eta1 = 0., eta2 = 0., noise = 0, shup = False):
        P = self.polarimeterMatrix(eta1 = eta1, eta2 = eta2, shup = shup)
        I = P.dot(self.stokes)[0]
        if noise > 0:
            I = I * (1 + np.random.uniform(0, noise))
        return I

    
    # ----------------------
    def polarimeterScan(self, eta1 = [0, 360, 5], eta2 = [0, 360, 5], noise = 0, plot = True, shup = False, **kwargs):
        """
        Runs a polarimeter scan. Pass angular range for the retarder and analyzer (eta1, eta2) as [start, stop, step].
        Returns a dict with all data. This dict is also stored in the .polarimeter_scan attribute.

        Note: If you want to run ONLY the analyzer or ONLY the retarder, pass eta1 = [fixed angle] or eta2 = [fixed angle].

        Will update the class so that one can run the analyzer without a retarder present at all (and vice versa) at a later stage.

        Keyword arguments:
            figsize         tuple
            vmin, vmax      numbers
            ticks           number
            ticks_minor     number
            title           string
        """
        if not type(eta1) is list and type(eta2) is list:
            print(Fore.RED + "myPolarimeter.polarimeter.polarimeterScan(): eta1 and eta2 must be lists."); return {}
        if not ((len(eta1) == 3 or len(eta1) == 1) and (len(eta2) == 3 or len(eta2) == 1)):
            print(Fore.RED + "myPolarimeter.polarimeter.polarimeterScan(): eta1 and eta2 must be lists of length 1 or 3."); return {}
        if len(eta1) == 3:
            angle1 = np.linspace(eta1[0], eta1[1], int(abs((eta1[1]-eta1[0])/eta1[2])+1))
        else:
            angle1 = np.array(eta1)
        if len(eta2) == 3:
            angle2 = np.linspace(eta2[0], eta2[1], int(abs((eta2[1]-eta2[0])/eta2[2])+1))
        else:
            angle2 = np.array(eta2)
        imap = np.ones([len(angle1),len(angle2)]) * np.NaN
        for i1, a1 in enumerate(angle1):
            for i2, a2 in enumerate(angle2):
                imap[i1][i2] = self.polarimeterIntensity(eta1 = a1, eta2 = a2, noise = noise, shup = True)
        self.polarimeter_scan.update({'energy': self.energy})
        self.polarimeter_scan.update({'stokes': self.stokes})
        self.polarimeter_scan.update({'eta1': angle1})
        self.polarimeter_scan.update({'eta2': angle2})
        self.polarimeter_scan.update({'intensity': imap})
        #
        if plot:
            ticks = kwargs.get('ticks', 45)
            ticks_minor = kwargs.get('ticks_minor', 0)
            title = kwargs.get('title', 'none')
            if len(angle1) > 1 and len(angle2) > 1:
                figsize = kwargs.get('figsize', (3,3))
                vmin = kwargs.get('vmin', None)
                vmax = kwargs.get('vmax', None)
                _ = self.plotPolarimeterScan(scan = None, ax = None, figsize = figsize, vmin = vmin, vmax = vmax, 
                                                ticks = ticks, ticks_minor = ticks_minor, title = title)
            else:
                figsize = kwargs.get('figsize', (5,3))
                fig, ax = plt.subplots(figsize = figsize)
                if len(angle2) > 1:
                    ax.plot(angle2, imap[0])
                    ax.set_xlabel('$\eta_{2}$', size = 12)
                else:
                    ax.plot(angle1, self.polarimeter_scan['intensity'].transpose()[0])
                    ax.set_xlabel('$\eta_{1}$', size = 12)
        return deepcopy(self.polarimeter_scan)
    
    # ------------------------
    def plotPolarimeterScan(self, scan = None, ax = None, figsize = (3,3), vmin = None, vmax = None, ticks = 45, ticks_minor = 0,
                            title = 'none', cbar = False):
        """
        Plots a calculated polarimeter scan.
        """
        if type(scan) is type(None): scan = self.polarimeter_scan
        energy = scan.get('energy')
        stokes = scan.get('stokes')
        eta1 = scan.get('eta1')
        eta2 = scan.get('eta2')
        intensity = scan.get('intensity')
        if type(ax) is type(None):
            fig, ax = plt.subplots(figsize = figsize)
        extent = [eta2[0], eta2[-1], eta1[-1], eta1[0]]
        im = ax.imshow(intensity, extent = extent, aspect = 'equal', vmin = vmin, vmax = vmax)
        if cbar:
            plt.colorbar(im, ax = ax)
        ax.invert_yaxis()
        ax.xaxis.set_major_locator(ticker.MultipleLocator(ticks))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ticks))
        if ticks_minor > 0:
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(ticks_minor))
        if title == 'none':
            title = f"{self.energy} eV, {mys.stokes2string(self.stokes, 3)}"
            ax.set_title(title, size = 10)
        elif not title == '':
            ax.set_title(title)
        ax.set_xlabel('$\eta_{2}$', size = 12)
        ax.set_ylabel('$\eta_{1}$', size = 12, rotation = 0)
        return ax
    

    def analyzerScan(self, eta = [0,360,5], shup = False):
        """
        Analyzer scan as if there is no retarder present. Pass eta as scan range like [start, stop, step].
        Returns a dict containing the result (also stored also stored in the .analyzer_scan attribute).
        """
        result = {'energy': self.energy, 'stokes': self.stokes, 'eta': np.array([]), 'intensity': np.array([])}
        if not len(eta) == 3:
            print(Fore.RED + "myPolarimeter.polarimeter.analyzerScan(): eta must be a list [start, stop, step]."); return result
        if not eta[0] == eta[1] and not eta[2] > 0:
            print(Fore.RED + "myPolarimeter.polarimeter.analyzerScan(): eta must be a list [start, stop, step]."); return result
        angle = np.linspace(eta[0], eta[1], int(abs((eta[1]-eta[0])/eta[2])+1))
        result.update({'eta': angle})
        intensity = np.ones([len(angle)])*np.NaN
        result.update({'intensity': intensity})
        try:
            for i, a in enumerate(angle):
                intensity[i] = self.analyzerMatrix(angle = a, shup = True).dot(self.stokes)[0]
        except:
            result.update({'intensity': intensity})
            print(Fore.RED + "myPolarimeter.polarimeter.analyzerScan(): an unexpected error occured. probably due to sloppy coding."); return result
        result.update({'intensity': intensity})
        self.analyzer_scan = deepcopy(result)
        if not shup:
            print(f"Returning an analyzer scan in range {angle[0]}° to {angle[-1]}° for energy {self.energy} eV")
            print(f"and Stokes {self.stokes} as a dict. The scan is also stored in attribute .analyzer_scan.")
        return result
        


                

    

    # ---------------------------------------------------------------------

    def info(self):
        print(Fore.BLUE + "About the object:")
        print("  The optical data comes from:")
        for file in self.mirror_files:
            print(f"    {file}")
        print("Energy range:")
        print(f"  from {self.energies[0]} to {self.energies[-1]} eV ({len(self.energies)} values)")
        print("The current energy and polarization:")
        print(f"  .energy = {self.energy}")
        print(f"  .stokes = {self.stokes}")
        print("Last polarimeter scan:")
        if not len(self.polarimeter_scan['eta1']) > 0:
            print("  None.")
        else:
            print(f"  energy = {self.polarimeter_scan['energy']} eV")
            print(f"  stokes = {self.polarimeter_scan['stokes']}")
            print(f"  eta1: from {self.polarimeter_scan['eta1'][0]} to {self.polarimeter_scan['eta1'][-1]} ({len(self.polarimeter_scan['eta1'])} points)")
            print(f"  eta2: from {self.polarimeter_scan['eta2'][0]} to {self.polarimeter_scan['eta2'][-1]} ({len(self.polarimeter_scan['eta2'])} points)")
            print(f"  (accessable in .polarimeter_scan)")



    
    # ---------------------------------------------------------------------

    # ------------------------
    def _specialCaseIntensity(self, case = 1, angle = 0, noise = 0, shup = False):
        try:
            case = int(case)
        except:
            print(Fore.RED + "myPolarimeter.polarimeter._specialCaseIntensity(): Pass argument 'case' as 1, 2, 3, or 4 (integer)." + Fore.BLACK)
            return np.NaN
        try:
            angle = float(angle)
        except:
            print(Fore.RED + "myPolarimeter.polarimeter._specialCaseIntensity(): Pass argument 'angle' as a number (deg.)." + Fore.BLACK)
            return np.NaN
        if noise < 0: noise = 0
        angles = [0, 90, 45, -45]
        if case in [1,2,3,4]:
            intensity = self.polarimeterIntensity(angle1 = angle, angle2 = angle + angles[case-1], noise = noise, shup = True)
        else:
            print(Fore.RED + "myPolarimeter.polarimeter._specialCaseIntensity(): Pass argument 'case' as 1, 2, 3, or 4 (integer)." + Fore.BLACK)
            intensity = np.NaN
        if not shup:
            print("The intensity for special case {0} with angles \u03B71 = {1}° and \u03B72 = {2}° is {3}".format(case, angle, angle + angles[case -1], intensity))
        return intensity

    def specialCaseIntensity(self, case = 0, angle = 0, noise = 0, shup = False):
        """
        The intensity for the 4 special cases for a given angle. Case = 0 returns all 4 intensities.
        """
        if case == 0:
            I1 = self._specialCaseIntensity(case = 1, angle = angle, noise = noise, shup = shup)
            I2 = self._specialCaseIntensity(case = 2, angle = angle, noise = noise, shup = shup)
            I3 = self._specialCaseIntensity(case = 3, angle = angle, noise = noise, shup = shup)
            I4 = self._specialCaseIntensity(case = 4, angle = angle, noise = noise, shup = shup)
            return I1, I2, I3, I4
        else:
            return self._specialCaseIntensity(case = case, angle = 0, shup = shup)


    def fitSpecialCase(self, eta = None, intensity = None, p0 = [], shup = False):
        """
        p0 = [I, V, Psi]
        """
        if not (type(eta) is np.ndarray and type(intensity) is np.ndarray and (type(p0) is list or type(p0) is tuple)):
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCase(): Arguments 'eta' and iintensity' must be arrays")
            print("and argument 'p0' must be a list [I, V, Psi]" + Fore.BLACK)
            return [np.NaN, np.Nan, np.NaN]
        if len(p0) == 0:
            if not shup: print("As argument 'p0' was passed as an empty list the start parameters wil be estimated.")
        elif not len(p0) == 3:
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCase(): Arguments 'p0' must be a list [I, V, Psi]." + Fore.BLACK)
            return [np.NaN, np.Nan, np.NaN]
        if not len(eta) == len(intensity):
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCase(): Arguments 'eta' and 'intensity' must be of the same length." + Fore.BLACK)
            return [np.NaN, np.Nan, np.NaN]
        #
        if len(p0) == 0:
            p0 = [None, None, None]
        pass1, pass2, pass3 = '', '', ''
        _V = (intensity.max() - intensity.min()) / (intensity.max() + intensity.min()) 
        _I = intensity.mean()
        max_indx = abs(max(intensity) - intensity).argmin()
        _Psi = eta[max_indx] - 90
        while _Psi < 0: _Psi += 180
        while _Psi > 180: _Psi -= 180
        _Psi = np.mod(_Psi, 180)
        if type(p0[0]) is type(None):
            p0[0] = _I; pass1 = ' (estimated)'
        if type(p0[1]) is type(None):
            p0[1] = _V; pass2 = ' (estimated)'
        if type(p0[2]) is type(None):
            p0[2] = _Psi; pass3 = ' (estimated)'
        if not shup:
            print("myPolarimeter.polarimeter.fitSpecialCase():")
            print(f"  Start parameters: I = {p0[0]}{pass1}, V = {p0[1]}{pass2}, Psi = {p0[2]}{pass3}")
        #
        bnds = ([0, 0, -np.inf], [np.inf, np.inf, np.inf])
        par, cov = curve_fit(specialCaseEquation, xdata = eta, ydata = intensity, p0 = p0, bounds = bnds)
        perr = np.sqrt(np.diag(cov))
        if not shup:
            #while par[2] < 0: par[2] += 180
            print(f"  Fit result: I = {par[0]}, V = {par[1]}, Psi = {par[2]}, error = {perr}\n")            
        return par, perr
    

    # ---------------------------------------------------------------------

    def specialCases(self, step = 0, **kwargs):
        """

        optinal args:
            energy (current energy)
            eta1_start (0), eta1_stop (360), eta1_step (5)
            p0_lp_1, p0_lp_1, p0_lp_1, p0_lp_1 ([])  start params for lp
            stokes ([1,0,0,1])
            p0_1, p0_2, p0_3, p0_4 ([]), start params for abitrary polarization (step 5 and onwards)
        """
        print(Fore.MAGENTA + "\nmyPolarimeter.polarimeter.specialCases():")
        print("  This method is not ready. It does not work well (yet) for polarizations with a helical component.\n" + Fore.BLACK)
        #
        shup = kwargs.get('shup', False)
        #
        if not shup:
            print(Fore.BLUE)
            print("1. Generate data for special case 1 - 4 using linear polarized light.")
            print("   Do that by passing argument 'step = 1'. Optinally, pass argument 'energy'.")
            print("2. Fit those curves.")
            print("   Do that by passing argument 'step = 2'. Optinally, pass start parameters for the fits.")
            print("3. Extract the optical parameters for the polarimeter.")
            print("   Do that by passing argument 'step = 3'.")
            print("4. Generate data for an arbitrary polarization.")
            print("   Do that by passing argument 'step = 4' and a Stokes vector (stokes = [s0,s1,s2,s3]).")
            print("5. Fit the data.")
            print("   Do that by passing argument 'step = 5'.")
            print("6. Extract polarization data (with the help of the results from steps 1-3).")
            print("   Do that by passing argument 'step = 6'.")
            print(Fore.BLACK)
        #
        if not type(step) is int:
            print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Argument 'step' must be an integer." + Fore.BLACK)
            return
        if not step in [0, 1, 2, 3, 4, 5, 6]:
            print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Invalid value for argument 'step'." + Fore.BLACK)
            return
        #
        self._specialCases_clearDict(step)
        if step == 0:
            return
        #
        if step == 1:
            energy = kwargs.get('energy', np.NaN)
            if not np.isnan(energy):
                if energy > 0:
                    print()
                    _ = self.setEnergy(energy = energy, shup = False)
                else:
                    print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Invalid value for argument 'energy'." + Fore.BLACK)
            _ = self.setStokes(stokes = [1,1,0,0], shup = True)
            self.special_cases.update({'energy': self.energy})
            self.special_cases.update({'stokes_lp': self.stokes})
            print(Fore.GREEN)
            print(f"Energy = {self.energy} eV")
            eta_start = kwargs.get('eta_start', 0)
            eta_stop  = kwargs.get('eta_stop', 360)
            if eta_start > eta_stop: eta_start, eta_stop = eta_stop, eta_start
            eta_step  = abs(kwargs.get('eta_step', 2.5))
            eta = np.linspace(eta_start, eta_stop, int((eta_stop - eta_start)/eta_step +1) )
            self.special_cases.update({'eta': eta})
            print(f"eta range: {eta_start} to {eta_stop} (step {eta_step})")
            case1_lp = np.copy(eta) * np.NaN
            case2_lp = np.copy(case1_lp)
            case3_lp = np.copy(case1_lp)
            case4_lp = np.copy(case1_lp)
            noise = kwargs.get('noise', 0)
            for i, e in enumerate(eta):
                case1_lp[i] = self._specialCaseIntensity(case = 1, angle = e, noise = noise, shup = True)
                case2_lp[i] = self._specialCaseIntensity(case = 2, angle = e, noise = noise, shup = True)
                case3_lp[i] = self._specialCaseIntensity(case = 3, angle = e, noise = noise, shup = True)
                case4_lp[i] = self._specialCaseIntensity(case = 4, angle = e, noise = noise, shup = True)
            self.special_cases.update({'case1_lp': case1_lp, 'case2_lp': case2_lp, 'case3_lp': case3_lp, 'case4_lp': case4_lp})
            fig, ax = plt.subplots(figsize = (5,3))
            ax.scatter(x = eta, y = case1_lp, label = 'case 1', s = 2)
            ax.scatter(x = eta, y = case2_lp, label = 'case 2', s = 2)
            ax.scatter(x = eta, y = case3_lp, label = 'case 3', s = 2)
            ax.scatter(x = eta, y = case4_lp, label = 'case 4', s = 2)
            ax.set_title(f"Linear polarization, {self.energy} eV")
            ax.set_xlabel("$\eta_{1}$")
            ax.set_ylabel('Intensity')
            ax.set_ylim(0, None)
            ax.legend()
            fig.tight_layout()
            ret = {'eta': eta, 'case1': case1_lp, 'case2': case2_lp, 'case3': case3_lp, 'case4': case4_lp}
            return ret


        #
        if step == 2:
            if len(self.special_cases.get('case1_lp', [])) == 0:
                print(Fore.RED + "\nCan not find data for linear polarization. Did you do step 1?" + Fore.BLACK)
                return
            #
            p0_lp_1 = kwargs.get('p0_lp_1', [])
            p0_lp_2 = kwargs.get('p0_lp_2', [])
            p0_lp_3 = kwargs.get('p0_lp_3', [])
            p0_lp_4 = kwargs.get('p0_lp_4', [])
            if not (len(p0_lp_1) == 0 or len(p0_lp_1) == 3) or not (len(p0_lp_2) == 0 or len(p0_lp_2) == 3) or not (len(p0_lp_3) == 0 or len(p0_lp_3) == 3) or not (len(p0_lp_4) == 0 or len(p0_lp_4) == 3):
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): start parameters for fitting must be either [] or [I,V,Psi]." + Fore.BLACK)
                return
            print()
            eta = self.special_cases.get('eta')
            case1_lp_par, case1_lp_perr = self.fitSpecialCase(
                eta = eta, intensity = self.special_cases.get('case1_lp'), p0 = p0_lp_1, shup = shup)
            case2_lp_par, case2_lp_perr = self.fitSpecialCase(
                eta = eta, intensity = self.special_cases.get('case2_lp'), p0 = p0_lp_2, shup = shup)
            case3_lp_par, case3_lp_perr = self.fitSpecialCase(
                eta = eta, intensity = self.special_cases.get('case3_lp'), p0 = p0_lp_3, shup = shup)
            case4_lp_par, case4_lp_perr = self.fitSpecialCase(
                eta = eta, intensity = self.special_cases.get('case4_lp'), p0 = p0_lp_4, shup = shup)
            while case1_lp_par[2] < 0: case1_lp_par[2] += 180
            while case1_lp_par[2] > 180: case1_lp_par[2] -= 180
            while case2_lp_par[2] < 0: case2_lp_par[2] += 180
            while case2_lp_par[2] > 180: case2_lp_par[2] -= 180
            while case3_lp_par[2] < 0: case3_lp_par[2] += 180
            while case3_lp_par[2] > 180: case3_lp_par[2] -= 180
            while case4_lp_par[2] < 0: case4_lp_par[2] += 180
            while case4_lp_par[2] > 180: case4_lp_par[2] -= 180
            
            self.special_cases.update({'case1_lp_par': case1_lp_par, 'case1_lp_perr': case1_lp_perr})
            self.special_cases.update({'case2_lp_par': case2_lp_par, 'case2_lp_perr': case2_lp_perr})
            self.special_cases.update({'case3_lp_par': case3_lp_par, 'case3_lp_perr': case3_lp_perr})
            self.special_cases.update({'case4_lp_par': case4_lp_par, 'case4_lp_perr': case4_lp_perr})
            case1_lp_fit = specialCaseEquation(eta, case1_lp_par[0], case1_lp_par[1], case1_lp_par[2])
            case2_lp_fit = specialCaseEquation(eta, case2_lp_par[0], case2_lp_par[1], case2_lp_par[2])
            case3_lp_fit = specialCaseEquation(eta, case3_lp_par[0], case3_lp_par[1], case3_lp_par[2])
            case4_lp_fit = specialCaseEquation(eta, case4_lp_par[0], case4_lp_par[1], case4_lp_par[2])
            self.special_cases.update({'case1_lp_fit': case1_lp_fit, 'case2_lp_fit': case2_lp_fit})
            self.special_cases.update({'case3_lp_fit': case3_lp_fit, 'case4_lp_fit': case4_lp_fit})
            fig, ax = plt.subplots(figsize = (7,4))
            ax.scatter(x = eta, y = self.special_cases['case1_lp'], label = 'case 1', s = 2, color = 'tab:blue')
            ax.scatter(x = eta, y = self.special_cases['case2_lp'], label = 'case 2', s = 2, color = 'tab:orange')
            ax.scatter(x = eta, y = self.special_cases['case3_lp'], label = 'case 3', s = 2, color = 'tab:green')
            ax.scatter(x = eta, y = self.special_cases['case4_lp'], label = 'case 4', s = 2, color = 'tab:red')
            ax.plot(eta, self.special_cases['case1_lp_fit'], label = 'case 1 fit', linewidth = 0.7, linestyle = '-', color = 'tab:blue')
            ax.plot(eta, self.special_cases['case2_lp_fit'], label = 'case 2 fit', linewidth = 0.7, linestyle = '-', color = 'tab:orange')
            ax.plot(eta, self.special_cases['case3_lp_fit'], label = 'case 3 fit', linewidth = 0.7, linestyle = '-', color = 'tab:green')
            ax.plot(eta, self.special_cases['case4_lp_fit'], label = 'case 4 fit', linewidth = 0.7, linestyle = '-', color = 'tab:red')
            ax.set_title(f"Linear polarization, {self.energy} eV")
            ax.set_xlabel("$\eta_{1}$")
            ax.set_ylabel('Intensity')
            ax.set_ylim(0, None)
            ax.legend()
            fig.tight_layout()
            #
            print(Fore.GREEN)
            print("Fit results:")
            for i, par in enumerate([case1_lp_par, case2_lp_par, case3_lp_par, case4_lp_par]):
                print(f"  I{i+1} = {par[0]:6.4f}, V{i+1} = {par[1]:6.4f}, PSI{i+1} = {par[2]}°")
            print(Fore.BLACK)
            #
            if not shup:
                print(Fore.BLUE + "\nProceed with step 3 if the fittings are okay, otherwise adjust the start ")
                print("parameters and run step 2 again. The parameters can be passed as p0_lp_1,")
                print("p0_lp_2, p0_lp_3, and/or p0_lp_4 and shall have the form [I, V, Psi].")
            #
            ret = {'eta': eta, 'case1_par': case1_lp_par, 'case2_par': case2_lp_par, 'case3_par': case3_lp_par, 'case4_par': case4_lp_par}
            ret.update({'case1_perr': case1_lp_perr, 'case2_perr': case2_lp_perr, 'case3_perr': case3_lp_perr, 'case4_perr': case4_lp_perr})
            ret.update({'case1_fit': case1_lp_fit, 'case2_fit': case2_lp_fit, 'case3_fit': case3_lp_fit, 'case4_fit': case4_lp_fit})
            return ret
        #
        if step == 3:
            par1 = self.special_cases.get('case1_lp_par', [])
            par2 = self.special_cases.get('case2_lp_par', [])
            par3 = self.special_cases.get('case3_lp_par', [])
            par4 = self.special_cases.get('case4_lp_par', [])
            if par1 == [] or par2 == [] or par3 == [] or par4 == []:
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Missing fitting parameters. Did you run step 2?" + Fore.BLACK)
                return
            #
            I1 = par1[0]/par1[0]; V1 = par1[1]; PSI1 = np.deg2rad(par1[2])
            I2 = par2[0]/par1[0]; V2 = par2[1]; PSI2 = np.deg2rad(par2[2])
            I3 = par3[0]/par1[0]; V3 = par3[1]; PSI3 = np.deg2rad(par3[2])
            I4 = par4[0]/par1[0]; V4 = par4[1]; PSI4 = np.deg2rad(par4[2])

            # eq. 6
            PL = np.sqrt( (V1**2 * I1**2 - V2**2 * I2**2) / (I1**2 - I2**2))
            # eqs. 9 & 10
            S1 = PL * np.sin(PSI1)
            S2 = PL * np.cos(PSI1)
            # eqs. 7 & 8
            cos2psi3 = (V1*I1 + V2*I2) / (I1 + I2) / PL
            cos2psi4 = (V1*I1 - V2*I2) / (I1 + I2) / PL
            # eq.15 (check cos2psi3 and cos2psi4)
            eta = self.special_cases['eta']
            indx90 = abs(eta - 90).argmin()
            indx0  = abs(eta - 0).argmin()
            print(f"(debug: indx90 = {indx90}, indx0 = {indx0}, 90 = {eta[indx90]}, 0 = {eta[indx0]})")
            I1_90 = self.special_cases['case1_lp'][indx90] /par1[0]
            I1_0  = self.special_cases['case1_lp'][indx0]  /par1[0]
            I2_90 = self.special_cases['case2_lp'][indx90] /par1[0]
            I2_0  = self.special_cases['case2_lp'][indx0]  /par1[0]
            print(f"(debug: I1_90 = {I1_90}, I1_0 = {I1_0}, I2_90 = {I2_90}, I2_0 = {I2_0})")
            eq15_left_side = (I2_90 - I2_0) / (I1_90 + I1_0)
            eq15_right_side_1 = (cos2psi3 - cos2psi4) / (cos2psi3 + cos2psi4)
            eq15_right_side_2 = (cos2psi4 - cos2psi3) / (cos2psi4 + cos2psi3)
            if abs(eq15_left_side - eq15_right_side_1) > (eq15_left_side - eq15_right_side_2):  # should be > but...
                cos2psi3, cos2psi4 = cos2psi4, cos2psi3
            #
            psi3 = 0.5 * np.arccos(cos2psi3)
            sin2psi3 = np.sin(2*psi3)
            # eq. 11
            cosd3 = cos2psi3 / (cos2psi4 * sin2psi3) * (-S1 + S2 * np.tan(PSI3) / (S2 + S1 * np.tan(PSI3)))
            if abs(cosd3) >= 1:
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases():")
                print("  The value for delta3 can not be calculated. It will be 'estimated'." + Fore.BLACK)
                if cosd3 >= 1: cosd3 = 0.99999
                else: cosd3 = -0.99999
            self.special_cases.update({'lp_PL': PL, 'lp_S1': S1, 'lp_S2': S2})
            self.special_cases.update({'cos2psi3': cos2psi3, 'cos2psi4': cos2psi4})
            self.special_cases.update({'cosd3': cosd3})
            #
            print(Fore.GREEN)
            print("\nAccording to the fitted parameters for case 1 and case 2 are")
            print(f"  PL = {PL:.4f}  (should be 1)")
            print(f"  S1 = {S1:.4f}  (should be 1)")
            print(f"  S2 = {S2:.4f}  (should be 0)")
            print("\nThe optical parameters derived from case 1 and 2 are:")
            print(f"  cos(2*psi3) = {cos2psi3:.5f}")
            print(f"  cos(2*psi4) = {cos2psi4:.5f}")
            print("and from case 3:")
            print(f"  cos(d3)     = {cosd3:.5f}")
            #
            if not shup:
                print(Fore.BLUE + "\nProceed with step 4.")
                print(Fore.BLACK)
            #
            ret = {'PL': PL, 'S1': S1, 'S2': S2, 'cos2psi3': cos2psi3, 'cos2psi4': cos2psi4, 'cosd3': cosd3}
            return ret
        #
        if step == 4:
            if np.isnan(self.special_cases.get('cos2psi3', np.NaN)) or np.isnan(self.special_cases.get('cos2psi4', np.NaN)):
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Missing data for the optical parameters")
                print("They are not needed for this step but to not mess with my sense of order, please fix it." + Fore.BLACK)
                return
            #
            stokes = kwargs.get('stokes', None)
            if type(stokes) is type(None):
                if not shup: print("\nAs no Stokes vector was passed we set it to inclined polarization.")
                stokes = [1,0,1,0]
            if not len(stokes) == 4:
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Argument 'stokes' has to be a list or array of size 4." + Fore.BLACK)
                return
            _ = self.setStokes(stokes)
            self.special_cases.update({'stokes': self.stokes})
            #
            eta = self.special_cases.get('eta', np.array([]))
            imap = np.zeros([len(eta), len(eta)]) * np.NaN
            noise = kwargs.get('noise', 0)
            for i1, e1 in enumerate(eta):
                for i2, e2 in enumerate(eta):
                    imap[i1][i2] = self.polarimeterIntensity(angle1 = e1, angle2 = e2, noise = noise, shup = True)
            self.special_cases.update({'imap': imap})
            
            case1 = np.zeros([len(eta)]) * np.NaN
            case2, case3, case4 = np.copy(case1), np.copy(case1), np.copy(case1)
            for i, e in enumerate(eta):
                case1[i] = self._specialCaseIntensity(case = 1, angle = e, noise = noise, shup = True)
                case2[i] = self._specialCaseIntensity(case = 2, angle = e, noise = noise, shup = True)
                case3[i] = self._specialCaseIntensity(case = 3, angle = e, noise = noise, shup = True)
                case4[i] = self._specialCaseIntensity(case = 4, angle = e, noise = noise, shup = True)
            self.special_cases.update({'case1': case1})
            self.special_cases.update({'case2': case2})
            self.special_cases.update({'case3': case3})
            self.special_cases.update({'case4': case4})
            #
            fig, ax = plt.subplots(ncols = 3, figsize = (14,4))
            ax[0] = mys.plotEllipse(stokes = [stokes], ax = ax[0], norm = True, color = True, alpha = True,
                                    lw = 1, legend = False, figsize = (3,3), title = '', stokesastitle = True, 
                                    stokesaslabel = False)
            extent = [eta[0], eta[-1], eta[-1], eta[0]]
            im = ax[1].imshow(imap, extent = extent, aspect = "equal", vmin = 0)
            ax[1].invert_yaxis()
            ax[1].set_ylabel("$\eta_{1}$")
            ax[1].set_xlabel("$\eta_{2}$")
            ax[2].scatter(x = eta, y = case1, label = 'case 1', s = 2)
            ax[2].scatter(x = eta, y = case2, label = 'case 2', s = 2)
            ax[2].scatter(x = eta, y = case3, label = 'case 3', s = 2)
            ax[2].scatter(x = eta, y = case4, label = 'case 4', s = 2)
            ax[2].set_xlabel("$\eta_{1}$")
            ax[2].set_ylabel('Intensity')
            ax[2].set_ylim(0, None)
            ax[2].legend()
            ax[1].scatter(x = eta, y = eta, color = 'tab:blue', s = 3)
            ax[1].scatter(x = eta + 90, y = eta, color = 'tab:orange', s = 3)
            ax[1].scatter(x = eta + 45, y = eta, color = 'tab:green', s = 3)
            ax[1].scatter(x = eta - 45, y = eta, color = 'tab:red', s = 3)
            ax[1].set_xlim(eta.min(), eta.max())
            ax[1].set_ylim(eta.min(), eta.max())
            fig.tight_layout()
            #
            if not shup: print(Fore.BLUE + "\nProceed with step 5 (fit special cases)." + Fore.BLACK)
            #
            ret = {'eta': eta, 'map': imap, 'case1': case1, 'case2': case2, 'case3': case3, 'case4': case4}
            return ret
        #
        if step == 5:
            case1 = self.special_cases.get('case1', None)
            case2 = self.special_cases.get('case2', None)
            case3 = self.special_cases.get('case3', None)
            case4 = self.special_cases.get('case4', None)
            if type(case1) is type(None) or  type(case2) is type(None) or  type(case3) is type(None) or  type(case4) is type(None):
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Could not find data for the special cases. Did you do step 4?" + Fore.BLACK)
                return
            #
            p0_1 = kwargs.get('p0_1', [])
            p0_2 = kwargs.get('p0_2', [])
            p0_3 = kwargs.get('p0_3', [])
            p0_4 = kwargs.get('p0_4', [])
            if not (len(p0_1) == 0 or len(p0_1) == 3) or not (len(p0_2) == 0 or len(p0_2) == 3) or not (len(p0_3) == 0 or len(p0_3) == 3) or not (len(p0_4) == 0 or len(p0_4) == 3):
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): start parameters for fitting must be either [] or [I,V,Psi]." + Fore.BLACK)
                return
            print()
            eta = self.special_cases.get('eta')
            case1_par, case1_perr = self.fitSpecialCase(eta = eta, intensity = case1, p0 = p0_1, shup = shup)
            case2_par, case2_perr = self.fitSpecialCase(eta = eta, intensity = case2, p0 = p0_2, shup = shup)
            case3_par, case3_perr = self.fitSpecialCase(eta = eta, intensity = case3, p0 = p0_3, shup = shup)
            case4_par, case4_perr = self.fitSpecialCase(eta = eta, intensity = case4, p0 = p0_4, shup = shup)
            while case1_par[2] < 0: case1_par[2] += 180
            while case1_par[2] > 180: case1_par[2] -= 180
            while case2_par[2] < 0: case2_par[2] += 180
            while case2_par[2] > 180: case2_par[2] -= 180
            while case3_par[2] < 0: case3_par[2] += 180
            while case3_par[2] > 180: case3_par[2] -= 180
            while case4_par[2] < 0: case4_par[2] += 180
            while case4_par[2] > 180: case4_par[2] -= 180
            self.special_cases.update({'case1_par': case1_par, 'case1_perr': case1_perr})
            self.special_cases.update({'case2_par': case2_par, 'case2_perr': case2_perr})
            self.special_cases.update({'case3_par': case3_par, 'case3_perr': case3_perr})
            self.special_cases.update({'case4_par': case4_par, 'case4_perr': case4_perr})
            case1_fit = specialCaseEquation(eta, case1_par[0], case1_par[1], case1_par[2])
            case2_fit = specialCaseEquation(eta, case2_par[0], case2_par[1], case2_par[2])
            case3_fit = specialCaseEquation(eta, case3_par[0], case3_par[1], case3_par[2])
            case4_fit = specialCaseEquation(eta, case4_par[0], case4_par[1], case4_par[2])
            self.special_cases.update({'case1_fit': case1_fit, 'case2_fit': case2_fit})
            self.special_cases.update({'case3_fit': case3_fit, 'case4_fit': case4_fit})
            fig, ax = plt.subplots(figsize = (7,4))
            ax.scatter(x = eta, y = self.special_cases['case1'], label = 'case 1', s = 2, color = 'tab:blue')
            ax.scatter(x = eta, y = self.special_cases['case2'], label = 'case 2', s = 2, color = 'tab:orange')
            ax.scatter(x = eta, y = self.special_cases['case3'], label = 'case 3', s = 2, color = 'tab:green')
            ax.scatter(x = eta, y = self.special_cases['case4'], label = 'case 4', s = 2, color = 'tab:red')
            ax.plot(eta, self.special_cases['case1_fit'], label = 'case 1 fit', linewidth = 1, linestyle = '-', color = 'tab:blue')
            ax.plot(eta, self.special_cases['case2_fit'], label = 'case 2 fit', linewidth = 1, linestyle = '-', color = 'tab:orange')
            ax.plot(eta, self.special_cases['case3_fit'], label = 'case 3 fit', linewidth = 1, linestyle = '-', color = 'tab:green')
            ax.plot(eta, self.special_cases['case4_fit'], label = 'case 4 fit', linewidth = 1, linestyle = '-', color = 'tab:red')
            ax.set_title(f"S = {self.special_cases['stokes']}, E = {self.energy} eV")
            ax.set_xlabel("$\eta_{1}$")
            ax.set_ylabel('Intensity')
            ax.set_ylim(0, None)
            ax.legend()
            fig.tight_layout()
            #
            print(Fore.GREEN)
            print("Fit results:")
            for i, par in enumerate([case1_par, case2_par, case3_par, case4_par]):
                print(f"  I{i+1} = {par[0]:6.4f}, V{i+1} = {par[1]:6.4f}, PSI{i+1} = {par[2]}°")
            print(Fore.BLACK)
            #
            if not shup:
                print(Fore.BLUE + "\nProceed with step 6 if the fittings are okay, otherwise adjust the start parameters and run step 5 again.")
                print("The parameters can be passed as p0_1, p0_2, p0_3, and/or p0_4 and shall have the form [I, V, Psi].")
            #
            ret = {'eta': eta, 'case1_par': case1_par, 'case2_par': case2_par, 'case3_par': case3_par, 'case4_par': case4_par}
            ret.update({'case1_perr': case1_perr, 'case2_perr': case2_perr, 'case3_perr': case3_perr, 'case4_perr': case4_perr})
            ret.update({'case1_fit': case1_fit, 'case2_fit': case2_fit, 'case3_fit': case3_fit, 'case4_fit': case4_fit})
            return ret

        #
        if step == 6:
            par1 = self.special_cases.get('case1_par', None)
            par2 = self.special_cases.get('case2_par', None)
            par3 = self.special_cases.get('case3_par', None)
            par4 = self.special_cases.get('case4_par', None)
            if type(par1) is type(None) or type(par2) is type(None) or type(par3) is type(None) or type(par4) is type(None):
                print(Fore.RED + "\nmyPolarimeter.polarimeter.specialCases(): Could not find fits for the special cases. Did you do step 5?" + Fore.BLACK)
                return
            I1 = par1[0]; V1 = par1[1]; PSI1 = np.deg2rad(par1[2])
            I2 = par2[0]; V2 = par2[1]; PSI2 = np.deg2rad(par2[2])
            I3 = par3[0]; V3 = par3[1]; PSI3 = np.deg2rad(par3[2])
            I4 = par4[0]; V4 = par4[1]; PSI4 = np.deg2rad(par4[2])

            #cos2psi3 = self.special_cases['cos2psi3']
            #cos2psi4 = self.special_cases['cos2psi4']
            #psi3 = 0.5 * np.arccos(cos2psi3)
            #sin2psi3 = np.sin(2*psi3)
            #cosd3 = self.special_cases['cosd3']
            #d3 = np.arccos(cosd3)
            # eq. 6
            PL = np.sqrt( (V1**2 * I1**2 - V2**2 * I2**2) / (I1**2 - I2**2))
            # eqs. 9 & 10
            S1 = PL * np.sin(PSI1)
            S2 = PL * np.cos(PSI1)
            # eqs. 7 & 8
            cos2psi3 = (V1*I1 + V2*I2) / (I1 + I2) / PL
            cos2psi4 = (V1*I1 - V2*I2) / (I1 + I2) / PL
            # eq.15 (check cos2psi3 and cos2psi4)
            eta = self.special_cases['eta']
            indx90 = abs(eta - 90).argmin()
            indx0  = abs(eta - 0).argmin()
            print(f"(debug: indx90 = {indx90}, indx0 = {indx0}, 90 = {eta[indx90]}, 0 = {eta[indx0]})")
            I1_90 = self.special_cases['case1_lp'][indx90]
            I1_0  = self.special_cases['case1_lp'][indx0]
            I2_90 = self.special_cases['case2_lp'][indx90]
            I2_0  = self.special_cases['case2_lp'][indx0]
            print(f"(debug: I1_90 = {I1_90}, I1_0 = {I1_0}, I2_90 = {I2_90}, I2_0 = {I2_0})")
            eq15_left_side = (I2_90 - I2_0) / (I1_90 + I1_0)
            eq15_right_side_1 = (cos2psi3 - cos2psi4) / (cos2psi3 + cos2psi4)
            eq15_right_side_2 = (cos2psi4 - cos2psi3) / (cos2psi4 + cos2psi3)
            if abs(eq15_left_side - eq15_right_side_1) > (eq15_left_side - eq15_right_side_2):  # should be > but...
                cos2psi3, cos2psi4 = cos2psi4, cos2psi3
            #
            print(Fore.GREEN)
            print("From cases 1 and 2:")
            print(f"  PL = {PL:7.4f}")
            print(f"  S1 = {S1:7.4f}")
            print(f"  S2 = {S2:7.4f}")
            print(f"  cos(2*psi3) = {cos2psi3:8.5f} ({self.special_cases['cos2psi3']:.5f})")
            print(f"  cos(2*psi4) = {cos2psi4:8.5f} ({self.special_cases['cos2psi4']:.5f})" + Fore.BLACK)
            if not 0.95 < (cos2psi3/self.special_cases['cos2psi3']) and (cos2psi3/self.special_cases['cos2psi3']) < 1.05:
                print(Fore.MAGENTA + "Warning: the value for cos(2*psi3) is not the same as for the LP case." + Fore.BLACK)
            if not 0.95 < (cos2psi4/self.special_cases['cos2psi4']) and (cos2psi4/self.special_cases['cos2psi4']) < 1.05:
                print(Fore.MAGENTA + "Warning: the value for cos(2*psi3) is not the same as for the LP case." + Fore.BLACK)
            #
            psi3 = 0.5 * np.arccos(cos2psi3)
            sin2psi3 = np.sin(2*psi3)
            cosd3 = self.special_cases['cosd3']
            #
            cosd3_ = cos2psi3 / (cos2psi4 * sin2psi3) * (S2 * np.tan(PSI3) - S1)/(S2 + S1 * np.tan(PSI3))
            print(f"  cos(d3)     = {cosd3_:8.5f} ({cosd3:.5f})")
            return
        
            # eq. 12
            S3 = (I4-I3)/(I4+I3) /(cos2psi4 * np.sin(2*psi3) * np.sin(d3))
            # eq. 13
            PL = I4 * V4 / (I3 + I4)
            PL = PL * 2 / np.sqrt(cosd3**2 * sin2psi3**2 * cos2psi4**2 + cos2psi3**2)
            # eq. 14 (S1/S2 = ...)
            S1_S2 = (np.cos(PSI3) + np.cos(PSI4)) / (np.sin(PSI3) + np.sin(PSI4))
            S1_S2 = 1/S1_S2     # ================= ???
            # PL^2 = S1^2 + S2^2:
            S1 = PL / np.sqrt(1 + 1 / S1_S2**2)
            S2 = PL / np.sqrt(1 + S1_S2**2)
            #
            S0 = np.sqrt(S1**2 + S2**2 + S3**2)
            #

            print("\nFrom case 3 and 4 (using the optical parameters from step 3):")
            print(f"  PL = {PL:7.4f}")
            print(f"  S0 = {S0:7.4f}, s0 = {S0/S0:7.4f}")
            print(f"  S1 = {S1:7.4f}, s1 = {S1/S0:7.4f}")
            print(f"  S2 = {S2:7.4f}, s2 = {S2/S0:7.4f}")
            print(f"  S3 = {S3:7.4f}, s3 = {S3/S0:7.4f}" + Fore.BLACK)
                
            self.special_cases.update({'S0': S0, 'S1': S1, 'S2': S2, 'S3': S3, 'PL': PL})
            norm_stokes = mys.normalizeStokes([S0,S1,S2,S3], shup = True)
            self.special_cases.update({'STOKES': norm_stokes})

            # Using case 1 and 2 (redundant)
            PL_ = np.sqrt( (V1**2 * I1**2 - V2**2 * I2**2) / (I1**2 - I2**2))
            S1_ = PL * np.sin(PSI1)
            S2_ = PL * np.cos(PSI1)

            print(Fore.GREEN)
            print("From case 1 and 2:")
            print(f"  PL = {PL_:7.4f}")
            print(f"  S1 = {S1_:7.4f}")
            print(f"  S2 = {S2_:7.4f}" + Fore.BLACK)

            if not shup:
                print(Fore.BLUE + "\n----------\n")
                print("Result:")
                print("  Optical parameters:")
                print(f"    cos(2*psi3) = {cos2psi3:8.5f}")
                print(f"    cos(2*psi4) = {cos2psi4:8.5f}")
                print(f"    cos(d3)     = {cosd3:8.5f}")
                print("  Polarization:")
                print(f"    PL = {PL:.4f}")
                print(f"    s1 = {norm_stokes[1]:7.4f}")
                print(f"    s2 = {norm_stokes[2]:7.4f}")
                print(f"    s3 = {norm_stokes[3]:7.4f}")
                print(Fore.BLACK)

            _ = mys.plotEllipse(stokes = [norm_stokes, self.special_cases['stokes']], ax = None, norm = True, color = True, alpha = True,
               lw = 1, legend = False, figsize = (3,3), title = 'Result', stokesastitle = False, stokesaslabel = False)
            
            if not shup:
                print(Fore.BLUE)
                print("(Done. The results are stored in dict .special_cases)")
                print(Fore.BLACK)
            
            ret = {'info': "Still working on the return from case 6..."}
            return ret

    def _specialCases_clearDict(self, step):
        """Helper"""
        if step <= 0 or step > 6:
            self.special_cases = {}
            return
        if step in [6]:
            # do nothing (for the moment)
            return
        if step in [1,2,3,4,5]:
            for poop in ['STOKES', 'S0', 'S1', 'S2', 'S3', 'PL']: _ = self.special_cases.pop(poop, None)
        if step in [1,2,3,4]:
            for poop in ['case1_par',  'case2_par',  'case3_par',  'case4_par']: _ = self.special_cases.pop(poop, None)
            for poop in ['case1_perr', 'case2_perr', 'case3_perr', 'case4_perr']: _ = self.special_cases.pop(poop, None)
            for poop in ['case1_fit',  'case2_fit',  'case3_fit',  'case4_fit']: _ = self.special_cases.pop(poop, None)
        if step in [1,2,3]:
            for poop in ['imap', 'case1', 'case2', 'case3', 'case4']: _ = self.special_cases.pop(poop, None)
        if step in [1,2]:
            for poop in ['lp_PL', 'lp_S1', 'lp_S2', 'cos2psi3', 'cos2psi4']: _ = self.special_cases.pop(poop, None)
        if step in [1]:
            for poop in ['case1_lp_par', 'case2_lp_par', 'case3_lp_par', 'case4_lp_par']: _ = self.special_cases.pop(poop, None)
            for poop in ['case1_lp_perr', 'case2_lp_perr', 'case3_lp_perr', 'case4_lp_perr']: _ = self.special_cases.pop(poop, None)
            for poop in ['case1_lp_fit', 'case2_lp_fit', 'case3_lp_fit', 'case4_lp_fit']: _ = self.special_cases.pop(poop, None)
        return
            
            
            



            







            


        #
         


    def fitSpecialCases_(self, eta = [], intensity = [], p0 = [], shup = False, plot = True, DEBUG = True):
        """
        eta:        a list of angle arrays, i.e. [eta1, eta2, ...]
        intensity:  a list of intensities, i.e. [int1, int2, ...]
        p0:         a list of start parameters, i.e. [[I1,V1,Psi1], [I2,V2,Psi2], ...]
        """
        print(Fore.RED)
        print("Note: this is under development and it does not work yet. Ironing out issues.")
        print(Fore.BLACK)

        nan = np.NaN
        result = {'S0': nan, 'S1': nan, 'S2': nan, 'S3': nan, 'pl': nan, 'cos2psi3': nan, 'cos2psi4': nan, 'cosd3': nan, 'pars': nan}

        if not(type(eta) is list and type(intensity) is list and type(p0) is list):
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCases(): Arguments 'eta', 'intensity', and 'p0' must be lists." + Fore.BLACK)
            return result
        if not (len(eta) == 4 and len(intensity) == 4 and len(p0) == 4):
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCases(): Arguments 'eta', 'intensity', and 'p0' must be lists of length 4." + Fore.BLACK)
            return result
        ok = True
        for i in range(4):
            if not type(eta[i]) is np.ndarray: ok = False
            if not type(intensity[i]) is np.ndarray: ok = False
            if not type(p0[i]) is list: ok = False
        if not ok:
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCases(): The elements in 'eta', 'intensity', and 'p0' must be arrays, arrays, and lists." + Fore.BLACK)
            return result
        ok = True
        for i in range(4):
            if not len(eta[i]) == len(intensity[i]): ok = False
            if not(len(p0[i]) == 0 or len(p0[i]) == 3): ok = False
        if not ok:
            print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCases(): ..." + Fore.BLACK)
            return result
        #
        if plot:
            fig, ax = plt.subplots(figsize = (8,4))
            colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
            for i in [0,1,2,3]:
                ax.plot(eta[i], intensity[i], linewidth = 0.6, color = colors[i], label = f"case {i+1}")
            ax.legend()
            ax.set_xlabel('$\eta_{1}$'); ax.set_ylabel('intensity')
            ax.set_ylim(0, None)
        #
        pars, perrs = [], []
        for i in range(4):
            par, perr = self.fitSpecialCase(eta[i], intensity[i], p0[i], shup = shup)
            pars.append(par)
            perrs.append(perr)
        #
        result.update({'pars': pars, 'perrs': perrs})
        for perr in perrs:
            if np.inf in perr:
                print(Fore.RED + "myPolarimeter.polarimeter.fitSpecialCases(): Bad fittings." + Fore.BLACK)
                return result


        if plot:
            for i in [0,1,2,3]:
                ax.plot(eta[0], specialCaseEquation(eta[i], pars[i][0], pars[i][1], pars[i][2]), colors[i], linestyle = '-.', label = f"fit {i+1}")
            ax.legend()
        if not shup:
            print(F"myPolarimeter.polarimeter.fitSpecialCases():")
            print("  Fit results:")
            for i in [0,1,2,3]:
                #pars[i][2] = np.mod(pars[i][2], 180)
                print(f"  I{i+1}, V{i+1}, Psi{i+1} = {pars[i][0], pars[i][1], pars[i][2]}, errors = {perrs[i]}")
            print()

        #
        I1 = pars[0][0]; V1 = pars[0][1]; PSI1 = np.deg2rad(pars[0][2])
        I2 = pars[1][0]; V2 = pars[1][1]; PSI2 = np.deg2rad(pars[1][2])
        I3 = pars[2][0]; V3 = pars[2][1]; PSI3 = np.deg2rad(pars[2][2])
        I4 = pars[3][0]; V4 = pars[3][1]; PSI4 = np.deg2rad(pars[3][2])
        if DEBUG:
            print(Fore.MAGENTA + f"  Debug:\n    I1 = {I1}, V1 = {V1}, PSI1 = {PSI1}, {np.rad2deg(PSI1)}")
            print(f"    I2 = {I2}, V2 = {V2}, PSI2 = {PSI2}, {np.rad2deg(PSI2)}")
            print(f"    I3 = {I3}, V3 = {V3}, PSI3 = {PSI3}, {np.rad2deg(PSI3)}")
            print(f"    I4 = {I4}, V4 = {V4}, PSI4 = {PSI4}, {np.rad2deg(PSI4)}" + Fore.BLACK)

        PL = np.sqrt( (V1**2 * I1**2 - V2**2 * I2**2) / (I1**2 - I2**2))
        if DEBUG:
            print(Fore.MAGENTA + f"\n    PL = {PL}" + Fore.BLACK)

        cos2psi3_1 = (V1*I1 + V2*I2) / (I1 + I2) / PL
        cos2psi3_2 = (V1*I1 - V2*I2) / (I1 + I2) / PL
        cos2psi4_1 = (V1*I1 - V2*I2) / (I1 + I2) / PL
        cos2psi4_2 = (V1*I1 + V2*I2) / (I1 + I2) / PL
        if DEBUG:
            print(Fore.MAGENTA + f"\n    cos2psi3_1 = {cos2psi3_1}, cos2psi3_2 = {cos2psi3_2}")
            print(f"    cos2psi4_1 = {cos2psi4_1}, cos2psi4_2 = {cos2psi4_2}" + Fore.BLACK)

        S1 = PL * np.sin(PSI1)
        S2 = PL * np.cos(PSI1)
        if DEBUG:
            print(Fore.MAGENTA + f"\n    S1 = {S1}")
            print(f"    S2 = {S2}" + Fore.BLACK)

        # - extra
        if DEBUG:
            ratioS1S2 = (np.cos(PSI3) + np.cos(PSI4)) / (np.sin(PSI3) + np.sin(PSI4))
            print(Fore.MAGENTA + f"\n    ratioS1S2 = {ratioS1S2}")
            print(f"    check w. S1 and S2: S1/S2 = {S1/S2}" + Fore.BLACK)

        # eq15
        indx1_90, indx1_0 = abs(eta[0] - 90).argmin(), abs(eta[0] - 0).argmin()
        indx2_90, indx2_0 = abs(eta[1] - 90).argmin(), abs(eta[1] - 0).argmin()
        if DEBUG:
            print(Fore.MAGENTA + f"\n    eta90 = {eta[0][indx1_90]}, eta0 = {eta[0][indx1_0]}")
            print(f"    eta90 = {eta[1][indx2_90]}, eta0 = {eta[1][indx2_0]}" + Fore.BLACK)
        eq15LS = (intensity[1][indx2_90] - intensity[1][indx2_0]) / (intensity[0][indx1_90] + intensity[0][indx1_0]) 
        eq15RS1 = (cos2psi3_1 - cos2psi4_1) / (cos2psi3_1 + cos2psi4_1)
        eq15RS2 = (cos2psi3_2 - cos2psi4_2) / (cos2psi3_2 + cos2psi4_2)
        if DEBUG:
            print(Fore.MAGENTA + f"\n    eq15LS = {eq15LS}, eq15RS1 = {eq15RS1}, eq15RS2 = {eq15RS2}" + Fore.BLACK)
        
        if abs(eq15LS - eq15RS1) > abs(eq15LS - eq15RS2):
            cos2psi3 = cos2psi3_1
            cos2psi4 = cos2psi4_1
        else:
            cos2psi3 = cos2psi3_2
            cos2psi4 = cos2psi4_2
        if DEBUG:
            print(Fore.MAGENTA + f"\n    cos2psi3 = {cos2psi3}, cos2psi4 = {cos2psi4}" + Fore.BLACK)
        
        psi3 = 0.5 * np.arccos(cos2psi3)
        if DEBUG: 
            print(Fore.MAGENTA + f"\n    psi3 = {psi3}, {np.rad2deg(psi3)}°" + Fore.BLACK)
        cosd3 = cos2psi3 / (cos2psi4 * np.sin(2*psi3)) * (-S1 + S2 * np.tan(PSI3) / (S2 + S1 * np.tan(PSI3)))
        if DEBUG:
            print(Fore.MAGENTA + f"\n    cosd3 = {cosd3}" + Fore.BLACK)
        d3 = np.arccos(cosd3)
        if DEBUG:
            print(Fore.MAGENTA + f"\n    d3 = {d3}, {np.rad2deg(d3)}°" + Fore.BLACK)
        S3 = 1/(cos2psi4 * np.sin(2*psi3) * np.sin(d3)) * (I4-I3)/(I4+I3)
        if DEBUG:
            print(Fore.MAGENTA + f"\n    S3 = {S3}°" + Fore.BLACK)

        S0 = np.sqrt(S1**2 + S2**2 + S3**2)
        if DEBUG:
            print(Fore.MAGENTA + f"\n    S0 = {S0}" + Fore.BLACK)

        result.update({'S0': S1, 'S1': S1, 'S2': S2, 'S3': S3})
        result.update({'PL': PL, 'cos2psi3': cos2psi3, 'cos2psi4': cos2psi4, 'cosd3': cosd3})
        result.update({'pars': pars, 'perrs': perrs})

        if not shup:
            print('\n  Result:')
            print(f"  S1 = {S1}, {S1/S0}")
            print(f"  S2 = {S2}, {S2/S0}")
            print(f"  S3 = {S3}, {S3/S0}")


        return result








    

def specialCaseEquation(eta, I, V, Psi):
    """
    intensity = I * (1 + V * sin(2*eta + Psi)
    """
    return I * (1 + V * np.sin(2 * np.deg2rad(eta)+ np.deg2rad(Psi)))
                    

def polarimeterEquation(S0, S1, S2, S3, phi_3, delta_3, phi_4, eta_1, eta_2):

    from numpy import sin, cos

    return (
    ((sin(2*eta_1)*sin(2*eta_2)+cos(2*eta_1)*cos(2*eta_2))*cos(2*phi_3)*cos(2*phi_4)+1)*S0 +
    (((cos(delta_3)*sin(2*eta_1)**2*cos(2*eta_2)-cos(delta_3)*cos(2*eta_1)*sin(2*eta_1)*sin(2*eta_2))*
    sin(2*phi_3)+cos(2*eta_1)*sin(2*eta_1)*sin(2*eta_2)+cos(2*eta_1)**2*cos(2*eta_2))*cos(2*phi_4)+cos(2*eta_1)*cos(2*phi_3))*S1+
    (((cos(delta_3)*cos(2*eta_1)**2*sin(2*eta_2)-cos(delta_3)*cos(2*eta_1)*sin(2*eta_1)*cos(2*eta_2))*
    sin(2*phi_3)+sin(2*eta_1)**2*sin(2*eta_2)+cos(2*eta_1)*sin(2*eta_1)*cos(2*eta_2))*cos(2*phi_4)+sin(2*eta_1)*cos(2*phi_3))*S2+
    (sin(delta_3)*cos(2*eta_1)*sin(2*eta_2)-sin(delta_3)*sin(2*eta_1)*cos(2*eta_2))*
    sin(2*phi_3)*cos(2*phi_4)*S3
    )

def polarimeterpsi(rs = 1, rp = 1):
    """
    This is psi in cos(2*psi) that appeares in Walan's simplified analysis (cos(2*psi3) and cos(2*psi4)).
    """
    return np.arctan(rp/rs)











# ========================================================================= Plot polarimeter data
        
def plotPolInt(eta1 = None, eta2 = None, intensity = None, ax = None, **kwargs):
    """
    Plots (1) polarimeter line scan(s) or (2) a polarimeter intensity maps. Returns a pyplot ax.

    (1) Pass eta1 as a (n,) list/array and intensity as a (n,) or (N,n) list/array. n is the
        length of eta1 and N is the number of curves.
    
    (2) Pass eta1 as a (n,) list array, eta2 as a (m,) list/array and intensity as a (n,m) list array.

    Pass ax as a pyplot subplot ax, or one will be created.

    Optional arguments:
        figsize     tuple       (1),(2)     (5,3),(3,3)
        linewidth   number      (1)         0.75
        xlabel      string      (1),(2)     eta1
        ylabel      string      (1),(2)     intensity, eta2
        title       string      (1),(2)
        legend      bool        (1)         True for many scans, False for one scan
        vmin        number      (2)         None
        vmax        number      (2)         None
        aspect      string      (2)         equal
        cmap        string      (2)         viridis
        colorbar    bool        (2)         False
        xticks      number      (1),(2)
        yticks      number      (1)
        xlim        tuple       (1),(2)
        ylim        tuple       (1),(2)
        scatter     bool        (1)         False
        marker      string      (1)         x
        s           number      (1)         10
    """
    if kwargs.get('help', False):
        help(plotPolInt)
        return
    #
    if type(eta1) is list or type(eta1) is np.ndarray: eta1 = np.array(eta1)
    else: eta1 = np.zeros([])
    if type(eta1) is list or type(eta1) is np.ndarray: eta2 = np.array(eta2)
    else: eta2 = np.zeros([])
    # ------------
    # figure out what the user is trying to plot (2d intensity or line scan(s))
    if len(np.shape(eta1)) == 1 and len(np.shape(eta2)) == 0 and len(np.shape(intensity)) == 1:
        kind_of_plot = "line_scan"
    elif len(np.shape(eta1)) == 1 and len(np.shape(eta2)) == 0 and len(np.shape(intensity)) == 2:
        kind_of_plot = "line_scans"
    elif len(np.shape(eta1)) == 1 and len(np.shape(eta2)) == 1 and len(np.shape(intensity)) == 2:
        kind_of_plot = "intensity_map"
    else:
        print(Fore.RED + "plotPolInt():")
        print("  I could not, via the passed arguments, figure out what you want to plot.")
        print("  The options are line scans or a 2d intensity map. For line scan(s) pass 'eta1'")
        print("  and 'intensity'. 'intensity' and have shape (n,) or (N, n). 'eta1' has shape")
        print("  (n,). For a 2d intensity map pass 'eta1', 'eta2', and 'intensity' with shapes")
        print("  (n,), (m,), and (n,m)." + Fore.BLACK)
        return ax
    # check if the dimensions are ok
    if kind_of_plot == "line_scan" and not len(eta1) == len(intensity):
        print(Fore.RED + "plotPolInt():")
        print("  Arguments 'eta1' and 'intensity' must have the same length.")
        print(f"  len(eta1) = {len(eta1)} and len(intensity) = {len(intensity)}" + Fore.BLACK)
        return ax
    if kind_of_plot == "line_scans" and not len(eta1) == np.shape(intensity)[1]:
        print(Fore.RED + "plotPolInt():")
        print("  The length of 'eta1' does not match the number of points in 'intensity'.")
        print(f"  'eta1' has {len(eta1)} points and the curves in 'intensity' have {np.shape(intensity)[1]} points." + Fore.BLACK)
        return ax
    if kind_of_plot == "intensity_map":
        _mpshp = np.shape(intensity)
        if not (len(eta1) == _mpshp[0] and len(eta2) == _mpshp[1]):
            print(Fore.RED + "plotPolInt():")
            print("  The lengths of 'eta1' and 'eta2' does not match the size of 'intensity'.")
            print(f"  'eta1' and 'eta2' have lengths {len(eta1)} and {len(eta2)} but the corresponding numbers")
            print(f"  for 'intensity' are {_mpshp[0]} and {_mpshp[1]}.")
            return ax
    # simplify
    if kind_of_plot == "line_scan":
        intensity = np.array([intensity])
        kind_of_plot = 'line_scans'
    # special cases?
    if kind_of_plot == "line_scans" and np.shape(intensity)[0] == 4:
        if kwargs.get('special_cases', True):
            kind_of_plot = "special_cases"
    # ------------

    if not hasattr(ax, 'plot'):
        if kind_of_plot == 'intensity_map': figsize = (3,3)
        else: figsize = (5,3)
        figsize = kwargs.get('figsize', figsize)
        fig, ax = plt.subplots(figsize = figsize)
    
    def _getticks():
        """Helper"""
        xticks, yticks, ticks = abs(kwargs.get('xticks', 0)), abs(kwargs.get('yticks', 0)), abs(kwargs.get('ticks', 0))
        if ticks > 0: return ticks, ticks
        else: return xticks, yticks
    

    if kind_of_plot in ["line_scans", "special_cases"]:
        linewidth = kwargs.get('linewidth', 0.75)
        scatter = kwargs.get('scatter', False)
        marker = kwargs.get('marker', 'x')
        s = kwargs.get('s', 10)
        for i, intensity_curve in enumerate(intensity):
            if kind_of_plot == "special_cases": label = f"case {i+1}"
            else: label = f"scan {i+1}"
            if scatter:
                ax.scatter(x = eta1, y = intensity_curve, marker = marker, s = s, label = label)
            else:
                ax.plot(eta1, intensity_curve, linewidth = linewidth, label = label)
        ax.set_xlabel( kwargs.get('xlabel', "$\eta_1$") )
        ax.set_ylabel( kwargs.get('ylabel', 'intensity') )
        title = kwargs.get('title', 'None')
        if not title == 'None': ax.set_title(title)
        if len(intensity) == 1 and kwargs.get('legend', False): ax.legend()
        else:
            if kwargs.get('legend', True): ax.legend()
        xticks, yticks = _getticks()
        if xticks > 0: ax.xaxis.set_major_locator(ticker.MultipleLocator(xticks))
        ax.set_xlim(kwargs.get('xlim', ax.get_xlim()))
        ax.set_ylim(kwargs.get('ylim', ax.get_ylim()))
        plt.tight_layout()
        return ax
    

    elif kind_of_plot == "intensity_map":
        im = ax.imshow(intensity, 
                       extent = [eta2[0], eta2[-1], eta1[-1], eta2[0]], 
                       aspect = kwargs.get('aspect', "equal"), 
                       cmap = kwargs.get('cmap', 'viridis'), 
                       vmin = kwargs.get('vmin', None), 
                       vmax = kwargs.get('vmax', None))
        ax.invert_yaxis()
        ax.set_xlim(kwargs.get('xlim', ax.get_xlim()))
        ax.set_ylim(kwargs.get('ylim', ax.get_ylim()))
        xticks, yticks = _getticks()
        if xticks > 0: ax.xaxis.set_major_locator(ticker.MultipleLocator(xticks))
        if yticks > 0: ax.yaxis.set_major_locator(ticker.MultipleLocator(yticks))
        ax.set_xlabel( kwargs.get('xlabel', "$\eta_2$") )
        ax.set_ylabel( kwargs.get('ylabel', "$\eta_1$") )
        if kwargs.get('colorbar', False): plt.colorbar(im, ax = ax)
        title = kwargs.get('title', 'None')
        if not title == 'None': ax.set_title(title)
        plt.tight_layout()
        return ax
    

    




    


    
        

    
    







    
    






        
