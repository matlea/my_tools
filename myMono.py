"""
Compiling main methods from old scripts.

Grating diffraction methods: 
    cpgm_e2a(energy, cff = 2.25)   opt: order = 1, line_density = 800
    cpgm_a2e(mirror, grating)           order = 1, line_density = 800
    nim_e2a(energy)                     order = 1, line_density = 4000, mirror = 15
    nim_a2e(grating)                    order = 1, line_density = 4000, mirror = 15

Wavelength and energy methods:
    wl2e(wavelength)    wavelength (m) to energy (eV)
    e2wl(energy)        energy (eV) to wavelength (m)

Other methods:
    blaze(energy,...)   to visualize the blaze angle's effect.


Versions:
    23.08.31    Moved blaze() to myMonoExtra(). Keeping it here as well but not updating it.
    23.08.26    Added method blaze(). Visualizes the blaze angle.
    23.07.29    As a part of a scripts-wide update: added global variable default_shup.
    23.07.18    Bux-fixes and touch-ups. some issues with diffreaction order for NIM so disabled.
                Added wl2e() and e2wl() to convert between wavelenth and energy.
    23.07.17    Added methods for NIM.
    23.07.16    Copied and re-wrote the methods for cPGM.

"""

__version__ = "23.08.31"
__author__  = "Mats Leandersson"



import numpy as np
from numpy import pi
from colorama import Fore
import matplotlib.pyplot as plt



# ==============================================================================================
print(Fore.LIGHTWHITE_EX + f"--- myMono, version {__version__}" + Fore.BLACK)
# ==============================================================================================





default_shup = False

# =============================================================================================
# ============================================================================================= C P G M
# =============================================================================================


# --------------------------------------------------------------------------------------------- E N E R G Y   T O   A N G L E S

def e2a(energy = None, cff = 2.25, order = 1, line_density = 800, shup = None, return_dict = False):
    """
    Calculates and returns the mirror and grating angles for energy (eV), cff, and line_density (l/mm) for cPGM.
    Arguments:
        energy, number > 0
        cff, number > 1
        order, integer
        line_density, number > 0
    Returns:
        alpha, beta, and theta if return_dict = False, otherwise a dict with all info
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    return cpgm_e2a(energy = energy, cff = cff, order = order, line_density = line_density, shup = shup, return_dict = return_dict)


def cpgm_e2a(energy = None, cff = 2.25, order = 1, line_density = 800, shup = None, return_dict = False):
    """
    Calculates and returns the mirror and grating angles for energy (eV), cff, and line_density (l/mm) for cPGM.
    Arguments:
        energy, number > 0
        cff, number > 1
        order, integer
        line_density, number > 0
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    ld = line_density * 1e3
    #
    Alpha = np.NaN; Beta = np.NaN; Theta = np.NaN
    #
    try:
        energy = float(energy); cff = float(cff); line_density = float(line_density)
    except:
        print(Fore.RED + "Arguments 'energy', 'cff', and 'line_density' must all be numbers (energy > 0, cff > 1, line_density > 0)."); return Alpha, Beta
    if not energy > 0 and cff > 1 and line_density > 0:
        print(Fore.RED + "Arguments 'energy', 'cff', and 'line_density' must all be numbers (energy > 0, cff > 1, line_density > 0)."); return Alpha, Beta
    #
    #
    if not order == 1:
        if not shup: print(Fore.BLUE + "Note: I have not included diffraction order in this method yet.\n" + Fore.BLACK)
    #
    #
    try:
        Lambda = e2wl(energy = energy, shup = True)
        alpha = np.arccos( Lambda * ld / (cff**2 - 1) * ( np.sqrt( cff**2 + (cff**2 -1)**2 / Lambda**2 / ld**2 )-1 ))
        beta = np.arcsin( cff * np.sin(alpha) )
        theta = 0.5  * (alpha + beta)
        Alpha = 90 - np.rad2deg(alpha)
        Beta  = np.rad2deg(beta) - 90
        Theta = np.rad2deg(theta)
    except:
        print(Fore.RED + "Something went wrong. What did you do???" + Fore.BLACK)
        Alpha = np.NaN; Beta = np.NaN; Theta = np.NaN

    if not shup:
        print(Fore.BLUE + f"energy = {energy:.3f} eV, cff = {cff:.3f}, line density = {line_density} l/mm, order = {order}:")
        print(f"\u03B1 = {Alpha:8.4f}°")
        print(f"\u03B2 = {Beta:8.4f}°")
        print(f"\u03B8 = {Theta:8.4f}°")
    #
    if return_dict:
            return {'alpha': Alpha, 'beta': Beta, 'theta': Theta,
                    'mirror': Theta, 'grating': 90+Beta,
                    'energy': energy, 'cff': cff, 'line_density': line_density, 'order': order}
    else: return Alpha, Beta, Theta




# --------------------------------------------------------------------------------------------- A N G L E S   T O   E N E R G Y

def a2e(mirror = None, grating = None, order = 1, line_density = 800, shup = None, return_dict = False):
    """"""
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    return cpgm_a2e(mirror = mirror, grating = grating, order = order, line_density = line_density, shup = shup, return_dict = return_dict)

def cpgm_a2e(mirror = None, grating = None, order = 1, line_density = 800, shup = None, return_dict = False):
    """
    Calculates the energy and cff value for mirror and grating angles (in deg.)
    The angles are defined from the horizon.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    Energy, Cff = np.NaN, np.NaN
    #
    if type(mirror) is type(None) and type(grating) is type(None):
        print(Fore.RED + "Pass arguments 'mirror' and 'grating' as angles in degrees (>0°)." + Fore.BLACK); return Energy, Cff
    try:
        mirror = float(mirror)
        grating = float(grating)
    except:
        print(Fore.RED + "Arguments 'mirror' and 'grating' must be numbers." + Fore.BLACK); return Energy, Cff
    if not order == int(order):
        print(Fore.RED + "Argument 'order' must be an integer." + Fore.BLACK); return Energy, Cff
    if order == 0:
        Energy = np.inf; Cff = 1; return Energy, Cff        # is this true? check it up.
    try:
        line_density = float(line_density)
    except:
        print(Fore.RED + "Argument 'line_density' must be an number > 0 (in l/mm)." + Fore.BLACK); return Energy, Cff
    if line_density <= 0:
        print(Fore.RED + "Argument 'line_density' must be an number > 0 (in l/mm)." + Fore.BLACK); return Energy, Cff
    #
    a = np.deg2rad(2*mirror - grating)
    b = np.deg2rad(grating)
    ld = line_density * 1e3
    #
    if a == b:
        Energy = np.inf; Cff = 1
        if not shup:
            print(Fore.BLUE + f"Zero order reflection. \u03B1 = \u03B2 = {np.rad2deg(a):.4f}°" + Fore.BLACK)
        return Energy, Cff
    elif np.sin(a) == 0:
        Energy = np.NaN; Cff == np.inf
        if not shup:
            print(Fore.BLUE + f"Diffraction impossible. \u03B1 = {np.rad2deg(a):.4f}° and \u03B2 = {np.rad2deg(b):.4f}°" + Fore.BLACK)
        return Energy, Cff
    elif (np.cos(a) - np.cos(b)) == 0:    # this is the same as the first condition, right? why did I do this back then? check it up.
        Energy = np.inf; Cff = 1
        if not shup:
            print(Fore.BLUE + f"Zero order reflection. \u03B1 = \u03B2 = {np.rad2deg(a):.4f}°" + Fore.BLACK)
        return Energy, Cff
    else:
        try:
            Lambda = (1/ld*order) * (np.cos(a) - np.cos(b))
            Energy = wl2e(wavelength = Lambda, shup = True)
            Cff = np.sin(b) / np.sin(a)
        except:
            Energy, Cff = np.NaN, np.NaN
            print(Fore.RED + f"Can not calculate the energy and cff for \u03B1 = \u03B2 = {np.rad2deg(a):.4f}°." + Fore.BLACK); return Energy, Cff
    #
    if Energy <= 0:
        print(Fore.RED + f"Negative energy ({Energy})." + Fore.BLACK)
        Energy = np.NaN
    if Cff < 1:
        print(Fore.RED + f"cff is below 1 ({Cff})." + Fore.BLACK)
    #
    Alpha = 90 - np.rad2deg(a)
    Beta = np.rad2deg(b) - 90
    Theta = np.rad2deg((a+b)/2)
    if not shup:
        print(Fore.BLUE + f"mirror  = {mirror:8.4f}°, grating = {grating:8.4f}°, line density = {line_density} l/mm:")
        #print(f"(\u03B1 = {Alpha:.4f}°, \u03B2 = {Beta:.4f}°, \u03B8 = {Theta:.4f}°)")
        print(f"Energy = {Energy:6.2f} eV")
        print(f"cff    = {Cff:6.2f}" + Fore.BLACK)
    #
    if return_dict:
        return {'energy': Energy, 'cff': Cff,
                'mirror': mirror, 'grating': grating,
                'line_density': line_density, 'order': order,
                'alpha': Alpha, 'beta': Beta, 'theta': Theta}
    else:
        return Energy, Cff





# =============================================================================================
# ============================================================================================= N I M
# =============================================================================================



# --------------------------------------------------------------------------------------------- E N E R G Y   T O   A N G L E S

def nim_e2a(energy = None, order = 1, line_density = 4000, mirror = 15, shup = None, return_dict = False):
    """
    Calculates and returns the mirror and grating angles for energy (eV), cff, and line_density (l/mm) for NIM.
    Arguments:
        energy, number > 0, eV
        order, integer, default 1
        line_density, number > 0, l/mm, default 4000
        mirror, number > 0, deg, default 15
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    h = 4.135667e-15
    c = 299792458
    ld = line_density * 1e3
    #
    Alpha = np.NaN; Beta = np.NaN; Theta = np.NaN
    #
    try:
        energy = float(energy); order = float(order); line_density = float(line_density); mirror = float(mirror)
    except:
        print(Fore.RED + "Arguments 'energy', 'order', 'line_density', and 'mirror' must all be numbers (energy > 0, order != 0, line_density > 0, mirror > 0)." + Fore.BLACK); return Beta
    if not energy > 0 and line_density > 0 and mirror > 0:
        print(Fore.RED + "Arguments 'energy', 'line_density', and 'mirror' must all be > 0." + Fore.BLACK); return Beta
    #
    if not order == int(order):
        print(Fore.RED + "Argument 'order' must be an integer." + Fore.BLACK); return Beta
    if order == 0:
        if not shup: print(Fore.CYAN + "Can not calculate the energy (order = 0)." + Fore.BLACK); return np.NaN
    #
    if not order == 1:
        print(Fore.MAGENTA + "Note: orders != 1 is not yet accounted for. Setting order = 1.\n" + Fore.BLACK)
        order = 1 
    #
    Theta = mirror
    theta = np.deg2rad(Theta)
    beta = np.arccos( order * ld * h * c / np.sqrt(2.) / np.sqrt(1. + np.cos(2. * theta)) / energy ) + theta
    beta = beta - np.pi/2
    alpha = 2. * theta - beta
    Alpha = np.rad2deg(alpha)
    Beta = np.rad2deg(beta)
    #
    #a = np.sin(2*theta)
    #b = -np.cos(2*theta)
    #k = np.sign(a) * np.sqrt(a**2 + b**2)
    #psi = np.arctan(-b/a)
    #A = k * np.cos(psi)
    #B = 1 - k * np.sin(psi)
    #K = np.sign(A) * np.sqrt(A**2 + B**2)
    #PSI = np.arctan(-B/A)
    #beta = np.arccos(order * h * c * ld / K / energy) - PSI
    #Beta = np.rad2deg(beta)
    #
    if not shup:
        print(Fore.BLUE + f"energy = {energy:.2f} eV, order = {order}, line density = {line_density:.1f} l/mm, mirror = {mirror}°:")
        print(f"\u03B1  = {Alpha:8.4f}°")
        print(f"\u03B2  = {Beta:8.4f}°")
        print(f"(grating angle = {mirror - Beta:.4f}°)" + Fore.BLACK)
    #
    if return_dict:
        return {'alpha': Alpha, 'beta': Beta,
                'grating': mirror - Beta, 'mirror': mirror,
                'energy': energy, 'line_density': line_density, 'order': order}
    else:
        return Alpha, Beta




# --------------------------------------------------------------------------------------------- A N G L E S   T O   E N E R G Y

def nim_a2e(grating = None, order = 1, line_density = 4000, mirror = 15, shup = None, return_dict = False):
    """
    Grating = 0° means beta = alpha = mirror (default 15°)

    Arguments:
        grating number  deg.    bla bla, non-optinal
        order   integer         diffraction order, default 1
        mirror  number  deg.    the incidence angle of the mirror, measured from the normal, default 15°
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    Energy = np.NaN
    #
    if type(grating) is type(None):
        print(Fore.RED + "The argument 'grating' must be a number (deg.)."); return Energy
    try:
        grating = float(grating)
    except:
        print(Fore.RED + "The argument 'grating' must be a number (deg.)."); return Energy
    if grating == 0:
        print(Fore.RED + "Zero order reflection."); return np.inf
    try:
        order = float(order)
    except:
        print(Fore.RED + "The argument 'order' must be an integer."); return Energy
    if not int(order) == order:
        print(Fore.RED + "The argument 'order' must be an integer."); return Energy
    order = int(order)
    if order == 0:
        if not shup: print(Fore.CYAN + "DIffraction order: 0" + Fore.BLACK)
        return Energy
    #
    if not order == 1:
        print(Fore.MAGENTA + "Note: orders != 1 is not yet accounted for. Setting order = 1.\n" + Fore.BLACK)
        order = 1
    #
    try:
        line_density = float(line_density)
    except:
        print(Fore.RED + "The argument 'line_density' must be a number > 0 (l/mm)."); return Energy
    ld = line_density * 1e3
    if not line_density > 0:
        print(Fore.RED + "The argument 'line_density' must be a number > 0 (l/mm)."); return Energy
    try:
        mirror = float(mirror)
    except:
        print(Fore.RED + "The argument 'mirror' must be a number > 0 (deg.)."); return Energy
    if not mirror > 0:
        print(Fore.RED + "The argument 'mirror' must be a number > 0 (deg.)."); return Energy
    #
    beta = mirror - grating
    alpha = mirror + grating
    Lambda = 1/order * 1/ld * (np.sin(np.deg2rad(alpha)) + np.sin(np.deg2rad(-beta)))
    Energy = wl2e(wavelength = Lambda, shup = True)
    if not shup:
        print(Fore.BLUE + f"mirror = {mirror:.4f}°, grating = {grating:.4f}°, order = {order}, line density = {line_density:.1f} l/mm:")
        #print(f"\u03B1 = {alpha:6.4f}°")
        #print(f"\u03B2 = {beta:6.4f}°")
        print(f"E = {Energy:.3f} eV" + Fore.BLACK)
    #
    if return_dict:
        return {'energy': Energy,
                'grating': grating, 'line_density': line_density, 'mirror': mirror, 'order': order,
                'alpha': alpha, 'beta': beta}
    else:
        return Energy






# ============================================================================================= wavelength
# =============================================================================================    <-> 
# =============================================================================================  energy


def e2wl(energy = None, shup = None):
    """
    Energy-to-wavelength converter. Argument 'energy' in eV, return in m.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    try: energy = abs(float(energy))
    except:
        print(Fore.RED + "Pass argument 'energy' as a number (eV, >0)." + Fore.BLACK)
        return np.NaN
    if energy == 0:
        print(Fore.BLUE + "Zero energy means infinite wavelength." + Fore.BLACK)
        return np.inf
    h = 4.135667e-15; c = 299792458
    wavelength = h * c / energy
    if not shup: print(Fore.BLUE + f"E = {energy} eV: \u03BB = {wavelength*1e9} nm" + Fore.BLACK)
    return wavelength


def wl2e(wavelength = None, shup = None):
    """
    Wavelength-to-energy converter. Argument 'wavelength' in m, return in eV.
    """
    if type(shup) is type(None):
        global default_shup
        shup = default_shup
    #
    try: wavelength = abs(float(wavelength))
    except:
        print(Fore.RED + "Pass argument 'wavelength' as a number." + Fore.BLACK)
        return np.NaN
    h = 4.135667e-15; c = 299792458
    energy = h * c / wavelength
    if not shup: print(Fore.BLUE + f"\u03BB = {wavelength*1e9} nm: E = {energy:.2f} eV" + Fore.BLACK)
    return energy
    





# ============================================================================================= 
# =============================================================================================    blaze 
# =============================================================================================


def blaze(energy = 100, cff = 2.25, blaze_angle = 2.17, line_density = 800, shup = False, plot = True, title = 1):
    """
    Draw the diffraction and blaze condition. Returns an pyplot axis and a dict with values, 
    or just a dict if argument plot is False.
    """
    print(Fore.BLUE + ("Note: this method is moved to myMonoExtra(), where any updates will or might be implemented.)\n")+ Fore.BLACK)
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
    mono_config = cpgm_e2a(energy = energy, cff = cff, order = 1, line_density = line_density, shup = True, return_dict = True)
    alpha = mono_config['alpha']
    beta = mono_config['beta']
    alpha_r, beta_r = np.deg2rad(alpha), np.deg2rad(beta)
    #
    alpha_r, beta_r, blaze_angle_r = np.deg2rad(alpha), np.deg2rad(beta), np.deg2rad(blaze_angle)
    d = 1/line_density # line spacing
    h = d * np.tan(blaze_angle_r) # groove height
    beta_reflected = 2 * blaze_angle - alpha
    beta_reflected_r = np.deg2rad(beta_reflected)
    #
    q = {'energy': energy, 'cff': cff, 'line_density': line_density, 'blaze_angle': blaze_angle, 'alpha': alpha, 'beta': beta}
    q.update({'beta_on_blaze': beta_reflected, 'mirror': mono_config['mirror'], 'grating': mono_config['grating']})
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
        print(f"mirror:         {mono_config['mirror']:9.5f}°")
        print(f"grating:        {mono_config['grating']:9.5f}°")
        print(f"reflected beam: {q['reflected_beam']:9.5f}° (from diffracted)")
    #
    if plot: return ax, q
    else: return q

    
