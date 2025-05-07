#print('mySource module, Mats Leandersson
# edit: December, 2019
# edit: changed to the common 'shup' instead of 'P' in method source (May, 2022)
# edit: added (for convenience) e2wl() and wl2e(). Also added some minor stuff
#       to source() and kValues() (which also was renamed to that from KValues()).

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from colorama import Fore

# ---------------------

epu = 'bloch'
acc = 'r1'

# ---------------------

Plank = 6.62606957E-34
SpeedOfLight = 299792458
ElectronCharge = 1.60217657e-19
pi = np.pi


def source(E=20., shup = True):
    """
    Returns an array [EmmitanceY, nx, dEE, BetaX, ...] where the elements in order are:

    0   EmittanceX
    1   EmittanceY
    2   nx
    3   dEE
    4   BetaX
    5   BetaY 
    6   SigmaX
    7   SigmaY
    8   SigmaXp
    9   SigmaYp
    10  UnduP
    11  UnduL
    12  int(UnduL/ UnduP), 
    13  PhotonSigma
    14  PhotonSigmaP
    15  SizeX
    16  SizeY
    17  DivergenceX
    18  DivergenceY


    """
    
    global epu
    global acc

    # constants - The Undulator
    if epu.lower() == 'bloch':
        UnduL = 2.8754
        UnduP = 0.084
    elif epu.lower() == 'finest':
        UnduL = 2.4752
        UnduP = 0.0952
    else:
        UnduL = np.nan
        UnduP = np.nan
    
    # constants - The Ring
    if acc.lower() == 'r1':
        EmittanceX = 6.00E-09      # m*rad
        EmittanceY = 6.00E-11      # m*rad
        nx=0.05
        dEE = 0.00074
        BetaX = 5.66               # m
        BetaY = 2.85               # m
    else:
        EmittanceX = np.nan; EmittanceY = np.nan 
        nx = np.nan; dEE = np.nan; BetaX = np.nan; BetaY = np.nan
    
    SigmaX = np.sqrt(EmittanceX * BetaX + (nx * dEE)**2)
    SigmaY = np.sqrt(EmittanceY * BetaY)
	
    SigmaXp = EmittanceX / SigmaX
    SigmaYp = EmittanceY / SigmaY
    
    Wavelength = Plank * SpeedOfLight / ElectronCharge / E
    PhotonSigma = np.sqrt( Wavelength * UnduL)/(2 * pi *np.sqrt(2) ) * 2
    PhotonSigmaP = np.sqrt( Wavelength / (2 * UnduL))

    SizeX = np.sqrt(SigmaX**2 + PhotonSigma**2)
    SizeY = np.sqrt(SigmaY**2 + PhotonSigma**2)
    DivergenceX = np.sqrt(SigmaXp**2 + PhotonSigmaP**2)
    DivergenceY = np.sqrt(SigmaYp**2 + PhotonSigmaP**2)
    
    if not shup:
        print( Fore.BLUE + 'Electron beam')
        print( '')
        print( 'Emittance x: {0} m*rad'.format(EmittanceX))
        print( 'Emittance y: {0} m*rad'.format(EmittanceY))
        print()
        print( 'nx:      {0:4.2f}'.format(nx))
        print( 'dEE:     {0:7.5f}'.format(dEE))
        print( 'Beta x:  {0:4.2f}'.format(BetaX))
        print( 'Beta y:  {0:4.2f}'.format(BetaY))
        print()
        print( 'Sigma x:  {0:.6e}'.format(SigmaX))
        print( 'Sigma y:  {0:.6e}'.format(SigmaY))
        print( 'Sigma xp: {0:.6e}'.format(SigmaXp))
        print( 'Sigma yp: {0:.6e}'.format(SigmaYp))
        print()
        print( 'Undulator')
        print()
        print( 'Undulator period: {0} m'.format(UnduP))
        print( 'Undulator length: {0} m'.format(UnduL))
        print( 'Num. of periods:  {0}'.format(int(UnduL/ UnduP)))
        print()
        print( 'Photon beam')
        print()
        print( 'Sigma:  {0:.6e}'.format(PhotonSigma))
        print( 'Sigmap: {0:.6e}'.format(PhotonSigmaP))
        print()
        print( 'Size x: {:.6e} mm'.format(1000.*SizeX) )
        print( 'Size y: {:.6e}'.format(1000.*SizeY) )
        print( 'Div. x: {:.6e} mrad'.format(1000.*DivergenceX) )
        print( 'Div. y: {:.6e}'.format(1000.*DivergenceY) + Fore.RESET)
    
    else:
        ret = [EmittanceX, EmittanceY, nx, dEE, BetaX, BetaY, 
                SigmaX, SigmaY, SigmaXp, SigmaYp, UnduP, UnduL, int(UnduL/ UnduP), 
                PhotonSigma, PhotonSigmaP, SizeX, SizeY, DivergenceX, DivergenceY]
        return np.array(ret)


def kValues(filename='', epuperiod=0.084, bxcutoff=0.5, bycutoff=0.5, shup = False):
    """
    
    filename:   a 3-column text file without header: z (mm), Bx (T), By (T)
    epuperiod:  (m)
    bxcutoff:   minimum intensity of the peaks to consider (as part of maximum)
    bycutoff:   -"-

    """
    e = 1.602e-19 # C
    me = 9.10938356e-31 #kg
    c = 299792458 # m/s

    try:
        data = np.loadtxt(filename)
    except:
        print(Fore.RED + 'mySource.Kvalues(): The data file should contain three columns: z, Bx, By.' + Fore.RESET)
        return
    
    fig = plt.figure(figsize=(14,4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.set_title('Bx'); ax2.set_title('By')
    ax1.set_ylabel('B (T)'); ax2.set_ylabel('B (T)')
    ax1.set_xlabel('z (mm)'); ax2.set_xlabel('z (mm)')
    ax1.plot(data[:,0],data[:,1])
    ax2.plot(data[:,0],data[:,2])
    
    peaks_ix, peaks_intx = find_peaks( np.abs(data[:,1]), np.max(np.abs(data[:,1])*bxcutoff))
    peaks_iy, peaks_inty = find_peaks( np.abs(data[:,2]), np.max(np.abs(data[:,2])*bycutoff))
    
    for ip in peaks_ix:
        ax1.scatter(data[:,0][ip], data[:,1][ip], color='k')
    for ip in peaks_iy:
        ax2.scatter(data[:,0][ip], data[:,2][ip], color='k')
    
    B0x = np.sum(peaks_intx['peak_heights'])/(len(peaks_intx['peak_heights'])+1)
    B0y = np.sum(peaks_inty['peak_heights'])/(len(peaks_inty['peak_heights'])+1)
    
    ax1.plot( [data[:,0][0],data[:,0][-1]], [B0x,B0x], linestyle=':', color='tab:orange' )
    ax1.plot( [data[:,0][0],data[:,0][-1]], [-B0x,-B0x], linestyle=':', color='tab:orange' )
    ax2.plot( [data[:,0][0],data[:,0][-1]], [B0y,B0y], linestyle=':', color='tab:orange' )
    ax2.plot( [data[:,0][0],data[:,0][-1]], [-B0y,-B0y], linestyle=':', color='tab:orange' )
    
    if not shup: print(Fore.BLUE + 'B0x, B0y = {0:5.3f}, {0:5.3f}'.format(B0x,B0y) + Fore.RESET)
    
    Kx = e * B0x * epuperiod / (2 * np.pi * me * c)
    Ky = e * B0y * epuperiod / (2 * np.pi * me * c)
    
    if not shup: print(Fore.BLUE + 'Kx, Ky   = {0:5.3f}, {0:5.3f}'.format(Kx,Ky) + Fore.RESET)

    return {"Kx": Kx, "Ky": Ky}
    

    


def e2wl(energy = 20.):
    """
    Energy (eV) to wavelength (m).
    """
    try: e = abs(float(energy))
    except:
        print(Fore.RED + "mySource.e2wl(): The attribute energy must be a scalar (eV)." + Fore.RESET)
        return np.NaN
    return Plank * SpeedOfLight / ElectronCharge / e

def wl2e(wavelength = 635e-9):
    """
    Wavelength (m) to energy (eV).
    """
    try: wl = abs(float(wavelength))
    except:
        print(Fore.RED + "mySource.wl2e(): The attribute wavelength must be a scalar (m)." + Fore.RESET)
        return np.NaN
    return Plank * SpeedOfLight / ElectronCharge / wl

