from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from ray_csv.read import *
from stokes_mueller.funcs import mirror_matrix as mirrormatrix
from stokes_mueller.funcs import normalize_stokes as norm_stokes
from stokes_mueller.funcs import rotation_matrix as rotmat

DATA = []
DATA_incidences = []
DATA_energies = []


def setIncidence(inc=15):
    # The incidence angle is defined as the angle between the surface
    # and the ray inpinging on the first mirror. This function then
    # returns the corresponding incidence angles for the rest of the
    # mirrors as POLINC.
    # The second returned array, iPOLINC, contain the indeces of those
    # angles in the by Ray generated data file.
    global DATA_incidences
    POLINC = np.array([inc,2.*inc,inc,45.])
    iPOLINC = []
    for i in range(len(POLINC)):
        iPOLINC.append( np.abs(POLINC[i]-DATA_incidences).argmin() )
    return POLINC, iPOLINC


def polarimeter(Energy=20, R1=0, R2=0, iPOLINC=[0,0,0,0]):
    # This function returns the polarimeter matrix for the retarder
    # being at angle R1 and the analyzer being at angle R2, for
    # a given energy.
    global DATA_energies
    global DATA
    iE = np.abs(DATA_energies - Energy).argmin()
    T = [ 1., DATA[1][iPOLINC[0]][iE], DATA[2][iPOLINC[0]][iE], DATA[3][iPOLINC[0]][iE] ]
    M1 = mirrormatrix(T)
    T = [ 1., DATA[1][iPOLINC[1]][iE], DATA[2][iPOLINC[1]][iE], DATA[3][iPOLINC[1]][iE] ]
    M2 = mirrormatrix(T)
    M3 = np.copy(M1)
    RET = M3.dot(M2.dot(M1))
    T = [ 1., DATA[1][iPOLINC[3]][iE], DATA[1][iPOLINC[3]][iE], DATA[3][iPOLINC[3]][iE] ]
    ANA = mirrormatrix(T)
    Rp = rotmat(R1)
    Rn = rotmat(-R1)
    rRET = Rn.dot(RET.dot(Rp))
    Rp = rotmat(R2)
    Rn = rotmat(-R2)
    rANA = Rn.dot(ANA.dot(Rp))
    POL = rANA.dot(rRET)
    return POL

def runPolarimeter(Energy=20., Sin=[1., 1., 0., 0.], iPOLINC=[15,30,15,45]):
    NU = np.linspace(0.,360.,361)
    I1 = np.zeros(len(NU)); I2 = np.zeros(len(NU))
    I3 = np.zeros(len(NU)); I4 = np.zeros(len(NU))
    for i in range(len(NU)):
        m = polarimeter(Energy=Energy, R1=NU[i], R2=NU[i], iPOLINC = iPOLINC)
        I1[i] = m.dot(Sin)[0]
        m = polarimeter(Energy=Energy, R1=NU[i], R2=NU[i]+90, iPOLINC = iPOLINC)
        I2[i] = m.dot(Sin)[0]
        m = polarimeter(Energy=Energy, R1=NU[i], R2=NU[i]+45, iPOLINC = iPOLINC)
        I3[i] = m.dot(Sin)[0]
        m = polarimeter(Energy=Energy, R1=NU[i], R2=NU[i]-45, iPOLINC = iPOLINC)
        I4[i] = m.dot(Sin)[0]
    return NU, [I1,I2,I3,I4]




# ===============================

def fitfunc(x,I,V,T):
    x=np.deg2rad(x)
    return I*(1+V*np.sin(2*x+T))

def dofit(x,y,Par,Bnd):
    fp, covar = curve_fit(f=fitfunc, xdata=x, ydata=y, p0=Par) #, bounds=Bnd)
    return fp, covar

def fittedcurves(NU,P):
    I=[]
    for i in range(4):
        y=np.zeros(len(NU))
        for j in range(len(NU)):
            y[j] = fitfunc(NU[j], P[i][0], P[i][1], P[i][2])
        I.append(y)
    I = np.array(I)
    return I

# =========================================

def getImean(I):
    return [I[0].mean(), I[1].mean(), I[2].mean(), I[3].mean()]

def getVisibility(I):
    V1 = (I[0].max() - I[0].min()) / (2.*I[0].mean())
    V2 = (I[1].max() - I[1].min()) / (2.*I[1].mean())
    V3 = (I[2].max() - I[2].min()) / (2.*I[2].mean())
    V4 = (I[3].max() - I[3].min()) / (2.*I[3].mean())
    return [V1,V2,V3,V4]


def extractParameters(NU=[], I=(), PRINT=False, inPSI=[0.,0.,0.,0.], ALARM=False):
    RES = [ np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan ]
    
    # Mean intensities, and visibilities
    Imean = getImean(I)
    V = getVisibility(I)
    if PRINT:
        print('Mean intensity:    ', np.round(Imean,4))
        print('Visibility:        ', np.round(V,4))
    
    Vcheck = [False,False,False,False]
    vislim = 0.0001
    for i in range(len(Vcheck)):
        if np.abs(V[i])>0.0001:
            Vcheck[i] = True
        else:   
            if PRINT: print ('* Curve I' + str(i+1) + ' visibility = 0, or below ' +str(vislim))
    
    # Get phases (PSI) from fittings
    inPSI = np.deg2rad(inPSI)
    cov = []
    FPAR = []
    PSI = np.zeros((len(Imean)))
    for i in range(len(Imean)):
        if Vcheck[i]:
            PAR = [Imean[i]        , V[i]        ,  inPSI[i] ]
            LB =  [0.00*I[i].max() , 0.00        , -2*np.pi  ]
            UB =  [2.00*I[i].max() , 2.00        ,  2*np.pi  ]
            
            fPAR, covar = dofit(NU, I[i], PAR, (LB,UB)) 
            if np.max(covar)>1:
                if ALARM: print('* Bad fit for I' + str(i) + '.')
            PSI[i] = fPAR[2]
        else:
            if i==0: 
                PSI[0] = np.pi/2
                if PRINT: print('* Setting PSI1 to ' + str(PSI[0]))
            if i==1: 
                PSI[1] = -np.pi/2         
                if PRINT: print('* Setting PSI2 to ' + str(PSI[1]))
            if i>1:
                PSI[i] = np.pi/2
                if PRINT: print('* PSI' + str(i) + ' unknown. (look into this...)')
        FPAR.append([Imean[i], V[i], PSI[i]])
        
    if PRINT:
        print('Phases (PSI):      ', np.round(PSI,4))
        print('                   ', np.round(np.rad2deg(PSI),2), 'deg.')
    
    
    RES[0] = PSI[0]; RES[1] = PSI[1]; RES[2] = PSI[2]; RES[3] = PSI[3]
    
    # Linear polarization
    tmp1 = (V[0]**2 * Imean[0]**2  -  V[1]**2 * Imean[1]**2)
    tmp2 = (Imean[0]**2 - Imean[1]**2)
    if tmp1>0 and tmp2>0:
        PL =  np.sqrt( tmp1 / tmp2 )
    else:
        PL = 0.
    if np.isnan(PL):
        PL = 0.
    if PRINT: print('PL:                ', np.round(PL,4))
    
    RES[4] = PL
    
    # Stokes 1 and 2 
    S1 = PL * np.sin(PSI[0])
    S2 = PL * np.cos(PSI[0])
    S3 = 1.   # for coding-technical purposes
    if PRINT: 
        print( 'S1:                ', np.round(S1,4) )
        print( 'S2:                ', np.round(S2,4) )

    # cos_2psi3 > cos_2psi4 ?
    if not Vcheck[0] or not Vcheck[1]:
        sign = 1.
    else:
        I1_0 = I[0][ abs(NU-0.).argmin() ]
        I1_90 = I[0][ abs(NU-90.).argmin() ]
        I2_0 = I[1][ abs(NU-0.).argmin() ]
        I2_90 = I[1][ abs(NU-90.).argmin() ]
        Q = (I2_90 - I2_0) / (I1_90 - I1_0)
        sign = np.sign(Q)
    if np.isnan(sign):
        sign = 1.
    if PRINT: 
        print('(Sign:             ', sign, ')')
    
    # psi3 and psi4
    if not PL==0:
        cos_2psi3 = (1/PL) * (V[0] * Imean[0]  + sign *  V[1] * Imean[1]) / (Imean[0] + Imean[1])
        cos_2psi4 = (1/PL) * (V[0] * Imean[0]  - sign *  V[1] * Imean[1]) / (Imean[0] + Imean[1])
        psi3 = np.arccos(cos_2psi3)/2
        psi4 = np.arccos(cos_2psi4)/2
        if PRINT: 
            print( 'psi3:              ', np.round(psi3,4) )
            print( '                   ', np.round(np.rad2deg(psi3),2), 'deg.' )
            print( 'psi4:              ', np.round(psi4,4) ) 
            print( '                   ', np.round(np.rad2deg(psi4),2), 'deg.' )
    else:
        S3 = 1
        psi3, psi4 = np.nan,np.nan
        if ALARM: 
            print('* S1 and S0 are zero. Can not detrmine psi3, psi4, and delta3. Setting S3=1.')
            print('  (TO DO: determine if S3 is +1 or -1.)')
            print( 'psi3:               ?')
            print( 'psi4:               ?') 
    
    RES[5] = psi3; RES[6] = psi4
    
    delta3 = np.nan
    if not (np.isnan(psi3) or np.isnan(psi4)):
        # delta3
        cos23 = np.cos(2*psi3)
        cos24 = np.cos(2*psi4)
        sin23 = np.sin(2*psi3)
        cosd3 = cos23/(cos24*sin23) * (-S1 + S2*np.tan(PSI[2]))/(S2 + S1*np.tan(PSI[2]))
        delta3 = np.arccos(cosd3)
        if PRINT: 
            print( 'delta3:            ', np.round(delta3,4) )
            print( '                   ', np.round(np.rad2deg(delta3),4), 'deg' )
    
        # S3
        sind3 = np.sin(delta3)
        S3 = 1./(cos24*sin23*sind3)*((Imean[3]-Imean[2])/(Imean[3]+Imean[2]))
    else:
        S3=1.
        delta3 = np.nan
        if ALARM: print( 'delta3:             ?')
        
    RES[7] = delta3

    if PRINT: 
        print( 'S3:                ', np.round(S3,4) )
        
    # Normalized Stokes
    s0 = np.sqrt( S1**2 + S2**2 + S3**2)
    SN = [1., S1/s0, S2/s0, S3/s0]
    if PRINT: 
        print( 'Normalized Stokes: ', np.round(SN,4) )
    
    RES = [ PSI[0], PSI[1], PSI[2], PSI[3], PL, psi3, psi4, delta3, SN ]
    return RES


def goodness(psi3=0, psi4=0, delta3=0):
    # i.e. the form factor
    return np.abs( np.sin(2.*psi3) * np.cos(2.*psi4) * np.sin(delta3))