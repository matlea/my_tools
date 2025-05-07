#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:59:10 2020

@author: matsleandersson

How to:

    Use M = setupmirrors(alpha, height, separation, gap2) to obtain a variable containing the mirror setup.

    Use M = shiftmirrors(x, M) to shift the mirror setup up- or downstream (x>0 downstream).

    Use L = trace(dA, M) to trace a beam through the polarimeter. By default, the focus is on the
            centre of the first mirror (unless shiftmirror() was used) and dA is the total beamsize
            at the refocusing mirror.
    
    Use F = plottrace(...) to plot the result.

    Plus some other stuff...

"""

import numpy as np
import matplotlib.pyplot as plt

# ======================== Four-mirror setup


def intersect(L1,L2):
    """
    NOT A USER FUNCTION.
    Get the intersecting point between two lines, and the angle (deg.) between them.
    """
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
    a = np.rad2deg(a)
    return xi, yi, a

def getangle(Line):
    """
    NOT A USER FUNCTION.
    Get a lines angle.
    """
    a = np.arctan2((Line[3]-Line[1]),(Line[2]-Line[0]))
    a = np.rad2deg(a)
    return a

def shiftmirrors(x,M):
    """
    Shift a mirror setup (positive in beam direction).
    Use on M from M = setupmirrors(...)

    Argument M is an array describing the four mirrors and the diode.
    shiftmirrors(...) returns a similar array.
    """
    ret = []
    for i in range(len(M)):
        m = [M[i][0]+x, M[i][1], M[i][2]+x, M[i][3]]
        ret.append(m)
    return ret




def trace_L(da=0., M=[]):
    """
    NOT A USER FUNCTION.
    """
    L = []
    for i in range(len(M)+1): L.append([])
        
    L[0] = [-1000, 1000*np.sin(da), 0, 0]
    
    for i in range(len(M)):
        L1 = L[i]
        x1,y1, a1 = intersect(M[i],L1)
        L1 = [L1[0], L1[1], x1, y1]
        L[i] = L1
        a2 = a1 + getangle(M[i])
        a2r = np.deg2rad(a2)
        L2 = [L[i][2], L[i][3], L[i][2] + 1000*np.cos(a2r), L[i][3] + 1000*np.sin(a2r)]
        L[i+1] = L2
    del L[-1]
    
    return L


def trace(da = 0., M=[]):
    """
    Trace the center beam and the upper and lower edges of the beam.
    Use on M from M = setupmirrors(...)
    """
    Lu = trace_L(da/2., M)
    L = trace_L(0.00, M)
    Ld = trace_L(-da/2., M)
    return [Lu,L,Ld]




def plottrace(ID=505, L=[], M=[], Scale=False, Title='', fs=(10,4)):
    """
    Plots the mirror setup and the 'beam(s)'.
    Use with L and M from L = trace(...) and M = setupmirrors(...).
    If Scale is True, the figure is zoomed-in and the aspect ratio is 2:1.

    """
    Fig = plt.figure(ID, figsize=fs )
    Ax = plt.subplot(aspect='auto')
    Ax.set_title(Title)
    Ax.set_xlabel('horizontal')
    Ax.set_ylabel("vertical")
    
    def plotlines(l,c):
        for i in range(len(l)):
            Line = l[i]
            Ax.plot([Line[0],Line[2]], [Line[1],Line[3]], marker = '', color=c)
    
    plotlines(L[0],'r')
    plotlines(L[1],'g')
    plotlines(L[2],'b')
    plotlines(M,'k')
    
    x1,x2 = Ax.get_xlim()
    dy = int(abs(x2-x1)*(fs[1]/fs[0])/2)
    Ax.set_ylim(-dy,dy)
    if Scale:
        Ax.set_xlim(M[0][0]-20,M[-1][2]+5)
        x1,x2 = Ax.get_xlim()
        dy = int(abs(x2-x1)*(fs[1]/fs[0])/2)
        Ax.set_ylim(-dy,dy)
    
    return Fig,Ax



def footprints(L=[],PR=False):
    """
    Returns an array of footpront sizes on the mirrors and diode.
    Use on L from L = trace(...)
    """
    LU = L[0]; LL = L[2]
    sh = np.shape(LU)
    FP=[]
    for i in range(sh[0]):
        x1 = LU[i][2]; y1 = LU[i][3]
        x2 = LL[i][2]; y2 = LL[i][3]
        fp = np.sqrt( (x1-x2)**2 + (y1-y2)**2)
        FP.append(fp)
    if PR:
        for i,fp in enumerate(FP):
            print("M{0}:\t{1:5.1f}".format(i,fp))
    return np.array(FP)





def setupmirrors(alpha=15, height=5, separation=45, gap1=0, gap2=15):
    """
    Setup the polarimeter.
    Input arguments:
        alpha:      incidence angle (1st mirror)
        height:     verticasl distance between 1st and 2nd mirror center.
        separation: distance betweeb the last edge on M3 and the first edge on M4.
        gap1:       old variable, not valid any more.
        gap2:       distance between M4 centre and diode.
    """

    # Fix for old variables:
    if not gap1 == 0:
        separation = separation - gap1
        print("(using old variable set but that should be okay...)")
    a = np.deg2rad(alpha)
    ML = 2.*height / np.cos(a) / np.tan(2.*a)
    
    M1x1 = - 0.5*ML*np.cos(a); 
    M1y1 = - 0.5*ML*np.sin(a)
    M1x2 = + 0.5*ML*np.cos(a);
    M1y2 = + 0.5*ML*np.sin(a)
    M1 = [M1x1, M1y1, M1x2, M1y2]
    
    M2x1 = M1x2 - height/np.tan(2.*a)
    M2y1 = height
    M2x2 = M1x2 + height/np.tan(2.*a)
    M2y2 = height
    M2 = [M2x1, M2y1, M2x2, M2y2]
    
    M3x1 = M1x2
    M3y1 = M1y2
    M3x2 = M3x1 + ML*np.cos(a)
    M3y2 = M1y1 
    M3 = [M3x1, M3y1, M3x2, M3y2]
    
    M4x1 = M3x2 + separation + gap1   # Note that gap1 should be zero, but if it's not, see above.
    M4x2 = M4x1 + (M3y1-M3y2)
    M4y1 = M3y2 
    M4y2 = M4y1 + (M3y1-M3y2)
    M4 = [M4x1, M4y1, M4x2, M4y2]
    
    Dx1 = (M4x2+M4x1)/2 -5.
    Dx2 = Dx1 + 10.
    Dy1 = M4y2 + gap2
    Dy2 = Dy1
    D = [Dx1, Dy1, Dx2, Dy2]
    
    return [M1,M2,M3,M4,D]






def mirrordist(M=[]):
    def f(v):return "{:8.4f}".format(v)
    def centre(M):
        x0 = (M[0]+M[2])/2
        y0 = (M[1]+M[3])/2
        return x0,y0
    def dist(Ma,Mb):
        Mxa, Mya = centre(Ma)
        Mxb, Myb = centre(Mb)
        d = np.sqrt( (Mxb-Mxa)**2 + (Myb-Mya)**2)
        return d
    print('M1 - M2' + f(dist(M[0],M[1])))
    print('M2 - M3' + f(dist(M[1],M[2])))
    print('M3 - M4' + f(dist(M[2],M[3])))
    

def mirrorsizes(M=[], pr=False):
    """
    Returns the sizes of the 4 mirrors (as an array).
    Note that this is NOT the actual sizes of the mirrors but the
    mirror holders.
    Use on M from M = setupmirrors(...)
    """
    if pr:
        print('Sizes: (M5=diode)')
    S = []
    for i,m in enumerate(M):
        s = np.sqrt( (m[2]-m[0])**2 + (m[3]-m[1])**2 )
        S.append(s)
        if pr:
            print("M{0}:\t{1:6.2f}".format(i+1,s))
    S = np.array(S)
    return S

def mirrorcentres(M=[], pr=False):
    """
    Returns the centre positions of the 4 mirrors (as an array of two arrays).
    Use on M from M = setupmirrors(...)
    """
    if pr:
        print('Centres: (M5=diode)')
    Px, Py = [], []
    for i,m in enumerate(M):
        x0 = (m[0]+m[2])/2
        y0 = (m[1]+m[3])/2
        Px.append(x0)
        Py.append(y0)
        if pr:
            print("M{0}:\tx = {1:7.2f}\ty = {2:7.2f}".format(i+1, x0, y0))
    P = np.array([np.array(Px), np.array(Py)])
    return P




 














def printsettings(BC,RA,RH,SEP,AG,DG,M):
    # BC  -  Beam converbence
    # RA  -  Retardet inc. angle
    # AA  -  Analyzer inc. angle
    # RH  -  Retarder height
    # AG  -  Analyzer gap
    # DG  -  Diode gap 
    # M   -  Mirror
    
    def frm(v):return "{:7.2f}".format(v)
    print('Convergence          ' + frm(BC) + ' (Define beam size.)' )
    print('Retarder inc. angle  ' + frm(RA) + ' (Incidence angle on the 1st (and 3rd) mirror.)' )
    print('Analyzer inc. angle  ' + frm(45.0) + ' (Incidence angle on the last (4th) mirror.)' )
    print('Retarder height      ' + frm(RH) + ' (Height = vertical distance between 1st and 2nd reflection.)')
    print('Separation           ' + frm(SEP) + ' (The distance needed for two piezos (2x15mm) + something.)')
    print('Analyzer gap         ' + frm(AG) + ' (The distance from the edge of the 3rd mirror to the edge of the 4th mirror.)')
    print('Diode gap            ' + frm(DG) + ' (Distance between the 4th mirror and the diode surface.)')
    print('Total size (x)       ' + frm(M[-2][2] - M[0][0]))
