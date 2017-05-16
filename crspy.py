"""
crspy.py
collection of python functions

(c) 2104 christopher.robertson.sherwood@gmail.com

"""

import numpy as np
from scipy.optimize import fsolve
   
def crspy_loaded():
    """
    Quick way to verify this module is loaded.
    
    boolean = crspy_loaded()
    """
    print("crspy is accessible.")
    return True

def fcat( fname, s, echo=True ):
    """
    Append s to file fname
    fcat( fname, s [, echo=True/False] )
    
    if echo, also print to screen
    """
    fid = open( fname, 'a' )  
    fid.write(s)
    fid.close()
    if echo :
        print (s)
    return

def pcoord(x, y):
    """
    Convert x, y to polar coordinates r, az (geographic convention)
    r,az = pcoord(x, y)
    """
    r  = np.sqrt( x**2 + y**2 )
    az = np.degrees(np.arctan2(x, y) )
    # az[where(az<0.)[0]] += 360.
    az = (az+360.)%360.
    return r, az

def xycoord(r, az):
    """
    Convert r, az [degrees, geographic convention] to rectangular coordinates
    x,y = xycoord(r, az)
    """
    x = r * np.sin(np.radians(az))
    y = r * np.cos(np.radians(az))
    return x, y

def qkhfs( w, h ):
    """
    Quick iterative calculation of kh in gravity-wave dispersion relationship
    kh = qkhfs(w, h )
    
    Input
        w - angular wave frequency = 2*pi/T where T = wave period [1/s]
        h - water depth [m]
    Returns
        kh - wavenumber * depth [ ]

    Orbital velocities from kh are accurate to 3e-12 !

    RL Soulsby (2006) \"Simplified calculation of wave orbital velocities\"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    """
    g = 9.81
    x = w**2.0 *h/g
    y = np.sqrt(x) * (x<1.) + x *(x>=1.)
    # is this faster than a loop?
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    kh = y
    return kh

def cwave( w, h):
    """
    Calculate wave phase celerity
    c = cwave(w, h)
    """
    g = 9.81
    kh = qkhfs( w, h )
    c = g*w*np.tanh(kh)
    return c

def fetch_limited( depth, fetch, u10 ):
    """
    Calculate fetch-limited wave height and period
    hs, t = fetch_limited( depth, fetch, u10 )

    Input
        depth - water depth [m]
        fetch - fetch [m]
        u10   - wind speed at 10 m [m/s]

    Returns
        hs    - significant wave heigh [m]
        t     - wave period [s]

    US Army Corps of Engineers, 1984, Shore Protection Manual, eqn. 3-39
    """
    ua=0.71*(u10**1.23)
    uasq=ua*ua
    g=9.81

    fact1=0.53*(((g*depth)/uasq)**0.75)
    fact2=0.00565*(((g*fetch)/uasq)**0.50)
    hs=(uasq*(0.283)*np.tanh(fact1)*np.tanh((fact2/np.tanh(fact1))))/g

    fact3=0.833*(((g*depth)/uasq)**0.375)
    fact4=0.0379*(((g*fetch)/uasq)**0.333)
    t=(ua*7.54*np.tanh(fact3)*np.tanh((fact4/np.tanh(fact3))))/g
    return hs, t

def earth_dist(lon1,lat1,lon2,lat2):
    """
    Calculate great circle distance between two lon/lat points
    d = earth_dist(lon1,lat1,lon2,lat2)
    Uses Haversine formula
    """
    radius = 6371.*1000.
    dlat = np.radians(lat2-lat1)
    dlon = np.radians(lon2-lon1)
    a = np.sin(dlat/2.)**2 + np.cos(np.radians(lat1))*np.cos(np.radians(lat2))*np.sin(dlon/2.)**2.
    d = 2.*radius*np.arcsin( np.sqrt(a) )
    return d

def rouseh(z,za,ws,ustr,h):
    """
    Rouse profile decreasing to zero at h
    Cz = rouseh(z, za, ws, ustr, h)
    """ 
    vk = 0.41 # von Karman's constant
    Cz = (((h-z)/z) * za/(h-za))**(-ws /(vk*ustr))
    return Cz
    
def rouse(z,za,ws,ustr):
    """
    Rouse profile
    Cz = rouseh(z, za, ws, ustr)
    """ 
    vk = 0.41 # von Karman's constant
    Cz = (za/z)**(-ws /(vk*ustr))
    return Cz
    
    
def m94( ubr, wr, ucr, zr, phiwc, kN, iverbose=False ):
    """
    M94 - Grant-Madsen model from Madsen(1994)
    ustrc, ustrr, ustrm, dwc, fwc, zoa =
        m94( ubr, wr, ucr, zr, phiwc, kN, iverbose )

    Input:
        ubr = rep. wave-orbital velocity amplitude outside wbl [m/s]
        wr = rep. angular wave frequency = 2pi/T [rad/s]
        ucr = current velocity at height zr [m/s]
        zr = reference height for current velocity [m]
        phiwc = angle between currents and waves at zr (radians)
        kN = bottom roughness height (e.q. Nikuradse k) [m]
        iverbose = True/False; when True, extra output
    Returned as tuple
        ustrc  = current friction velocity         u*c [m/s]
        ustrr  = w-c combined friction velocity    u*r [m/s]
        ustrwm = wave max. friction velocity      u*wm [m/s]
        dwc = wave boundary layer thickness [m]
        fwc = wave friction factor [ ]
        zoa = apparent bottom roughness [m]
        
    Chris Sherwood, USGS
    November 2005: Removed when waves == 0
    July 2005: Removed bug found by JCW and RPS
    March 2014: Ported from Matlab to Python
    """

    vk = 0.41

    # ...junk return values
    ustrc = np.nan
    ustrwm = np.nan
    ustrr = np.nan
    fwc = np.nan
    zoa = kN/30.
    zoa = zoa
    dwc = kN

    # ...some data checks
    if( wr <= 0. ):
        print('WARNING: Bad value for frequency in M94: wr={0}\n'.format(wr))
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa
	
    if( ubr < 0. ):
        print('WARNING: Bad value for orbital vel. in M94: ub={0}\n'.format(ubr))
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa

    if( kN < 0. ):
        print('WARNING: Weird value for roughness in M94: kN={0}\n'.format(kN))
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa
	
    if( (zr<zoa or zr<0.05) and iverbose == True):
        print('WARNING: Low value for ref. level in M94: zr={0}\n'.format(zr))	

    zo = kN/30.
    if(ubr <= 0.01):
        
        if(ucr <= 0.01):
            # ...no waves or currents
            ustrc = 0.
            ustrwm = 0.
            ustrr = 0.
            return ustrc, ustrr, ustrwm, dwc, fwc, zoa
        # ...no waves
#        print('nowaves')
        ustrc = ucr * vk / np.log(zr/zo) 
        ustrwm = 0.
        ustrr = ustrc
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa
    
    cosphiwc =  np.abs(np.cos(phiwc))
    rmu = 0.
    Cmu = 1.
    cukw = Cmu*ubr/(kN*wr)
#    print(Cmu[0], cukw)
    fwc = fwc94( Cmu, cukw ) #Eqn. 32 or 33
    ustrwm2 = 0.5*fwc*ubr*ubr                 #Eqn. 29
    ustrr2 = Cmu*ustrwm2                   #Eqn. 26
    ustrr = np.sqrt( ustrr2 )
    dwci = kN
    if (cukw >= 8.):
        dwci= 2.*vk*ustrr/wr
    lnzr = np.log(zr/dwc)
    lndw = np.log(dwc/zo)
    lnln = lnzr/lndw
    bigsqr = (-1.+np.sqrt(1+ ((4.*vk*lndw)/(lnzr*lnzr))*ucr/ustrr ))
    ustrc = 0.5*ustrr*lnln*bigsqr

    
    def get_ustrs(_rmu):
        Cmu = np.sqrt(1.+2.*rmu*cosphiwc+rmu*rmu)
        cukw = Cmu*ubr/(kN*wr)
        fwc = fwc94(Cmu, cukw)
        ustrwm2 = 0.5*fwc*ubr*ubr
        ustrwm = np.sqrt(ustrwm2)
        ustrr2 = Cmu*ustrwm2
        ustrr = np.sqrt( ustrr2 )
        dwc = kN
        if ((Cmu*ubr/(kN*wr))>= 8.):
            dwc = 2.*vk*ustrr/wr #Eqn.36
        lnzr = np.log( zr/dwci )    
        lndw = np.log( dwci/zo )
        lnln = lnzr/lndw
        bigsqr = (-1.+np.sqrt(1+ ((4.*vk*lndw)/(lnzr*lnzr))*ucr/ustrr ))
        ustrc = 0.5*ustrr*lnln*bigsqr
        return ustrc, ustrr, ustrwm, dwc, fwc
    
    def solve_gm(_rmu):
        ustrc, ustrr, ustrwm, dwc, fwc = get_ustrs(_rmu)
        return ustrc*ustrc/np.square(ustrwm) - _rmu

    rmu = fsolve(solve_gm, ustrc*ustrc/ustrwm2)[0]
    ustrc, ustrr, ustrwm, dwc, fwc = get_ustrs(rmu)
    
    
    zoa = np.exp( np.log(dwc)-(ustrc/ustrr)*np.log(dwc/zo) ) #Eqn. 11

#    if(iverbose==True):
#        print("M94 nit=",nit)
#        for i in range(nit):
#            print('i={0} fwc={1} dwc={2} u*c={3} u*wm={4} u*r={5}'\
#	    .format(i,fwci[i],dwci[i],ustrci[i],np.sqrt(ustrwm2[i]),np.sqrt(ustrr2[i])))

    return ustrc, ustrr, ustrwm, dwc, fwc, zoa


def fwc94( cmu, cukw ):
    """
    fwc94 - Wave-current friction factor
    fwc = fwc94( cmu, cukw )
    Equations 32 and 33 in Madsen, 1994
 
    csherwood@usgs.gov 4 March 2014
    """

    fwc = np.nan #meaningless (small) return value
    if( cukw <= 0. ):
        print('ERROR: Cmu*ubr/(kN*wr) too small in fwc94: {0}\n'.format( cukw ))
        return fwc
    
    if( cukw < 0.2 ):
        fwc = np.exp( 7.02*0.2**(-0.078) - 8.82 )
        print('WARNING: Cmu*ubr/(kN*wr) very small in fwc94: {0}\n'.format( cukw ))
    if( (cukw >= 0.2) and (cukw <= 100.) ):
        fwc = cmu*np.exp( 7.02*cukw**(-0.078)-8.82 )
    elif( (cukw > 100.) and (cukw <= 10000.) ):
        fwc = cmu*np.exp( 5.61*cukw**(-0.109)-7.30 )
    elif( cukw > 10000.):
        fwc = cmu*np.exp( 5.61*10000.**(-0.109)-7.30 )
        print('WARNING: Cmu*ubr/(kN*wr) very large in fwc94: {0}\n'.format(cukw))

    return fwc
