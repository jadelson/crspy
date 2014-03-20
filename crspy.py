"""
crspy.py
collection of python functions

(c) 2104 christopher.robertson.sherwood@gmail.com

"""
from pylab import *
from numpy.lib.scimath import *

def fcat( fname, s, echo=True ):
    """
    fcat( fname, s [, echo=True/False] )
    Append s to file fname
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
    r  = sqrt( x**2 + y**2 )
    az=degrees( arctan2(x, y) )
    # az[where(az<0.)[0]] += 360.
    az = (az+360.)%360.
    return r, az

def xycoord(r, az):
    """
    Convert r, az [degrees, geographic convention] to rectangular coordinates
    x,y = xycoord(r, az)
    """
    x = r * sin(radians(az))
    y = r * cos(radians(az))
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
    y = sqrt(x) * (x<1.) + x *(x>=1.)
    # is this faster than a loop?
    t = tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = tanh( y )
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
    c = g*w*tanh(kh)
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
    hs=(uasq*(0.283)*tanh(fact1)*tanh((fact2/tanh(fact1))))/g

    fact3=0.833*(((g*depth)/uasq)**0.375)
    fact4=0.0379*(((g*fetch)/uasq)**0.333)
    t=(ua*7.54*tanh(fact3)*tanh((fact4/tanh(fact3))))/g
    return hs, t

def earth_dist(lon1,lat1,lon2,lat2):
    """
    Calculate great circle distance between two lon/lat points
    d = earth_dist(lon1,lat1,lon2,lat2)
    Uses Haversine formula
    """
    radius = 6371.*1000.
    dlat = radians(lat2-lat1)
    dlon = radians(lon2-lon1)
    a = sin(dlat/2.)**2 + cos(radians(lat1))*cos(radians(lat2))*sin(dlon/2.)**2.
    d = 2.*radius*math.asin( sqrt(a) )
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
    MAXIT = 20
    vk = 0.41
    rmu=zeros((MAXIT,1))
    Cmu=zeros((MAXIT,1))
    fwci=zeros((MAXIT,1))
    dwci=zeros((MAXIT,1))
    ustrwm2=zeros((MAXIT,1))
    ustrr2=zeros((MAXIT,1))
    ustrci=zeros((MAXIT,1))

    # ...junk return values
    ustrc = 99.99
    ustrwm = 99.99
    ustrr = 99.99
    fwc = .4
    zoa = kN/30.
    zoa = zoa
    dwc = kN

    # ...some data checks
    if( wr <= 0. ):
        print 'WARNING: Bad value for frequency in M94: wr={0}\n'.format(wr)
	return ustrc, ustrr, ustrwm, dwc, fwc, zoa
	
    if( ubr < 0. ):
        print 'WARNING: Bad value for orbital vel. in M94: ub={0}\n'.format(ubr)
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa

    if( kN < 0. ):
	print 'WARNING: Weird value for roughness in M94: kN={0}\n'.format(kN)
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa
	
    if( (zr<zoa or zr<0.05) and iverbose == True):
	print 'WARNING: Low value for ref. level in M94: zr={0}\n'.format(zr)	

    zo = kN/30.
    if(ubr <= 0.01):
        if(ucr <= 0.01):
            # ...no waves or currents
            ustrc = 0.
            ustrwm = 0.
            ustrr = 0.
            return ustrc, ustrr, ustrwm, dwc, fwc, zoa
        # ...no waves
        ustrc = ucr * vk / log(zr/zo) 
        ustrwm = 0.
        ustrr = ustrc
        return ustrc, ustrr, ustrwm, dwc, fwc, zoa
  
    cosphiwc =  abs(cos(phiwc))
    rmu[0] = 0.
    Cmu[0] = 1.
    cukw = Cmu[0]*ubr/(kN*wr)
    print Cmu[0], cukw
    fwci[0] = fwc94( Cmu[0], cukw ) #Eqn. 32 or 33
    ustrwm2[0]= 0.5*fwci[0]*ubr*ubr                 #Eqn. 29
    ustrr2[0] = Cmu[0]*ustrwm2[0]                   #Eqn. 26
    ustrr = sqrt( ustrr2[0] )
    dwci[0] = kN
    if (cukw >= 8.):
        dwci[0]= 2.*vk*ustrr/wr
    lnzr = log(zr/dwci[0])
    lndw = log(dwci[0]/zo)
    lnln = lnzr/lndw
    bigsqr = (-1.+sqrt(1+ ((4.*vk*lndw)/(lnzr*lnzr))*ucr/ustrr ))
    ustrci[0] = 0.5*ustrr*lnln*bigsqr
    nit = 1

    for i in range(1,MAXIT):      
        rmu[i] = ustrci[i-1]*ustrci[i-1]/ustrwm2[i-1]
        Cmu[i] = sqrt(1.+2.*rmu[i]*cosphiwc+rmu[i]*rmu[i])#Eqn 27
        cukw = Cmu[i]*ubr/(kN*wr)
        fwci[i] = fwc94( Cmu[i], cukw )   #Eqn. 32 or 33
        ustrwm2[i]= 0.5*fwci[i]*ubr*ubr                   #Eqn. 29
        ustrr2[i] = Cmu[i]*ustrwm2[i]                     #Eqn. 26
        ustrr = sqrt( ustrr2[i] )
        dwci[i] = kN
        if ((Cmu[i]*ubr/(kN*wr))>= 8.):
            dwci[i]= 2.*vk*ustrr/wr #Eqn.36
        lnzr = log( zr/dwci[i] )
        lndw = log( dwci[i]/zo )
        lnln = lnzr/lndw
        bigsqr = (-1.+sqrt(1+ ((4.*vk*lndw)/(lnzr*lnzr))*ucr/ustrr ))
        ustrci[i] = 0.5*ustrr*lnln*bigsqr                  #Eqn. 38
        diffw = abs( (fwci[i]-fwci[i-1])/fwci[i] )
        # print i,diffw
        if(diffw < 0.0005):
            break
        nit = nit+1
        ustrwm = sqrt( ustrwm2[nit] )
        ustrc = ustrci[nit]
        ustrr = sqrt( ustrr2[nit] )

    zoa = exp( log(dwci[nit])-(ustrc/ustrr)*log(dwci[nit]/zo) ) #Eqn. 11
    fwc = fwci[nit]
    dwc = dwci[nit]
    if(iverbose==True):
        print "M94 nit=",nit
        for i in range(nit):
            print \
	    'i={0} fwc={1} dwc={2} u*c={3} u*wm={4} u*r={5}'\
	    .format(i,fwci[i],dwci[i],ustrci[i],sqrt(ustrwm2[i]),sqrt(ustrr2[i]))

    return ustrc, ustrr, ustrwm, dwc, fwc, zoa


def fwc94( cmu, cukw ):
    """
    fwc94 - Wave-current friction factor
    fwc = fwc94( cmu, cukw )
    Equations 32 and 33 in Madsen, 1994
 
    csherwood@usgs.gov 4 March 2014
    """
    fwc = 0.00999 #meaningless (small) return value
    if( cukw <= 0. ):
        print 'ERROR: cukw too small in fwc94: {0}\n'.format( cukw )
        return fwc
    
    if( cukw < 0.2 ):
        fwc = exp( 7.02*0.2**(-0.078) - 8.82 )
        print 'WARNING: cukw very small in fwc94: {0}\n'.format( cukw )
    if( (cukw >= 0.2) and (cukw <= 100.) ):
        fwc = cmu*exp( 7.02*cukw**(-0.078)-8.82 )
    elif( (cukw > 100.) and (cukw <= 10000.) ):
        fwc = cmu*exp( 5.61*cukw**(-0.109)-7.30 )
    elif( cukw > 10000.):
        fwc = cmu*exp( 5.61*10000.**(-0.109)-7.30 )
        print 'WARNING: cukw very large in fwc94: {0}\n'.format(cukw)

    return fwc
