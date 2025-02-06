#!/usr/bin/env python
#
# Typical System: /usr/bin/env python
# LBT MODS runtime: /lbt/mods_runtime/anaconda/bin/python
#
# modsView - View a MODS target acquisition (.acq) or imaging (.img) script
# 
# Reads the contents of a MODS target acquisition (.acq) or the target
# acquisition blocks of a MODS imaging (.img) script, and displays the
# results on a ds9 window, overlaying the MODS focal plane on a
# digitized sky survey image of the target field.  It also finds and
# overlays catalog stars and their R magnitudes.
#
# It tests to see if the guide star is inside the guide patrol field
# given the instrument rotator angle.  If there are any telescope
# offsets as part of the target acquisition, it will test the location
# of the guide star after the offset to verify that it is still inside
# the MODS guide patrol field.
#
# Uses the Python ds9 module to interact with a named DS9 window.  The
# window is launched as needed, turning off all of the IRAF imtool
# pipes so a person running IRAF will not step on this display (and
# vis-versa).
#
# Creates a temporary modsView.reg file containing the DS9 region file
# used to draw the MODS focal plane on the screen.
#
# It can also be used to create a PNG-format finder chart by using the
# --finder option.
#
# Dependencies:
#   Requires your system have SAOImage DS9 installed and running
#   https://sites.google.com/cfa.harvard.edu/saoimageds9
#
#   Starting with version 3 we are using the SAMP messaging protocol
#   to interact with ds9 as pyds9 is no longer supported.  We are adopting
#   the astropy.samp implementatino.
#
# Distribution:
#   The primary distribution is now on GitHub
#     github.com/rwpogge/modsView
#
# Author:
#   R. Pogge, OSU Astronomy Department
#   pogge.1@osu.edu
#   2012 May 4
#
# Modification History:
#   2012 May 04 - first beta release, does only guide-star checking
#                 and *very* minimal syntax checks. [rwp/osu]
#   2012 May 06 - first experiments with interactive mode, added a
#                 number of functions [rwp/osu]
#   2012 May 20 - Added probe shadow option (--shadow) and plot OBJNAME 
#                 label or suppress with --nolabel [rwp/osu]
#   2012 Oct 01 - Experimental interactive version, call it 1.0.0...
#   2012 Nov 08 - Added ability to use B magnitudes in catalogs instead
#                 of just R, and introduced the --find option to give an
#                 R_mag-sorted list of candidate guide stars to pick from
#                 interactively (type number at prompt) [rwp/osu]
#   2012 Nov 28 - added NOMAD1 catalog support as default [rwp/osu]
#   2012 Dec 03 - restored --rotate function [rwp/osu]
#   2012 Dec 18 - fixed bug in display of the instrument rotator sweep
#                 and long-slit positions after offsets [rwp/osu]
#   2012 Dec 19 - minor patch, exception catching for python 2.7
#                 (raises exception of you try to close a catalog pane that
#                 is already closed - python 2.6 doesn't care) [rwp/osu]
#   2013 Jan 7  - minor bug if nonsensical min/max magnitude ranges set
#                 [rwp/osu]
#   2013 Mar 19 - bugs in the offseting logic [rwp/osu]
#   2013 Apr 11 - Allow for no + sign on Decs in MMS files (not strictly
#                 permitted, but apparent can happen), and fixed a
#                 previously undetected bug in guide probe shaddow
#                 rendering under RA/Dec offsets [rwp/osu]
#   2014 Feb 23 - adjusted guide patrol field size for changes in 
#                 maximum Y stage travel limits with the WFS hotspot
#                 offset.  Found conflicting limits in different parts of
#                 the program resulting in the patrol field being drawn
#                 too large, and tests failing to alert a star just outside
#                 the effective guide patrol field. [rwp/osu]
#   ================================================================
#   2014 Apr 28 - Start of binocular MODS hooks, including support for the
#                 different MODS1 and MODS2 AGw parameters [rwp/osu]
#   2015 Nov 20 - First v2.0 release for MODS1 or MODS2 use [rwp/osu]
#   2016 Oct 11 - MODS1 and 2 have the same AGw WFS configuration as of
#                 October 2016, so old MODS1 offset hotspot has been removed,
#                 and --mods1/mods2 is no longer required [rwp/osu]
#   2016 Oct 15 - Minor mods for pyds9 vs ds9 back compatibility [rwp/osu]
#   2018 May 22 - Support for experimental SNS masks [rwp/osu]
#   2018 Jul 22 - Patches for Python 3 & MacOS operation, first release
#                 using GitHub [rwp/osu]
#   2018 Sep 05 - fixed input/raw_input problem P2/3 issue [rwp/osu]
#   2019 Nov 24 - Updated AGw patrol field coordinates [rwp/osu]
#   2022 Nov 11 - Updated for changes in XPA with ds9 version 8.x [rwp/osu]
#   ================================================================
#   2025 Feb 05 - start of v3 SAMP development [rwp/osu]
#
#---------------------------------------------------------------------------

import sys
import os
import math
import getopt
import subprocess
import shlex
import time

from pathlib import Path

from operator import itemgetter

# input vs raw_input for Python 3/2 compatibility

try:
    input = raw_input
except NameError:
    pass

# Version number and date, update as needed

versNum  = '3.0.0'
versDate = '2025-02-05'

# Some useful global defaults (mostly so we can report them in usage)

lbtScale = 0.600 # LBT focal plane scale in mm/arcsec
minRMag = 16.5   # Typical guide star R magnitude limits for the MODS AGw unit
maxRMag = 12.0
minBMag = 16.5   # Typical guide star B magnitude limits for the MODS AGw unit
maxBMag = 11.0
fsFoV    = 5.5   # Radius of the focal station field-of-view in arcminutes
defCatalog = 'nomad' # default is the NOMAD1 Catalog, other option is ub1=USNO-B1
defAGwFilt = 'R'     # default AGw guide camera filter is 'R'
defServer  = 'stsci' # default is the STScI image server
defSurvey  = 'all'   # default is the STScI composite survey image catalog

#---------------------------------------------------------------------------
#
# internal functions and classes
#

from astropy.samp import SAMPIntegratedClient, SAMPHubError
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

# DS9 class

class DS9():
    '''
    DS9 interaction class to replace using pyds9 with SAMP interaction
    
    '''
    
    def __init__(self,ds9ID):
        '''
        Initialize a DS9 class instance

        Parameters
        ----------
        ds9ID : string
            name of the ds9 instance we will create.

        Raises
        ------
        ValueError
            raised if bad inputs are provided.
        RuntimeError
            raised if unrecoverable errors occur..

        Returns
        -------
        None.

        '''
        
        if len(ds9ID) == 0:
            raise ValueError("required argument ds9ID missing")
        
        self.ds9ID = ds9ID
        
        # instantiate a SAMPIntegratedClient instance
        
        self.ds9 = SAMPIntegratedClient(name=f"myDS9_{ds9ID}",description="myDS9 instance")
        self.connected = self.ds9.is_connected
        
        self.clientID = None
        self.haveDS9 = False
        
        try:
            self.ds9.connect()
            self.connected = self.ds9.is_connected
            self.clientID = self.getID()
            if self.clientID:
                self.haveDS9 = True
            else:
                subprocess.Popen(shlex.split(f"ds9 -samp -fifo none -port none -unix none -title {self.ds9ID}"))
                time.sleep(5)
                self.clientID = self.getID()
                if not self.clientID:
                    raise RuntimeError(f"Failed to start ds9 for {self.ds9ID}, aborting")
                else:
                    self.connected = self.ds9.is_connected
                    
        except SAMPHubError:
            subprocess.Popen(shlex.split(f"ds9 -samp -fifo none -port none -unix none -title {self.ds9ID}"))
            time.sleep(5)
            try:
                self.ds9.connect()
                self.clientID = self.getID()
                if not self.clientID:
                    raise RuntimeError(f"Failed to start ds9 for {self.ds9ID}, aborting")
                else:
                    self.haveDS9 = True
                    self.connected = self.ds9.is_connected
            except Exception as exp:
                raise RuntimeError(f"Cannot connect to a samp hub or ds9: {exp}")
            
        except Exception as exp:
            raise RuntimeError(f"DS9 startup failed with exception: {exp}")


    def getID(self):
        '''
        Get the samp hub client ID corresponding to the ds9 instance we need

        Returns
        -------
        clientID : string
            SAMP client ID of the ds9 instance of interest

        Description
        -----------
        Uses the astropy.samp get_registered_client() and get_methdata()
        methods to find the named ds9 window of interest.  The name
        attached to a ds9 process with the -title command-line argument
        is in the samp.name metadata parameter.  We require an exact
        case match.
        
        '''
        if not self.ds9.is_connected:
            return None
        for clientID in self.ds9.get_registered_clients():
            cliMD = self.ds9.get_metadata(clientID)
            if "samp.name" in cliMD:
                if cliMD["samp.name"] == self.ds9ID:
                    self.clientID = clientID
                    return self.clientID
        return None
    
    
    def set(self,cmdStr):
        '''
        Send a commmand to the ds9 instance

        Parameters
        ----------
        cmdStr : string
            SAOImage ds9 set command to send

        Raises
        ------
        ValueError
            if bad information is provided
        RuntimeError
            if there are unrecoverable runtime errors

        Returns
        -------
        None.

        '''
        
        if len(cmdStr) == 0:
            return        
        
        if self.haveDS9:
            try:
                self.ds9.ecall_and_wait(self.clientID,"ds9.set","10",cmd=cmdStr)
            except Exception as exp:
                raise ValueError(f"ds9 set command {cmdStr} returned error: {exp}")
        else:
            raise RuntimeError(f"named ds9 instance {self.ds9ID} not connected")

    
    def get(self,cmdStr):
        '''
        Execute a ds9 get directive to the ds9 instance

        Parameters
        ----------
        cmdStr : string
            SAOImage ds9 get command to send.

        Raises
        ------
        RuntimeError
            if there are unrecoverable errors.

        Returns
        -------
        dictionary
            SAMP dictionary with the return from the ds9.get command.

        See Also
        --------
        getCursKey() for a function that uses this for cursor interaction
        
        '''
        
        if len(cmdStr) == 0:
            return None
        
        if self.haveDS9:
            try:
                sampRet = self.ds9.ecall_and_wait(self.clientID,"ds9.get","0",cmd=cmdStr)
                if "samp.result" in sampRet:
                    if "value" in sampRet["samp.result"]:
                        return sampRet["samp.result"]["value"]
                    else:
                        return sampRet["samp.result"]
                else:
                    return sampRet
            except Exception as exp:
                raise RuntimeError(f"ds9 get command {cmdStr} returned error: {exp}")
        else:
            raise RuntimeError(f"named ds9 instances {self.ds9ID} not connected")

            
    def getCursKey(self,prompt=None,coords="image"):
        '''
        Put an interactive cursor on the image and wait for a key press,
        returning coordinates and the key pressed.

        Parameters
        ----------
        prompt : string, optional
            Prompt to show on stdout at the start. The default is None.
        coords : string, optional
            coordinates to return, must be one of "image", "fk5", "wcs fk5", or 
            "wcs galactic", The default is "image".

        Raises
        ------
        RuntimeError
            if unrecoverable errors occur in execution.

        Returns
        -------
        key : string
            The key hit (also space if spacebar, etc.)
        cursX : float
            the X coordinate of the cursor when key was hit
        cursY : float
            the Y coordinate of the cursor when key was hit

        Description
        -----------
        A wrapper for "iexam key coordinate image", or substitute
        one of (fk5,wcs fk5,wsc galactic) if alternatives given.
         * "image" returns pixel coordinates (default)
         * "fk5" returns ra/dec in decimal hours/degrees
        Returns the key hit, does not respond to mouse click.
        
        '''
        
        if not prompt:
            print(f"Put cursor on the {self.ds9ID} display image and hit any key")
        else:
            print(f"{prompt}")
        
        try:
            cursData = self.get(f"iexam key coordinate {coords}")
        except Exception as exp:
            raise RuntimeError(f"iexam returned error: {exp}")
        
        cursBits = cursData.split(' ')
        if len(cursBits) == 3:
            return cursBits[0],float(cursBits[1]),float(cursBits[2])
        else:
            return cursData
    

#------------------
#
# replacements for coordinate conversions using astropy
#

def rdToStd(ra,dec,ra0,dec0):
    '''
    compute offset from (ra0,dec0) to (ra,dec) in arcsec

    Parameters
    ----------
    ra : float
        destination right ascension in decimal hours
    dec : float
        destination declination in decimal degrees
    ra0 : float
        origin right ascension in decimal hours
    dec0 : float
        origin declination in decimal degrees

    Returns
    -------
    xi : float
        RA offset from ra0 to ra in arcseconds on the standard plane
    eta : float
        Dec offset from dec0 to dec in arcseconds on the standard plane.

    See Also
    --------
    `stdToRD` for the reverse offset
    
    '''
    rd0 = SkyCoord(ra0*u.hour,dec0*u.deg,frame="icrs")
    rd = SkyCoord(ra*u.hour,dec*u.deg,frame="icrs")
    xi, eta = rd0.spherical_offsets_to(rd)
    return xi.to(u.arcsec).value,eta.to(u.arcsec).value


def stdToRD(xi,eta,ra0,dec0):
    '''
    convert standard coordinates (xi,eta) on the standard plane to RA/Dec

    Parameters
    ----------
    xi : float
        RA offset from ra0 on the standard plane in decimal arcseconds.
    eta : float
        Dec offset from dec0 on the standard plane in decimal arcseconds.
    ra0 : float
        Right ascension of the reference (tangent) point in decimal hours
    dec0 : float
        Declination of the reference (tangent) point in decimal degrees

    Returns
    -------
    ra
        Right ascension of the offset position in decimal hours.
    dec
        Declination of the offset position in decimal degrees.

    See Also
    --------
    `rdToStd` for the reverse offset
    
    '''
    rd0 = SkyCoord(ra0*u.hour,dec0*u.deg,frame="icrs")
    newRD = rd0.spherical_offsets_by(xi*u.arcsec,eta*u.arcsec)
    return newRD.ra.value/15.0, newRD.dec.value


def dec2sex(decAng,precision=2,sign=False):
    '''
    Convert decimal angle into a sexigesmial string

    Parameters
    ----------
    decAng : float
        angle, either degrees or hours.
    precision : int, optional
        precision of the seconds part. The default is 2.
    sign : boolean, optional
        include + sign if decAng>0. The default is False (no + sign).

    Returns
    -------
    string
        Sexigesimal string representation of decAng with a colon (:)
        separator.
            
    See Also
    --------
    sex2dec

    '''    
    return Angle(decAng,unit=u.deg).to_string(sep=":",precision=precision,pad=True,alwayssign=sign)


def sex2dec(sexStr,precision=8):
    '''
    Convert a sexigesimal string to decimal

    Parameters
    ----------
    sexStr : string
        sexigesimal representation of an angle (see Notes)
    precision : integer, optional
        precision of the conversion. The default is 8 figures beyond the decimal point.

    Returns
    -------
    float
        decimal representation of the sexigesimal string
            
    Notes
    -----
    If the string is in 12h13m18.23s format, the "h" will make sure the
    conversion preserves units of hours.  Correctly recognizes strings
    that don't use the colon (:) separator, so 12:13:14.15 = 12d13m14.15s and 12:13:14.15 = 12h13m14.5s preserve
    correct hours or degrees units on the conversion to decimal
                
    See Also
    --------
    dec2sex

    '''
    if "h" in sexStr:
        units = u.hour
    else:
        units = u.deg
    return float(Angle(sexStr,unit=units).to_string(decimal=True,precision=precision))

#---------------------------------------------------------------------------
#
# inTriangle() - Test to see if (x,y) is inside a triangle
#
# Inputs:
#   x1,y1,x2,y2,x3,y3 = Cartesian coordinates of the triangle vertices
#   x,y = Cartesian coordinates of the test point.
#
# Returns:
#   True if (x,y) is inside the triangle, False if outside.
#
# Description:
#   Solves this classic geometry problem by transforming the
#   coordinates of the triangle vertices and the test point into the
#   barycentric coordinate system of the triangle.
#
#   A very lucid description of the problem is in the Wikipedia
#   article on the barycentric coordinate system:
#
#     en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 2
#

def inTriangle(x1,y1,x2,y2,x3,y3,x,y):
    '''
    Test to determine if (x,y) are inside a triangle

    Parameters
    ----------
    x1 : float
        x coordinate of triangle vertex 1.
    y1 : float
        y coordinate of triangle vertex 1.
    x2 : float
        x coordinate of triangle vertex 2.
    y2 : float
        y coordinate of triangle vertex 2.
    x3 : float
        x coordinate of triangle vertex 3.
    y3 : float
        y coordinate of triangle vertex 3.
    x : float
        x coordinate to test.
    y : float
        y coordinate to test.

    Returns
    -------
    bool
        True if (x,y) are inside the triangle, False if they are outside

    Description
    -----------
    Solves this classic geometry problem by transforming the
    coordinates of the triangle vertices and the test point into the
    barycentric coordinate system of the triangle.

    A very lucid description of the problem is in the Wikipedia
    article on the barycentric coordinate system:
      
    en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
    
    '''

    n1 = (y2-y3)*(x -x3) + (x3-x2)*(y -y3)
    d1 = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
    n2 = (y3-y1)*(x -x3) + (x1-x3)*(y -y3)
    d2 = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)

    if d1 != 0.0:
        lam1 = n1/d1
    else:
        lam1 = 0.0

    if d2 !=0.0:
        lam2 = n2/d2
    else:
        lam2 = 0.0

    lam3 = 1 - lam1 - lam2

    if (lam1<0.0 or lam2<0.0 or lam3<0.0):
        return False
    else:
        return True


def inBox(x,y,rect,rotAng):
    '''
    Test to determine if (x,y) are inside a rectangular box including rotation

    Parameters
    ----------
    x : float
        x coordinate to test.
    y : float
        y coordinate to test.
    rect : float list
        rectange with center (xc,yc) and size (dx,dy).
    rotAng : float
        rotation angle for the box about (0,0).

    Returns
    -------
    bool
        True if (x,y) are in the box, False if they are outside.

    Description
    -----------
    Creats a box from the center (xc,yc) and size (dx,dy) but allows
    rotation of the box about the (0,0) origin of the (xc,yc) coordinates.
    It then divides the box diagonally into 2 triangles and uses the
    `inTriangle()` function to make the test.
    
    See Also
    --------
    inTriangle() and rotXY()
    
    '''
    
    (xc,yc,dx,dy) = rect
    
    # Coordinates of the center and rectangle vertices after rotation
    
    (xr0,yr0)=rotXY(xc,yc,-posAng)
    (xr1,yr1)=rotXY(xc-dx/2,yc-dy/2,-rotAng)
    (xr2,yr2)=rotXY(xc+dx/2,yc-dy/2,-rotAng)
    (xr3,yr3)=rotXY(xc+dx/2,yc+dy/2,-rotAng)
    (xr4,yr4)=rotXY(xc-dx/2,yc+dy/2,-rotAng)

    # Slice the box diagonally into two triangles and test if
    # (x,y) are in one or other triangle
    
    if inTriangle(xr1,yr1,xr2,yr2,xr3,yr3,x,y) or inTriangle(xr1,yr1,xr3,yr3,xr4,yr4,x,y):
        return True
    else:
        return False

    
#---------------------------------------------------------------------------
#
# rotXY - rotate XY coordinates
#
# Inputs:
#     x,y = float positions relative to origin
#  rotAng = rotation angle in degrees
#
# Returns:
#    xr,yr = (x,y) in the rotated frame
#
# Description:
#   Convenience function to evaluate the standard Cartesian 2D
#   coordinate system rotation.  Note that for applying this to
#   astronomical standard coordinates (xi,eta), the helicity of rotAng
#   has the opposite sign (e.g., compute xi,eta when rotating by
#   celestial position angle posAng, use rotAng=-posAng).
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 2
#

def rotXY(x,y,rotAng):
    '''
    rotate (x,y) relative to the origin (0,0)

    Parameters
    ----------
    x : float
        x position relative to x=0.
    y : float
        y position relative to y=0
    rotAng : float
        rotation angle in degrees about (0,0).

    Returns
    -------
    xr : float
        rotated x coordinate.
    yr : float
        rotate y coordinate.
        
    Description
    -----------
    Convenience function to evaluate the standard Cartesian 2D
    coordinate system rotation.  Note that for applying this to
    astronomical standard coordinates (xi,eta), the helicity of rotAng
    has the opposite sign of celestial position angle (e.g., to compute 
    xi,eta when rotating by position angle posAng, use rotAng=-posAng).    

    '''
    
    if rotAng==0.0 or rotAng==-0.0:
        return x, y
    else:
        sinPA = math.sin(math.radians(rotAng))
        cosPA = math.cos(math.radians(rotAng))
        xr = x*cosPA - y*sinPA
        yr = x*sinPA + y*cosPA
        return xr, yr


def parseMMS(file):
    '''
    Parse a MODS Mask Specification (MMS) file

    Parameters
    ----------
    file : string
        name of the MMS file to open and parse.

    Returns
    -------
    ra : float
        Right ascencsion in decimal hours.
    dec : float
        Declination in decimal degrees.
    swid : float
        slit width in arcseconds.
    slen : float
        slit lengthin arcseconds.
    srot : float
        slit rotation angle in decimal degrees.

    Description
    -----------
    Opens and reads in the contents of a MODS mask design (MMS) file,
    and extracts the coordinates and dimensions of the slits,
    returning them as useful floating-point representations.  The
    trick lies in reading the MMS file's coding for slit RA/Dec
    coordinates:

        RA format: INS.TARGnnn.ALPHA 203448.232
        Dec format: INS.TARGnnn.DELTA -001415.123
    
    '''

    ra = []
    dec = []
    swid = []
    slen = []
    srot = []
    F = open(file, 'r')
    M = F.readlines()[::]
    F.close()
    numRef = 0 
    for i in range(len(M)):
        inStr = M[i].strip()
        if not inStr.startswith('#') and len(inStr)>0: # ignore comments and blank lines
            (param,datum) = inStr.split()
            (slitType,targID,field) = param.split('.')
            if slitType == 'INS':
                if targID=='RSLIT' and field=='NUMBER':
                    numRef = int(datum)
                    firstTarget = 100+numRef+1  # first non-reference slit ID number
                elif targID[0:4] == 'TARG':
                    targNum = int(targID[4:])
                    if targNum >= firstTarget:
                        if field == 'ALPHA':
                            raStr = datum
                            sexStr = '%s:%s:%s' % (raStr[0:2],raStr[2:4],raStr[4:]) 
                            ra.append(sex2dec(sexStr))
                        elif field == 'DELTA':
                            decStr = datum
                            if decStr[0] == '-' or decStr[0]=='+':
                                sexStr = '%s:%s:%s' % (decStr[0:3],decStr[3:5],decStr[5:])
                            else:
                                sexStr = '%s:%s:%s' % (decStr[0:2],decStr[2:4],decStr[4:])
                            dec.append(sex2dec(sexStr))
                        elif field == 'WID':
                            swid.append(float(datum))
                        elif field == 'LEN':
                            slen.append(float(datum))
                        elif field == 'ROT':
                            srot.append(float(datum))

    return ra, dec, swid, slen, srot


def findStar(ra,dec,catRAd,catDec,radius):
    '''
    find a star in a star catalog given RA/Dec coordinates within search radius

    Parameters
    ----------
    ra : float
        right ascension in decimal degrees.
    dec : float
        declination in decimal degrees.
    catRAd : float list
        star catalog RA list.
    catDec : float list
        star catalog Dec list.
    radius : float
        search radius for matches in arcseconds

    Returns
    -------
    index : integer
        array index of the star in the catalog closest to (ra,dec).

    '''
    if len(catRAd)==0:
        return -1
    for i in range(len(catRAd)):
        dist = math.sqrt(math.pow(ra-catRAd[i],2) + math.pow(dec-catDec[i],2))
        if i==0:
            dmin = dist
            iFound = 0
        else:
            if dist<dmin:
                iFound = i
                dmin = dist
        
    darcs = 3600.0*dmin
    if darcs < radius: 
        return iFound
    else:
        return -1


def loadCat(catFile):
    '''
    load a star catalog into working arrays

    Parameters
    ----------
    catFile : string
        Name of the star catalog to load.

    Returns
    -------
    numStars : integer
        number of stars in the catalog.
    catID : string
        catalog identifier (e.g., USNO-B1.0).
    catRAd : float list
        catalog star RAs in decial degrees.
    catDec : float list
        catalog star declinations in decimal degrees.
    catBmag : float list
        catalog star B magnitudes.
    catRmag : float list
        catalog star R magnitudes.
    catName : string
        catalog star names.

    Description
    -----------
    Knows how to read star catalog downloaded from a catalog
    server by the ds9 instance.
    
    '''
    SC=open(catFile,'r')
    catLines = SC.readlines()[::]
    SC.close()
    catRAd  = []
    catDec  = []
    catRmag = []
    catBmag = []
    catName = []
    catID = "None"
    numStars = len(catLines)-1
    for i in range(len(catLines)):
        catBits = catLines[i].split('\t') # split on tabs
        if i==0: # first line is the header, get the element count
            numItems = len(catBits)
            catID = catBits[2]
        else:
            catRAd.append(float(catBits[0]))
            catDec.append(float(catBits[1]))
            catName.append(catBits[2])
            if catID.upper() == 'NOMAD1':
                try:
                    catRmag.append(float(catBits[15]))
                except ValueError:
                    catRmag.append(99.99)
                try:
                    catBmag.append(float(catBits[11]))
                except ValueError:
                    catBmag.append(99.99)
            else:
                catRmag.append(float(catBits[numItems-2]))
                try:
                    catBmag.append(float(catBits[numItems-3]))
                except ValueError:
                    catBmag.append(99.99)  # placeholder if B is blank in catalog

    return numStars,catID,catRAd,catDec,catBmag,catRmag,catName

#---------------------------------------------------------------------------
#
# drawMODS - Create the MODS instrument overlay as a DS9 region file
#
# Inputs:
#   objName  - [string] Object Name
#   target   - [tuple] target (RA,Dec) in decimal hours and degrees
#   posAng   - [float] mask celestial position angle in decimal degrees
#   gstar    - [tuple] guide star (RA,Dec) in decimal h/deg or None
#   slitMask - [string] slit mask ID
#   offRD    - [tuple] RADEC offset (dRA,dDec) in decimal arcseconds
#   offXY    - [tuple] DETXY offset (dX,dY) in decimal arcseconds
#   mmsFile  - [string] name of a MODS Mask Specification (mms) file
#   gprobe   - [bool] show the guide probe shadow
#   boxSize  - [int] size of the image display box in arcminutes 
#   

def drawMODS(objName,target,posAng,gstar,slitMask,offRD,offXY,mmsFile,gprobe,boxSize):
    '''
    Create the MODS instrument overlay as a DS9 regions file

    Parameters
    ----------
    objName : string
        object name.
    target : float list
        object RA/Dec in decimal hours/degrees.
    posAng : float
        mask celestial position angle in decimal degrees, +=North to East.
    gstar : float list
        guide star RA/Dec in decimal hours/degrees.
    slitMask : string
        slit mask ID.
    offRD : float list
        offset position RA/Dec in decimal arcseconds (RADEC).
    offXY : float lost
        offset position in detector XY decimal arcseconds (DETXY).
    mmsFile : string
        filename of a MODS Mask Specification (MMS) file
    gprobe : boolean
        if True include the approximate guide probe shadow outline.
    boxSize : int
        Size of the image display box in decimal arcminutes.

    Returns
    -------
    None.

    '''

    # MODS science field at PA=0 (center and width/height in arcsec).

    (sciX,sciY) = (0,0)
    (sciW,sciH) = (360,360)

    # MODS guide patrol field at PA=0 (center and width/height in arcsec)

    (gpfX,gpfY) = (0,-150.0)   
    (gpfW,gpfH) = (290,300)   
    (xr0,yr0) = rotXY(gpfX,gpfY,-posAng)
    gsBox = (gpfX,gpfY,gpfW,gpfH)

    # Nominal guide probe shadow region (center & dimensions) relative to
    # the guide star position in arcsec in (xi,eta) coordinates.

    (gpsX,gpsY) = (35,0)
    (gpsW,gpsH) = (150,85)

    # Guide probe carrier arm shadow
        
    (armX,armY) = (76,-105)
    (armW,armH) = (80,125)

    # Guide Field of View position & dimensions relative to the guide star.
    # Both MODS have a 40x40-arcsec FoV

    (gcX,gcY) = (0.042,0.404)   # really a slight offset (0.042,0.404)
    (gcW,gcH) = (40,40)

    # Centers of the facility long-slit mask segments in arcsec
  
    slitCen = [-126,-63,0,63,126]

    # Process Arguments

    # Target RA/Dec

    (targRA,targDec) = target
    targRAd = 15.0*targRA

    # Guide Star?  If not, turn off guide star plotting features

    if gstar==None:
        hasGS = False
        gprobe = False
    else:
        hasGS = True
        (gsRA,gsDec) = gstar
        gsRAd = 15.0*gsRA
        (gsXi,gsEta) = rdToStd(gsRA,gsDec,targRA,targDec)
        if inBox(gsXi,gsEta,gsBox,posAng):
            gsValid = True
        else:
            gsValid = False
            print('\n  ** WARNING: The guide star is OUTSIDE the guide patrol field.')

    # RA/Dec Offset

    (offsetRA,offsetDec)=offRD
    dRAd = (offsetRA/math.cos(math.radians(targDec)))/3600.0
    dRA  = dRAd/15.0
    dDec = offsetDec/3600.0

    # DETXY Offset

    (offsetX,offsetY)=offXY
    (dX,dY) = rotXY(offsetX,offsetY,-posAng)

    # Mask ID for slit overlay?

    if slitMask==None:
        showSlit = False
    else:
        if slitMask.upper().startswith('LS'):
            showSlit = True
        else:
            showSlit = False

    # MMS File for MOS mask overlay?

    if mmsFile==None:
        showMMS = False
    else:
        showMMS = True

    #
    # Build the regions file
    #
    # We have to create an external regions file to draw our boxes on
    # the image.  Why?  Because if you send individual regions
    # commands directly (e.g., regions command {box ...})  it reverts
    # to physical coordinates not WCS, and screws up, even if you
    # explicitly include 'fk5;' in the command string (proper syntax
    # ala the manual).  A blind-spot of ds9 and XPA...
    #

    regFile = str(Path.cwd() / "modsView.reg")
    if os.path.isfile(regFile):
        os.remove(regFile)
    RF = open(regFile,'w')

    RF.write('#\n# modsView regions file\n#\nfk5\n')

    # Draw the 'aim point', the original preset coordinates before any
    # offsets are applied, as a yellow circle

    RF.write(f'circle {targRAd}d {targDec}d 3.0\" # color=yellow width=2\n')

    # RA/Dec of the instrument aim point after all offsets are applied

    (instRA,instDec) = stdToRD(sciX+dX,sciY+dY,targRA+dRA,targDec+dDec)
    instRAd = 15.0*instRA

    # The circle shows the full sweep of the MODS science and guide patrol fields

    RF.write(f'circle {instRAd}d {instDec}d {fsFoV-0.0333}\' # color=red\n')
    RF.write(f'circle {instRAd}d {instDec}d {fsFoV}\' # color=cyan\n')
    RF.write(f'circle {instRAd}d {instDec}d {fsFoV+0.0333}\' # color=red\n')
    
    # RA/Dec of the science field center after all offsets are applied

    (sciRA,sciDec) = stdToRD(sciX+dX,sciY+dY,targRA+dRA,targDec+dDec)
    sciRAd = 15.0*sciRA
    
    # Draw the MODS science field

    RF.write(f'box {sciRAd}d {sciDec}d 6\' 6\' {posAng} # width=2 color=green\n')
    
    # RA/Dec of the guide patrol field after all offsets are applied

    (gpfRA,gpfDec) = stdToRD(xr0+dX,yr0+dY,targRA+dRA,targDec+dDec)
    gpfRAd = 15.0*gpfRA

    # Draw the MODS AGw Guide Patrol Field

    RF.write(f'box {gpfRAd}d {gpfDec}d {gpfW/60.0}\' {gpfH/60.0}\' {posAng} # color=cyan width=2\n')
    
    # Show the long-slit locations if using a facility mask

    if showSlit:
        if slitMask.upper() == 'LS60X5':
            slitWid=5.0
            slitLen=60.0
            RF.write(f'box {sciRAd}d {sciDec}d {slitWid/60.0}\' {slitLen/60.0}\' {posAng}\n')
            
        elif slitMask.upper() == 'LS10X0.8SNS':
            slitWid = 0.8
            slitLen = 10.0
            ys = -120.0
            (xsr,ysr) = rotXY(0,ys,-posAng)
            (rs,ds) = stdToRD(xsr,ysr,sciRA,sciDec)
            rsd = 15.0*rs
            RF.write(f'box {rsd}d {ds}d {slitWid/60.0}\' {slitLen/60.0}\' {posAng}\n')
            
        else:
            bits = slitMask.upper().split('X')
            slitWid = float(bits[len(bits)-1])
            slitLen=60.0
            for ys in slitCen:
                (xsr,ysr) = rotXY(0,ys,-posAng)
                (rs,ds) = stdToRD(xsr,ysr,sciRA,sciDec)
                rsd = 15.0*rs
                RF.write(f'box {rsd}d {ds}d {slitWid/60.0}\' {slitLen/60.0}\' {posAng}\n')

   # If using a MOS mask and showMMS, display the MOS mask slitlets

    if showMMS:
        (mmsRA,mmsDec,mmsWid,mmsLen,mmsRot) = parseMMS(mmsFile)
        for i in range(len(mmsRA)): 
            if mmsWid[i]==mmsLen[i]: # assume square = alignment box, color magenta
                RF.write(f'box {15.0*mmsRA[i]}d {mmsDec[i]}d {mmsWid[i]:.1f}\" {mmsLen[i]:.1f}\" {posAng+mmsRot[i]:.2f} # color=magenta width=2\n')
            else:
                RF.write(f'box {15.0*mmsRA[i]}d {mmsDec[i]}d {mmsWid[i]:.1f}\" {mmsLen[i]:.1f}\" {posAng+mmsRot[i]:.2f} # color=green width=2\n')

    # Put in an orientation compass.  Arms are 40-arcsec long

    xCompass = -4.5*boxSize/10.0
    yCompass =  4.0*boxSize/10.0
    (cRA,cDec) = stdToRD(60.0*xCompass,60.0*yCompass,targRA,targDec)
    cRAd=15.0*cRA
    RF.write(f'compass {cRAd}d {cDec}d 40\" # compass=fk5 {{N}} {{E}} 1 1 color=yellow\n')
    
    # Put in the target name from the OBJNAME parameter in the ACQ file.

    if showLabel and len(objName) > 0:
        xName = 0.0
        yName = 4.75*boxSize/10.0
        (cRA,cDec) = stdToRD(60.0*xName,60.0*yName,targRA,targDec)
        cRAd=15.0*cRA
        RF.write(f'text {cRAd}d {cDec}d # text={{objName}} color=yellow font=\'helvetica 14 normal roman\'\n')
        
    # If requested, show the nominal guide probe shadow.  Note that
    # the probe stays fixed on the guide star, so we don't follow any
    # offsets.  We only do this if the guide star is actually in the
    # patrol field

    if hasGS and (gprobe and gsValid):
        # Pickoff shadow, including sensor cable
        (dXgps,dYgps) = rotXY(gpsX,gpsY,-posAng)
        (rsh,dsh) = stdToRD(gsXi+dXgps,gsEta+dYgps,targRA,targDec) 
        rshd = 15.0*rsh
        RF.write(f'box {rshd}d {dsh}d {gpsW/60.0}\' {gpsH/60.0}\' {posAng} # color=yellow width=2\n')
        # Pickoff actuator arm
        (dXarm,dYarm) = rotXY(armX,armY,-posAng)
        (rsh,dsh) = stdToRD(gsXi+dXarm,gsEta+dYarm,targRA,targDec)
        rshd = 15.0*rsh
        RF.write(f'box {rshd}d {dsh}d {armW/60.0}\' {armH/60.0}\' {posAng} # color=yellow width=2\n')
        # Guide camera FoV
        (dXgc,dYgc) = rotXY(gcX,gcY,-posAng)
        (rsh,dsh) = stdToRD(gsXi+dXgc,gsEta+dYgc,targRA,targDec)
        rshd = 15.0*rsh
        RF.write(f'box {rshd}d {dsh}d {gcW/60.0}\' {gcH/60.0}\' {posAng} # color=yellow width=2\n')
        # WFS pickoff FoV ('hot spot')
        RF.write(f'box {gsRAd}d {gsDec}d 8\" 8\" {posAng} # color=yellow\n')
        
    # Draw a heavy cyan circle around the guide star if valid, red if
    # it is invalid

    if hasGS:
        if gsValid:
            RF.write(f'circle {gsRAd}d {gsDec}d 10\" # width=3 color=cyan\n')
        else:
            RF.write(f'circle {gsRAd}d {gsDec}d 10\" # width=3 color=red\n')

    # Close the regions file as our work here is done.  The calling
    # routine is responsible for loading the regions file onto the ds9
    # display

    RF.close()


def printUsage():
    '''
    print the usage message

    Returns
    -------
    None.

    '''
    print('\nUsage: modsView [options] modsScript [fitsFile]')
    print('\nWhere:')
    print('  modsScript is a MODS .acq or .img script')
    print('  fitsFile  optional: use this FITS image with a WCS instead of DSS')
    print('\nOptions:')
    print(' --mms mmsFile overlay slits from an MMS multi-object mask file')
    print(' --shadow    overlay the guide probe pickoff shadow region (default: no shadow)')
    print(' --finder    create a PNG finder chart')
    print(' --grid      overlay celestial coordinate grid (default: no grid)')
    print(' --rotate    rotate to fixed-MODS orientation (default: N=up/E=left)')
    print(' --noalign   do not align the DSS image to N=up/E=left, default: align')
    print(' --size s    change the size of the image to s arcmin (default: 12 arcmin)')
    print(' --cat catID use catalog catID, options: nomad, ub1 or ua2 (default: %s)' % (defCatalog))
    print(' --nocat     do not overlay catalog stars')
    print(' --keepcat   do not delete star catalog working files (default: delete catalogs)')
    print(' --find      Print a list of candidate guide stars to select from')
    print(' --minmag x  specify the catalog faint magnitude limit. default: %.1f' % (minRMag))
    print(' --maxmag x  specify the catalog bright magnitude limit, default: %.1f' % (maxRMag))
    print(' --server x  image server to use, must be one of stsci or eso, default: %s' % (defServer))
    print(' --survey x  sky survey to use, must be valid for server')
    print('               defaults: stsci=all, eso=DSS2-Red')
    print(' --nolabel   do not label the image with OBJNAME')
    print(' --nodisp    only print the analysis and quit without displaying in ds9')
    print(' --kill      kill any delinquent/hidden modsView ds9 window and exit')
    print('  -V         print version info and exit')
    print('\nSee the ds9 manual for server/survey options')
    print('')

#===========================================================================
#
# main program
#

# Other runtime defaults

useDSS = True         # Use the Digitized Sky Survey as the image source
showField = True      # Display the target field in DS9 by default
makeFinder = False    # Do not make a finder unless asked to
alignWCS = True       # Align the image N=up/E=left by default
alignMODS = False     # Align MODS view with the sky (False to align sky to MODS)
showCat = True        # Overlay catalog stars by default
keepCat = False       # Delete the catalog file when done (True to keep)
dssServer = defServer # Default image server
skySurvey = defSurvey # Default sky survey to use (must be valid)
catFile = 'modsView.cat' # Placeholder star catalog file
starCat = defCatalog  # Default optical star catalog
catFilt = defAGwFilt  # Default AGw Guide Camera Filter (options: 'R' or 'B')
boxSize = 12          # Size of the sky box in arcminutes (square)
showSlit = False      # Don't show a slit unless there is one to show
showGrid = False      # Do not overlay a coordinate grid
killDS9  = False      # If True, kill the DS9 window if up and exit
fitsFile = 'none'
interact = False      # Non-interactive by default
catRadius = 10.0      # catalog star to cursor search radius in arcsec
showShadow = False    # Do not draw the nominal pickoff shadow
objName = None        # Blank object name by default
showLabel = True      # Show the objName label by default
showMMS = False       # Show MOS mask slitlets
mmsFile = None        # Name of the MMS file, if showMMS=True
findGStars = False    # Search for candidate guide stars if True, default: False

# Typical limiting magnitudes for guide stars. The manual suggests
# 12.0-16.5, this offers some margin as the USNO photometry is not
# always that great

minMag = minRMag
maxMag = maxRMag

# Parse the command-line arguments using GNU-style getopt

try:
    opts, files = getopt.gnu_getopt(sys.argv[1:],'igflrs:Vc:m:',
                                    ['interact','grid','finder','version','catalog=',
                                     'longslit','noalign','cat=','size=','rotate',
                                     'nocat','minmag=','min=','maxmag=','max=',
                                     'server=','survey=','kill','nodisp','shadow',
                                     'nolabel','mask=','mms=','keepcat','agwfilt=','find',
                                     'mods1','mods2','mods='])
except getopt.GetoptError as err:
    print('\n** ERROR: %s' % (err))
    printUsage()
    sys.exit(2)

if len(opts)==0 and len(files)==0:
    printUsage()
    sys.exit(1)

for opt, arg in opts:
    if opt in ('-g','--grid'):
        showGrid = True
    elif opt in ('-f','--finder'):
        makeFinder = True
    elif opt in ('-i','--interact'):
        interact = True
    elif opt in ('-r','--rotate'):
        alignMODS = True
    elif opt in ('-l','--longslit'):
        showSlit = True
    elif opt in ('-s','--size'):
        boxSize = float(arg)
    elif opt in ('-m','--mask','--mms'):
        mmsFile = arg
        showMMS = True
    elif opt in ('--agwfilt'):
        catFilt = arg
        if catFilt.lower() == 'r':
            catFilt = 'R'
        elif catFilt.lower() == 'b':
            catFilt = 'B'
        else:
            print('\n**ERROR: Invalid AGw Filter option %s, must be R or B' % (catFilt))
            sys.exit(1)
    elif opt in ('--nocat'):
        showCat = False
    elif opt in ('--keepcat'):
        keepCat = True
    elif opt in ('--minmag','--min'):
        minMag = float(arg)
    elif opt in ('--maxmag','--max'):
        maxMag = float(arg)
    elif opt in ('--noalign'):
        alignWCS = False
    elif opt in ('--nodisp'):
        showField = False
    elif opt in ('--nolabel'):
        showLabel = False
    elif opt in ('--shadow'):
        showShadow = True
    elif opt in ('--find'):
        findGStars = True
    elif opt in ('--server'):
        dssServer = arg
        if dssServer.lower() == 'stsci' or dssServer.lower() == 'dss':
            skySurvey = 'DSS2-red' # DSS second epoch
            # skySurvey = 'all'   # best of all available surveys for the field
        elif dssServer.lower() == 'eso':
            skySurvey = 'DSS2RED' # DSS second epoch
        else:
            print('\n**ERROR: Unrecognized image server option %s, must be stsci or eso\n' % (dssServer))
            sys.exit(1)
    elif opt in ('--survey'):
        skySurvey = arg
    elif opt in ('-V','--version'):
        print('modsView.py v%s [%s]' % (versNum,versDate))
        sys.exit(0)
    elif opt in ('-c','--catalog','--cat'):
        showCat = True
        starCat = arg
        if not starCat.lower() in ('ua2','ub1','nomad'):
            print('\n** ERROR: unrecognized star catalog %s, must be nomad, ub1, or ua2\n' % (starCat))
            sys.exit(1)
    # elif opt in ('--kill'):
    #     if isds9up('modsView'):
    #         disp = ds9.DS9('modsView')
    #         disp.set('exit')
    #         print('\nmodsView ds9 window killed\n')
    #     else:
    #         print('\nNo modsView ds9 window is running, nothing to kill.\n')
        sys.exit(0)

# MODS science field at PA=0 (center and width/height in arcsec)

(sciX,sciY) = (0,0)
(sciW,sciH) = (360,360)

# MODS guide patrol field at PA=0 (center and width/height in arcsec)

(gpfX,gpfY) = (0,-150)
(gpfW,gpfH) = (290,300)  

gsBox = (gpfX,gpfY,gpfW,gpfH)

# Nominal guide probe shadow region (center & dimensions) relative to
# the guide star position in arcsec in (xi,eta) coordinates.

(gpsX,gpsY) = (35,0)
(gpsW,gpsH) = (150,85)

# Guide probe carrier arm shadow
        
(armX,armY) = (76,-105)
(armW,armH) = (80,125)

# Guide Field of View position & dimensions relative to the guide star.
# Both MODS have a 40x40-arcsec FoV

(gcX,gcY) = (0.0,0.0)
(gcW,gcH) = (40.0,40.0)


# Centers of the facility long-slit mask segments

slitCen = [-126,-63,0,63,126]

# Take care of file arguments

numFiles = len(files)
       
if numFiles<1:
    print('\n**ERROR: insufficient command-line arguments')
    printUsage()
    sys.exit(1)
elif numFiles>2:
    print('\n**ERROR: too many command-line arguments')
    printUsage()
    sys.exit(1)

if numFiles==2:
    inFile = files[0]
    fitsFile = files[1]
    useDSS = False
else:
    inFile = files[0]
    useDSS = True

# Get the name of the MODS script to view

if not os.path.isfile(inFile):
    print('\n**ERROR: Could not find MODS script file %s\n' % (inFile))
    sys.exit(1)

(fileRoot,fileExt) = os.path.splitext(inFile)

# If using an MMS file, make sure it exists

if showMMS and not os.path.isfile(mmsFile):
    print('\n**ERROR: Could not open MMS file %s\n' % (mmsFile))
    sys.exit(1)

# Check the magnitude limits and make sure they are logical

if maxMag > minMag:
    print('\n**ERROR: Invalid magnitude limits - min=%.1f max=%.1f' % (minMag,maxMag))
    print('         maxMag (bright) must be less than minMag (faint)!\n')
    sys.exit(1)

# This program works with .acq and .img files, not .obs files

if fileExt == '.obs':
    print('\n**ERROR: %s is a MODS observing script.' % (inFile))
    print('         modsView only works with MODS .acq or .img scripts\n')
    sys.exit(1)

# Open the file, scan it in, and extract the coordinates of the
# target and guide star, and the instrument position angle

MF=open(inFile,'r')
fileLines = MF.readlines()[::]
MF.close()

hasTarget = False
hasGStar  = False
hasPosAng = False
hasSlitMask = False
slitMask = None
hasOffset = False
hasOffsetXY = False
offsetRA = 0.0
offsetDec = 0.0
offsetX = 0.0
offsetY = 0.0
offsetType = 'none'
agwFilt = 'none'
acqCamera = None
acqFilter = None
acqExpTime = 0.0

for i in range(len(fileLines)):
    testStr = fileLines[i].strip()
    strBits = testStr.split()
    if not testStr.startswith('#') and len(strBits)>0:
        if strBits[0].upper() == 'OBJCOORDS':
            if len(strBits)<3:
                print('\n**ERROR: OBJCOORDS has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            hasTarget = True
            tRAStr = strBits[1]
            targRA = sex2dec(strBits[1])
            targRAd = 15.0*targRA
            tDecStr = strBits[2]
            targDec = sex2dec(strBits[2])
        elif strBits[0].upper() == 'GUICOORDS':
            if len(strBits)<3:
                print('\n**ERROR: GUICOORDS has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            hasGStar = True
            gRAStr = strBits[1]
            gsRA = sex2dec(strBits[1])
            gsRAd = 15.0*gsRA
            gDecStr = strBits[2]
            gsDec = sex2dec(strBits[2])
        elif strBits[0].upper() == 'POSANGLE':
            if len(strBits)<2:
                print('\n**ERROR: POSANGLE has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            hasPosAng = True
            posAng = float(strBits[1])
        elif strBits[0].upper() == 'ROTATOR':
            if len(strBits)<3:
                print('\n**ERROR: ROTATOR has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            hasPosAng = True
            posAng = float(strBits[1])
            print('\n**NOTE: ROTATOR is deprecated, please use POSANGLE xy.z')
        elif strBits[0].upper() == 'OFFSET':
            if len(strBits)<4:
                print('\n**ERROR: OFFSET has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            hasOffset = True
            offsetRA = float(strBits[1])
            offsetDec = float(strBits[2])
            offsetType = strBits[3].lower()
        elif strBits[0].upper() == 'OFFSETXY':
            if len(strBits)<4:
                print('\n**ERROR: OFFSETXY has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            hasOffsetXY = True
            offsetX = float(strBits[1])
            offsetY = float(strBits[2])
            offsetType = strBits[3].lower()
        elif strBits[0].upper() == 'OBJNAME':
            objName = testStr.split(None,1)[1]
        elif strBits[0].upper() == 'SLITMASK':
            if len(strBits)<2:
                print('**ERROR: SLITMASK has insufficient arguments')
                print('         %s' % (testStr))
                sys.exit(1)
            slitMask = strBits[1]
            hasSlitMask = True
        elif strBits[0].upper() == 'AGWFILT':
            agwFilt = testStr.split(None,1)[1]
            if agwFilt.upper() == 'B_BESSEL':
                catFilt = 'B'
            else:
                catFilt = 'R' # effective bandpass of Clear et al.
        elif strBits[0].upper() == 'ACQCAMERA':
            acqCamera = strBits[1]
        elif strBits[0].upper() == 'ACQFILTER':
            acqFilter = strBits[1]
        elif strBits[0].upper() == 'ACQEXPTIME':
            acqExpTime = float(strBits[1])

# See if we are missing anything in the script.  Behavior depends on
# whether or not we are in interactive mode, searching for guide
# stars, etc. (this is a bit complicated)

if not interact:
    if not findGStars and not hasGStar:
        print('\n*** ERROR: GUICOORDS not found in the script')
        print('\nIs this a properly formatted MODS .acq or .img script?\n')
        print('To search for a guide star, add the --find option.\n')
        sys.exit(1)

    if not hasTarget or not hasPosAng or not hasSlitMask:
        print('\n**ERROR: One or more required elements are missing from the script:')
        if not hasTarget:
            print('         OBJCOORDS not found')
        if not hasPosAng:
            print('         POSANGLE not found')
        if not hasSlitMask:
            print('         SLITMASK not found')
        print('\nIs this a properly formatted MODS .acq or .img script?\n')
        sys.exit(1)
else: 
    if not hasTarget or not hasSlitMask:
        print('\n**ERROR: One or more required elements are missing from the script')
        print('         for interactive mode operation:')
        if not hasTarget:
            print('         OBJCOORDS not found')
        if not hasSlitMask:
            print('         SLITMASK not found')
        print('\nIs this a properly formatted MODS .acq or .img script?\n')
        sys.exit(1)

# Star catalog magnitude to use for filtering

if starCat.lower()=='nomad' or starCat.lower()=='ua2':
    if catFilt == 'R':
        starMag = 'Rmag'  # NOMAD1 or USNOA2.0 R magnitude field
    else:
        starMag = 'Bmag'  # NOMAD1 or USNOA2.0 B magnitude field
else:
    if catFilt == 'R':
        starMag = 'R2mag' # USNOB1.0 R magnitude field
    else:
        starMag = 'B2mag' # USNOB1.0 B magnitude field

# Print a summary:

print('\nMODS %s Script: %s' % (fileExt,inFile))
print('\nSummary:')
print('      Object: %s' % (objName))
print('      Coords: %s %s' % (tRAStr,tDecStr))

if hasPosAng:
    print('  Rotator PA: %.1f deg' % (posAng))
else:
    print('  Rotator PA not given, assuming 0 degrees, select in interactive mode')
    posAng = 0.0

if hasGStar:
    print('  Guide Star: %s %s' % (gRAStr,gDecStr))
else:
    if findGStars:
        print('  Guide Star: NONE - search requested')
    else:
        print('  Guide Star: NONE')

print('  AGw Filter: %s' % (agwFilt))

if hasSlitMask:
    print('   Slit Mask: %s' % (slitMask))
    if slitMask.upper().startswith('LS'):
        showSlit = True
    else:
        showSlit = False
else:
    print('   Slit Mask: None')
    showSlit = False

print(' Acquisition: %s Camera, %s Filter, Exp=%.1f sec' % (acqCamera,acqFilter,acqExpTime))

if hasOffset:
    print('RADEC Offset: dRA=%.2f arcsec dDec=%.2f arcsec type=%s' % (offsetRA,offsetDec,offsetType))
        
elif hasOffsetXY:
    print('DETXY Offset: dX=%.2f arcsec dY=%.2f arcsec type=%s' % (offsetX,offsetY,offsetType))
    
# Compute the position of the guide star in standard coordinates
# (xi,eta) relative to the target at the science field center

if hasGStar:
    print('\nGuide Star Check:')
    (gsXi,gsEta) = rdToStd(gsRA,gsDec,targRA,targDec)
    (xr0,yr0) = rotXY(gpfX,gpfY,-posAng)

# Compute the location of the center and vertices of the guide patrol
# field rotated to pa.  Note the sign convention: pa is flipped from the
# input because of the way xi/eta and celestial PA are defined.

    gsValid = False
    if inBox(gsXi,gsEta,gsBox,posAng):
        print('  The guide star is inside the MODS patrol field after the telescope preset.')
        gsValid = True
    else:
        print('  *** ERROR: The guide star is OUTSIDE the patrol field after the telescope preset.')
        gsValid = False

    # If an offset was requested, test at the offset position 

    dX = 0.0
    dY = 0.0
    dRA = 0.0
    dRAd = 0.0
    dDec = 0.0
    if hasOffsetXY:
        gsOffBox = (gpfX+offsetX,gpfY+offsetY,gpfW,gpfH)
        if inBox(gsXi,gsEta,gsOffBox,posAng):
            print('  The guide star is inside the MODS patrol field after OFFSETXY.')
        else:
            print('  *** ERROR: The guide star is OUTSIDE the patrol field after OFFSETXY.')
        (dX,dY) = rotXY(offsetX,offsetY,-posAng)

    elif hasOffset:
        # Convert offset in arcsec to offset in decimal forms to use later... 
        dRAd = (offsetRA/math.cos(math.radians(targDec)))/3600.0
        dRA  = dRAd/15.0
        dDec = offsetDec/3600.0
        # Test the guide star by offsetting the guide star coordinates in the opposite 
        # direction of the requested telescope RADEC offset
        if inBox(gsXi-offsetRA,gsEta-offsetDec,gsBox,posAng):
            print('  The guide star is inside the MODS patrol field after telescope OFFSET.')
        else:
            print('  *** ERROR: The guide star is OUTSIDE the patrol field after telescope OFFSET.')

    else:
        gsOffBox = gsBox


# Now try to draw the fool thing in ds9, unless they gave --nodisp

if not showField:
    sys.exit(0)

# Setup the DS9 instance to make it look distinctive

disp = DS9('modsView')
disp.set('width 800')
disp.set('height 800')
disp.set('frame clear all')
#disp.set('view layout vertical') # look very different from generic ds9
disp.set('view image no')        # redundant pixel readout
disp.set('view colorbar no')     # don't need color bar

# If using a DSS image, retrieve and display it, otherwise display the
# FITS file given on the command line (provisional, other ways to do
# this

print('\nDisplaying MODS sky view:')
if useDSS:
    imgServer = f"dss{dssServer}"
    print(f" Downloading DSS image from the {imgServer} server...")
    disp.set(f"{imgServer} size {boxSize:.2f} {boxSize:.2f} arcmin")
    disp.set(f"{imgServer} survey {skySurvey}")
    try:
        disp.set(f"{imgServer} coord {tRAStr} {tDecStr}")
    except Exception as exp:
        print(f"*** ERROR: Could not connect to the image server, aborting: {exp}")
        sys.exit(1)

    disp.set(f"{imgServer} close")
    disp.set("scale mode 99.5")
else:
    print(f"  Displaying FITS image {fitsFile}")
    disp.set(f"file fits {fitsFile}")
    disp.set("scale mode zscale")

# Clear any regions, reset display and zoom

disp.set("regions delete all")
disp.set("scale linear")

# Grid and WCS alignment options...

if showGrid:
    disp.set('grid yes')
else:
    disp.set('grid no')
if alignWCS:
    disp.set('wcs align yes')
else:
    disp.set('wcs align no')

if alignMODS:
    disp.set(f"rotate {-posAng:.2f}")

# Finally, zoom the view to fit the display

disp.set('zoom to fit')

# Generate instrument overlay regions file

regFile = str(Path.cwd() / "modsView.reg")

if hasGStar:
    drawMODS(objName,(targRA,targDec),posAng,(gsRA,gsDec),slitMask,(offsetRA,offsetDec), \
                 (offsetX,offsetY),mmsFile,showShadow,boxSize)
else:
    drawMODS(objName,(targRA,targDec),posAng,None,slitMask,(offsetRA,offsetDec), \
                 (offsetX,offsetY),mmsFile,showShadow,boxSize)

# and load it onto the display

disp.set(f"regions file {regFile}")

# Overlay catalog stars with magnitudes if requested

numStars = 0
catUp = False
if showCat:
    ds9cmd = 'catalog %s' % starCat
    disp.set(f"catalog {starCat}")
    disp.set(f"catalog filter ${starMag}<{minMag:.2f}&&${starMag}>{maxMag:.2f}")
    disp.set(f"catalog symbol text ${starMag}")
    catFile = str(Path.cwd() / f"{fileRoot}_{starCat}.cat")
    
    disp.set(f"catalog export tsv {catFile}")
    catUp = True

    # Open and read the catalog into working vectors, for now only
    # store RA, Dec, Rmag, and the catalog star names

    catID = 'None'
    (numStars,catID,catRAd,catDec,catBmag,catRmag,catName) = loadCat(catFile)

    print(f'  Plotting {catID} catalog stars with their {catFilt} magnitudes...')
    print(f'  Magnitude Range: {minMag:.1f} < {catFilt} < {maxMag:.1f}')

    # Since we've loaded a catalog, if we have a guide star seee if it
    # is in the catalog and print info about it
    
    if hasGStar:
        iFound = findStar(gsRAd,gsDec,catRAd,catDec,catRadius)
        if iFound < 0:
            print('  NB: No star found in the %s catalog within the magnitude limits' % (catID))
        else:
            raStar = catRAd[iFound]/15.0
            decStar = catDec[iFound]
            print('\nGuide Star Catalog Info:')
            print('  Star ID: %s %s' % (catID,catName[iFound]))
            print('   Coords: %s %s' % (dec2sex(raStar),dec2sex(decStar)))
            print('     Phot: R=%.2f B=%.2f' % (catRmag[iFound],catBmag[iFound]))

    # If asked to find candidate guide stars, print the list here and
    # ask them to pick one.

    if findGStars:
        if not interact:
            disp.set('catalog close')

        print(f'\nSearching the {catID} catalog excerpt for candidate guide stars...')
        starList = []
        starTable = []
        numFound = 0
        for i in range(numStars):
            raStar = catRAd[i]/15.0
            decStar = catDec[i]
            (xiStar,etaStar) = rdToStd(raStar,decStar,targRA,targDec)
            # Test position after preset
            psetValid = inBox(xiStar,etaStar,gsBox,posAng) 
            # Test position after a requested offset
            offXYValid = True
            if hasOffsetXY:
                gsOffBox = (gpfX+offsetX,gpfY+offsetY,gpfW,gpfH)
                offXYValid = inBox(xiStar,etaStar,gsOffBox,posAng)
            offRDValid = True
            if hasOffset:
                offRDValid = inBox(xiStar-offsetRA,etaStar-offsetDec,gsBox,posAng)
            # A valid guide star must work for preset and offset combination
            if psetValid and offXYValid and offRDValid:
                starList.append(i)
                starTable.append((i+1,catName[i],dec2sex(raStar),dec2sex(decStar),catRmag[i],catBmag[i]))
                numFound += 1

        # Did we find any candidates?  If so, sort the table on R mag
        # and ask them to pick one.

        if numFound == 0:
            print('*** No valid guide stars found within the guide patrol field!')
            hasPicked = True

        else:
            starTable.sort(key=itemgetter(4))
            print(f'\nFound {numFound} candidates:')
            print('  Star     ID           RA          Dec        R     B')
            print('  -------------------------------------------------------')
            for stuff in starTable:
                (i,ID,raStr,decStr,rMag,bMag) = stuff
                if bMag == 99.99:
                    print(f'  {i:3d}  {ID} {raStr} {decStr} {rMag:.2f}')
                else:
                    print(f'  {i:3d}  {ID} {raStr} {decStr} {rMag:.2f} {bMag:.2f}')
            print('  -------------------------------------------------------')
            hasPicked = False

        while not hasPicked:
            pStr = 'Select a guide star (0 to abort): '
            s = input(pStr)
            try:
                iPick = int(s)-1
                if iPick == -1:
                    print('\nGuide Star Selection aborted...')
                    hasGStar = False
                    hasPicked = True
                elif iPick in starList:
                    print('\nGuide Star Selection:')
                    print(f'\n  GUINAME {catID} {catName[iPick]}')
                    raStar = catRAd[iPick]/15.0
                    decStar = catDec[iPick]
                    if decStar>0:
                        print(f'  GUICOORDS {dec2sex(raStr)} +{dec2sex(decStr)} # R={catRmag[iPick]:.2f} B={catBmag[iPick]:.2f}')
                    else:
                        print(f'  GUICOORDS {dec2sex(raStr)} {dec2sex(decStr)} # R={catRmag[iPick]:.2f} B={catBmag[iPick]:.2f}')
                    gsRA = raStar
                    gsDec = decStar
                    hasGStar = True
                    hasPicked = True
                else:
                    print('*** Star %d is not in the candidate list ***' % (iPick+1))
                    hasPicked = False
            except ValueError:
                if len(s) > 0:
                    print(f'{s} is not a valid number, try again...')
                hasPicked = False

        if hasGStar:
            disp.set('regions delete all')
            drawMODS(objName,(targRA,targDec),posAng,(gsRA,gsDec),slitMask,(offsetRA,offsetDec), \
                         (offsetX,offsetY),mmsFile,showShadow,boxSize)
            disp.set(f"regions file {regFile}")

# Make an PNG finder chart if requested

if makeFinder:
    fcFile = str(Path.cwd() / f"{fileRoot}.png")
    disp.set(f"saveimage png {fcFile}")
    print(f"\nCreated finder chart {fileRoot}.png")

# Cleanup: remove delinquent catalog files, close open catalog tools,
#          etc.

if os.path.isfile(catFile):
    if keepCat:
        print(f"Star catalog saved in {catFile}")
    else:
        os.remove(catFile)

if showCat:
    try:
        disp.set('catalog close')
    except:
        pass

# All Done!

print('\n')
sys.exit(0)

