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
#   Requires your system have SAOImage DS9 and the XPA utilities
#   installed and in your default path (hea-www.harvard.edu/saord/ds9/)
#
#   Requires the ds9 Python module for DS9 interaction
#     github.com/ericmandel/pyds9
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
#
#---------------------------------------------------------------------------

import sys
import os
import readline
import math
import getopt
import subprocess
import shlex
from time import sleep
from operator import itemgetter

# pyds9 has moved from SAO to github and the module name changed

try:
    import pyds9 as ds9
except:
    import ds9

# input vs raw_input for Python 3/2 compatibility

try:
    input = raw_input
except NameError:
    pass

# Version number and date, update as needed

versNum  = '2.2.1'
versDate = '2022-11-11'

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
# sex2dec - Sexagesimal to decimal conversion
#
# Inputs:
#   inStr = string with the sexagesimal number in +/-dd:mm:ss.ss format
#
# Returns:
#   Decimal representation of the sexagesmial string
#
# Description:
#   Converts a signed or unsigned sexagesimal string in dd:mm:ss.sss
#   notation into a decimal equivalent.  It is agnostic about units
#   (e.g., degrees or hours).  Does not test validity on any range.
#   It correctly handles the conversion for the -00:mm:ss.ss case.
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 2
#

def sex2dec(inStr):
    bits = inStr.split(':')
    dsign = 1.0
    major = float(bits[0])
    minor = float(bits[1])
    sec = float(bits[2])
    if major==0.0 and inStr.startswith('-'):
        dsign = -1.0
    if major < 0.0 :
        decimal = major - minor/60.0 - sec/3600.0
    else :
        decimal = major + minor/60.0 + sec/3600.0
    decimal *= dsign
    return decimal

#---------------------------------------------------------------------------
# 
# dec2sex - Decimal to Sexagesimal to conversion
#
# Inputs:
#   angle = decimal angle to convert
#
# Returns:
#   string with the sexagesimal number in +/-dd:mm:ss.ss format
#
# Description:
#   Converts a floating-point decimal 'angle' into a sexagesimal
#   string in dd:mm:ss.sss notation.  It is agnostic about units
#   (e.g., degrees or hours).  Does not test validity on any range.
#   The final 'seconds' part is rounded to the nearest 0.01 seconds.
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 7
#

def dec2sex(angle):
    arg = math.fabs(angle)
    dd = int(arg)
    temp = 60.0*(arg - float(dd))
    mm = int(temp)
    ss = 3600.0*(arg - (float(dd)+(float(mm)/60.0)))
    if angle < 0:
        sexStr = '-%02d:%02d:%05.2f' % (dd,mm,ss)
    else:
        sexStr = '%02d:%02d:%05.2f' % (dd,mm,ss)
    return sexStr

#---------------------------------------------------------------------------
#
# rdToStd - convert celestial coordinates (ra,dec) to standard coordinates
#           (xi,eta) on the tanget plane.
#
# Inputs:
#   ra,dec = target RA and Dec in decimal hours and degrees, respectively
#   ra0,dec0 = reference (tangent) point RA and Dec in decimal hrs/deg
#
# Returns: 
#   xi,eta = standard coordinates in arcseconds
#  
# Description:
#   Given the celestial coordinates (ra,dec) of an object, and the
#   celestial coordinates of the field center (ra0,dec0), compute
#   standard coordinates (xi,eta) of the object in the tangent plane
#   to the celestial sphere.  This geometric problem is described in
#   W.W. Smart, Textbook on Spherical Astronomy, Chapter XII, sections 
#   160 and 161.
#
# See also stdToRD()
#
# R. Pogge, OSU Astronomy Dept
# pogge.1@osu.edu
# 2012 May 3
#

def rdToStd(ra,dec,ra0,dec0):
    dra = math.radians(15.0*(ra-ra0))
    cosd0 = math.cos(math.radians(dec0))
    sind0 = math.sin(math.radians(dec0))
    tand  = math.tan(math.radians(dec))
    cosdra = math.cos(dra)
    denom = sind0*tand + cosd0*cosdra
    xi = 3600.0 * math.degrees(math.sin(dra)/denom)
    eta = 3600.0 * math.degrees((cosd0*tand - sind0*cosdra)/denom)
    return xi, eta

#---------------------------------------------------------------------------
#
# stdToRD - convert standard coordinates (xi,eta) on the tanget plane to
#           celestial coordinates (ra,dec)
#
# Inputs:
#   xi,eta   = object standard coordinates in arcseconds
#   ra0,dec0 = reference (tangent) point RA and Dec in decimal hrs/deg
#
# Returns: 
#   ra,dec = celestial coordinates in decimal hours and degrees, respetively
#  
# Description:
#   Given standard coordinates (xi,eta) of an object and the celestial
#   coordinates of the field center (ra0,dec0), compute the celestial
#   coordinates (ra,dec) of the object on the celestial sphere.  This
#   geometric problem is described in W.W. Smart, Textbook on Spherical
#   Astronomy, Chapter XII, sections 160 and 161.
#
# See also rdToStd()
#
# R. Pogge, OSU Astronomy Dept
# pogge.1@osu.edu
# 2012 May 3
#

def stdToRD(xi,eta,ra0,dec0):
    xi_r  = math.radians(xi/3600.0)
    eta_r = math.radians(eta/3600.0)
    ra0_r = math.radians(15.0*ra0)

    cosd0 = math.cos(math.radians(dec0))
    tand0 = math.tan(math.radians(dec0))

    denom = 1.0 - eta_r*tand0

    dra = math.atan(xi_r/(cosd0*denom))

    ra = math.degrees(ra0_r+dra)/15.0
    dec = math.degrees(math.atan((math.cos(dra)*(eta_r+tand0))/denom))

    return ra, dec

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

#---------------------------------------------------------------------------
#
# inBox - is (x,y) inside a given rectangle at a given field rotation?
#
# Inputs:
#   x,y = coordinates of the test point
#   rect = rectangle (unrotated frame)
#   rotAng = rotation angle about (0,0)
#
# Breaks the box into two triangles and then tests the point against
# these triangles using inTriangle()
#

def inBox(x,y,rect,rotAng):
    (xc,yc,dx,dy) = rect;
    # Coordinates of the center and rectangle vertices after rotation
    (xr0,yr0)=rotXY(xc,yc,-posAng)
    (xr1,yr1)=rotXY(xc-dx/2,yc-dy/2,-rotAng)
    (xr2,yr2)=rotXY(xc+dx/2,yc-dy/2,-rotAng)
    (xr3,yr3)=rotXY(xc+dx/2,yc+dy/2,-rotAng)
    (xr4,yr4)=rotXY(xc-dx/2,yc+dy/2,-rotAng)

    # Slice the guide patrol field diagonally into two triangles.  The
    # guide star must be inside one or other triangle
    
    if inTriangle(xr1,yr1,xr2,yr2,xr3,yr3,x,y) or \
       inTriangle(xr1,yr1,xr3,yr3,xr4,yr4,x,y):
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
    if rotAng==0.0 or rotAng==-0.0:
        return x, y
    else:
        sinPA = math.sin(math.radians(rotAng))
        cosPA = math.cos(math.radians(rotAng))
        xr = x*cosPA - y*sinPA
        yr = x*sinPA + y*cosPA
        return xr, yr

#----------------------------------------------------------------
#
# parseMMS - parse the contents of a MODS Mask design (MMS) file
#
# Input: 
#   file [str] = name of the MMS file to open and parse
#
# Returns:
#   ra,dec = RA, Dec in decimal hours/degrees
#   wid,len,rot = slit width and length in decimal arcseconds
#   
# Description:
#   Opens and reads in the contents of a MODS mask design (MMS) file,
#   and extracts the coordinates and dimensions of the slits,
#   returning them as useful floating-point representations.  The
#   trick lies in reading the MMS file's coding for slit RA/Dec
#   coordinates:
#
#     RA format: INS.TARGnnn.ALPHA 203448.232
#    Dec format: INS.TARGnnn.DELTA -001415.123
# 
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 21
#
# 2012 Dec 16 - Added slitlet rotation angles [rwp/osu]
# 2013 Apr 11 - Allow for (technically incorrect) lack of + sign
#               declination entries. [rwp/osu]
#

def parseMMS(file): 
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
  numSlits = len(ra)
  return ra, dec, swid, slen, srot


#---------------------------------------------------------------------------
#
# findStar - find a star in a star catalog given RA, Dec and a search radius
#
# inputs:
#   ra, dec = RA/Dec in decimal degrees of the test point
#   catRAd,catDec = arrays with the catalog star RA/Dec in decimal degrees
#                   (how it comes from the DS9 catalog servers)
#   radius = search radius in arcseconds.  Star must be this close to
#            the cursor to be 'found'
#
# returns 0..N-1, the index of the star found, or -1 if no star found
#

def findStar(ra,dec,catRAd,catDec,radius):
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

#---------------------------------------------------------------------------
#
# isds9up - see if a named DS9 window is up
#
# Inputs:
#   ds9ID = ID ('title') of a DS9 display window
#
# Returns:
#   True if ds9ID is running, False otherwise.
#
# Description:
#   Uses the shell's xpaaccess method to see if the named ds9 window
#   is up and running.
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 3
#

def isds9up(ds9ID):
    test = subprocess.Popen(['xpaaccess','-n',ds9ID],
                            stdout=subprocess.PIPE).communicate()[0]
    if int(test):
        return True
    else:
        return False


#---------------------------------------------------------------------------
#
# startDS9 - launch a named ds9 window
#
# Inputs:
#   ds9ID = ID ('title') of a DS9 display window to open
#
# Description:
#   Launches a named DS9 instance, making sure all of the IRAF
#   imtool pipes are suppresed so that IRAF won't interfere
#   with it (and vis-vers).  It sleeps for 2 seconds to allow
#   the tool to open.  This may have to be increased on slower
#   or more loaded systems.
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#   2012 May 3
#

def startDS9(ds9ID):
    cmdStr = 'ds9 -fifo none -port none -unix none -title %s' % (ds9ID)
    args = shlex.split(cmdStr)
    subprocess.Popen(args)
    sleep(2)

#---------------------------------------------------------------------------
#
# loadCat - load in a star catalog into working arrays
#
# inputs:
#   catFile - name of the catalog file
#
# returns:
#   numStars, catID, catRAd, catDec, catBmag, catRmag, catName
# where:
#     numStars = [scalar] number of stars in the catalog
#     catID = [string] catalog ID (e.g., USNO-B1.0)
#     catRAd = [vector] RA in decimal degrees (float)
#     catDec = [vector] Dec in decimal degrees (float)
#     catBmag = [vector] B magnitude (float)
#     catRmag = [vector] R magnitude (float)
#     catName = [vector] star ID (string)
#

def loadCat(catFile):
    SC=open(catFile,'r')
    catLines = SC.readlines()[::]
    SC.close()
    catRAd  = []
    catDec  = []
    catRmag = []
    catBmag = []
    catName = []
    catID = 'None'
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

    # Where are we working? (required for getting paths right for DS9 later)

    myDir = os.getcwd()

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

    regFile = os.path.join(myDir,'modsView.reg')
    if os.path.isfile(regFile):
        os.remove(regFile)
    RF = open(regFile,'w')

    RF.write('#\n# modsView regions file\n#\nfk5\n')

    # Draw the 'aim point', the original preset coordinates before any
    # offsets are applied, as a yellow circle

    regCmd = 'circle %fd %fd 3.0\" # color=yellow width=2\n' % (targRAd,targDec)
    RF.write(regCmd)

    # RA/Dec of the instrument aim point after all offsets are applied

    (instRA,instDec) = stdToRD(sciX+dX,sciY+dY,targRA+dRA,targDec+dDec)
    instRAd = 15.0*instRA

    # The circle shows the full sweep of the MODS science and guide patrol fields

    regCmd = 'circle %fd %fd %f\' # color=red\n' % (instRAd,instDec,fsFoV-0.0333)
    RF.write(regCmd)
    regCmd = 'circle %fd %fd %f\' # color=cyan\n' % (instRAd,instDec,fsFoV)
    RF.write(regCmd)
    regCmd = 'circle %fd %fd %f\' # color=red\n' % (instRAd,instDec,fsFoV+0.033)
    RF.write(regCmd)
    
    # RA/Dec of the science field center after all offsets are applied

    (sciRA,sciDec) = stdToRD(sciX+dX,sciY+dY,targRA+dRA,targDec+dDec)
    sciRAd = 15.0*sciRA
    
    # Draw the MODS science field

    regCmd = 'box %fd %fd 6\' 6\' %f # width=2 color=green\n' % (sciRAd,sciDec,posAng)
    RF.write(regCmd)

    # RA/Dec of the guide patrol field after all offsets are applied

    (gpfRA,gpfDec) = stdToRD(xr0+dX,yr0+dY,targRA+dRA,targDec+dDec)
    gpfRAd = 15.0*gpfRA

    # Draw the MODS AGw Guide Patrol Field

    regCmd = 'box %fd %fd %f\' %f\' %f # color=cyan width=2\n' % (gpfRAd,gpfDec,(gpfW/60.0),(gpfH/60.0),posAng)
    RF.write(regCmd)

    # Show the long-slit locations if using a facility mask

    if showSlit:
        if slitMask.upper() == 'LS60X5':
            slitWid=5.0
            slitLen=60.0
            regCmd = 'box %fd %fd %f\' %f\' %f\n' % \
                (sciRAd,sciDec,(slitWid/60.0),(slitLen/60.0),posAng)
            RF.write(regCmd)

        elif slitMask.upper() == 'LS10X0.8SNS':
            slitWid = 0.8
            slitLen = 10.0
            ys = -120.0
            (xsr,ysr) = rotXY(0,ys,-posAng)
            (rs,ds) = stdToRD(xsr,ysr,sciRA,sciDec)
            rsd = 15.0*rs
            regCmd = 'box %fd %fd %f\' %f\' %f\n' % \
                (rsd,ds,(slitWid/60.0),(slitLen/60.0),posAng)
            RF.write(regCmd)

        else:
            bits = slitMask.upper().split('X')
            slitWid = float(bits[len(bits)-1])
            slitLen=60.0
            for ys in slitCen:
                (xsr,ysr) = rotXY(0,ys,-posAng)
                (rs,ds) = stdToRD(xsr,ysr,sciRA,sciDec)
                rsd = 15.0*rs
                regCmd = 'box %fd %fd %f\' %f\' %f\n' % \
                    (rsd,ds,(slitWid/60.0),(slitLen/60.0),posAng)
                RF.write(regCmd)

   # If using a MOS mask and showMMS, display the MOS mask slitlets

    if showMMS:
        (mmsRA,mmsDec,mmsWid,mmsLen,mmsRot) = parseMMS(mmsFile)
        for i in range(len(mmsRA)): 
            if mmsWid[i]==mmsLen[i]: # assume square = alignment box, color magenta
                regCmd = 'box %fd %fd %.1f\" %.1f\" %.2f # color=magenta width=2\n' % \
                    (15.0*mmsRA[i],mmsDec[i],mmsWid[i],mmsLen[i],posAng+mmsRot[i])
            else:
                regCmd = 'box %fd %fd %.1f\" %.1f\" %.2f # color=green\n' % \
                    (15.0*mmsRA[i],mmsDec[i],mmsWid[i],mmsLen[i],posAng+mmsRot[i])
            RF.write(regCmd)

    # Put in an orientation compass.  Arms are 40-arcsec long

    xCompass = -4.5*boxSize/10.0
    yCompass =  4.0*boxSize/10.0
    (cRA,cDec) = stdToRD(60.0*xCompass,60.0*yCompass,targRA,targDec)
    cRAd=15.0*cRA
    regCmd = 'compass %fd %fd %f\" # compass=fk5 {N} {E} 1 1 color=yellow\n' % (cRAd,cDec,40)
    RF.write(regCmd)

    # Put in the target name from the OBJNAME parameter in the ACQ file.

    if showLabel and len(objName) > 0:
        xName = 0.0
        yName = 4.75*boxSize/10.0
        (cRA,cDec) = stdToRD(60.0*xName,60.0*yName,targRA,targDec)
        cRAd=15.0*cRA
        regCmd = 'text %fd %fd # text={%s} color=yellow font=\'helvetica 14 normal roman\'\n' \
            % (cRAd,cDec,objName)
        RF.write(regCmd)

    # If requested, show the nominal guide probe shadow.  Note that
    # the probe stays fixed on the guide star, so we don't follow any
    # offsets.  We only do this if the guide star is actually in the
    # patrol field

    if hasGS and (gprobe and gsValid):
        # Pickoff shadow, including sensor cable
        (dXgps,dYgps) = rotXY(gpsX,gpsY,-posAng)
        (rsh,dsh) = stdToRD(gsXi+dXgps,gsEta+dYgps,targRA,targDec) 
        rshd = 15.0*rsh
        regCmd = 'box %fd %fd %f\' %f\' %f # color=yellow width=2\n' % \
            (rshd,dsh,(gpsW/60.0),(gpsH/60.0),posAng)
        RF.write(regCmd)
        # Pickoff actuator arm
        (dXarm,dYarm) = rotXY(armX,armY,-posAng)
        (rsh,dsh) = stdToRD(gsXi+dXarm,gsEta+dYarm,targRA,targDec)
        rshd = 15.0*rsh
        regCmd = 'box %fd %fd %f\' %f\' %f # color=yellow width=2\n' % \
            (rshd,dsh,(armW/60.0),(armH/60.0),posAng)
        RF.write(regCmd)
        # Guide camera FoV
        (dXgc,dYgc) = rotXY(gcX,gcY,-posAng)
        (rsh,dsh) = stdToRD(gsXi+dXgc,gsEta+dYgc,targRA,targDec)
        rshd = 15.0*rsh
        regCmd = 'box %fd %fd %f\' %f\' %f # color=yellow\n' % \
            (rshd,dsh,(gcW/60.0),(gcH/60.0),posAng)
        RF.write(regCmd)
        # WFS pickoff FoV ('hot spot')
        regCmd = 'box %fd %fd 8\" 8\" %f # color=yellow\n' % (gsRAd,gsDec,posAng)
        RF.write(regCmd)
    
    # Draw a heavy cyan circle around the guide star if valid, red if
    # it is invalid

    if hasGS:
        if gsValid:
            regCmd = 'circle %fd %fd 10\" # width=3 color=cyan\n' % (gsRAd,gsDec)
        else:
            regCmd = 'circle %fd %fd 10\" # width=3 color=red\n' % (gsRAd,gsDec)
        RF.write(regCmd)

    # Close the regions file as our work here is done.  The calling
    # routine is responsible for loading the regions file onto the ds9
    # display

    RF.close()


#---------------------------------------------------------------------------
#
# printUsage() - print the usage message
#

def printUsage():
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
    elif opt in ('--kill'):
        if isds9up('modsView'):
            disp = ds9.DS9('modsView')
            disp.set('exit')
            print('\nmodsView ds9 window killed\n')
        else:
            print('\nNo modsView ds9 window is running, nothing to kill.\n')
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

# Where are we working? (required for getting paths right for DS9 later)

myDir = os.getcwd()

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

if not isds9up('modsView'):
    startDS9('modsView')

# Setup the DS9 instance to make it look distinctive

disp = ds9.DS9('modsView')
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
    imgServer = 'dss'+dssServer
    print('  Downloading DSS image from the server...')
    ds9cmd = '%s size %.2f %.2f arcmin' % (imgServer,boxSize,boxSize)
    disp.set(ds9cmd)
    ds9cmd = '%s survey %s' % (imgServer,skySurvey)
    disp.set(ds9cmd)
    ds9cmd = f'{imgServer} coord {tRAStr} {tDecStr}'
    try:
        disp.set(ds9cmd)
    except:
        print('  *** ERROR: Could not connect to the image server, aborting')
        sys.exit(1)

    ds9cmd = '%s close' % (imgServer)
    disp.set(ds9cmd)
    disp.set('scale mode 99.5')
#    disp.set('scale mode minmax')
else:
    print('  Displaying FITS image %s...' % (fitsFile))
    ds9cmd = 'file fits %s' % (fitsFile)
    disp.set(ds9cmd)
    disp.set('scale mode zscale')

# Clear any regions, reset display and zoom

disp.set('regions delete all')
disp.set('scale linear')

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
    ds9cmd = 'rotate %.2f' % (-posAng)
    disp.set(ds9cmd)

# Finally, zoom the view to fit the display

disp.set('zoom to fit')

# Generate instrument overlay regions file

regFile = os.path.join(myDir,'modsView.reg')

if hasGStar:
    drawMODS(objName,(targRA,targDec),posAng,(gsRA,gsDec),slitMask,(offsetRA,offsetDec), \
                 (offsetX,offsetY),mmsFile,showShadow,boxSize)
else:
    drawMODS(objName,(targRA,targDec),posAng,None,slitMask,(offsetRA,offsetDec), \
                 (offsetX,offsetY),mmsFile,showShadow,boxSize)


# and load it onto the display

ds9cmd = 'regions file %s' % (regFile)
disp.set(ds9cmd)

# Overlay catalog stars with magnitudes if requested

numStars = 0
catUp = False
if showCat:
    ds9cmd = 'catalog %s' % starCat
    disp.set(ds9cmd)
    ds9cmd = 'catalog filter $%s<%.2f&&$%s>%.2f' % (starMag,minMag,starMag,maxMag)
    disp.set(ds9cmd)
    ds9cmd = 'catalog symbol text $%s' % (starMag)
    disp.set(ds9cmd)
    tmpFile = fileRoot+'_'+starCat+'.cat'
    catFile = os.path.join(myDir,tmpFile)
    #ds9cmd = 'catalog save tsv %s' % (catFile)
    ds9cmd = 'catalog export tsv %s' % (catFile)
    disp.set(ds9cmd)
    #disp.set('catalog close')
    catUp = True

    # Open and read the catalog into working vectors, for now only
    # store RA, Dec, Rmag, and the catalog star names

    catID = 'None'
    (numStars,catID,catRAd,catDec,catBmag,catRmag,catName) = loadCat(catFile)

    print('  Plotting %s catalog stars with their %s magnitudes...' % (catID,catFilt))
    print('  Magnitude Range: %.1f < %s < %.1f' % (minMag,catFilt,maxMag))

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

        print('\nSearching the %s catalog excerpt for candidate guide stars...' % (catID))
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
            print('\nFound %d candidates:' % (numFound))
            print('  Star     ID           RA          Dec        R     B')
            print('  -------------------------------------------------------')
            for stuff in starTable:
                (i,ID,raStr,decStr,rMag,bMag) = stuff
                if bMag == 99.99:
                    print('  %3d  %s %s %s %.2f' % (i,ID,raStr,decStr,rMag))
                else:
                    print('  %3d  %s %s %s %.2f %.2f' % (i,ID,raStr,decStr,rMag,bMag))
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
                    print('\n  GUINAME %s %s' % (catID,catName[iPick]))
                    raStar = catRAd[iPick]/15.0
                    decStar = catDec[iPick]
                    if decStar>0:
                        print('  GUICOORDS %s +%s # R=%.2f B=%.2f' % (dec2sex(raStar),dec2sex(decStar),catRmag[iPick],catBmag[iPick]))
                    else:
                        print('  GUICOORDS %s %s # R=%.2f B=%.2f' % (dec2sex(raStar),dec2sex(decStar),catRmag[iPick],catBmag[iPick]))
                    gsRA = raStar
                    gsDec = decStar
                    hasGStar = True
                    hasPicked = True
                else:
                    print('*** Star %d is not in the candidate list ***' % (iPick+1))
                    hasPicked = False
            except ValueError:
                if len(s) > 0:
                    print('%s is not a valid number, try again...' % (s))
                hasPicked = False

        if hasGStar:
            disp.set('regions delete all')
            drawMODS(objName,(targRA,targDec),posAng,(gsRA,gsDec),slitMask,(offsetRA,offsetDec), \
                         (offsetX,offsetY),mmsFile,showShadow,boxSize)
            ds9cmd = 'regions file %s' % (regFile)
            disp.set(ds9cmd)

#----------------------------------------
#
# Experimental interactive cursor mode
#

if interact:
    cursorLive = True
    print('Entering interactive cursor mode...')
    while cursorLive:
        cursData = ''
        gotCursor = False
        try:
            cursData = disp.get('imexam key coordinate fk5')
        except:
            # (type, value, tb) = sys.exc_info()
            # print('** ds9 cursor glitch - cursData=%s - %s (%s)' % (cursData,type,value.message))
            gotCursor = False

        else:
            if len(cursData) > 0:
                gotCursor = True
            else:
                gotCursor = False
                
        if not gotCursor:
            continue

        # We have imexam output to process

        (cursKey, cursRAStr, cursDecStr) = cursData.split()
        cursRAd = float(cursRAStr)
        cursDec = float(cursDecStr)
        
        # Process Cursor Keys

        # q - quit interactive mode

        if cursKey.lower() == 'q':
            cursorLive = False
            continue

        # */i - get catalog star info
        #       * = mark nearest star as a guide star
        #       i = just print info about the nearest star

        elif cursKey.lower() == 'asterisk' or cursKey.lower()=='i':
            if numStars==0:
                print('  ** No stars in the catalog...')
            else:
                iFound = findStar(cursRAd,cursDec,catRAd,catDec,catRadius)
                if iFound < 0:
                    print('  No %s star was found within %d-arcsec of the cursor' % (catID,int(catRadius)))
                    print('  NOTE: Only stars with %.1f<%s<%.1f are plotted.' % (minMag,catFilt,maxMag))
                    print('        Use the m key to change the magnitude limits.')
                else:
                    radStar = catRAd[iFound]
                    raStar = radStar/15.0
                    decStar = catDec[iFound]
                    (xiStar,etaStar) = rdToStd(raStar,decStar,targRA,targDec)
                    #
                    # * key = select a guide star
                    #
                    if cursKey.lower() == 'asterisk':
                        print('\nGuide Star Selection:')
                        print('\n  GUINAME %s %s' % (catID,catName[iFound]))
                        if decStar>0:
                            print('  GUICOORDS %s +%s # R=%.2f B=%.2f\n' % (dec2sex(raStar),dec2sex(decStar),catRmag[iFound],catBmag[iFound]))
                        else:
                            print('  GUICOORDS %s %s # R=%.2f\n' % (dec2sex(raStar),dec2sex(decStar),catRmag[iFound],catBmag[iFound]))
                        gsRA = raStar
                        gsDec = decStar
                        hasGStar = True
                    else:
                        if decStar>0:
                            print('  \nStar nearest cursor: %s %s +%s Rmag=%.2f Bmag=%.2f' % \
                                (catName[iFound],dec2sex(raStar),dec2sex(decStar),catRmag[iFound],catBmag[iFound]))
                        else:
                            print('  \nStar nearest cursor: %s %s %s Rmag=%.2f Bmag=%.2f' % \
                                (catName[iFound],dec2sex(raStar),dec2sex(decStar),catRmag[iFound],catBmag[iFound]))
                    if inBox(xiStar,etaStar,gsBox,posAng):
                        ds9cmd = 'regions command {circle %fd %fd 10\' # width=3 color=cyan}' % (radStar,decStar)
                        print('  This star is inside the post-preset MODS AGw guide field')
                    else:
                        ds9cmd = 'regions command {circle %fd %fd 10\' # width=3 color=red}' % (radStar,decStar)
                        if cursKey.lower() == 'asterisk':
                            print('  **WARNING: This star is OUTSIDE the post-preset MODS AGw guide field')
                        else:
                            print('  NOTE: This star is OUTSIDE the post-preset MODS AGw guide field')
                    disp.set(ds9cmd)

        # r - redraw/refresh the display

        elif cursKey.lower() == 'r':
            disp.set('regions delete all')
            if hasGStar:
                drawMODS(objName,(targRA,targDec),posAng,(gsRA,gsDec),slitMask,(offsetRA,offsetDec), \
                             (offsetX,offsetY),mmsFile,showShadow,boxSize)
            else:
                drawMODS(objName,(targRA,targDec),posAng,None,slitMask,(offsetRA,offsetDec), \
                             (offsetX,offsetY),mmsFile,showShadow,boxSize)
            ds9cmd = 'regions file %s' % (regFile)
            disp.set(ds9cmd)

        # g - toggle the WCS coordinate grid on/off

        elif cursKey.lower() == 'g':
            if showGrid:
                disp.set('grid no')
                showGrid = False
            else:
                disp.set('grid yes')
                showGrid = True

        # w - toggle WCS alignment on/off 
        #     On : N=up/E=left
        #     Off: use native pixel orientation)

        elif cursKey.lower() == 'w':
            if alignWCS:
                disp.set('wcs align no')
                alignWCS = False
            else:
                disp.set('wcs align yes')
                alignWCS = True

        # w - toggle display of catalog stars on/off

        elif cursKey.lower() == 'c':
            if showCat and catUp:
                catUp = False
                disp.set('catalog hide')
            elif showCat and not catUp:
                catUp = True
                disp.set('catalog show')
            else:
                print('  No star catalog to toggle.')

        # m - change catalog star display min/max magnitude limits

        elif cursKey.lower() == 'm':
            if not showCat:
                print('No star catalog is loaded')
            else:
                print('\nCurrent Magnitude Limits: %.1f to %.1f' % (minMag,maxMag))
                pStr = '  New faint limit  [%.1f] > ' % (minMag)
                s = input(pStr)
                try:
                    newMin = float(s)
                    minMag = newMin
                except ValueError:
                    if len(s)> 0:
                        print('%s is not a valid number, still using faint limit = %.1f mag' % (s,minMag))
                pStr = '  New bright limit [%.1f] > ' % (maxMag)
                s = input(pStr)
                try:
                    newMax = float(s)
                    maxMag = newMax
                except ValueError: 
                    if len(s)> 0:
                        print('%s is not a valid number, still using bright limit = %.1f mag' % (s,maxMag))
                   
                # Re-generate and re-load the catalog

                print('Re-plotting %s stars within limits %.1f<%s<%.1f...' % (catID,minMag,catFilt,maxMag))
                ds9cmd = 'catalog filter $%s<%.2f&&$%s>%.2f' % (starMag,minMag,starMag,maxMag)
                disp.set(ds9cmd)
                ds9cmd = 'catalog save tsv %s' % (catFile)
                disp.set(ds9cmd)
                (numStars,catID,catRAd,catDec,catBmag,catRmag,catName) = loadCat(catFile)

        # p - change the slit mask position angle

        elif cursKey.lower() == 'p':
            print('\nChange Slit Mask Position Angle:\n')
            pStr = '  New Slit PA [%.2f] > ' % (posAng)
            s = input(pStr)
            try:
                newPA = float(s)
                if newPA>-270.0 and newPA<270.0:
                    posAng = newPA
                    disp.set('regions delete all')
                    if hasGStar:
                        drawMODS(objName,(targRA,targDec),posAng,(gsRA,gsDec),slitMask,(offsetRA,offsetDec), \
                                     (offsetX,offsetY),mmsFile,showShadow,boxSize)
                    else:
                        drawMODS(objName,(targRA,targDec),posAng,None,slitMask,(offsetRA,offsetDec), \
                                     (offsetX,offsetY),mmsFile,showShadow,boxSize)

                    ds9cmd = 'regions file %s' % (regFile)
                    disp.set(ds9cmd)

                else:
                    print('Invalid PA %s, must be -270..270 degrees' % (s))

            except ValueError:
                if len(s)> 0:
                    print('%s is not a valid number, still PA=%.2f deg' % (s,posAng))


        # v - print version info (engineering)

        elif cursKey.lower() == 'v':
            print('\nmodsView.py v%s [%s]\n' % (versNum,versDate))
            print('pyds9 v%s' % (disp.__version__))

        # ? - print list of cursor keys

        elif cursKey == 'question':
            print('\nInteractive Cursor Keys:')
            print('  q = QUIT interactive mode')
            print('  r = redraw the MODS focal plane overlay')
            print('  i = get info about the catalog star nearest cursor')
            print('  * = select the catalog star nearest the cursor as the guide star')
            print('  c = toggle the star catalog overlay on/off')
            print('  m = change the star catalog magnitude limits')
            print('  p = change the rotator position angle')
            print('  g = toggle the WCS coordinate grid on/off')
            print('  w = toggle image WCS alignment on/off')
            print('  v = show version information about modsView')
            print('  ? = show a list of the cursor keys')

        # Bottom of cursKey check
    # Bottom of interactive cursor loop

    print('\nExiting interactive Mode.  Cleanup in progress.')

# Make an PNG finder chart if requested

if makeFinder:
    tmpFile = fileRoot+'.png'
    fcFile = os.path.join(myDir,tmpFile)
    ds9cmd = 'saveimage png %s' % (fcFile)
    disp.set(ds9cmd)
    print('\nCreated finder chart %s' % (tmpFile))

# Cleanup: remove delinquent catalog files, close open catalog tools,
#          etc.

if os.path.isfile(catFile):
    if keepCat:
        print('Star catalog saved in %s' % (catFile))
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

