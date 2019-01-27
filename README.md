# modsView
MODS Target Visualization and Guide Star Selection Tool

## Overview

modsView is a Python program that uses DS9 to view and verify the
target and guide star information in a MODS target acquisition (.acq)
or imaging (.img) script.  It displays a digitized sky survey image of
the field in DS9 and overlays the outlines of the MODS focal plane
(science field, AGw patrol field, and any slit apertures), as well as
marking the guide star and overlaying NOMAD catalog stars labeled with
their magnitudes.

While by default it uses digitized sky survey (DSS) images, it can
also accept a user FITS image with a valid world coordinate system
(WCS) in its header as the finder image.

modsView will test to make sure the guide star is inside the MODS AGw
unit's guide patrol field.  If any offset commands follow the
telescope preset in the .acq or .img files, it also will test the
post-offset position of the guide star to make sure it does not move
out of the patrol field.  The final MODS focal plane position drawn
over the sky image will be at the offset position, allowing users to
verify that they have the signs of blind offsets correct.

If a standard slit mask is used (long-slit or spectrophotometric
slit), it will also overlay a view of the slit.  The --mms option is
used to allow display of multi-slit mask apertures given an MMS file.

The --shadow option will overlay a depiction of the region of the
focal plane shadowed by the AGw guide pickoff probe, and show the
effective fields-of-view of the guide camera and wavefront sensor
cameras.

The --finder option creates a PNG-format finding chart suitable for
use with Phase 2 program submission.

Other command-line options allow users to use different catalogs or
image servers, add a grid, and make other minor changes.

## Documentation

An illustrated user's guide and downloads of the latest version are
available online at [www.astronomy.ohio-state.edu/MODS/ObsTools/modsView/](www.astronomy.ohio-state.edu/MODS/ObsTools/modsView), but will migrate off that platform in the future.

## System Requirements

modsView is written in Python and has been tested on the Linux (CentOS
5.x and 6.x) and Mac OS X (v10.7 Lion and v10.13 High Sierra)
operating systems.  Python 2.3 or later is required, and the programs
have been used primarily with Python versions 2.3 through 2.7 to date.

modsView requires the following packages:

 * [pyds9](https://github.com/ericmandel/pyds9)
 * [ds9](http://ds9.si.edu/site/Download.html) - version 7.2 or later
 * [xpa](http://ds9.si.edu/site/XPA.html) - version 2.1.14 or later

### Python 3 Compatibility

The current version has been verified to run on Python 3 (v3.6.5,
anaconda distro), using the latest releases of pyds9. Still finding
bugs, so your mileage may vary (but it's getting better).

## Installation

You have two options: Personal or Public installation

For a personal installation, you have two ways to install the files to use:

1. Keep the modsView Python programs in place, and put the modsView directory in your default execution path.
2. Copy the executable Python programs into your ~/bin/ or  ~/programs/ directory where you put executables in your default execution path.

For a public installation, e.g., on a central disk to share one copy
of the package among many users, we suggest logging in as root and
unpacking the tarball in your usual place for public add-on programs,
e.g., /usr/local/, and then installing the executables in, e.g.,
/usr/local/bin/.


## MacOS Notes

If running on MacOS (e.g, v10.13 HighSierra), we have had good results
with the X-Quartz version of ds9 but not with the Aqua version.  This
is due to us having to launch ds9 with a named pipe, and once launched
pyds9 looks for ds9 in your default path which appears to have some
effect on xpa interaction.  You also need to make sure all the xpa
executables are the same path as ds9.


## Release Notes

These notes are chronological from most recent to oldest
releases. Older copies of modsView are available below, but we
strongly recommend that you only download the latest version. This is
the version that is installed at LBT.

### Version 2

This is the current operational version of modsView to support MODS1+MODS2
Binocular operations

<dl>
<dt>Version 2.1.6 (2018 Sep 05) - Python 3 fixes
<dd>Fixed bug in Python 3 compatibility (raw_input instead of input). Caught by Johnny Greco, tested
    on iMac running MacOS HighSierra and anaconda python 3.6.5

<dt>Verison 2.1.4 (2018 May 22) - SNS Masks
<dd>Support for experimental SNS masks, will not impact normal users

<dt>Version 2.1.3 (2016 Dec 04) - Cosmetic Patch
<dd>This is a patch for version 2.1 that addresses a minor cosmetic
issue in one print statement that elicited comment from some users.

<dt>Version 2.1.2 (2016 Nov 15) - min/maxmag Patch
<dd>This is a patch for version 2.1 that addresses a problem whereby
users could not change the minimum (faint) and maximum (bright)
catalog star magnitudes for display with --minmag/--minmax.

<dt>Version 2.1.1 (2016 Oct 15) - ds9 Patch
<dd>This is a patch for version 2.1 that addresses a problem
encountered by users with an older version of pyds9. While we urge
people to update to the latest version of pyds9, this patch will
bridge the gap.

<dt>Version 2.1 (2016 Oct 12) - MAJOR UPDATE
<dd>Version 2 is a major update to support MODS1 and MODS2 binocular
    operation starting with LBT Semester 2016B when the MODS1 WFS 
    pickoff optic was updated so that it is the same as in MOD2.  
    As a result, there is now no practical distinction between MODS1 and 
    MODS2 for modsView, but, you should check all pre-2016B MODS1 
    scripts with modsView to verify that your guide star selections are
    still valid with the reconfigured MODS1 guider (MODS2 scripts will
    not have this problem).
</dl>

### Previous Versions

These are notes on Version 1 which was pre-binocular MODS operation.

<dl>
<dt>Version 2.0 (2015 Nov 20) - Initial MODS2 support
<dd>Version 2.0 added support for MODS2, and correctly drew the MODS1
    and MODS2 guide-probe positions for the different WFS pickoff optics
    then in MODS1 and MODS2.  After we updated the MODS1 WFS pickoff
    optics to be the same as MODS2 in October 2016, this version was
    retired, and must not be used.

<dt>Version 1.6 (2015 Jan 22) - new features
<dd>Added additional interaction to the --find feature for guide star
    selection.  It now requests you confirm the new guide star choice
    or allow you to repeat the selection or abort.  This makes
    modsView more convenient for initial guide star selection, and at
    the telescope.  Added some additional script syntax checking - for
    example it now requires that OBJNAME and GUINAME have an argument
    (at least "None") otherwise it aborts.  A blank OBJNAME/GUINAME
    will cause an acq script to fail at the telescope.

<dd>NOTE: Version 1.6 is the last release of modsView before
   Version 2.0.  It only supports MODS1. Use Version 2 to pick
   guide stars for MODS1 and MODS2.

<p>
<dt>Version 1.4.3 (2014 Feb 23) - bug patch
<dd>Bug found in plotting the guide patrol field that would allow
    selection of stars that are outside the AGw stage Y-axis limits
    after an offset of ~13mm to put the guide star into the WFS
    hotspot.  Also found a place where the dimensions of the guide
    patrol field had been entered twice and different.  This bug patch
    Addresses issues encountered at the telescope during engineering
    tests that selected guide stars as far as possible away from the
    optical axis on UTC 2014 Feb 23.

<dt>Version 1.4.2 (2013 Apr 11) - bug patch
<dd>Allows for the case of MMS files with no + signs on declinations.
    While this technically not permitted in well-formatted MMS files,
    it can still happen, especially in MMS files that might be altered
    by hand.  We now test for and allow this case.  Also fixed a
    hitherto undetected bug in the rendering of the guide probe shadow
    under RA/Dec offsets.

<dt>Version 1.4.1 (2013 Mar 19) - bug patch
<dd>Fixed a bug in the offsetting logic that would only have been
    apparent in certain RADEC offset combinations at high declinations.

<dt>Version 1.3.4 (2013 Jan 7) - incremental bug patch
<dd>Fixed a bug in --minmag/--maxmag whereby it was possible to give
    modsView nonsensical magnitude limits (bright limit fainter than
    the faint limit).  This condition is now trapped at startup and the 
    crash-point it causes in findStar() patched.

<dt>Version 1.3 (2012 Dec 19) - feature update and bug patch
<dd>Fixed bugs in instrument view after offsets.  Added the ability to
    display tilted MOS slitlets with --mms and restored the --rotate
    feature.  A patch added Dec 19 addresses a minor problem with 
    python 2.7 with the --find option.

<dt>Version 1.2 (2012 Dec 03) - major release
<dd>Added numerous (mostly transparent) enhancements to
    execution.  Added --find to print a candidate list of guide
    stars to select from, and the NOMAD1 catalog is now the
    default star catalog.  The AGWFILT value in the acq script
    is used to select B magnitudes to display on the finder if
    AGWFILT=B_Bessel.  The summary report has been expanded to
    give more details, including star catalog data (R and B mags)
    for the guide star if available.  Also fixed a number of bugs,
    cleaned up some ds9 setup, and updated the webpages.

<dt>Version 0.6.6 (2012 Jun 10) - incremental bug patch 
<dd>Fixes a bug in the display of the guide probe shadow when a blind
    RADEC offset is requested in the acquisition script.

<dt>Version 0.6 (2012 May 21)
<dd>Introduces MMS file overlay and user FITS image support, customization of
    the DS9 window appearance, and quick fixes for some ugly bugs that snuck 
    into v0.5.  This is the last beta testing release before v1.0 (the full 
    public version) is launched.

<dt>Version 0.5 (2012 May 20)
<dd>Second beta release.  Adds the --shadow option to overlay the
    outline of the nominal guide probe shadow and OBJNAME labeling.

<dt>Version 0.4 (2012 May 4)
<dd>First beta release of modsView.

</dl>

## Acknowledging modsView

If you used modsView in your research work, we ask that you follow emerging
[software citation principles](https://doi.org/10.7717/peerj-cs.86) being adopted by the astronomical community
to ensure the proper citation of software in scientific publications. 

[![DOI](https://zenodo.org/badge/142055392.svg)](https://zenodo.org/badge/latestdoi/142055392)

modsView was developed for the MODS1 and MODS2 instruments for the Large Binocular Telescope, which were built with with
major support provided by grants from the U.S. National Science Foundation's Division of Astronomical Sciences Advanced 
Technologies and Instrumention (AST-9987045), the NSF/NOAO TSIP Program, and matching funds provided by the Ohio State
University Office of Research and the Ohio Board of Regents. Additional support was provided by NSF Grant AST-1108693.

## Credits & Future Revisions

modsView was based on an earlier, nasty bit of Perl code that made a
DS9 regions file you could then load by hand (or via xpaset from the
shell).  It was developed mostly for planning purposes early in the
MODS project, but not developed further until we started MODS on-sky
operations.  An intermediate version named luciView was developed at
OSU for checking LUCI scripts, as syntax and guide-star errors in LUCI
scripts was a major problem early on in OSU/RC queue operations.

The need for a similar capability for MODS led to reviving modsView
and rewriting it in Python, building off our experiences with the
modsDisp program deployed for a quasi-realtime display of raw MODS
images at the telescope, and our experiences (good and bad) developing
the lbtView program for selecting guide stars for MODS and other LBT
instruments using the ESO skycat package.  The success with modsView
led us to halt all development of lbtView, which was phased out
completely in 2013.
