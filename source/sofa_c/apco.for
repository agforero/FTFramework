      SUBROUTINE iau_APCO ( DATE1, DATE2, EBPV, EHP, X, Y, S,
     :                      THETA, ELONG, PHI, HM, XP, YP, SP,
     :                      REFA, REFB, ASTROM )
*+
*  - - - - - - - - -
*   i a u _ A P C O
*  - - - - - - - - -
*
*  For a terrestrial observer, prepare star-independent astrometry
*  parameters for transformations between ICRS and observed coordinates.
*  The caller supplies the Earth ephemeris, the Earth rotation
*  information and the refraction constants as well as the site
*  coordinates.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1    d       TDB as a 2-part...
*     DATE2    d       ...Julian Date (Note 1)
*     EBPV     d(3,2)  Earth barycentric pos/vel (au, au/day, Note 2)
*     EHP      d(3)    Earth heliocentric position (au, Note 2)
*     X,Y      d       CIP X,Y (components of unit vector)
*     S        d       the CIO locator s (radians)
*     THETA    d       Earth rotation angle (radians)
*     ELONG    d       longitude (radians, east +ve, Note 3)
*     PHI      d       latitude (geodetic, radians, Note 3)
*     HM       d       height above ellipsoid (m, geodetic, Note 3)
*     XP,YP    d       polar motion coordinates (radians, Note 4)
*     SP       d       the TIO locator s' (radians, Note 4)
*     REFA     d       refraction constant A (radians, Note 5)
*     REFB     d       refraction constant B (radians, Note 5)
*
*  Returned:
*     ASTROM   d(30)   star-independent astrometry parameters:
*               (1)      PM time interval (SSB, Julian years)
*               (2-4)    SSB to observer (vector, au)
*               (5-7)    Sun to observer (unit vector)
*               (8)      distance from Sun to observer (au)
*               (9-11)   v: barycentric observer velocity (vector, c)
*               (12)     sqrt(1-|v|^2): reciprocal of Lorenz factor
*               (13-21)  bias-precession-nutation matrix
*               (22)     longitude + s' (radians)
*               (23)     polar motion xp wrt local meridian (radians)
*               (24)     polar motion yp wrt local meridian (radians)
*               (25)     sine of geodetic latitude
*               (26)     cosine of geodetic latitude
*               (27)     magnitude of diurnal aberration vector
*               (28)     "local" Earth rotation angle (radians)
*               (29)     refraction constant A (radians)
*               (30)     refraction constant B (radians)
*
*  Notes:
*
*  1) The TDB date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TDB)=2450123.7 could be expressed in any of these ways, among
*     others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in cases
*     where the loss of several decimal digits of resolution is
*     acceptable.  The J2000 method is best matched to the way the
*     argument is handled internally and will deliver the optimum
*     resolution.  The MJD method and the date & time methods are both
*     good compromises between resolution and convenience.  For most
*     applications of this routine the choice will not be at all
*     critical.
*
*     TT can be used instead of TDB without any significant impact on
*     accuracy.
*
*  2) The vectors EB, EH, and all the ASTROM vectors, are with respect
*     to BCRS axes.
*
*  3) The geographical coordinates are with respect to the WGS84
*     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN
*     CONVENTION:  the longitude required by the present routine is
*     right-handed, i.e. east-positive, in accordance with geographical
*     convention.
*
*  4) XP and YP are the coordinates (in radians) of the Celestial
*     Intermediate Pole with respect to the International Terrestrial
*     Reference System (see IERS Conventions), measured along the
*     meridians 0 and 90 deg west respectively.  SP is the TIO locator
*     s', in radians, which positions the Terrestrial Intermediate
*     Origin on the equator.  For many applications, XP, YP and
*     (especially) SP can be set to zero.
*
*     Internally, the polar motion is stored in a form rotated onto the
*     local meridian.
*
*  5) The refraction constants REFA and REFB are for use in a
*     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
*     (i.e. refracted) zenith distance and dZ is the amount of
*     refraction.
*
*  6) It is advisable to take great care with units, as even unlikely
*     values of the input parameters are accepted and processed in
*     accordance with the models used.
*
*  7) In cases where the caller does not wish to provide the Earth
*     Ephemeris, the Earth rotation information and refraction
*     constants, the routine iau_APCO13 can be used instead of the
*     present routine.  This starts from UTC and weather readings etc.
*     and computes suitable values using other SOFA routines.
*
*  8) This is one of several routines that inserts into the ASTROM
*     array star-independent parameters needed for the chain of
*     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
*
*     The various routines support different classes of observer and
*     portions of the transformation chain:
*
*            routines           observer        transformation
*
*        iau_APCG iau_APCG13    geocentric      ICRS <-> GCRS
*        iau_APCI iau_APCI13    terrestrial     ICRS <-> CIRS
*        iau_APCO iau_APCO13    terrestrial     ICRS <-> observed
*        iau_APCS iau_APCS13    space           ICRS <-> GCRS
*        iau_APER iau_APER13    terrestrial     update Earth rotation
*        iau_APIO iau_APIO13    terrestrial     CIRS <-> observed
*
*     Those with names ending in "13" use contemporary SOFA models to
*     compute the various ephemerides.  The others accept ephemerides
*     supplied by the caller.
*
*     The transformation from ICRS to GCRS covers space motion,
*     parallax, light deflection, and aberration.  From GCRS to CIRS
*     comprises frame bias and precession-nutation.  From CIRS to
*     observed takes account of Earth rotation, polar motion, diurnal
*     aberration and parallax (unless subsumed into the ICRS <-> GCRS
*     transformation), and atmospheric refraction.
*
*  9) The context array ASTROM produced by this routine is used by
*     iau_ATIOQ, iau_ATOIQ, iau_ATCIQ*, and iau_ATICQ*.
*
*  Called:
*     iau_APER     astrometry parameters: update ERA
*     iau_C2IXYS   celestial-to-intermediate matrix, given X,Y and s
*     iau_PVTOB    position/velocity of terrestrial station
*     iau_TRXPV    product of transpose of r-matrix and pv-vector
*     iau_APCS     astrometry parameters, ICRS-GCRS, space observer
*     iau_CR       copy r-matrix
*
*  This revision:   2017 March 12
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION DATE1, DATE2, EBPV(3,2), EHP(3), X, Y, S,
     :                 THETA, ELONG, PHI, HM, XP, YP, SP, REFA, REFB,
     :                 ASTROM(30)

      DOUBLE PRECISION SL, CL, R(3,3), PVC(3,2), PV(3,2)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Longitude with adjustment for TIO locator s'.
      ASTROM(22) = ELONG + SP

*  Polar motion, rotated onto the local meridian.
      SL = SIN(ASTROM(22))
      CL = COS(ASTROM(22))
      ASTROM(23) = XP*CL - YP*SL
      ASTROM(24) = XP*SL + YP*CL

*  Functions of latitude.
      ASTROM(25) = SIN(PHI)
      ASTROM(26) = COS(PHI)

*  Refraction constants.
      ASTROM(29) = REFA
      ASTROM(30) = REFB

*  Local Earth rotation angle.
      CALL iau_APER ( THETA, ASTROM )

*  Disable the (redundant) diurnal aberration step.
      ASTROM(27) = 0D0

*  CIO based BPN matrix.
      CALL iau_C2IXYS ( X, Y, S, R )

*  Observer's geocentric position and velocity (m, m/s, CIRS).
      CALL iau_PVTOB ( ELONG, PHI, HM, XP, YP, SP, THETA, PVC )

*  Rotate into GCRS.
      CALL iau_TRXPV ( R, PVC, PV )

*  ICRS <-> GCRS parameters.
      CALL iau_APCS ( DATE1, DATE2, PV, EBPV, EHP, ASTROM )

*  Store the CIO based BPN matrix.
      CALL iau_CR ( R, ASTROM(13) )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2020
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
*  CONDITIONS WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The names of all routines in your derived work shall not
*        include the prefix "iau" or "sofa" or trivial modifications
*        thereof such as changes of case.
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  5. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  6. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  In any published work or commercial product which uses the SOFA
*  software directly, acknowledgement (see www.iausofa.org) is
*  appreciated.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*-----------------------------------------------------------------------

      END
