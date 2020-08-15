      SUBROUTINE iau_APCO13 ( UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP,
     :                        PHPA, TC, RH, WL, ASTROM, EO, J )
*+
*  - - - - - - - - - - -
*   i a u _ A P C O 1 3
*  - - - - - - - - - - -
*
*  For a terrestrial observer, prepare star-independent astrometry
*  parameters for transformations between ICRS and observed coordinates.
*  The caller supplies UTC, site coordinates, ambient air conditions and
*  observing wavelength, and SOFA models are used to obtain the Earth
*  ephemeris, CIP/CIO and refraction constants.
*
*  The parameters produced by this routine are required in the parallax,
*  light deflection, aberration, and bias-precession-nutation parts of
*  the ICRS/CIRS transformations.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     UTC1     d      UTC as a 2-part...
*     UTC2     d      ...quasi Julian Date (Notes 1,2)
*     DUT1     d      UT1-UTC (seconds, Note 3)
*     ELONG    d      longitude (radians, east +ve, Note 4)
*     PHI      d      latitude (geodetic, radians, Note 4)
*     HM       d      height above ellipsoid (m, geodetic, Notes 4,6)
*     XP,YP    d      polar motion coordinates (radians, Note 5)
*     PHPA     d      pressure at the observer (hPa = mB, Note 6)
*     TC       d      ambient temperature at the observer (deg C)
*     RH       d      relative humidity at the observer (range 0-1)
*     WL       d      wavelength (micrometers, Note 7)
*
*  Returned:
*     ASTROM   d(30)  star-independent astrometry parameters:
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
*     EO       d      equation of the origins (ERA-GST)
*     J        i      status: +1 = dubious year (Note 2)
*                              0 = OK
*                             -1 = unacceptable date
*
*  Notes:
*
*  1)  UTC1+UTC2 is quasi Julian Date (see Note 2), apportioned in any
*      convenient way between the two arguments, for example where UTC1
*      is the Julian Day Number and UTC2 is the fraction of a day.
*
*      However, JD cannot unambiguously represent UTC during a leap
*      second unless special measures are taken.  The convention in the
*      present routine is that the JD day represents UTC days whether
*      the length is 86399, 86400 or 86401 SI seconds.
*
*      Applications should use the routine iau_DTF2D to convert from
*      calendar date and time of day into 2-part quasi Julian Date, as
*      it implements the leap-second-ambiguity convention just
*      described.
*
*  2)  The warning status "dubious year" flags UTCs that predate the
*      introduction of the time scale or that are too far in the
*      future to be trusted.  See iau_DAT for further details.
*
*  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
*      one second at the end of each positive UTC leap second,
*      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
*      practice is under review, and in the future UT1-UTC may grow
*      essentially without limit.
*
*  4)  The geographical coordinates are with respect to the WGS84
*      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
*      longitude required by the present routine is east-positive
*      (i.e. right-handed), in accordance with geographical convention.
*
*  5)  The polar motion XP,YP can be obtained from IERS bulletins.  The
*      values are the coordinates (in radians) of the Celestial
*      Intermediate Pole with respect to the International Terrestrial
*      Reference System (see IERS Conventions 2003), measured along the
*      meridians 0 and 90 deg west respectively.  For many applications,
*      XP and YP can be set to zero.
*
*      Internally, the polar motion is stored in a form rotated onto
*      the local meridian.
*
*  6)  If hm, the height above the ellipsoid of the observing station
*      in meters, is not known but phpa, the pressure in hPa (=mB), is
*      available, an adequate estimate of hm can be obtained from the
*      expression
*
*            hm = -29.3 * tsl * log ( phpa / 1013.25 );
*
*      where tsl is the approximate sea-level air temperature in K
*      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
*      52).  Similarly, if the pressure phpa is not known, it can be
*      estimated from the height of the observing station, hm, as
*      follows:
*
*            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
*
*      Note, however, that the refraction is nearly proportional to
*      the pressure and that an accurate phpa value is important for
*      precise work.
*
*  7)  The argument WL specifies the observing wavelength in
*      micrometers.  The transition from optical to radio is assumed to
*      occur at 100 micrometers (about 3000 GHz).
*
*  8)  It is advisable to take great care with units, as even unlikely
*      values of the input parameters are accepted and processed in
*      accordance with the models used.
*
*  9)  In cases where the caller wishes to supply his own Earth
*      ephemeris, Earth rotation information and refraction constants,
*      the routine iau_APCO can be used instead of the present routine.
*
*  10) This is one of several routines that inserts into the ASTROM
*      array star-independent parameters needed for the chain of
*      astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
*
*      The various routines support different classes of observer and
*      portions of the transformation chain:
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
*      Those with names ending in "13" use contemporary SOFA models to
*      compute the various ephemerides.  The others accept ephemerides
*      supplied by the caller.
*
*      The transformation from ICRS to GCRS covers space motion,
*      parallax, light deflection, and aberration.  From GCRS to CIRS
*      comprises frame bias and precession-nutation.  From CIRS to
*      observed takes account of Earth rotation, polar motion, diurnal
*      aberration and parallax (unless subsumed into the ICRS <-> GCRS
*      transformation), and atmospheric refraction.
*
*  11) The context array ASTROM produced by this routine is used by
*      iau_ATIOQ, iau_ATOIQ, iau_ATCIQ* and iau_ATICQ*.
*
*  Called:
*     iau_UTCTAI   UTC to TAI
*     iau_TAITT    TAI to TT
*     iau_UTCUT1   UTC to UT1
*     iau_EPV00    Earth position and velocity
*     iau_PNM06A   classical NPB matrix, IAU 2006/2000A
*     iau_BPN2XY   extract CIP X,Y coordinates from NPB matrix
*     iau_S06      the CIO locator s, given X,Y, IAU 2006
*     iau_ERA00    Earth rotation angle, IAU 2000
*     iau_SP00     the TIO locator s', IERS 2000
*     iau_REFCO    refraction constants for given ambient conditions
*     iau_APCO     astrometry parameters, ICRS-observed
*     iau_EORS     equation of the origins, given NPB matrix and s
*
*  This revision:   2013 December 5
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP,
     :                 PHPA, TC, RH, WL, ASTROM(30), EO
      INTEGER J

      INTEGER JS, JW
      DOUBLE PRECISION TAI1, TAI2, TT1, TT2, UT11, UT12,
     :                 EHPV(3,2), EBPV(3,2), R(3,3), X, Y, S, THETA,
     :                 SP, REFA, REFB

      DOUBLE PRECISION iau_S06, iau_ERA00, iau_SP00, iau_EORS

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  UTC to other time scales.
      CALL iau_UTCTAI ( UTC1, UTC2, TAI1, TAI2, JS )
      IF ( JS.LT.0 ) GO TO 9
      CALL iau_TAITT ( TAI1, TAI2, TT1, TT2, JS )
      CALL iau_UTCUT1 ( UTC1, UTC2, DUT1, UT11, UT12, JS )
      IF ( JS.LT.0 ) GO TO 9

*  Earth barycentric & heliocentric position/velocity (au, au/d).
      CALL iau_EPV00 ( TT1, TT2, EHPV, EBPV, JW )

*  Form the equinox based BPN matrix, IAU 2006/2000A.
      CALL iau_PNM06A ( TT1, TT2, R )

*  Extract CIP X,Y.
      CALL iau_BPN2XY ( R, X, Y )

*  Obtain CIO locator s.
      S = iau_S06 ( TT1, TT2, X, Y )

*  Earth rotation angle.
      THETA = iau_ERA00 ( UT11, UT12 )

*  TIO locator s'.
      SP = iau_SP00 ( TT1, TT2 )

*  Refraction constants A and B.
      CALL iau_REFCO ( PHPA, TC, RH, WL, REFA, REFB )

*  Compute the star-independent astrometry parameters.
      CALL iau_APCO ( TT1, TT2, EBPV, EHPV, X, Y, S, THETA,
     :                ELONG, PHI, HM, XP, YP, SP, REFA, REFB,
     :                ASTROM )

*  Equation of the origins.
      EO = iau_EORS ( R, S )

*  Return the status.
 9    CONTINUE
      J = JS

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
