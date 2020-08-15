      SUBROUTINE iau_ATIOQ ( RI, DI, ASTROM, AOB, ZOB, HOB, DOB, ROB )
*+
*  - - - - - - - - - -
*   i a u _ A T I O Q
*  - - - - - - - - - -
*
*  Quick CIRS to observed place transformation.
*
*  Use of this routine is appropriate when efficiency is important and
*  where many star positions are all to be transformed for one date.
*  The star-independent astrometry parameters can be obtained by
*  calling iau_APIO[13] or iau_APCO[13].
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     RI       d      CIRS right ascension
*     DI       d      CIRS declination
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
*
*  Returned:
*     AOB      d      observed azimuth (radians: N=0,E=90)
*     ZOB      d      observed zenith distance (radians)
*     HOB      d      observed hour angle (radians)
*     DOB      d      observed declination (CIO-based, radians)
*     ROB      d      observed right ascension (CIO-based, radians)
*
*  Notes:
*
*  1) This routine returns zenith distance rather than altitude in
*     order to reflect the fact that no allowance is made for
*     depression of the horizon.
*
*  2) The accuracy of the result is limited by the corrections for
*     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
*     Providing the meteorological parameters are known accurately and
*     there are no gross local effects, the predicted observed
*     coordinates should be within 0.05 arcsec (optical) or 1 arcsec
*     (radio) for a zenith distance of less than 70 degrees, better
*     than 30 arcsec (optical or radio) at 85 degrees and better than
*     20 arcmin (optical) or 30 arcmin (radio) at the horizon.
*
*     Without refraction, the complementary routines iau_ATIOQ and
*     iau_ATOIQ are self-consistent to better than 1 microarcsecond all
*     over the celestial sphere.  With refraction included, consistency
*     falls off at high zenith distances, but is still better than
*     0.05 arcsec at 85 degrees.
*
*  3) It is advisable to take great care with units, as even unlikely
*     values of the input parameters are accepted and processed in
*     accordance with the models used.
*
*  4) The CIRS RA,Dec is obtained from a star catalog mean place by
*     allowing for space motion, parallax, the Sun's gravitational lens
*     effect, annual aberration and precession-nutation.  For star
*     positions in the ICRS, these effects can be applied by means of
*     the iau_ATCI13 (etc.) routines.  Starting from classical "mean
*     place" systems, additional transformations will be needed first.
*
*  5) "Observed" Az,El means the position that would be seen by a
*     perfect geodetically aligned theodolite.  This is obtained from
*     the CIRS RA,Dec by allowing for Earth orientation and diurnal
*     aberration, rotating from equator to horizon coordinates, and then
*     adjusting for refraction.  The HA,Dec is obtained by rotating back
*     into equatorial coordinates, and is the position that would be
*     seen by a perfect equatorial with its polar axis aligned to the
*     Earth's axis of rotation.  Finally, the RA is obtained by
*     subtracting the HA from the local ERA.
*
*  6) The star-independent CIRS-to-observed-place parameters in ASTROM
*     may be computed with iau_APIO[13] or iau_APCO[13].  If nothing has
*     changed significantly except the time, iau_APER[13] may be used
*     to perform the requisite adjustment to the ASTROM array.
*
*  Called:
*     iau_S2C      spherical coordinates to unit vector
*     iau_C2S      p-vector to spherical
*     iau_ANP      normalize angle into range 0 to 2pi
*
*  This revision:   2017 July 10
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RI, DI, ASTROM(30), AOB, ZOB, HOB, DOB, ROB

*  Minimum sine and cosine of altitude for refraction purposes
      DOUBLE PRECISION SELMIN, CELMIN
      PARAMETER ( SELMIN = 0.05D0, CELMIN = 1D-6 )

      DOUBLE PRECISION V(3), X, Y, Z, XHD, YHD, ZHD, F,
     :                 XHDT, YHDT, ZHDT, XAET, YAET, ZAET, AZOBS,
     :                 R, TZ, W, DEL, COSDEL, XAEO, YAEO, ZAEO,
     :                 ZDOBS, HMOBS, DCOBS, RAOBS
      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  CIRS RA,Dec to Cartesian -HA,Dec.
      CALL iau_S2C ( RI-ASTROM(28), DI, V )
      X = V(1)
      Y = V(2)
      Z = V(3)

*  Polar motion.
      XHD = X + ASTROM(23)*Z
      YHD = Y - ASTROM(24)*Z
      ZHD = Z - ASTROM(23)*X + ASTROM(24)*Y

*  Diurnal aberration.
      F = ( 1D0 - ASTROM(27)*YHD )
      XHDT = F * XHD
      YHDT = F * ( YHD + ASTROM(27) )
      ZHDT = F * ZHD

*  Cartesian -HA,Dec to Cartesian Az,El (S=0,E=90).
      XAET = ASTROM(25)*XHDT - ASTROM(26)*ZHDT
      YAET = YHDT
      ZAET = ASTROM(26)*XHDT + ASTROM(25)*ZHDT

*  Azimuth (N=0,E=90).
      IF ( XAET.NE.0D0 .OR. YAET.NE.0D0 ) THEN
         AZOBS = ATAN2 ( YAET, -XAET )
      ELSE
         AZOBS = 0D0
      END IF

*  ----------
*  Refraction
*  ----------

*  Cosine and sine of altitude, with precautions.
      R = MAX ( SQRT ( XAET*XAET + YAET*YAET ), CELMIN)
      Z = MAX ( ZAET, SELMIN )

*  A*tan(z)+B*tan^3(z) model, with Newton-Raphson correction.
      TZ = R/Z
      W = ASTROM(30)*TZ*TZ
      DEL = ( ASTROM(29) + W ) * TZ /
     :      ( 1D0 + ( ASTROM(29) + 3D0*W ) / ( Z*Z ) )

*  Apply the change, giving observed vector.
      COSDEL = 1D0 - DEL*DEL/2D0
      F = COSDEL - DEL*Z/R
      XAEO = XAET*F
      YAEO = YAET*F
      ZAEO = COSDEL*ZAET + DEL*R

*  Observed ZD.
      ZDOBS = ATAN2 ( SQRT ( XAEO*XAEO + YAEO*YAEO ), ZAEO )

*  Az/El vector to HA,Dec vector (both right-handed).
      V(1) = ASTROM(25)*XAEO + ASTROM(26)*ZAEO
      V(2) = YAEO
      V(3) = - ASTROM(26)*XAEO + ASTROM(25)*ZAEO

*  To spherical -HA,Dec.
      CALL iau_C2S ( V, HMOBS, DCOBS )

*  Right ascension (with respect to CIO).
      RAOBS = ASTROM(28) + HMOBS

*  Return the results.
      AOB = iau_ANP(AZOBS)
      ZOB = ZDOBS
      HOB = -HMOBS
      DOB = DCOBS
      ROB = iau_ANP(RAOBS)

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
