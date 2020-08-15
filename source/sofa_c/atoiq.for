      SUBROUTINE iau_ATOIQ ( TYPE, OB1, OB2, ASTROM, RI, DI )
*+
*  - - - - - - - - - -
*   i a u _ A T O I Q
*  - - - - - - - - - -
*
*  Quick observed place to CIRS, given the star-independent astrometry
*  parameters.
*
*  Use of this routine is appropriate when efficiency is important and
*  where many star positions are all to be transformed for one date.
*  The star-independent astrometry parameters can be obtained by calling
*  iau_APIO[13] or iau_APCO[13].
*
*  Status:  support routine.
*
*  Given:
*     TYPE     c*(*)  type of coordinates: 'R', 'H' or 'A' (Note 2)
*     OB1      d      observed Az, HA or RA (radians; Az is N=0,E=90)
*     OB2      d      observed ZD or Dec (radians)
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
*     RI       d      CIRS right ascension (CIO-based, radians)
*     DI       d      CIRS declination (radians)
*
*  Notes:
*
*  1) "Observed" Az,El means the position that would be seen by a
*     perfect geodetically aligned theodolite.  This is related to the
*     observed HA,Dec via the standard rotation, using the geodetic
*     latitude (corrected for polar motion), while the observed HA and
*     RA are related simply through the Earth rotation angle and the
*     site longitude.  "Observed" RA,Dec or HA,Dec thus means the
*     position that would be seen by a perfect equatorial with its polar
*     axis aligned to the Earth's axis of rotation.  By removing from
*     the observed place the effects of atmospheric refraction and
*     diurnal aberration, the CIRS RA,Dec is obtained.
*
*  2) Only the first character of the type argument is significant.
*     'R' or 'r' indicates that OB1 and OB2 are the observed right
*     ascension and declination;  'H' or 'h' indicates that they are
*     hour angle (west +ve) and declination;  anything else ('A' or
*     'a' is recommended) indicates that OB1 and OB2 are azimuth (north
*     zero, east 90 deg) and zenith distance.  (Zenith distance is used
*     rather than altitude in order to reflect the fact that no
*     allowance is made for depression of the horizon.)
*
*  3) The accuracy of the result is limited by the corrections for
*     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
*     Providing the meteorological parameters are known accurately and
*     there are no gross local effects, the predicted observed
*     coordinates should be within 0D05 arcsec (optical) or 1 arcsec
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
*  4) It is advisable to take great care with units, as even unlikely
*     values of the input parameters are accepted and processed in
*     accordance with the models used.
*
*  5) The star-independent astrometry parameters in ASTROM may be
*     computed with iau_APIO13 (or iau_APIO).  If nothing has changed
*     significantly except the time, iau_APER13 (or iau_APER) may be
*     used to perform the requisite adjustment to the ASTROM array.
*
*  Called:
*     iau_S2C      spherical coordinates to unit vector
*     iau_C2S      p-vector to spherical
*     iau_ANP      normalize angle into range 0 to 2pi
*
*  This revision:   2013 August 3
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) TYPE
      DOUBLE PRECISION OB1, OB2, ASTROM(30), RI, DI

      CHARACTER C
      DOUBLE PRECISION C1, C2, SPHI, CPHI, CE, XAEO, YAEO, ZAEO, V(3),
     :                 XMHDO, YMHDO, ZMHDO, AZ, SZ, ZDO, REFA, REFB,
     :                 TZ, DREF, ZDT, XAET, YAET, ZAET,
     :                 XMHDA, YMHDA, ZMHDA, F, XHD, YHD, ZHD,
     :                 XPL, YPL, W, HMA

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Coordinate type.
      C = TYPE(:1)

*  Coordinates.
      C1 = OB1
      C2 = OB2

*  Sin, cos of latitude.
      SPHI = ASTROM(25)
      CPHI = ASTROM(26)

*  Standardize coordinate type.
      IF ( C.EQ.'r' .OR. C.EQ.'R' ) THEN
         C = 'R'
      ELSE IF ( C.EQ.'h' .OR. C.EQ.'H' ) THEN
         C = 'H'
      ELSE
         C = 'A'
      END IF

*  If Az,ZD, convert to Cartesian (S=0,E=90).
      IF ( C.EQ.'A' ) THEN
         CE = SIN(C2)
         XAEO = - COS(C1) * CE
         YAEO = SIN(C1) * CE
         ZAEO = COS(C2)
      ELSE

*     If RA,Dec, convert to HA,Dec.
         IF ( C.EQ.'R' ) C1 = ASTROM(28) - C1

*     To Cartesian -HA,DeC.
         CALL iau_S2C ( -C1, C2, V )
         XMHDO = V(1)
         YMHDO = V(2)
         ZMHDO = V(3)

*     To Cartesian Az,El (S=0,E=90).
         XAEO = SPHI*XMHDO - CPHI*ZMHDO
         YAEO = YMHDO
         ZAEO = CPHI*XMHDO + SPHI*ZMHDO
      END IF

*  Azimuth (S=0,E=90).
      IF ( XAEO.NE.0D0 .OR. YAEO.NE.0D0 ) THEN
         AZ = ATAN2(YAEO,XAEO)
      ELSE
         AZ = 0D0
      END IF

*  Sine of observed ZD, and observed ZD.
      SZ = SQRT ( XAEO*XAEO + YAEO*YAEO )
      ZDO = ATAN2 ( SZ, ZAEO )

*
* Refraction
* ----------

*  Fast algorithm using two constant model.
      REFA = ASTROM(29)
      REFB = ASTROM(30)
      TZ = SZ / ZAEO
      DREF = ( REFA + REFB*TZ*TZ ) * TZ
      ZDT = ZDO + DREF

*  To Cartesian Az,ZD.
      CE = SIN(ZDT)
      XAET = COS(AZ) * CE
      YAET = SIN(AZ) * CE
      ZAET = COS(ZDT)

*  Cartesian Az,ZD to Cartesian -HA,DeC.
      XMHDA = SPHI*XAET + CPHI*ZAET
      YMHDA = YAET
      ZMHDA = - CPHI*XAET + SPHI*ZAET

*  Diurnal aberration.
      F = ( 1D0 + ASTROM(27)*YMHDA )
      XHD = F * XMHDA
      YHD = F * ( YMHDA - ASTROM(27) )
      ZHD = F * ZMHDA

*  Polar motion.
      XPL = ASTROM(23)
      YPL = ASTROM(24)
      W = XPL*XHD - YPL*YHD + ZHD
      V(1) = XHD - XPL*W
      V(2) = YHD + YPL*W
      V(3) = W - ( XPL*XPL + YPL*YPL ) * ZHD

*  To spherical -HA,DeC.
      CALL iau_C2S ( V, HMA, DI )

*  Right ascension.
      RI = iau_ANP ( ASTROM(28) + HMA )

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
