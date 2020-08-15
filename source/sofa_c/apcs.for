      SUBROUTINE iau_APCS ( DATE1, DATE2, PV, EBPV, EHP, ASTROM )
*+
*  - - - - - - - - -
*   i a u _ A P C S
*  - - - - - - - - -
*
*  For an observer whose geocentric position and velocity are known,
*  prepare star-independent astrometry parameters for transformations
*  between ICRS and GCRS.  The Earth ephemeris is supplied by the
*  caller.
*
*  The parameters produced by this routine are required in the space
*  motion, parallax, light deflection and aberration parts of the
*  astrometric transformation chain.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1    d       TDB as a 2-part...
*     DATE2    d       ...Julian Date (Note 1)
*     PV       d(3,2)  observer's geocentric pos/vel (m, m/s)
*     EBPV     d(3,2)  Earth barycentric position/velocity (au, au/day)
*     EHP      d(3)    Earth heliocentric position (au)
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
*               (22)     unchanged
*               (23)     unchanged
*               (24)     unchanged
*               (25)     unchanged
*               (26)     unchanged
*               (27)     unchanged
*               (28)     unchanged
*               (29)     unchanged
*               (30)     unchanged
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
*  2) All the vectors are with respect to BCRS axes.
*
*  3) Providing separate arguments for (i) the observer's geocentric
*     position and velocity and (ii) the Earth ephemeris is done for
*     convenience in the geocentric, terrestrial and Earth orbit cases.
*     For deep space applications it maybe more convenient to specify
*     zero geocentric position and velocity and to supply the
*     observer's position and velocity information directly instead of
*     with respect to the Earth.  However, note the different units:
*     m and m/s for the geocentric vectors, au and au/day for the
*     heliocentric and barycentric vectors.
*
*  4) In cases where the caller does not wish to provide the Earth
*     ephemeris, the routine iau_APCS13 can be used instead of the
*     present routine.  This computes the Earth ephemeris using the
*     SOFA routine iau_EPV00.
*
*  5) This is one of several routines that inserts into the ASTROM
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
*  6) The context array ASTROM produced by this routine is used by
*     iau_ATCIQ* and iau_ATICQ*.
*
*  Called:
*     iau_CP       copy p-vector
*     iau_PM       modulus of p-vector
*     iau_PN       decompose p-vector into modulus and direction
*     iau_IR       initialize r-matrix to identity
*
*  This revision:   2017 March 16
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION DATE1, DATE2, PV(3,2), EBPV(3,2), EHP(3),
     :                 ASTROM(30)

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian year
      DOUBLE PRECISION DJY
      PARAMETER ( DJY = 365.25D0 )

*  Seconds per day
      DOUBLE PRECISION DAYSEC
      PARAMETER ( DAYSEC = 86400D0 )

*  Speed of light (m/s)
      DOUBLE PRECISION CMPS
      PARAMETER ( CMPS = 299792458D0 )

*  Astronomical unit (m, IAU 2012)
      DOUBLE PRECISION AUM
      PARAMETER ( AUM = 149597870.7D3 )

*  Light time for 1 au (s)
      DOUBLE PRECISION AULT
      PARAMETER ( AULT = AUM/CMPS )

*  au/d to m/s
      DOUBLE PRECISION AUDMS
      PARAMETER ( AUDMS = AUM/DAYSEC )

*  Light time for 1 au (day)
      DOUBLE PRECISION CR
      PARAMETER ( CR = AULT/DAYSEC )

      INTEGER I
      DOUBLE PRECISION DP, DV, PB(3), VB(3), PH(3), V2, W

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Time since reference epoch, years (for proper motion calculation).
      ASTROM(1) = ( ( DATE1 - DJ00 ) + DATE2 ) / DJY

*  Adjust Earth ephemeris to observer.
      DO 1 I=1,3
         DP = PV(I,1) / AUM
         DV = PV(I,2) / AUDMS
         PB(I) = EBPV(I,1) + DP
         VB(I) = EBPV(I,2) + DV
         PH(I) = EHP(I) + DP
 1    CONTINUE

*  Barycentric position of observer (au).
      CALL iau_CP ( PB, ASTROM(2) )

*  Heliocentric direction and distance (unit vector and au).
      CALL iau_PN ( PH, ASTROM(8), ASTROM(5) )

*  Barycentric vel. in units of c, and reciprocal of Lorenz factor.
      V2 = 0D0
      DO 3 I=1,3
         W = VB(I) * CR
         ASTROM(8+I) = W
         V2 = V2 + W*W
 3    CONTINUE
      ASTROM(12) = SQRT ( 1D0 - V2 )

*  Reset the NPB matrix.
      CALL iau_IR ( ASTROM(13) )

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
