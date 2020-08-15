      SUBROUTINE iau_ATCIQN ( RC, DC, PR, PD, PX, RV, ASTROM, N, B,
     :                        RI, DI )
*+
*  - - - - - - - - - - -
*   i a u _ A T C I Q N
*  - - - - - - - - - - -
*
*  Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
*  star-independent astrometry parameters plus a list of light-
*  deflecting bodies.
*
*  Use of this routine is appropriate when efficiency is important and
*  where many star positions are to be transformed for one date.  The
*  star-independent parameters can be obtained by calling one of the
*  routines iau_APCI[13], iau_APCG[13], iau_APCO[13] or iau_APCS[13].
*
*  If the only light-deflecting body to be taken into account is the
*  Sun, the iau_ATCIQ routine can be used instead.  If in addition the
*  parallax and proper motions are zero, the iau_ATCIQZ routine can be
*  used.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     RC,DC    d      ICRS RA,Dec at J2000.0 (radians, Note 1)
*     PR       d      RA proper motion (radians/year; Note 2)
*     PD       d      Dec proper motion (radians/year)
*     PX       d      parallax (arcsec)
*     RV       d      radial velocity (km/s, +ve if receding)
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
*     N        i       number of bodies (Note 3)
*     B        d(8,N)  data for each of the NB bodies (Notes 3,4):
*               (1,I)    mass of the body (solar masses, Note 5)
*               (2,I)    deflection limiter (Note 6)
*               (3-5,I)  barycentric position of the body (au)
*               (6-8,I)  barycentric velocity of the body (au/day)
*
*  Returned:
*     RI,DI    d      CIRS RA,Dec (radians)
*
*  Notes:
*
*  1) Star data for an epoch other than J2000.0 (for example from the
*     Hipparcos catalog, which has an epoch of J1991.25) will require a
*     preliminary call to iau_PMSAFE before use.
*
*  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
*
*  3) The array B contains N entries, one for each body to be
*     considered.  If N = 0, no gravitational light deflection will be
*     applied, not even for the Sun.
*
*  4) The array B should include an entry for the Sun as well as for any
*     planet or other body to be taken into account.  The entries should
*     be in the order in which the light passes the body.
*
*  5) In the entry in the B array for body I, the mass parameter B(1,I)
*     can, as required, be adjusted in order to allow for such effects
*     as quadrupole field.
*
*  6) The deflection limiter parameter B(2,I) is phi^2/2, where phi is
*     the angular separation (in radians) between star and body at which
*     limiting is applied.  As phi shrinks below the chosen threshold,
*     the deflection is artificially reduced, reaching zero for phi = 0.
*     Example values suitable for a terrestrial observer, together with
*     masses, are as follows:
*
*        body I     B(1,I)         B(2,I)
*
*        Sun        1D0            6D-6
*        Jupiter    0.00095435D0   3D-9
*        Saturn     0.00028574D0   3D-10
*
*  7) For efficiency, validation of the B array is omitted.  The
*     supplied masses must be greater than zero, the position and
*     velocity vectors must be right, and the deflection limiter
*     greater than zero.
*
*  Called:
*     iau_PMPX     proper motion and parallax
*     iau_LDN      light deflection by n bodies
*     iau_AB       stellar aberration
*     iau_RXP      product of r-matrix and pv-vector
*     iau_C2S      p-vector to spherical
*     iau_ANP      normalize angle into range 0 to 2pi
*
*  This revision:   2013 September 30
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RC, DC, PR, PD, PX, RV, ASTROM(30)
      INTEGER N
      DOUBLE PRECISION B(8,N), RI, DI

      DOUBLE PRECISION PCO(3), PNAT(3), PPR(3), PI(3), W

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Proper motion and parallax, giving BCRS coordinate direction.
      CALL iau_PMPX ( RC, DC, PR, PD, PX, RV, ASTROM(1), ASTROM(2),
     :                PCO )

*  Light deflection, giving BCRS natural direction.
      CALL iau_LDN ( N, B, ASTROM(2), PCO, PNAT )

*  Aberration, giving GCRS proper direction.
      CALL iau_AB (PNAT, ASTROM(9), ASTROM(8), ASTROM(12), PPR )

*  Bias-precession-nutation, giving CIRS proper direction.
      CALL iau_RXP ( ASTROM(13), PPR, PI )

*  CIRS RA,Dec.
      CALL iau_C2S ( PI, W, DI )
      RI = iau_ANP ( W )

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
