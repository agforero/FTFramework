      SUBROUTINE iau_ATCIQZ ( RC, DC, ASTROM, RI, DI )
*+
*  - - - - - - - - - - -
*   i a u _ A T C I Q Z
*  - - - - - - - - - - -
*
*  Quick ICRS to CIRS transformation, given precomputed star-independent
*  astrometry parameters, and assuming zero parallax and proper motion.
*
*  Use of this routine is appropriate when efficiency is important and
*  where many star positions are to be transformed for one date.  The
*  star-independent parameters can be obtained by calling one of the
*  routines iau_APCI[13], iau_APCG[13], iau_APCO[13] or iau_APCS[13].
*
*  The corresponding routine for the case of non-zero parallax and
*  proper motion is iau_ATCIQ.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     RC,DC    d      ICRS astrometric RA,Dec (radians)
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
*     RI,DI    d      CIRS RA,Dec (radians)
*
*  Note:
*
*     All the vectors are with respect to BCRS axes.
*
*  References:
*
*     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
*     the Astronomical Almanac, 3rd ed., University Science Books
*     (2013).
*
*     Klioner, Sergei A., "A practical relativistic model for micro-
*     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
*
*  Called:
*     iau_S2C      spherical coordinates to unit vector
*     iau_LDSUN    light deflection due to Sun
*     iau_AB       stellar aberration
*     iau_RXP      product of r-matrix and p-vector
*     iau_C2S      p-vector to spherical
*     iau_ANP      normalize angle into range +/- pi
*
*  This revision:   2013 August 31
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RC, DC, ASTROM(30), RI, DI

      DOUBLE PRECISION PCO(3), PNAT(3), PPR(3), PI(3), W

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  BCRS coordinate direction (unit vector).
      CALL iau_S2C ( RC, DC, PCO )

*  Light deflection by the Sun, giving BCRS natural direction.
      CALL iau_LDSUN ( PCO, ASTROM(5), ASTROM(8), PNAT )

*  Aberration, giving GCRS proper direction.
      CALL iau_AB ( PNAT, ASTROM(9), ASTROM(8), ASTROM(12), PPR )

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
