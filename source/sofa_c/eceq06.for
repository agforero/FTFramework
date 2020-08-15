      SUBROUTINE iau_ECEQ06 ( DATE1, DATE2, DL, DB, DR, DD )
*+
*  - - - - - - - - - - -
*   i a u _ E C E Q 0 6
*  - - - - - - - - - - -
*
*  Transformation from ecliptic coordinates (mean equinox and ecliptic
*  of date) to ICRS RA,Dec, using the IAU 2006 precession model.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2  d      TT as a 2-part Julian Date (Note 1)
*     DL,DB        d      ecliptic longitude and latitude (radians)
*
*  Returned:
*     DR,DD        d      ICRS right ascension and declination (radians)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) No assumptions are made about whether the coordinates represent
*     starlight and embody astrometric effects such as parallax or
*     aberration.
*
*  3) The transformation is approximately that from ecliptic longitude
*     and latitude (mean equinox and ecliptic of date) to mean J2000.0
*     right ascension and declination, with only frame bias (always less
*     than 25 mas) to disturb this classical picture.
*
*  Called:
*     iau_S2C      spherical coordinates to unit vector
*     iau_ECM06    J2000.0 to ecliptic rotation matrix, IAU 2006
*     iau_TRXP     product of transpose of r-matrix and p-vector
*     iau_C2S      unit vector to spherical coordinates
*     iau_ANP      normalize angle into range 0 to 2pi
*     iau_ANPM     normalize angle into range +/- pi
*
*  This revision:  2016 February 9
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DL, DB, DR, DD

      DOUBLE PRECISION RM(3,3), V1(3), V2(3), A, B

      DOUBLE PRECISION iau_ANP, iau_ANPM

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


*  Spherical to Cartesian.
      CALL iau_S2C ( DL, DB, V1 )

*  Rotation matrix, ICRS equatorial to ecliptic.
      CALL iau_ECM06 ( DATE1, DATE2, RM )

*  The transformation from ecliptic to ICRS.
      CALL iau_TRXP ( RM, V1, V2 )

*  Cartesian to spherical.
      CALL iau_C2S ( V2, A, B )

*  Express in conventional ranges.
      DR = iau_ANP ( A )
      DD = iau_ANPM ( B )

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
