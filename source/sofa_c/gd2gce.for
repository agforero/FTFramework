      SUBROUTINE iau_GD2GCE ( A, F, ELONG, PHI, HEIGHT, XYZ, J )
*+
*  - - - - - - - - - - -
*   i a u _ G D 2 G C E
*  - - - - - - - - - - -
*
*  Transform geodetic coordinates to geocentric for a reference
*  ellipsoid of specified form.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     A         d     equatorial radius (Notes 1,4)
*     F         d     flattening (Notes 2,4)
*     ELONG     d     longitude (radians, east +ve)
*     PHI       d     latitude (geodetic, radians, Note 4)
*     HEIGHT    d     height above ellipsoid (geodetic, Notes 3,4)
*
*  Returned:
*     XYZ       d(3)  geocentric vector (Note 3)
*     J         i     status:  0 = OK
*                             -1 = illegal case (Note 4)
*
*  Notes:
*
*  1) The equatorial radius, A, can be in any units, but meters is
*     the conventional choice.
*
*  2) The flattening, F, is (for the Earth) a value around 0.00335,
*     i.e. around 1/298.
*
*  3) The equatorial radius, A, and the height, HEIGHT, must be
*     given in the same units, and determine the units of the
*     returned geocentric vector, XYZ.
*
*  4) No validation is performed on individual arguments.  The error
*     status J=-1 protects against (unrealistic) cases that would lead
*     to arithmetic exceptions.  If an error occurs, XYZ is unchanged.
*
*  5) The inverse transformation is performed in the routine iau_GC2GDE.
*
*  6) The transformation for a standard ellipsoid (such as WGS84) can
*     more conveniently be performed by calling iau_GD2GC, which uses a
*     numerical code (1 for WGS84) to identify the required A and F
*     values.
*
*  References:
*
*     Green, R.M., Spherical Astronomy, Cambridge University Press,
*     (1985) Section 4.5, p96.
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992),
*     Section 4.22, p202.
*
*  This revision:  2009 November 2
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION A, F, ELONG, PHI, HEIGHT, XYZ(3)
      INTEGER J

      DOUBLE PRECISION SP, CP, W, D, AC, AS, R

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


*  Functions of geodetic latitude.
      SP = SIN(PHI)
      CP = COS(PHI)
      W = 1D0-F
      W = W*W
      D = CP*CP + W*SP*SP
      IF ( D .GT. 0D0 ) THEN
         AC = A / SQRT(D)
         AS = W * AC

*     Geocentric vector.
         R = ( AC + HEIGHT ) * CP
         XYZ(1) = R * COS(ELONG)
         XYZ(2) = R * SIN(ELONG)
         XYZ(3) = ( AS + HEIGHT ) * SP

*     Success.
         J = 0
      ELSE

*     Fail.
         J = -1
      END IF

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
