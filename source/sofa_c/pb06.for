      SUBROUTINE iau_PB06 ( DATE1, DATE2, BZETA, BZ, BTHETA )
*+
*  - - - - - - - - -
*   i a u _ P B 0 6
*  - - - - - - - - -
*
*  This routine forms three Euler angles which implement general
*  precession from epoch J2000.0, using the IAU 2006 model.  Frame
*  bias (the offset between ICRS and mean J2000.0) is included.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2   d    TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     BZETA         d    1st rotation: radians clockwise around z
*     BZ            d    3rd rotation: radians clockwise around z
*     BTHETA        d    2nd rotation: radians counterclockwise around y
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the arguments DATE1 and DATE2.  For
*     example, JD(TT)=2450123.7 could be expressed in any of these
*     ways, among others:
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
*  2) The traditional accumulated precession angles zeta_A, z_A, theta_A
*     cannot be obtained in the usual way, namely through polynomial
*     expressions, because of the frame bias.  The latter means that two
*     of the angles undergo rapid changes near this date.  They are
*     instead the results of decomposing the precession-bias matrix
*     obtained by using the Fukushima-Williams method, which does not
*     suffer from the problem.  The decomposition returns values which
*     can be used in the conventional formulation and which include
*     frame bias.
*
*  3) The three angles are returned in the conventional order, which
*     is not the same as the order of the corresponding Euler rotations.
*     The precession-bias matrix is R_3(-z) x R_2(+theta) x R_3(-zeta).
*
*  4) Should zeta_A, z_A, theta_A angles be required that do not contain
*     frame bias, they are available by calling the SOFA routine
*     iau_P06E.
*
*  Called:
*     iau_PMAT06   PB matrix, IAU 2006
*     iau_RZ       rotate around Z-axis
*
*  This revision:  2020 May 27
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, BZETA, BZ, BTHETA

      DOUBLE PRECISION R(3,3), Y, X

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Precession matrix via Fukushima-Williams angles.
      CALL iau_PMAT06 ( DATE1, DATE2, R )

*  Solve for z, choosing the +/- pi alternative.
      Y = R(2,3)
      X = -R(1,3)
      IF ( X .LT. 0D0 ) THEN
         Y = -Y
         X = -X
      END IF
      IF ( X.NE.0D0 .OR. Y.NE.0D0 ) THEN
         BZ = - ATAN2 ( Y, X )
      ELSE
         BZ = 0D0
      END IF

*  Derotate it out of the matrix.
      CALL iau_RZ ( BZ, R )

*  Solve for the remaining two angles.
      Y = R(1,3)
      X = R(3,3)
      IF ( X.NE.0D0 .OR. Y.NE.0D0 ) THEN
         BTHETA = - ATAN2 ( Y, X )
      ELSE
         BTHETA = 0D0
      END IF

      Y = -R(2,1)
      X = R(2,2)
      IF ( X.NE.0D0 .OR. Y.NE.0D0 ) THEN
         BZETA = - ATAN2 ( Y, X )
      ELSE
         BZETA = 0D0
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
