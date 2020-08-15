      SUBROUTINE iau_PAP ( A, B, THETA )
*+
*  - - - - - - - -
*   i a u _ P A P
*  - - - - - - - -
*
*  Position-angle from two p-vectors.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     A        d(3)      direction of reference point
*     B        d(3)      direction of point whose PA is required
*
*  Returned:
*     THETA    d         position angle of B with respect to A (radians)
*
*  Notes:
*
*  1) The result is the position angle, in radians, of direction B with
*     respect to direction A.  It is in the range -pi to +pi.  The sense
*     is such that if B is a small distance "north" of A the position
*     angle is approximately zero, and if B is a small distance "east" of
*     A the position angle is approximately +pi/2.
*
*  2) A and B need not be unit vectors.
*
*  3) Zero is returned if the two directions are the same or if either
*     vector is null.
*
*  4) If A is at a pole, the result is ill-defined.
*
*  Called:
*     iau_PN       decompose p-vector into modulus and direction
*     iau_PM       modulus of p-vector
*     iau_PXP      vector product of two p-vectors
*     iau_PMP      p-vector minus p-vector
*     iau_PDP      scalar product of two p-vectors
*
*  This revision:  2006 November 13
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION A(3), B(3), THETA

      DOUBLE PRECISION AM, AU(3), BM, ST, CT, XA, YA, ZA, ETA(3),
     :                 XI(3), A2B(3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Modulus and direction of the A vector.
      CALL iau_PN ( A, AM, AU )

*  Modulus of the B vector.
      CALL iau_PM ( B, BM )

*  Deal with the case of a null vector.
      IF ( AM.EQ.0D0 .OR. BM.EQ.0D0 ) THEN
         ST = 0D0
         CT = 1D0
      ELSE

*     The "north" axis tangential from A (arbitrary length).
         XA = A(1)
         YA = A(2)
         ZA = A(3)
         ETA(1) = - XA * ZA
         ETA(2) = - YA * ZA
         ETA(3) = XA*XA + YA*YA

*     The "east" axis tangential from A (same length).
         CALL iau_PXP ( ETA, AU, XI )

*     The vector from A to B.
         CALL iau_PMP ( B, A, A2B )

*     Resolve into components along the north and east axes.
         CALL iau_PDP ( A2B, XI, ST )
         CALL iau_PDP ( A2B, ETA, CT )

*     Deal with degenerate cases.
         IF ( ST.EQ.0D0 .AND. CT.EQ.0D0 ) CT = 1D0

      END IF

*  Position angle.
      THETA = ATAN2(ST,CT)

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
