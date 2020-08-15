      SUBROUTINE iau_FK54Z ( R2000, D2000, BEPOCH,
     :                       R1950, D1950, DR1950, DD1950 )
*+
*  - - - - - - - - - -
*   i a u _ F K 5 4 Z
*  - - - - - - - - - -
*
*  Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
*  proper motion in FK5 and parallax.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     R2000,D2000      d     J2000.0 FK5 RA,Dec (rad)
*     BEPOCH           d     Besselian epoch (e.g. 1950D0)
*
*  Returned:
*     R1950,D1950      d     B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH
*     DR1950,DD1950    d     B1950.0 FK4 proper motions (rad/trop.yr)
*
*  Notes:
*
*  1) In contrast to the iau_FK524 routine, here the FK5 proper motions,
*     the parallax and the radial velocity are presumed zero.
*
*  2) This routine converts a star position from the IAU 1976 FK5
*     (Fricke) system to the former FK4 (Bessel-Newcomb) system, for
*     cases such as distant radio sources where it is presumed there is
*     zero parallax and no proper motion.  Because of the E-terms of
*     aberration, such objects have (in general) non-zero proper motion
*     in FK4, and the present routine returns those fictitious proper
*     motions.
*
*  3) Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
*     Conversions involving other equinoxes would require additional
*     treatment for precession.
*
*  4) The position returned by this routine is in the B1950.0 FK4
*     reference system but at Besselian epoch BEPOCH.  For comparison
*     with catalogs the BEPOCH argument will frequently be 1950D0.  (In
*     this context the distinction between Besselian and Julian epoch is
*     insignificant.)
*
*  5) The RA component of the returned (fictitious) proper motion is
*     dRA/dt rather than cos(Dec)*dRA/dt.
*
*  Called:
*     iau_ANP      normalize angle into range 0 to 2pi
*     iau_C2S      p-vector to spherical
*     iau_FK524    FK4 to FK5
*     iau_S2C      spherical to p-vector
*
*  This revision:   2018 January 11
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R2000, D2000, BEPOCH,
     :                 R1950, D1950, DR1950, DD1950

      DOUBLE PRECISION R, D, PR, PD, PX, RV, P(3), W, V(3)
      INTEGER I

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  FK5 equinox J2000.0 to FK4 equinox B1950.0.
      CALL iau_FK524 ( R2000, D2000, 0D0, 0D0, 0D0, 0D0,
     :                 R, D, PR, PD, PX, RV )

*  Spherical to Cartesian.
      CALL iau_S2C ( R, D, P )

*  Fictitious proper motion (radians per year).
      V(1) = - PR*P(2) - PD*COS(R)*SIN(D)
      V(2) =   PR*P(1) - PD*SIN(R)*SIN(D)
      V(3) =             PD*COS(D)

*  Apply the motion.
      W = BEPOCH - 1950D0
      DO 1 I=1,3
         P(I) = P(I) + W*V(I)
 1    CONTINUE

*  Cartesian to spherical.
      CALL iau_C2S ( P, W, D1950 )
      R1950 = iau_ANP ( W )

*  Fictitious proper motion.
      DR1950 = PR
      DD1950 = PD

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
