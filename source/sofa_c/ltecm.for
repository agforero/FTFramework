      SUBROUTINE iau_LTECM ( EPJ, RM )
*+
*  - - - - - - - - - -
*   i a u _ L T E C M
*  - - - - - - - - - -
*
*  ICRS equatorial to ecliptic rotation matrix, long-term.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     EPJ       d        Julian epoch (TT)
*
*  Returned:
*     RM        d(3,3)   ICRS to ecliptic rotation matrix
*
*  Notes:
*
*  1) The matrix is in the sense
*
*        E_ep = RM x P_ICRS,
*
*     where P_ICRS is a vector with respect to ICRS right ascension
*     and declination axes and E_ep is the same vector with respect to
*     the (inertial) ecliptic and equinox of epoch EPJ.
*
*  2) P_ICRS is a free vector, merely a direction, typically of unit
*     magnitude, and not bound to any particular spatial origin, such as
*     the Earth, Sun or SSB.  No assumptions are made about whether it
*     represents starlight and embodies astrometric effects such as
*     parallax or aberration.  The transformation is approximately that
*     between mean J2000.0 right ascension and declination and ecliptic
*     longitude and latitude, with only frame bias (always less than
*     25 mas) to disturb this classical picture.
*
*  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
*     agrees with the IAU 2006 precession at J2000.0 and stays within
*     100 microarcseconds during the 20th and 21st centuries.  It is
*     accurate to a few arcseconds throughout the historical period,
*     worsening to a few tenths of a degree at the end of the
*     +/- 200,000 year time span.
*
*  Called:
*     iau_LTPEQU   equator pole, long term
*     iau_LTPECL   ecliptic pole, long term
*     iau_PXP      vector product
*     iau_PN       normalize vector
*
*  References:
*
*     Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
*     expressions, valid for long time intervals, Astron.Astrophys. 534,
*     A22
*
*     Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
*     expressions, valid for long time intervals (Corrigendum),
*     Astron.Astrophys. 541, C1
*
*  This revision:  2015 December 6
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION EPJ, RM(3,3)

*  Arcseconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER ( AS2R = 4.848136811095359935899141D-6 )

*  Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33)
      DOUBLE PRECISION DX, DE, DR
      PARAMETER ( DX = -0.016617D0 * AS2R,
     :            DE = -0.0068192D0 * AS2R,
     :            DR = -0.0146D0 * AS2R )

      DOUBLE PRECISION P(3), Z(3), W(3), S, X(3), Y(3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


*  Equator pole.
      CALL iau_LTPEQU ( EPJ, P )

*  Ecliptic pole (bottom row of equatorial to ecliptic matrix).
      CALL iau_LTPECL ( EPJ, Z )

*  Equinox (top row of matrix).
      CALL iau_PXP ( P, Z, W )
      CALL iau_PN ( W, S, X )

*  Middle row of matrix.
      CALL iau_PXP ( Z, X, Y )

*  Combine with frame bias.
      RM(1,1) =   X(1)    - X(2)*DR + X(3)*DX
      RM(1,2) =   X(1)*DR + X(2)    + X(3)*DE
      RM(1,3) = - X(1)*DX - X(2)*DE + X(3)
      RM(2,1) =   Y(1)    - Y(2)*DR + Y(3)*DX
      RM(2,2) =   Y(1)*DR + Y(2)    + Y(3)*DE
      RM(2,3) = - Y(1)*DX - Y(2)*DE + Y(3)
      RM(3,1) =   Z(1)    - Z(2)*DR + Z(3)*DX
      RM(3,2) =   Z(1)*DR + Z(2)    + Z(3)*DE
      RM(3,3) = - Z(1)*DX - Z(2)*DE + Z(3)

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
