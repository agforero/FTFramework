      SUBROUTINE iau_FK5HIP ( R5H, S5H )
*+
*  - - - - - - - - - - -
*   i a u _ F K 5 H I P
*  - - - - - - - - - - -
*
*  FK5 to Hipparcos rotation and spin.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Returned:
*     R5H     d(3,3)   r-matrix: FK5 rotation wrt Hipparcos (Note 2)
*     S5H     d(3)     r-vector: FK5 spin wrt Hipparcos (Note 3)
*
*  Notes:
*
*  1) This routine models the FK5 to Hipparcos transformation as a
*     pure rotation and spin;  zonal errors in the FK5 catalogue are
*     not taken into account.
*
*  2) The r-matrix R5H operates in the sense:
*
*           P_Hipparcos = R5H x P_FK5
*
*     where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is
*     the equivalent Hipparcos p-vector.
*
*  3) The r-vector S5H represents the time derivative of the FK5 to
*     Hipparcos rotation.  The units are radians per year (Julian,
*     TDB).
*
*  Called:
*     iau_RV2M     r-vector to r-matrix
*
*  Reference:
*
*     F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
*
*  This revision:  2017 October 12
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R5H(3,3), S5H(3)

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  FK5 to Hipparcos orientation and spin (radians, radians/year)
      DOUBLE PRECISION EPX, EPY, EPZ
      DOUBLE PRECISION OMX, OMY, OMZ

      PARAMETER ( EPX = -19.9D-3 * DAS2R,
     :            EPY =  -9.1D-3 * DAS2R,
     :            EPZ = +22.9D-3 * DAS2R )

      PARAMETER ( OMX = -0.30D-3 * DAS2R,
     :            OMY = +0.60D-3 * DAS2R,
     :            OMZ = +0.70D-3 * DAS2R )

      DOUBLE PRECISION V(3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  FK5 to Hipparcos orientation expressed as an r-vector.
      V(1) = EPX
      V(2) = EPY
      V(3) = EPZ

*  Re-express as an r-matrix.
      CALL iau_RV2M ( V, R5H )

*  Hipparcos wrt FK5 spin expressed as an r-vector.
      S5H(1) = OMX
      S5H(2) = OMY
      S5H(3) = OMZ

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
