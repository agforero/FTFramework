      SUBROUTINE iau_TCGTT ( TCG1, TCG2, TT1, TT2, J )
*+
*  - - - - - - - - - -
*   i a u _ T C G T T
*  - - - - - - - - - -
*
*  Time scale transformation:  Geocentric Coordinate Time, TCG, to
*  Terrestrial Time, TT.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical.
*
*  Given:
*     TCG1,TCG2    d      TCG as a 2-part Julian Date
*
*  Returned:
*     TT1,TT2      d      TT as a 2-part Julian Date
*     J            i      status:  0 = OK
*
*  Note:
*
*     TCG1+TCG2 is Julian Date, apportioned in any convenient way
*     between the two arguments, for example where TCG1 is the Julian
*     Day Number and TCG2 is the fraction of a day.  The returned
*     TT1,TT2 follow suit.
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
*     IERS Technical Note No. 32, BKG (2004)
*
*     IAU 2000 Resolution B1.9
*
*  This revision:  2019 June 20
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION TCG1, TCG2, TT1, TT2
      INTEGER J

*  JD for MJD 0
      DOUBLE PRECISION DJM0
      PARAMETER (DJM0 = 2400000.5D0 )

*  1977 Jan 1 00:00:32.184 TT, as MJD
      DOUBLE PRECISION T77T
      PARAMETER ( T77T = 43144.0003725D0 )

*  L_G = 1 - dTT/dTCG
      DOUBLE PRECISION ELG
      PARAMETER ( ELG = 6.969290134D-10 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Result, safeguarding precision.
      IF ( ABS(TCG1).GT.ABS(TCG2) ) THEN
         TT1 = TCG1
         TT2 = TCG2 - ( ( TCG1-DJM0 ) + ( TCG2-T77T ) ) * ELG
      ELSE
         TT1 = TCG1 - ( ( TCG2-DJM0 ) + ( TCG1-T77T ) ) * ELG
         TT2 = TCG2
      END IF

*  Status (always OK).
      J = 0

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
