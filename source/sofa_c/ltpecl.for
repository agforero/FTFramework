      SUBROUTINE iau_LTPECL ( EPJ, VEC )
*+
*  - - - - - - - - - - -
*   i a u _ L T P E C L
*  - - - - - - - - - - -
*
*  Long-term precession of the ecliptic.
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
*     VEC       d(3)     ecliptic pole unit vector
*
*  Notes:
*
*  1) The returned vector is with respect to the J2000.0 mean equator
*     and equinox.
*
*  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
*     agrees with the IAU 2006 precession at J2000.0 and stays within
*     100 microarcseconds during the 20th and 21st centuries.  It is
*     accurate to a few arcseconds throughout the historical period,
*     worsening to a few tenths of a degree at the end of the
*     +/- 200,000 year time span.
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
*  This revision:  2016 February 9
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION EPJ, VEC(3)

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Obliquity at J2000.0 (radians).
      DOUBLE PRECISION EPS0
      PARAMETER ( EPS0 = 84381.406D0 * DAS2R )

*  Number of polynomial terms
      INTEGER NPOL
      PARAMETER ( NPOL = 4 )

*  Number of periodic terms
      INTEGER NPER
      PARAMETER ( NPER = 8 )

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION T, P, Q, W, A, S, C

*  Polynomial and periodic coefficients
      DOUBLE PRECISION PQPOL(NPOL,2), PQPER(5,NPER)

*  Polynomials
      DATA ((PQPOL(I,J),I=1,NPOL),J=1,2) /
     :
     :       +5851.607687    D0,
     :          -0.1189000   D0,
     :          -0.00028913  D0,
     :          +0.000000101 D0,
     :
     :       -1600.886300    D0,
     :          +1.1689818   D0,
     :          -0.00000020  D0,
     :          -0.000000437 D0 /

*  Periodics
       DATA ((PQPER(I,J),I=1,5),J=1,NPER) /
     :
     :   708.15D0,  -5486.751211D0,   -684.661560D0,
     :                667.666730D0,  -5523.863691D0,
     :  2309.00D0,    -17.127623D0,   2446.283880D0,
     :              -2354.886252D0,   -549.747450D0,
     :  1620.00D0,   -617.517403D0,    399.671049D0,
     :               -428.152441D0,   -310.998056D0,
     :   492.20D0,    413.442940D0,   -356.652376D0,
     :                376.202861D0,    421.535876D0,
     :  1183.00D0,     78.614193D0,   -186.387003D0,
     :                184.778874D0,    -36.776172D0,
     :   622.00D0,   -180.732815D0,   -316.800070D0,
     :                335.321713D0,   -145.278396D0,
     :   882.00D0,    -87.676083D0,    198.296701D0,
     :               -185.138669D0,    -34.744450D0,
     :   547.00D0,     46.140315D0,    101.135679D0,
     :               -120.972830D0,     22.885731D0 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


*  Centuries since J2000.
      T = (EPJ-2000D0)/100D0

*  Initialize P_A and Q_A accumulators.
      P = 0D0
      Q = 0D0

*  Periodic terms.
      W = D2PI*T
      DO 1 I=1,NPER
         A = W/PQPER(1,I)
         S = SIN(A)
         C = COS(A)
         P = P + C*PQPER(2,I) + S*PQPER(4,I)
         Q = Q + C*PQPER(3,I) + S*PQPER(5,I)
 1    CONTINUE

*  Polynomial terms.
      W = 1D0
      DO 2 I=1,NPOL
         P = P + PQPOL(I,1)*W
         Q = Q + PQPOL(I,2)*W
         W = W*T
 2    CONTINUE

*  P_A and Q_A (radians).
      P = P*DAS2R
      Q = Q*DAS2R

*  Form the ecliptic pole vector.
      W = SQRT(MAX(1D0-P*P-Q*Q,0D0))
      S = SIN(EPS0)
      C = COS(EPS0)
      VEC(1) = P
      VEC(2) = - Q*C - W*S
      VEC(3) = - Q*S + W*C

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
