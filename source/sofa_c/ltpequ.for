      SUBROUTINE iau_LTPEQU ( EPJ, VEQ )
*+
*  - - - - - - - - - - -
*   i a u _ L T P E Q U
*  - - - - - - - - - - -
*
*  Long-term precession of the equator.
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
*     VEQ       d(3)     equator pole unit vector
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
      DOUBLE PRECISION EPJ, VEQ(3)

*  Arcseconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER ( AS2R = 4.848136811095359935899141D-6 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Number of polynomial terms
      INTEGER NPOL
      PARAMETER ( NPOL = 4 )

*  Number of periodic terms
      INTEGER NPER
      PARAMETER ( NPER = 14 )

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION T, X, Y, W, A, S, C

*  Polynomial and periodic coefficients
      DOUBLE PRECISION XYPOL(NPOL,2), XYPER(5,NPER)

*  Polynomials
      DATA ((XYPOL(I,J),I=1,NPOL),J=1,2) /
     :
     :       +5453.282155    D0,
     :          +0.4252841   D0,
     :          -0.00037173  D0,
     :          -0.000000152 D0,
     :
     :      -73750.930350    D0,
     :          -0.7675452   D0,
     :          -0.00018725  D0,
     :          +0.000000231 D0 /

*  Periodics
       DATA ((XYPER(I,J),I=1,5),J=1,NPER) /
     :
     :   256.75D0,   -819.940624D0,  75004.344875D0,
     :              81491.287984D0,   1558.515853D0,
     :   708.15D0,  -8444.676815D0,    624.033993D0,
     :                787.163481D0,   7774.939698D0,
     :   274.20D0,   2600.009459D0,   1251.136893D0,
     :               1251.296102D0,  -2219.534038D0,
     :   241.45D0,   2755.175630D0,  -1102.212834D0,
     :              -1257.950837D0,  -2523.969396D0,
     :  2309.00D0,   -167.659835D0,  -2660.664980D0,
     :              -2966.799730D0,    247.850422D0,
     :   492.20D0,    871.855056D0,    699.291817D0,
     :                639.744522D0,   -846.485643D0,
     :   396.10D0,     44.769698D0,    153.167220D0,
     :                131.600209D0,  -1393.124055D0,
     :   288.90D0,   -512.313065D0,   -950.865637D0,
     :               -445.040117D0,    368.526116D0,
     :   231.10D0,   -819.415595D0,    499.754645D0,
     :                584.522874D0,    749.045012D0,
     :  1610.00D0,   -538.071099D0,   -145.188210D0,
     :                -89.756563D0,    444.704518D0,
     :   620.00D0,   -189.793622D0,    558.116553D0,
     :                524.429630D0,    235.934465D0,
     :   157.87D0,   -402.922932D0,    -23.923029D0,
     :                -13.549067D0,    374.049623D0,
     :   220.30D0,    179.516345D0,   -165.405086D0,
     :               -210.157124D0,   -171.330180D0,
     :  1200.00D0,     -9.814756D0,      9.344131D0,
     :                -44.919798D0,    -22.899655D0 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


*  Centuries since J2000.
      T = (EPJ-2000D0)/100D0

*  Initialize X and Y accumulators.
      X = 0D0
      Y = 0D0

*  Periodic terms.
      W = D2PI*T
      DO 1 I=1,NPER
         A = W/XYPER(1,I)
         S = SIN(A)
         C = COS(A)
         X = X + C*XYPER(2,I) + S*XYPER(4,I)
         Y = Y + C*XYPER(3,I) + S*XYPER(5,I)
 1    CONTINUE

*  Polynomial terms.
      W = 1D0
      DO 2 I=1,NPOL
         X = X + XYPOL(I,1)*W
         Y = Y + XYPOL(I,2)*W
         W = W*T
 2    CONTINUE

*  X and Y (direction cosines).
      X = X*AS2R
      Y = Y*AS2R

*  Form the equator pole vector.
      VEQ(1) = X
      VEQ(2) = Y
      VEQ(3) = SQRT(MAX(1D0-X*X-Y*Y,0D0))

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
