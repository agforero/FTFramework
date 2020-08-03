FUNCTION hermite(npt, func, ier) RESULT(fn_val)

!       HERMITE INTEGRATION, I.E. EVALUATION OF THE INTEGRAL OF
!       exp(-x**2)*f(x) FROM -INFINITY TO +INFINITY

!       NPT     = INPUT, NO. OF POINTS AT WHICH f(x) IS TO BE EVALUATED.
!                 NPT MUST BE ONE OF 2, 4, 6, 8, 10, 12, 16, 20.
!       FUNC    = INPUT, NAME OF THE USER'S FUNCTION TO SUPPLY VALUES
!                 OF F(X).   THE FUNCTION MAY HAVE ANY NAME BUT ONLY THE
!                 ONE ARGUMENT, X.
!       IER     = OUTPUT, ERROR INDICATOR
!                   = 0 NO ERROR DETECTED
!                   = 1 ILLEGAL VALUE OF NPT

!       LATEST REVISION - 5 December 1999

!************************************************************************

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)  :: npt
INTEGER, INTENT(OUT) :: ier
REAL (dp)            :: fn_val

INTERFACE
  FUNCTION func(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fn_val
  END FUNCTION func
END INTERFACE

!     Local variables

REAL (dp) :: xpt(39) = (/ 0.707106781186548_dp, &
               0.524647623275290_dp, 1.650680123885785_dp, &
               0.436077411927617_dp, 1.335849074013697_dp, 2.350604973674492_dp, &
               0.381186990207322_dp, 1.157193712446780_dp, 1.981656756695843_dp, &
               2.930637420257244_dp, &
               0.342901327223705_dp, 1.036610829789514_dp, 1.756683649299882_dp, &
               2.532731674232790_dp, 3.436159118837738_dp, &
               0.314240376254359_dp, 0.947788391240164_dp, 1.597682635152605_dp, &
               2.279507080501060_dp, 3.020637025120890_dp, 3.889724897869782_dp, &
               0.27348104613815_dp,  0.82295144914466_dp, 1.38025853919888_dp, &
               1.95178799091625_dp,  2.54620215784748_dp, 3.17699916197996_dp, &
               3.86944790486012_dp,  4.68873893930582_dp, &
               0.2453407083009_dp,   0.7374737285454_dp, 1.2340762153953_dp, &
               1.7385377121166_dp,   2.2549740020893_dp, 2.7888060584281_dp, &
               3.3478545673832_dp,   3.9447640401156_dp, 4.6036824495507_dp, &
               5.3874808900112_dp /),                                    &
   wt(39) = (/ 8.862269254528D-1, &
               8.049140900055D-1, 8.131283544725D-2, &
               7.246295952244D-1, 1.570673203229D-1, 4.530009905509D-3, &
               6.611470125582D-1, 2.078023258149D-1, 1.707798300741D-2, &
               1.996040722114D-4, &
               6.108626337353D-1, 2.401386110823D-1, 3.387439445548D-2, &
               1.343645746781D-3, 7.640432855233D-6, &
               5.701352362625D-1, 2.604923102642D-1, 5.160798561588D-2, &
               3.905390584629D-3, 8.573687043588D-5, 2.658551684356D-7, &
               5.079294790166D-1, 2.806474585285D-1, 8.381004139899D-2, &
               1.288031153551D-2, 9.322840086242D-4, 2.711860092538D-5, &
               2.320980844865D-7, 2.654807474011D-10, &
               4.622436696006D-1, 2.866755053628D-1, 1.090172060200D-1, &
               2.481052088746D-2, 3.243773342238D-3, 2.283386360163D-4, &
               7.802556478532D-6, 1.086069370769D-7, 4.399340992273D-10, &
               2.229393645534D-13 /), zero = 0.0_dp, x
INTEGER  :: npts(8) = (/ 2, 4, 6, 8, 10, 12, 16, 20 /), ipos, i, nby2

!       CHECK FOR PERMISSIBLE VALUE OF NPT.

fn_val = zero
ipos = 1
DO i = 1, 8
  IF(npts(i).EQ.npt) GO TO 20
  ipos = ipos + (npts(i) + 1)/2
END DO
ier = 1
RETURN

!       EVALUATE SUM OF WT(I) * FUNC(X(I))

20 nby2 = (npt + 1)/2
DO i = 1, nby2
  x = xpt(ipos)
  fn_val = fn_val + wt(ipos) * (func(x) + func(-x))
  ipos = ipos + 1
END DO
ier = 0
RETURN

END FUNCTION hermite
