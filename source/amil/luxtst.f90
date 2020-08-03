PROGRAM luxtst
!         Exercise for the RANLUX Pseudorandom number generator.

USE luxury

IMPLICIT NONE

REAL    :: rvec(1000)
INTEGER :: i1, i2, i3, i4, li

!         check that we get the right numbers (machine-indep.)
WRITE (6,'(/A)') '  CALL RANLUX(RVEC,100)'
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers 101-105:', rvec(1:5)

WRITE (6,'(/A)') ' CALL RLUXGO(0,0,0,0)'
CALL rluxgo(0,0,0,0)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0,   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0, 101-105:', rvec(1:5)

WRITE (6,'(/A)') '   CALL RLUXGO(389,1,0,0)'
CALL rluxgo(389,1,0,0)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389,   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389, 101-105:', rvec(1:5)

WRITE (6,'(/A)') '  CALL RLUXGO(75,0,0,0)'
CALL rluxgo(75,0,0,0)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75,   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75, 101-105:', rvec(1:5)

WRITE (6,'(/A)') '  test restarting from the full vector'
CALL rluxut
WRITE (6,'(/A/(1X,5I14))') '  current RANLUX status saved:', isdext
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:', rvec(1:5)

WRITE (6,'(/A)') '   previous RANLUX status will be restored'
CALL rluxin
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:', rvec(1:5)

WRITE (6,'(/A)') '     test the restarting by skipping'
CALL rluxgo(4,7674985,0,0)
CALL rluxat(i1,i2,i3,i4)
WRITE (6,'(A,4I10)') '  RLUXAT values =', i1, i2, i3, i4
DO li = 1, 10
  CALL ranlux(rvec,1000)
END DO
CALL rluxat(i1,i2,i3,i4)
WRITE (6,'(A,4I10)') '  RLUXAT values =', i1, i2, i3, i4
CALL ranlux(rvec,200)
WRITE (6,'(A,2F10.6)') '  Next and 200th numbers are:', rvec(1), rvec(200)
CALL rluxgo(i1,i2,i3,i4)
CALL ranlux(rvec,200)
WRITE (6,'(A,2F10.6)') '  Next and 200th numbers are:', rvec(1), rvec(200)

WRITE (6,'(/A)') ' The following should provoke an error message'
CALL rluxgo(4,11111,31,0)
STOP

!   OUTPUT FROM THE ABOVE TEST PROGRAM SHOULD BE:
!   --------------------------------------------
!  CALL RANLUX(RVEC,100)
! RANLUX DEFAULT INITIALIZATION:    314159265
! RANLUX DEFAULT LUXURY LEVEL =   3      p = 223
! RANLUX default numbers   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX default numbers 101-105:
!           0.43156743  0.03774416  0.24897110  0.00147784  0.90274453

!  CALL RLUXGO(0,0,0,0)
! RANLUX LUXURY LEVEL SET BY RLUXGO : 0     P=  24
! RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
! RANLUX luxury level 0,   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX luxury level 0, 101-105:
!           0.41538775  0.05330932  0.58195311  0.91397446  0.67034441

!   CALL RLUXGO(389,1,0,0)
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS           1           0           0
! RANLUX luxury p=389,   1-  5:
!           0.94589490  0.47347850  0.95152789  0.42971975  0.09127384
! RANLUX luxury p=389, 101-105:
!           0.02618265  0.03775346  0.97274780  0.13302165  0.43126065

!  CALL RLUXGO(75,0,0,0)
! RANLUX P-VALUE SET BY RLUXGO TO:   75
! RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
! RANLUX luxury p= 75,   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX luxury p= 75, 101-105:
!           0.25600731  0.23443210  0.59164381  0.59035838  0.07011414

!  test restarting from the full vector

!  current RANLUX status saved:
!       16156027      16534309      15243811       2751687       6002207
!        7979506       1301976       4567313       4305996       5872599
!       12003090       2146823      12606367       4111505       5979640
!       12739666      10489318      14036909      11729352       8061448
!        7832659       6069758       3197719       1832730      75080216
! RANLUX numbers 1- 5:
!           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
! RANLUX numbers 101-105:
!           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233

!   previous RANLUX status will be restored
! FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:
!         16156027    16534309    15243811     2751687     6002207
!          7979506     1301976     4567313     4305996     5872599
!         12003090     2146823    12606367     4111505     5979640
!         12739666    10489318    14036909    11729352     8061448
!          7832659     6069758     3197719     1832730    75080216
! RANLUX P-VALUE SET BY RLUXIN TO:   75
! RANLUX numbers 1- 5:
!           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
! RANLUX numbers 101-105:
!           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233

!     test the restarting by skipping
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985           0           0
!  RLUXAT values =         4   7674985         0         0
!  RLUXAT values =         4   7674985    161840         0
!  Next and 200th numbers are:  0.019648  0.590586
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985      161840           0
!  Next and 200th numbers are:  0.019648  0.590586

! The following should provoke an error message
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS       11111          31           0
!  Error in RESTARTING with RLUXGO:
!  The values      11111         31          0 cannot occur at luxury level    4
END PROGRAM luxtst
