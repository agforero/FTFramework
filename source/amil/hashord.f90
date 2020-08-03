PROGRAM hashord

!===SKIP, THIS IS THE ORDERED HASH PROGRAM I WAS TALKING ABOUT
!   AND PROMISED YOU A WHILE BACK... WELL HERE IT IS, HAVE FUN!!
!   S. W. WHARTON 2/5/83

! Code converted using TO_F90 by Alan Miller
! Date: 2000-03-26  Time: 09:35:46
 
IMPLICIT NONE
INTEGER :: i, ncell, nentry, nprobe, npts, nstat, nswap, rndseq(10000),  &
           tswap, tprobe, tfound

INTERFACE
  SUBROUTINE orhash(ncell, ENTRY, mode, STATUS, nprobe, nswap)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: ncell
    INTEGER, INTENT(IN)   :: ENTRY
    INTEGER, INTENT(IN)   :: mode
    INTEGER, INTENT(OUT)  :: STATUS
    INTEGER, INTENT(OUT)  :: nprobe
    INTEGER, INTENT(OUT)  :: nswap
  END SUBROUTINE orhash
END INTERFACE

nentry = 450
ncell = 499
npts = 10000
CALL irand(rndseq, 10000, -30000, 30000)
tfound = 0
tprobe = 0
tswap = 0

!     ...ENTER NENTRY POINTS INTO THE ORDERED HASH TABLE
DO  i = 1, nentry
  CALL orhash(ncell, rndseq(i), 2, nstat, nprobe, nswap)
  IF (nstat /= 3) tfound = tfound + 1
  tswap = tswap + nswap
  tprobe = tprobe + nprobe
END DO
WRITE(6,600) nentry, tfound, tprobe, tswap
tfound = 0
tprobe = 0
tswap = 0

!     ...NOW LOOK UP NPTS POINTS IN THE HASH TABLE
!     ...THE ORDERED HASH IS MORE EFFICIENT FOR UNSUCCESSFUL SEARCHES
DO  i = 1, npts
  CALL orhash(ncell, rndseq(i), 1, nstat, nprobe, nswap)
  IF (nstat /= 3) tfound = tfound + 1
  tswap = tswap + nswap
  tprobe = tprobe + nprobe
END DO
WRITE(6,600) npts, tfound, tprobe, tswap
STOP
600 FORMAT(' N=', i10, '     NFOUND =', i10, '     TPROBES =', i10,  &
           '     TSWAPS =', i10)

CONTAINS


SUBROUTINE irand(randseq, n, low, high)
INTEGER, INTENT(OUT)  :: randseq(:)
INTEGER, INTENT(IN)   :: n, low, high

REAL    :: ran(n)

CALL RANDOM_NUMBER(ran)
randseq = low + (high + 1 - low) * ran

RETURN
END SUBROUTINE irand

END PROGRAM hashord



SUBROUTINE orhash(ncell, ENTRY, mode, STATUS, nprobe, nswap)
!=
!=    USE A HASH FUNCTION TO LOOK UP OR INSERT KEYS INTO A TABLE.
!=    ORHASH USES THE ORDERED HASHING ALGORITHM WITH PASS BITS BY AMBLE AND
!=    KNUTH, (1974) ORDERED HASH TABLES, COMPUTER JOURNAL 17(2):135-142.
!=
!=    ORDERED HASHING REDUCES THE NUMBER OF PROBES FOR UNSUCCESSFUL
!=    SEARCHING AND SO REDUCES THE SCANNER CPU TIME.
!=
!=    VARIABLE    DESCRIPTION
!=    --------    -----------
!=    NCELL       I - HASH TABLE LENGTH (PRIME NUMBER)
!=    NDIM        I - LENGTH OF EACH KEY
!=    ENTRY       I - NUMBER TO BE ENTERED OR PROBED FOR IN HASH TABLE
!=       NOTE: BECAUSE ORHASH MAY SWAP KEYS, ENTRY IS COPIED INTO KEY
!=             FOR SEARCHING.  NEGATIVE AND ZERO VALUES OF KEY ARE OK.
!=    MODE        I - IS GIVEN AS AN INDICATOR AS FOLLOWS...
!=          (1) - MEANS TO ONLY LOOK UP THE KEY IN TABLE
!=          (2) - MEANS LOOK UP AND ENTER KEY IF IT IS NOT YET IN TABLE
!=    STATUS      I - FINAL STATUS OF THE PROBE SEQUENCE AS FOLLOWS:
!=          (1) - THE KEY WAS ADDED TO THE HASH TABLE AS A NEW ENTRY
!=          (2) - THE KEY MATCHED AN ENTRY IN THE HASH TABLE
!=          (3) - THE KEY WAS NOT ENTERED NOR WAS A MATCH FOUND
!=    NPROBE      I - THE NUMBER OF PROBES NEEDED TO ENTER THE KEY.
!=    NSWAP       I - THE NUMBER OF SWAPS NEEDED TO ENTER KEY
!=
!=    COMMON BLOCKS:    (NONE)
!=    I/O DEVICE #'S:  INPUT - (NONE);  OUTPUT - (NONE)
!=

IMPLICIT NONE

INTEGER, INTENT(IN)   :: ncell
INTEGER, INTENT(IN)   :: ENTRY
INTEGER, INTENT(IN)   :: mode
INTEGER, INTENT(OUT)  :: STATUS
INTEGER, INTENT(OUT)  :: nprobe
INTEGER, INTENT(OUT)  :: nswap

INTEGER :: afp, i, ikey, iknt, inc, ka, key, keytab(5000), newknt

!     ...FOR ORHASH TO WORK PROPERLY, THE ARRAYS BINKNT AND PASSBT
!     ...MUST BE INITIALIZED TO ALL ZEROES AND FALSE, RESPECTIVELY.
! DATA binknt/005000*0/, passbt/005000*.false./
INTEGER, SAVE :: binknt(5000) = (/ (0, i=1,5000) /)
LOGICAL, SAVE :: passbt(5000) = (/ (.FALSE., i=1,5000) /)

key = ENTRY
newknt = 1
nprobe = 1
nswap = 0
STATUS = 3

!     ...COMPUTE ADDRESS OF FIRST PROBE (AFP) AND INCREMENT (INC).
afp = ABS(MOD(key,ncell))  + 1
inc = ABS(MOD(key,ncell-2))  + 1
ka = afp

!     ...BEGIN THE SEARCH LOOP TO FIND AN ENTRY IN THE ORDERED
!     ...HASH TABLE (KEYTAB) THAT MATCHES THE NEW KEY

!     ...EXIT SEARCH LOOP IF EMPTY LOCATION FOUND
10 IF( binknt(ka) == 0 ) GO TO 100

!        ...CHECK IF KEYTAB(KA) MATCHES THE KEY
IF( key /= keytab(ka) ) GO TO 20

!           ...A MATCH WAS FOUND, INCREMENT BIN COUNTER AND EXIT LOOP
binknt(ka) = binknt(ka) + 1
STATUS = 2
GO TO 999

!        ...EXIT SEARCH LOOP IF NO SMALLER ENTRIES FOLLOW KEYTAB(KA)
20 IF ( .NOT. passbt(ka) ) GO TO 100

!           ...EXIT SEARCH LOOP IF KEYTAB(KA) IS SMALLER THAN KEY
IF ( key > keytab(ka) ) GO TO 100

!              ...COMPUTE ADDRESS OF NEXT PROBE
ka = ka - inc
IF( ka < 1 ) ka = ka + ncell
nprobe = nprobe + 1
GO TO 10

!     ...BEGIN THE INSERTION LOOP TO ENTER THE NEW KEY IN THE
!     ...ORDERED HASH TABLE IN THE APPROPRIATE POSITION SUCH
!     ...THAT THE KEYS ARE ARRANGED IN DECREASING ORDER

!     ...TERMINATE HASHING IF NO ENTRIES ARE TO BE MADE (MODE=1)
100 IF( mode == 1 ) GO TO 999
110 IF( binknt(ka) /= 0 ) GO TO 120

!        ...LOAD KEY INTO EMPTY SPACE AT POSITION KA
binknt(ka) = newknt
keytab(ka) = key
STATUS = 1
GO TO 999

120 IF( key <= keytab(ka) ) GO TO 130

!        ...SWAP KEY WITH KEYTAB(KA) SO THAT THE KEYS ALONG A GIVEN
!        ...SEARCH SEQUENCE ARE ARRANGED IN DECREASING ORDER
ikey = key
key = keytab(ka)
keytab(ka) = ikey
iknt = newknt
newknt = binknt(ka)
binknt(ka) = iknt
nswap = nswap + 1
inc = ABS(MOD(key,ncell-2)) + 1
130   passbt(ka) = .true.

!     ...COMPUTE ADDRESS OF NEXT PROBE.
ka = ka - inc
IF( ka < 1 ) ka = ka + ncell
nprobe = nprobe + 1
GO TO 110

999 RETURN
END SUBROUTINE orhash
