MODULE common_Brent
! COMMON /brent/ LEN, keytab(73), kvtab(73)

IMPLICIT NONE
INTEGER, SAVE :: LEN, keytab(73), kvtab(73)

END MODULE common_Brent



SUBROUTINE hash(key, mode, ka, found)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-03-26  Time: 09:35:42
 
!===============BRENT'S HASHING ALGORITHM.
!               SEE CACM 16(2):105-109.

!       THIS ROUTINE WILL LOOK UP AND INSERT OR DELETE
!       INTEGER KEYS IN THE TABLE KEYTAB.  VALUES ASSOCIATED
!       WITH THE KEYS MAY BE STORED IN THE TABLE KVTAB. THE
!       ROUTINE IS DESIGNED TO BE EFFICIENT, EVEN IF THE TABLE
!       IS NEARLY FULL, PROVIDED MOST ENTRIES ARE LOOKED UP
!       SEVERAL TIMES.
!  PARAMETERS:
!  /BRENT/ COMMON BLOCK DESCRIBES:
!          LEN:    THE TABLE LENGTH.  MUST BE GIVEN AS A PRIME NUMBER
!                  BEFORE HASH IS CALLED.  NO CHECK IS MADE FOR PRIMALITY

!          KEYTAB: IS GIVEN AS A LEN-VECTOR FOR STORING THE KEYS.  IT MUST BE
!                  INITIALIZED TO ALL ZEROS BEFORE THE FIRST CALL TO HASH.

!          KVTAB:  IS USED BY BOTH HASH AND THE CALLING PROGRAM AS A LEN-VECTOR
!                  OF VALUES ASSOCIATED WITH KEYS IN THE HASH TABLE, KEYTAB.
!                  KVTAB COULD ALSO BE A POINTER VECTOR INITIALIZED TO THE
!                  INTEGERS 1,2,...LEN IF MULTIPLE VALUES ARE ASSOCIATED WITH
!                  EACH KEY.

!  (KEY, MODE, KA, FOUND):
!          KEY:    IS GIVEN AS A POSITIVE INTEGER KEY TO BE ENTERED,
!                  DELETED, OR PROBED FOR.
!                  THE VALUES 0 AND -1 ARE RESERVED FOR EMPTY SPACE
!                  OR DELETED ITEM RESPECTIVELY.

!          MODE:   IS GIVEN AS AN INDICATOR AS FOLLOWS...
!          MODE=1: MEANS LOOK UP ONLY.  IF THE KEY IS IN KEYTAB, THEN
!                  KA IS RETURNED AS THE SUBSCRIPT POINTER; OTHERWISE
!                  KA IS RETURNED 0.  THE VALUE ASSOCIATED WITH
!                  KEYTAB(KA) IS KVTAB(KA) WHEN KA IS NOT ZERO.
!          MODE=2: MEANS LOOK UP AND ENTER.  THIS IS USED TO BUILD THE TABLE.
!                  OTHERWISE THE KEY IS ENTERED AT KEYTAB(KA).
!                  KA IS RETURNED AS THE POINTER SUBSCRIPT IN EITHER CASE,
!                  OR KA=0 IF AN ENTRY IS ATTEMPTED AND THE TABLE IS FULL OR
!                  KA=0 IF KEY=0.  UPON RETURN, WHEN FOUND=.FALSE., THE CALLING
!                  PROGRAM MUST ENTER THE VALUE ASSOCIATED WITH KEY WHEN MODE=2.
!          MODE=3: MEANS LOOK UP AND DELETE.  IF THE KEY IS IN THE TABLE
!                  IT IS DELETED AND ITS FORMER ADDRESS KA IS RETURNED.
!                  IF THE KEY IS NOT THERE, KA=0 IS RETURNED.  (THE SUBROUTINE
!                  CAN BE SIMPLIFIED CONSIDERABLY IF KEYS NEVER DELETED).

!          KA:     IS RETURNED AS THE SUBSCRIPT OF KEYTAB SUCH THAT
!                  KEY=KEYTAB(KA) OR KA IS RETURNED AS ZERO IF KEY=0
!                  OR KEY IS NOT IN KEYTAB.

!          FOUND:  IS A LOGICAL FLAG RETURNED .TRUE. IF THE KEY WAS
!                  FOUND IN KEYTAB; ELSE FOUND IS RETURNED .FALSE. .

USE common_Brent
IMPLICIT NONE

INTEGER, INTENT(IN)   :: key
INTEGER, INTENT(IN)   :: mode
INTEGER, INTENT(OUT)  :: ka
LOGICAL, INTENT(OUT)  :: found

LOGICAL :: del
INTEGER :: ia, ic, iq, ir, is, ix, jq, jr, kt, len2

!---------THE COMMON STATEMENT BELOW MUST BE CHANGED TO MEET THE PROBLEM SIZE.
!         It has been replaced by saved variables in module common_brent
! COMMON /brent/ LEN, keytab(73), kvtab(73)

len2 = LEN-2
ic = -1
!--------COMPUTE ADDRESS OF FIRST PROBE(IR) AND INCREMENT(IQ).
!        ANY INDEPENDENT PSEUDO-RANDOM FUNCTIONS OF THE KEY MAY BE USED,
!        PROVIDED 0 < IQ < LEN AND 0 < IR <= LEN .
iq = MOD(ABS(key),len2) + 1
ir = MOD(ABS(key),LEN) + 1
ka = ir
!-----------LOOK IN THE TABLE.
20 kt = keytab(ka)
!----------CHECK FOR AN EMPTY SPACE, A DELETE ENTRY, OR A MATCH.
IF(kt == 0) GO TO 30
IF(kt == -1) GO TO 40
IF(kt == key) GO TO 60
ic = ic + 1
!---------COMPUTE ADDRESS OF NEXT PROBE.
ka = ka + iq
IF(ka > LEN) ka = ka - LEN
!----------SEE IF WHOLE TABLE HAS BEEN SEARCHED.
IF(ka /= ir) GO TO 20
!---------THE KEY IS NOT IN THE TABLE.
30 found = .false.
!---------RETURN WITH KA = 0 UNLES AN ENTRY HAS TO BE MADE.
IF(mode == 2 .AND. ic <= len2 .AND. key /= 0 .AND. key /= -1) GO TO 70
ka = 0
RETURN
!-----------A DELETED ENTRY HAS BEEN FOUND.
40 ia = ka
!----------COMPUTE ADDRESS OF NEXT PROBE.
50 ia = ia + iq
IF(ia > LEN) ia = ia - LEN
is = keytab(ia)
!----------CHECK FOR AN EMPTY SPACE OR A COMPLETE SCAN OF THE TABLE.
IF(is == 0 .OR. ia == ir) GO TO 30
!------------CHECK FOR A MISMATCH OR A DELETED ENTRY.
IF(is /= key .OR. is == -1) GO TO 50
!----------KEY FOUND.  MOVE IT AND THE ASSOCIAATED VALUE TO SAVE PROBES
!          ON THE NEXT SEARCH FOR THE SAME KEY.
kvtab(ka) = kvtab(ia)
keytab(ka) = is
keytab(ia) = -1
!----------THE KEY IS IN THE TABLE.
60 found = .true.
!----------DELETE IT IF MODE = 3
IF(mode == 3) keytab(ka) = -1
RETURN
!-----------LOOK FOR THE BEST WAY TO MAKE AN ENTRY.
70 IF(ic <= 0) GO TO 120
!----------SET DEL IF A DELETED ENTRY HAS BEEN FOUND.
del = kt /= 0
ia = ka
is = 0
!-----------COMPUTE THE MAXIMUM LENGTH TO SEARCH ALONG CURRENT CHAIN.
80 ix = ic - is
!-----------COMPUTE INCREMENT JQ FOR CURRENT CHAIN.
jq = MOD(ABS(keytab(ir)),len2) + 1
jr = ir
!----------LOOK ALONG THE CHAIN.
90 jr = jr + jq
IF(jr > LEN) jr = jr - LEN
kt = keytab(jr)
!----------CHECK FOR A HOLE (AN EMPTY SPACE OR DELETED ENTRY).
IF(kt == 0 .OR. kt == -1) GO TO 100
ix = ix - 1
IF(ix > 0) GO TO 90
GO TO 110
!----------SKIP IF THIS IS AN EMPTY SPACE AND A DELETED ENTRY HAS
!          ALREADY BEEN FOUND.
100 IF(del .AND. kt == 0) GO TO 110
!------------CHECK FOR A DELETED ENTRY.
IF(kt /= 0) del = .true.
!-----------SAVE LOCATION OF HOLE.
ia = jr
ka = ir
ic = ic - ix
!------------MOVE DOWN TO THE NEXT CHAIN.
110 is = is + 1
ir = ir + iq
IF(ir > LEN) ir = ir - LEN
!---------GO BACK IF A BETTER HOLE MIGHT STILL BE FOUND.
IF(ic > is) GO TO 80
!---------SKIP IF THERE IS NOTHING TO MOVE.
IF(ia == ka) GO TO 120
!---------MOVE AN OLD ENTRY AND ITS ASSOCIATED VALUE TO MAKE ROOM FOR
!         THE NEW ENTRY.
kvtab(ia) = kvtab(ka)
keytab(ia) = keytab(ka)
!---------ENTER THE NEW KEY, BUT NOT ITS ASSICIATED VALUE.
120 keytab(ka) = key

RETURN
END SUBROUTINE hash
