      SUBROUTINE iau_D2DTF ( SCALE, NDP, D1, D2, IY, IM, ID, IHMSF, J )
*+
*  - - - - - - - - - -
*   i a u _ D 2 D T F
*  - - - - - - - - - -
*
*  Format for output a 2-part Julian Date (or in the case of UTC a
*  quasi-JD form that includes special provision for leap seconds).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     SCALE      c*(*)  time scale ID (Note 1)
*     NDP        i      resolution (Note 2)
*     D1,D2      d      time as a 2-part Julian Date (Notes 3,4)
*
*  Returned:
*     IY,IM,ID   i      year, month, day in Gregorian calendar (Note 5)
*     IHMSF      i(4)   hours, minutes, seconds, fraction (Note 1)
*     J          i      status: +1 = dubious year (Note 5)
*                                0 = OK
*                               -1 = unacceptable date (Note 6)
*
*  Notes:
*
*  1) SCALE identifies the time scale.  Only the value 'UTC' (in upper
*     case) is significant, and enables handling of leap seconds (see
*     Note 4).
*
*  2) NDP is the number of decimal places in the seconds field, and can
*     have negative as well as positive values, such as:
*
*     NDP         resolution
*     -4            1 00 00
*     -3            0 10 00
*     -2            0 01 00
*     -1            0 00 10
*      0            0 00 01
*      1            0 00 00.1
*      2            0 00 00.01
*      3            0 00 00.001
*
*     The limits are platform dependent, but a safe range is -5 to +9.
*
*  3) D1+D2 is Julian Date, apportioned in any convenient way between
*     the two arguments, for example where D1 is the Julian Day Number
*     and D2 is the fraction of a day.  In the case of UTC, where the
*     use of JD is problematical, special conventions apply:  see the
*     next note.
*
*  4) JD cannot unambiguously represent UTC during a leap second unless
*     special measures are taken.  The SOFA internal convention is that
*     the quasi-JD day represents UTC days whether the length is 86399,
*     86400 or 86401 SI seconds.  In the 1960-1972 era there were
*     smaller jumps (in either direction) each time the linear UTC(TAI)
*     expression was changed, and these "mini-leaps" are also included
*     in the SOFA convention.
*
*  5) The warning status "dubious year" flags UTCs that predate the
*     introduction of the time scale or that are too far in the future
*     to be trusted.  See iau_DAT for further details.
*
*  6) For calendar conventions and limitations, see iau_CAL2JD.
*
*  Called:
*     iau_JD2CAL   JD to Gregorian calendar
*     iau_D2TF     decompose days to hms
*     iau_DAT      delta(AT) = TAI-UTC
*
*  This revision:  2014 February 15
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) SCALE
      INTEGER NDP
      DOUBLE PRECISION D1, D2
      INTEGER IY, IM, ID, IHMSF(4), J

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

      LOGICAL LEAP
      CHARACTER S
      INTEGER IY1, IM1, ID1, JS, IY2, IM2, ID2, IHMSF1(4), I
      DOUBLE PRECISION A1, B1, FD, DAT0, DAT12, W, DAT24, DLEAP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  The two-part JD.
      A1 = D1
      B1 = D2

*  Provisional calendar date.
      CALL iau_JD2CAL ( A1, B1, IY1, IM1, ID1, FD, JS )
      IF ( JS.NE.0 ) GO TO 9

*  Is this a leap second day?
      LEAP = .FALSE.
      IF ( SCALE.EQ.'UTC' ) THEN

*     TAI-UTC at 0h today.
         CALL iau_DAT ( IY1, IM1, ID1, 0D0, DAT0, JS )
         IF ( JS.LT.0 ) GO TO 9

*     TAI-UTC at 12h today (to detect drift).
         CALL iau_DAT ( IY1, IM1, ID1, 0.5D0, DAT12, JS )
         IF ( JS.LT.0 ) GO TO 9

*     TAI-UTC at 0h tomorrow (to detect jumps).
         CALL iau_JD2CAL ( A1+1.5D0, B1-FD, IY2, IM2, ID2, W, JS )
         IF ( JS.NE.0 ) GO TO 9
         CALL iau_DAT ( IY2, IM2, ID2, 0D0, DAT24, JS )
         IF ( JS.LT.0 ) GO TO 9

*     Any sudden change in TAI-UTC (seconds).
         DLEAP = DAT24 - ( 2D0 * DAT12 - DAT0 )

*     If leap second day, scale the fraction of a day into SI.
         LEAP = DLEAP.NE.0D0
         IF ( LEAP ) FD = FD + FD*DLEAP/D2S

      END IF

*  Provisional time of day.
      CALL iau_D2TF ( NDP, FD, S, IHMSF1 )

*  Has the (rounded) time gone past 24h?
      IF ( IHMSF1(1).GT.23 ) THEN

*     Yes.  We probably need tomorrow's calendar date.
         CALL iau_JD2CAL ( A1+1.5D0, B1-FD, IY2, IM2, ID2, W, JS )
         IF ( JS.LT.0 ) GO TO 9

*     Is today a leap second day?
         IF ( .NOT. LEAP ) THEN

*        No.  Use 0h tomorrow.
            IY1 = IY2
            IM1 = IM2
            ID1 = ID2
            IHMSF1(1) = 0
            IHMSF1(2) = 0
            IHMSF1(3) = 0

         ELSE

*        Yes.  Are we past the leap second itself?
            IF ( IHMSF1(3).GT.0 ) THEN

*           Yes.  Use tomorrow but allow for the leap second.
               IY1 = IY2
               IM1 = IM2
               ID1 = ID2
               IHMSF1(1) = 0
               IHMSF1(2) = 0
               IHMSF1(3) = 0

            ELSE

*           No.  Use 23 59 60... today.
               IHMSF1(1) = 23
               IHMSF1(2) = 59
               IHMSF1(3) = 60
            END IF

*        If rounding to 10s or coarser always go up to new day.
            IF ( NDP.LT.0 .AND. IHMSF1(3).EQ.60 ) THEN
               IY1 = IY2
               IM1 = IM2
               ID1 = ID2
               IHMSF1(1) = 0
               IHMSF1(2) = 0
               IHMSF1(3) = 0
            END IF
         END IF
      END IF

*  Results.
      IY = IY1
      IM = IM1
      ID = ID1
      DO 2 I=1,4
         IHMSF(I) = IHMSF1(I)
 2    CONTINUE

*  Status.
 9    CONTINUE
      J = JS

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
