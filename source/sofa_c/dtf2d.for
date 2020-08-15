      SUBROUTINE iau_DTF2D ( SCALE, IY, IM, ID, IHR, IMN, SEC,
     :                       D1, D2, J )
*+
*  - - - - - - - - - -
*   i a u _ D T F 2 D
*  - - - - - - - - - -
*
*  Encode date and time fields into 2-part Julian Date (or in the case
*  of UTC a quasi-JD form that includes special provision for leap
*  seconds).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     SCALE      c*(*)  time scale ID (Note 1)
*     IY,IM,ID   i      year, month, day in Gregorian calendar (Note 2)
*     IHR,IMN    i      hour, minute
*     SEC        d      seconds
*
*  Returned:
*     D1,D2      d      2-part Julian Date (Notes 3,4)
*     J          i      status: +3 = both of next two
*                               +2 = time is after end of day (Note 5)
*                               +1 = dubious year (Note 6)
*                                0 = OK
*                               -1 = bad year
*                               -2 = bad month
*                               -3 = bad day
*                               -4 = bad hour
*                               -5 = bad minute
*                               -6 = bad second (<0)
*
*  Notes:
*
*  1) SCALE identifies the time scale.  Only the value 'UTC' (in upper
*     case) is significant, and enables handling of leap seconds (see
*     Note 4).
*
*  2) For calendar conventions and limitations, see iau_CAL2JD.
*
*  3) The sum of the results, D1+D2, is Julian Date, where normally D1
*     is the Julian Day Number and D2 is the fraction of a day.  In the
*     case of UTC, where the use of JD is problematical, special
*     conventions apply:  see the next note.
*
*  4) JD cannot unambiguously represent UTC during a leap second unless
*     special measures are taken.  The SOFA internal convention is that
*     the quasi-JD day represents UTC days whether the length is 86399,
*     86400 or 86401 SI seconds.  In the 1960-1972 era there were
*     smaller jumps (in either direction) each time the linear UTC(TAI)
*     expression was changed, and these "mini-leaps" are also included
*     in the SOFA convention.
*
*  5) The warning status "time is after end of day" usually means that
*     the SEC argument is greater than 60D0.  However, in a day ending
*     in a leap second the limit changes to 61D0 (or 59D0 in the case of
*     a negative leap second).
*
*  6) The warning status "dubious year" flags UTCs that predate the
*     introduction of the time scale or that are too far in the future
*     to be trusted.  See iau_DAT for further details.
*
*  7) Only in the case of continuous and regular time scales (TAI, TT,
*     TCG, TCB and TDB) is the result D1+D2 a Julian Date, strictly
*     speaking.  In the other cases (UT1 and UTC) the result must be
*     used with circumspection;  in particular the difference between
*     two such results cannot be interpreted as a precise time
*     interval.
*
*  Called:
*     iau_CAL2JD   Gregorian calendar to JD
*     iau_DAT      delta(AT) = TAI-UTC
*     iau_JD2CAL   JD to Gregorian calendar
*
*  This revision:  2013 July 26
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) SCALE
      INTEGER IY, IM, ID, IHR, IMN
      DOUBLE PRECISION SEC, D1, D2
      INTEGER J

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

      INTEGER JS, IY2, IM2, ID2
      DOUBLE PRECISION DJ, W, DAY, SECLIM, DAT0, DAT12, DAT24,
     :                 DLEAP, TIME

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Today's Julian Day Number.
      CALL iau_CAL2JD ( IY, IM, ID, DJ, W, JS )
      IF ( JS.NE.0 ) GO TO 9
      DJ = DJ + W

*  Day length and final minute length in seconds (provisional).
      DAY = D2S
      SECLIM = 60D0

*  Deal with the UTC leap second case.
      IF ( SCALE.EQ.'UTC' ) THEN

*     TAI-UTC at 0h today.
         CALL iau_DAT ( IY, IM, ID, 0D0, DAT0, JS )
         IF ( JS.LT.0 ) GO TO 9

*     TAI-UTC at 12h today (to detect drift).
         CALL iau_DAT ( IY, IM, ID, 0.5D0, DAT12, JS )
         IF ( JS.LT.0 ) GO TO 9

*     TAI-UTC at 0h tomorrow (to detect jumps).
         CALL iau_JD2CAL ( DJ, 1.5D0, IY2, IM2, ID2, W, JS )
         IF ( JS.NE.0 ) GO TO 9
         CALL iau_DAT ( IY2, IM2, ID2, 0D0, DAT24, JS )
         IF ( JS.LT.0 ) GO TO 9

*     Any sudden change in TAI-UTC between today and tomorrow.
         DLEAP = DAT24 - ( 2D0 * DAT12 - DAT0 )

*     If leap second day, correct the day and final minute lengths.
         DAY = DAY + DLEAP
         IF ( IHR.EQ.23 .AND. IMN.EQ.59 ) SECLIM = SECLIM + DLEAP

*     End of UTC-specific actions.
      END IF

*  Validate the time.
      IF ( IHR.GE.0 .AND. IHR.LE.23 ) THEN
         IF ( IMN.GE.0 .AND. IMN.LE.59 ) THEN
            IF ( SEC.GE.0D0 ) THEN
               IF ( SEC.GE.SECLIM ) THEN
                  JS = JS + 2
               END IF
            ELSE
               JS = -6
            END IF
         ELSE
            JS = -5
         END IF
      ELSE
         JS = -4
      END IF
      IF ( JS.LT.0 ) GO TO 9

*  The time in days.
      TIME = (60D0*DBLE(60*IHR+IMN)+SEC) / DAY

*  Return the date and time.
      D1 = DJ
      D2 = TIME

*  Return the status.
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
