C$Procedure      SPKW13 ( Write SPK segment, type 13 )
 
      SUBROUTINE SPKW13 (  HANDLE,  BODY,    CENTER,  FRAME,
     .                     FIRST,   LAST,    SEGID,   DEGREE,
     .                     N,       STATES,  EPOCHS           )
 
C$ Abstract
C
C     Write a type 13 segment to an SPK file.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     NAIF_IDS
C     SPC
C     SPK
C     TIME
C
C$ Keywords
C
C     EPHEMERIS
C     FILES
C
C$ Declarations
 
      IMPLICIT NONE 
      
      INTEGER               HANDLE
      INTEGER               BODY
      INTEGER               CENTER
      CHARACTER*(*)         FRAME
      DOUBLE PRECISION      FIRST
      DOUBLE PRECISION      LAST
      CHARACTER*(*)         SEGID
      INTEGER               DEGREE
      INTEGER               N
      DOUBLE PRECISION      STATES ( 6, * )
      DOUBLE PRECISION      EPOCHS (    * )
 
      INTEGER               MAXDEG
      PARAMETER           ( MAXDEG = 27 )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of an SPK file open for writing.
C     BODY       I   NAIF code for an ephemeris object.
C     CENTER     I   NAIF code for center of motion of BODY.
C     FRAME      I   Reference frame name.
C     FIRST      I   Start time of interval covered by segment.
C     LAST       I   End time of interval covered by segment.
C     SEGID      I   Segment identifier.
C     DEGREE     I   Degree of interpolating polynomials.
C     N          I   Number of states.
C     STATES     I   Array of states.
C     EPOCHS     I   Array of epochs corresponding to states.
C     MAXDEG     P   Maximum allowed degree of interpolating polynomial.
C
C$ Detailed_Input
C
C     HANDLE         is the file handle of an SPK file that has been
C                    opened for writing.
C
C     BODY           is the NAIF integer code for an ephemeris object
C                    whose state relative to another body is described
C                    by the segment to be created.
C
C     CENTER         is the NAIF integer code for the center of motion
C                    of the object identified by BODY.
C
C     FRAME          is the NAIF name for a reference frame
C                    relative to which the state information for BODY
C                    is specified.
C
C     FIRST,
C     LAST           are, respectively, the start and stop times of
C                    the time interval over which the segment defines
C                    the state of BODY.
C
C     SEGID          is the segment identifier.  An SPK segment
C                    identifier may contain up to 40 characters.
C
C     DEGREE         is the degree of the Hermite polynomials used to
C                    interpolate the states.  All components of the
C                    state vectors are interpolated by polynomials of
C                    fixed degree.
C
C     N              is the number of states in the input state vector
C                    array.
C
C     STATES         contains a time-ordered array of geometric states
C                    ( x, y, z, dx/dt, dy/dt, dz/dt, in kilometers and
C                    kilometers per second ) of BODY relative to CENTER,
C                    specified relative to FRAME.
C
C     EPOCHS         is an array of epochs corresponding to the members
C                    of the state array.  The epochs are specified as
C                    seconds past J2000, TDB.
C
C$ Detailed_Output
C
C     None.  See $Particulars for a description of the effect of this
C     routine.
C
C$ Parameters
C
C     MAXDEG         is the maximum allowed degree of the interpolating
C                    polynomial.  If the value of MAXDEG is increased,
C                    the SPICELIB routine SPKPV must be changed
C                    accordingly.  In particular, the size of the
C                    record passed to SPKRnn and SPKEnn must be
C                    increased, and comments describing the record size
C                    must be changed.
C
C$ Exceptions
C
C     If any of the following exceptions occur, this routine will return
C     without creating a new segment.
C
C     1) If FRAME is not a recognized name, the error
C        SPICE(INVALIDREFFRAME) is signaled.
C
C     2) If the last non-blank character of SEGID occurs past index 40,
C        the error SPICE(SEGIDTOOLONG) is signaled.
C
C     3) If SEGID contains any nonprintable characters, the error
C        SPICE(NONPRINTABLECHARS) is signaled.
C
C     4) If DEGREE is not at least 1 or is greater than MAXDEG, the
C        error SPICE(INVALIDDEGREE) is signaled.
C
C     5) If DEGREE is not odd, the error SPICE(INVALIDDEGREE) is 
C        signaled.
C
C     6) If the number of states N is not at least (DEGREE+1)/2, 
C        the error SPICE(TOOFEWSTATES) will be signaled.
C
C     7) If FIRST is greater than or equal to LAST then the error
C        SPICE(BADDESCRTIMES) will be signaled.
C
C     8) If the elements of the array EPOCHS are not in strictly
C        increasing order, the error SPICE(TIMESOUTOFORDER) will be
C        signaled.
C
C     9) If the first epoch EPOCHS(1) is greater than FIRST, the error
C        SPICE(BADDESCRTIMES) will be signaled.
C
C    10) If the last epoch EPOCHS(N) is less than LAST, the error
C        SPICE(BADDESCRTIMES) will be signaled.
C
C
C$ Files
C
C     A new type 13 SPK segment is written to the SPK file attached
C     to HANDLE.
C
C$ Particulars
C
C     This routine writes an SPK type 13 data segment to the open SPK
C     file according to the format described in the type 13 section of
C     the SPK Required Reading. The SPK file must have been opened with
C     write access.
C
C$ Examples
C
C     Suppose that you have states and are prepared to produce
C     a segment of type 13 in an SPK file.
C
C     The following code fragment could be used to add the new segment
C     to a previously opened SPK file attached to HANDLE. The file must
C     have been opened with write access.
C
C        C
C        C     Create a segment identifier.
C        C
C                  SEGID = 'MY_SAMPLE_SPK_TYPE_13_SEGMENT'
C
C        C
C        C     Write the segment.
C        C
C              CALL SPKW13 (  HANDLE,  BODY,    CENTER,  FRAME,
C             .               FIRST,   LAST,    SEGID,   DEGREE,
C             .               N,       STATES,  EPOCHS          )
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 2.0.0, 23-DEC-2013 (NJB)
C
C        Increased MAXDEG to 27 for compatibility
C        with SPK type 21. Deleted declaration of 
C        unused parameter TYPIDX.
C        
C
C-    SPICELIB Version 1.0.0, 20-MAR-2000 (NJB)
C
C-&
 
C$ Index_Entries
C
C     write spk type_13 ephemeris data segment
C
C-&
 
  
C
C     SPICELIB functions
C
      INTEGER               LASTNB
      
      LOGICAL               EVEN
      LOGICAL               FAILED
      LOGICAL               RETURN
 
C
C     Local parameters
C
      INTEGER               SIDLEN
      PARAMETER           ( SIDLEN  =  40 )
 
      INTEGER               FPRINT
      PARAMETER           ( FPRINT  =  32 )
 
      INTEGER               LPRINT
      PARAMETER           ( LPRINT  = 126 )
 
C
C     Local variables
C
      INTEGER               DSCSIZ
      PARAMETER           ( DSCSIZ  =   5 )
           
      INTEGER               DTYPE
      PARAMETER           ( DTYPE   =  13 )

      INTEGER               DIRSIZ
      PARAMETER           ( DIRSIZ  = 100 )
 
      INTEGER               STATSZ
      PARAMETER           ( STATSZ  =   6 )
       
C
C     Local variables
C
      DOUBLE PRECISION      DESCR  ( DSCSIZ )
      DOUBLE PRECISION      MAXTIM
 
      INTEGER               CHRCOD
      INTEGER               I
      INTEGER               REFCOD
      INTEGER               WINSIZ
      
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'SPKW13' )
      END IF
      
C
C     Set the window size corresponding to the input degree.  This
C     size will be used in various places below.
C
      WINSIZ  =  ( DEGREE + 1 ) / 2
      
C
C     Get the NAIF integer code for the reference frame.
C
      CALL NAMFRM ( FRAME, REFCOD )
 
      IF ( REFCOD .EQ. 0 ) THEN
 
         CALL SETMSG ( 'The reference frame # is not supported.'   )
         CALL ERRCH  ( '#', FRAME                                  )
         CALL SIGERR ( 'SPICE(INVALIDREFFRAME)'                    )
         CALL CHKOUT ( 'SPKW13'                                    )
         RETURN
 
      END IF
 
C
C     Check to see if the segment identifier is too long.
C
      IF ( LASTNB(SEGID) .GT. SIDLEN ) THEN
 
         CALL SETMSG ( 'Segment identifier contains more than ' //
     .                 '40 characters.'                          )
         CALL SIGERR ( 'SPICE(SEGIDTOOLONG)'                     )
         CALL CHKOUT ( 'SPKW13'                                  )
         RETURN
 
      END IF
 
C
C     Now check that all the characters in the segment identifier
C     can be printed.
C
      DO I = 1, LASTNB(SEGID)
 
         CHRCOD = ICHAR( SEGID(I:I) )
 
         IF ( ( CHRCOD .LT. FPRINT ) .OR. ( CHRCOD .GT. LPRINT ) ) THEN
 
            CALL SETMSG ( 'The segment identifier contains '  //
     .                    'nonprintable characters'               )
            CALL SIGERR ( 'SPICE(NONPRINTABLECHARS)'              )
            CALL CHKOUT ( 'SPKW13'                                )
            RETURN
 
         END IF
 
      END DO
 
C
C     Make sure that the degree of the interpolating polynomials is
C     in range.
C
      IF (  ( DEGREE .LT. 1 ) .OR. ( DEGREE .GT. MAXDEG )  )  THEN
 
         CALL SETMSG ( 'The interpolating polynomials have degree #; '//
     .                 'the valid degree range is [1, #]'              )
         CALL ERRINT ( '#', DEGREE                                     )
         CALL ERRINT ( '#', MAXDEG                                     )
         CALL SIGERR ( 'SPICE(INVALIDDEGREE)'                          )
         CALL CHKOUT ( 'SPKW13'                                        )
         RETURN
 
      END IF
 
C
C     Make sure that the degree of the interpolating polynomials is odd.
C
      IF (  EVEN( DEGREE )  )  THEN
 
         CALL SETMSG ( 'The interpolating polynomials have degree #; '//
     .                 'for SPK type 13, the degree must be odd.'      )
         CALL ERRINT ( '#', DEGREE                                     )
         CALL SIGERR ( 'SPICE(INVALIDDEGREE)'                          )
         CALL CHKOUT ( 'SPKW13'                                        )
         RETURN
 
      END IF
 
C
C     Make sure that the number of states is sufficient to define a
C     polynomial whose degree is DEGREE.
C
      IF ( N .LT. WINSIZ ) THEN
 
         CALL SETMSG ( 'At least # states are required to define a '  //
     .                 'Hermite polynomial of degree #.  Number of '  //
     .                 'states supplied:  #'                           )
         CALL ERRINT ( '#', WINSIZ                                     )
         CALL ERRINT ( '#', DEGREE                                     )
         CALL ERRINT ( '#', N                                          )
         CALL SIGERR ( 'SPICE(TOOFEWSTATES)'                           )
         CALL CHKOUT ( 'SPKW13'                                        )
         RETURN
 
      END IF
 
C
C     The segment stop time should be greater then the begin time.
C
      IF ( FIRST .GE. LAST ) THEN
 
         CALL SETMSG ( 'The segment start time: # is greater then ' //
     .                 'the segment end time: #'                   )
         CALL ERRDP  ( '#', FIRST                                  )
         CALL ERRDP  ( '#', LAST                                   )
         CALL SIGERR ( 'SPICE(BADDESCRTIMES)'                      )
         CALL CHKOUT ( 'SPKW13'                                    )
         RETURN
 
      END IF
 
C
C     Make sure the epochs form a strictly increasing sequence.
C
      MAXTIM = EPOCHS(1)
 
      DO I = 2, N
 
         IF ( EPOCHS(I) .LE. MAXTIM ) THEN
 
            CALL SETMSG ( 'EPOCH # having index # is not greater '    //
     .                    'than its predecessor #.'                   )
            CALL ERRDP  ( '#',  EPOCHS(I)                             )
            CALL ERRINT ( '#',  I                                     )
            CALL ERRDP  ( '#',  EPOCHS(I-1)                           )
            CALL SIGERR ( 'SPICE(TIMESOUTOFORDER)'                    )
            CALL CHKOUT ( 'SPKW13'                                    )
            RETURN
         ELSE
            MAXTIM = EPOCHS(I)
         END IF
 
      END DO
 
C
C     Make sure that the span of the input epochs includes the interval
C     defined by the segment descriptor.
C
      IF (  EPOCHS(1) .GT. FIRST ) THEN
 
         CALL SETMSG ( 'Segment start time # precedes first epoch #.' )
         CALL ERRDP  ( '#',  FIRST                                    )
         CALL ERRDP  ( '#',  EPOCHS(1)                                )
         CALL SIGERR ( 'SPICE(BADDESCRTIMES)'                         )
         CALL CHKOUT ( 'SPKW13'                                       )
         RETURN
 
      ELSE IF (  EPOCHS(N) .LT. LAST ) THEN
 
         CALL SETMSG ( 'Segment end time # follows last epoch #.'     )
         CALL ERRDP  ( '#',  LAST                                     )
         CALL ERRDP  ( '#',  EPOCHS(N)                                )
         CALL SIGERR ( 'SPICE(BADDESCRTIMES)'                         )
         CALL CHKOUT ( 'SPKW13'                                       )
         RETURN
 
      END IF
 


C
C     If we made it this far, we're ready to start writing the segment.
C
C
C     Create the segment descriptor.
C
      CALL SPKPDS ( BODY, CENTER, FRAME, DTYPE, FIRST, LAST, DESCR )
 
C
C     Begin a new segment.
C
      CALL DAFBNA ( HANDLE, DESCR, SEGID )
 
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'SPKW13' )
         RETURN
      END IF
      
      
C
C     The type 13 segment structure is eloquently described by this
C     diagram from the SPK Required Reading:
C
C        +-----------------------+
C        | State 1               |
C        +-----------------------+
C        | State 2               |
C        +-----------------------+
C                    .
C                    .
C                    .
C        +-----------------------+
C        | State N               |
C        +-----------------------+
C        | Epoch 1               |
C        +-----------------------+
C        | Epoch 2               |
C        +-----------------------+
C                    .
C                    .
C                    .
C        +-----------------------+
C        | Epoch N               |
C        +-----------------------+
C        | Epoch 100             | (First directory)
C        +-----------------------+
C                    .
C                    .
C                    .
C        +-----------------------+
C        | Epoch ((N-1)/100)*100 | (Last directory)
C        +-----------------------+
C        | Window size - 1       |
C        +-----------------------+
C        | Number of states      |
C        +-----------------------+
C
C
 
      CALL DAFADA ( STATES,  N*STATSZ )
      CALL DAFADA ( EPOCHS,  N        )
 
      DO I = 1,   (N-1) / DIRSIZ
         CALL DAFADA (  EPOCHS( DIRSIZ*I ),  1  )
      END DO
 
      CALL DAFADA ( DBLE( WINSIZ - 1),  1 )
      CALL DAFADA ( DBLE( N         ),  1 )
 
C
C     As long as nothing went wrong, end the segment.
C
      IF ( .NOT. FAILED() ) THEN
         CALL DAFENA
      END IF
 
 
      CALL CHKOUT ( 'SPKW13' )
      RETURN
      END
