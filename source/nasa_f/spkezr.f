C$Procedure SPKEZR ( S/P Kernel, easier reader )

      SUBROUTINE SPKEZR ( TARG, ET, REF, ABCORR, OBS, STARG, LT )

C$ Abstract
C
C     Return the state (position and velocity) of a target body
C     relative to an observing body, optionally corrected for light
C     time (planetary aberration) and stellar aberration.
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
C     SPK
C     NAIF_IDS
C     FRAMES
C     TIME
C
C$ Keywords
C
C     EPHEMERIS
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE               'frmtyp.inc'
      INCLUDE               'zzctr.inc'

      CHARACTER*(*)         TARG
      DOUBLE PRECISION      ET
      CHARACTER*(*)         REF
      CHARACTER*(*)         ABCORR
      CHARACTER*(*)         OBS
      DOUBLE PRECISION      STARG    ( 6 )
      DOUBLE PRECISION      LT

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     TARG       I   Target body name.
C     ET         I   Observer epoch.
C     REF        I   Reference frame of output state vector.
C     ABCORR     I   Aberration correction flag.
C     OBS        I   Observing body name.
C     STARG      O   State of target.
C     LT         O   One way light time between observer and target.
C
C$ Detailed_Input
C
C     TARG        is the name of a target body. Optionally, you may
C                 supply the integer ID code for the object as
C                 an integer string. For example both 'MOON' and
C                 '301' are legitimate strings that indicate the
C                 moon is the target body.
C
C                 The target and observer define a state vector whose
C                 position component points from the observer to the
C                 target.
C
C     ET          is the ephemeris time, expressed as seconds past J2000
C                 TDB, at which the state of the target body relative to
C                 the observer is to be computed. ET refers to time at
C                 the observer's location.
C
C     REF         is the name of the reference frame relative to which
C                 the output state vector should be expressed. This may
C                 be any frame supported by the SPICE system, including
C                 built-in frames (documented in the Frames Required
C                 Reading) and frames defined by a loaded frame kernel
C                 (FK).
C
C                 When REF designates a non-inertial frame, the
C                 orientation of the frame is evaluated at an epoch
C                 dependent on the selected aberration correction.
C                 See the description of the output state vector STARG
C                 for details.
C
C     ABCORR      indicates the aberration corrections to be applied
C                 to the state of the target body to account for one-way
C                 light time and stellar aberration. See the discussion
C                 in the Particulars section for recommendations on
C                 how to choose aberration corrections.
C
C                 ABCORR may be any of the following:
C
C                    'NONE'     Apply no correction. Return the
C                               geometric state of the target body
C                               relative to the observer.
C
C                 The following values of ABCORR apply to the
C                 "reception" case in which photons depart from the
C                 target's location at the light-time corrected epoch
C                 ET-LT and *arrive* at the observer's location at ET:
C
C                    'LT'       Correct for one-way light time (also
C                               called "planetary aberration") using a
C                               Newtonian formulation. This correction
C                               yields the state of the target at the
C                               moment it emitted photons arriving at
C                               the observer at ET.
C
C                               The light time correction uses an
C                               iterative solution of the light time
C                               equation (see Particulars for details).
C                               The solution invoked by the 'LT' option
C                               uses one iteration.
C
C                    'LT+S'     Correct for one-way light time and
C                               stellar aberration using a Newtonian
C                               formulation. This option modifies the
C                               state obtained with the 'LT' option to
C                               account for the observer's velocity
C                               relative to the solar system
C                               barycenter. The result is the apparent
C                               state of the target---the position and
C                               velocity of the target as seen by the
C                               observer.
C
C                    'CN'       Converged Newtonian light time
C                               correction. In solving the light time
C                               equation, the 'CN' correction iterates
C                               until the solution converges (three
C                               iterations on all supported platforms).
C                               Whether the 'CN+S' solution is
C                               substantially more accurate than the
C                               'LT' solution depends on the geometry
C                               of the participating objects and on the
C                               accuracy of the input data. In all
C                               cases this routine will execute more
C                               slowly when a converged solution is
C                               computed. See the Particulars section
C                               below for a discussion of precision of
C                               light time corrections.
C
C                    'CN+S'     Converged Newtonian light time
C                               correction and stellar aberration
C                               correction.
C
C
C                 The following values of ABCORR apply to the
C                 "transmission" case in which photons *depart* from
C                 the observer's location at ET and arrive at the
C                 target's location at the light-time corrected epoch
C                 ET+LT:
C
C                    'XLT'      "Transmission" case:  correct for
C                               one-way light time using a Newtonian
C                               formulation. This correction yields the
C                               state of the target at the moment it
C                               receives photons emitted from the
C                               observer's location at ET.
C
C                    'XLT+S'    "Transmission" case:  correct for
C                               one-way light time and stellar
C                               aberration using a Newtonian
C                               formulation  This option modifies the
C                               state obtained with the 'XLT' option to
C                               account for the observer's velocity
C                               relative to the solar system
C                               barycenter. The position component of
C                               the computed target state indicates the
C                               direction that photons emitted from the
C                               observer's location must be "aimed" to
C                               hit the target.
C
C                    'XCN'      "Transmission" case:  converged
C                               Newtonian light time correction.
C
C                    'XCN+S'    "Transmission" case:  converged
C                               Newtonian light time correction and
C                               stellar aberration correction.
C
C
C                 Neither special nor general relativistic effects are
C                 accounted for in the aberration corrections applied
C                 by this routine.
C
C                 Case and blanks are not significant in the string
C                 ABCORR.
C
C     OBS         is the name of an observing body. Optionally, you
C                 may supply the ID code of the object as an integer
C                 string. For example, both 'EARTH' and '399' are
C                 legitimate strings to supply to indicate the
C                 observer is Earth.
C
C$ Detailed_Output
C
C     STARG       is a Cartesian state vector representing the position
C                 and velocity of the target body relative to the
C                 specified observer. STARG is corrected for the
C                 specified aberrations, and is expressed with respect
C                 to the reference frame specified by REF. The first
C                 three components of STARG represent the x-, y- and
C                 z-components of the target's position; the last three
C                 components form the corresponding velocity vector.
C
C                 The position component of STARG points from the
C                 observer's location at ET to the aberration-corrected
C                 location of the target. Note that the sense of the
C                 position vector is independent of the direction of
C                 radiation travel implied by the aberration
C                 correction.
C
C                 The velocity component of STARG is the derivative
C                 with respect to time of the position component of
C                 STARG.
C
C                 Units are always km and km/sec.
C
C                 Non-inertial frames are treated as follows: letting
C                 LTCENT be the one-way light time between the observer
C                 and the central body associated with the frame, the
C                 orientation of the frame is evaluated at ET-LTCENT,
C                 ET+LTCENT, or ET depending on whether the requested
C                 aberration correction is, respectively, for received
C                 radiation, transmitted radiation, or is omitted.
C                 LTCENT is computed using the method indicated by
C                 ABCORR.
C
C     LT          is the one-way light time between the observer and
C                 target in seconds. If the target state is corrected
C                 for aberrations, then LT is the one-way light time
C                 between the observer and the light time corrected
C                 target location.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If name of target or observer cannot be translated to its
C        NAIF ID code, the error SPICE(IDCODENOTFOUND) is signaled.
C
C     2) If the reference frame REF is not a recognized reference
C        frame the error 'SPICE(UNKNOWNFRAME)' is signaled.
C
C     3) If the loaded kernels provide insufficient data to
C        compute the requested state vector, the deficiency will
C        be diagnosed by a routine in the call tree of this routine.
C
C     4) If an error occurs while reading an SPK or other kernel file,
C        the error  will be diagnosed by a routine in the call tree
C        of this routine.
C
C$ Files
C
C     This routine computes states using SPK files that have been
C     loaded into the SPICE system, normally via the kernel loading
C     interface routine FURNSH. See the routine FURNSH and the SPK
C     and KERNEL Required Reading for further information on loading
C     (and unloading) kernels.
C
C     If the output state STARG is to be expressed relative to a
C     non-inertial frame, or if any of the ephemeris data used to
C     compute STARG are expressed relative to a non-inertial frame in
C     the SPK files providing those data, additional kernels may be
C     needed to enable the reference frame transformations required to
C     compute the state. Normally these additional kernels are PCK
C     files or frame kernels. Any such kernels must already be loaded
C     at the time this routine is called.
C
C$ Particulars
C
C     This routine is part of the user interface to the SPICE ephemeris
C     system. It allows you to retrieve state information for any
C     ephemeris object relative to any other in a reference frame that
C     is convenient for further computations.
C
C     This routine is identical in function to the routine SPKEZ except
C     that it allows you to refer to ephemeris objects by name (via a
C     character string).
C
C     Please refer to the Aberation Corrections Required Reading
C     (ABCORR.REQ) for detailed information describing the nature and
C     calculation of the applied corrections.
C
C$ Examples
C
C     1)  Load a planetary ephemeris SPK, then look up a
C         state of the MARS BARYCENTER relative to EARTH
C         in the J2000 frame with aberration correction LT+S.
C
C
C           PROGRAM SPKEZR_T
C           IMPLICIT NONE
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(32)        FRAME
C           CHARACTER*(32)        ABCORR
C           CHARACTER*(36)        OBS
C           CHARACTER*(36)        TARGET
C           CHARACTER*(36)        EPOCH
C
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      LT
C           DOUBLE PRECISION      STATE ( 6 )
C
C           INTEGER               I
C
C     C
C     C     Load a set of kernels: an SPK file, a PCK
C     C     file and a leapseconds file.
C     C
C           CALL FURNSH( 'naif0011.tls' )
C           CALL FURNSH( 'pck00010.tpc' )
C           CALL FURNSH( 'de430.bsp' )
C
C     C
C     C     Define parameters for a state lookup:
C     C
C     C     Return the state vector of Mars Barycenter (4) as seen
C     C     from Earth (399) in the J2000 frame  using aberration
C     C     correction LT+S (light time plus stellar aberration)
C     C     at the epoch JAN 1 2015 12:00:00.
C     C
C           TARGET   = 'MARS BARYCENTER'
C           EPOCH    = 'JAN 1 2015 12:00:00'
C           FRAME    = 'J2000'
C           ABCORR   = 'LT+S'
C           OBS      = 'EARTH'
C
C     C
C     C     Convert the epoch to ephemeris time.
C     C
C           CALL STR2ET( EPOCH, ET )
C
C     C
C     C     Look-up the state for the defined parameters.
C     C
C           CALL SPKEZR( TARGET, ET, FRAME, ABCORR, OBS, STATE, LT)
C
C     C
C     C     Output...
C     C
C           WRITE(*,*) 'The position of    : ', TARGET
C           WRITE(*,*) 'As observed from   : ', OBS
C           WRITE(*,*) 'In reference frame : ', FRAME
C           WRITE(*,*) 'At epoch           : ', EPOCH
C           WRITE(*,*) ' '
C
C     C
C     C     The first three entries of state contain the
C     C     X, Y, Z position components. The final three contain
C     C     the Vx, Vy, Vz velocity components.
C     C
C           WRITE(*,*) 'R (kilometers)     : '
C           WRITE(*,*) (STATE(I), I=1,3 )
C
C           WRITE(*,*) 'V (kilometers/sec) : '
C           WRITE(*,*) (STATE(I), I=4,6 )
C
C           WRITE(*,*) 'Light time (secs)  : ', LT
C
C           END
C
C   The program outputs:
C
C      The position of    : MARS BARYCENTER
C      As observed from   : EARTH
C      In reference frame : J2000
C      At epoch           : JAN 1 2015 12:00:00
C
C      R (kilometers)     :
C        229953013.74649832   -167125346.21158829   -78800343.963572651
C      V (kilometers/sec) :
C        35.380861440845095    28.653401530093195    12.861523564981141
C      Light time (secs)  :    983.97882466162321
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     SPK Required Reading.
C
C$ Author_and_Institution
C
C     C.H. Acton      (JPL)
C     B.V. Semenov    (JPL)
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 4.1.1, 19-JAN-2016 (EDW)
C
C       Example code replaced with a complete program and
C       the corresponding output.
C
C       Particulars updated to refer to Aberration Corrections
C       Required Reading document.
C
C-    SPICELIB Version 4.1.0, 03-JUL-2014 (NJB) (BVS)
C
C        Discussion of light time corrections was updated. Assertions
C        that converged light time corrections are unlikely to be
C        useful were removed.
C
C     Last update was 19-SEP-2013 (BVS)
C
C        Updated to save the input body names and ZZBODTRN state
C        counters and to do name-ID conversions only if the counters
C        have changed.
C
C-    SPICELIB Version 4.0.0, 27-DEC-2007 (NJB)
C
C        This routine was upgraded to more accurately compute
C        aberration-corrected velocity, and in particular, make it
C        more consistent with observer-target positions.
C
C        When light time corrections are used, the derivative of light
C        time with respect to time is now accounted for in the
C        computation of observer-target velocities. When the reference
C        frame associated with the output state is time-dependent, the
C        derivative of light time with respect to time is now accounted
C        for in the computation of the rate of change of orientation of
C        the reference frame.
C
C        When stellar aberration corrections are used, velocities
C        now reflect the rate of range of the stellar aberration
C        correction.
C
C-    SPICELIB Version 3.0.2, 20-OCT-2003 (EDW)
C
C        Added mention that LT returns in seconds.
C
C-    SPICELIB Version 3.0.1, 29-JUL-2003 (NJB) (CHA)
C
C        Various minor header changes were made to improve clarity.
C
C-    SPICELIB Version 3.0.0, 31-DEC-2001 (NJB)
C
C        Updated to handle aberration corrections for transmission
C        of radiation. Formerly, only the reception case was
C        supported. The header was revised and expanded to explain
C        the functionality of this routine in more detail.
C
C-    Spicelib Version 2.0.0, 21-FEB-1997 (WLT)
C
C        Extended the functionality of the routine. Users may
C        now entered the id code of an object as an ascii string
C        and the string will be converted to the corresponding
C        integer representation.
C
C-    Spicelib Version 1.1.0, 09-JUL-1996 (WLT)
C
C        Corrected the description of LT in the Detailed Output
C        section of the header.
C
C-    SPICELIB Version 1.0.0, 25-SEP-1995 (BVS)
C
C-&

C$ Index_Entries
C
C     using body names get target state relative to an observer
C     get state relative to observer corrected for aberrations
C     read ephemeris data
C     read trajectory data
C
C-&


C$ Revisions
C
C     None.
C
C-&

C
C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Saved body name length.
C
      INTEGER               MAXL
      PARAMETER           ( MAXL  = 36 )

C
C     Local variables
C
      INTEGER               TARGID
      INTEGER               OBSID

      LOGICAL               FOUND

C
C     Saved name/ID item declarations.
C
      INTEGER               SVCTR1 ( CTRSIZ )
      CHARACTER*(MAXL)      SVTARG
      INTEGER               SVTGID
      LOGICAL               SVFND1

      INTEGER               SVCTR2 ( CTRSIZ )
      CHARACTER*(MAXL)      SVOBSN
      INTEGER               SVOBSI
      LOGICAL               SVFND2

      LOGICAL               FIRST

C
C     Saved name/ID items.
C
      SAVE                  SVCTR1
      SAVE                  SVTARG
      SAVE                  SVTGID
      SAVE                  SVFND1

      SAVE                  SVCTR2
      SAVE                  SVOBSN
      SAVE                  SVOBSI
      SAVE                  SVFND2

      SAVE                  FIRST

C
C     Initial values.
C
      DATA                  FIRST   / .TRUE. /


C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'SPKEZR' )
      END IF

C
C     Initialization.
C
      IF ( FIRST ) THEN

C
C        Initialize counters.
C
         CALL ZZCTRUIN( SVCTR1 )
         CALL ZZCTRUIN( SVCTR2 )

         FIRST = .FALSE.

      END IF

C
C     Starting from translation of target name to its code
C
      CALL ZZBODS2C ( SVCTR1, SVTARG, SVTGID, SVFND1,
     .                TARG, TARGID, FOUND    )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'The target, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit. '
     .   //            'Alternatively you may call SPKEZ '
     .   //            'directly if you know the SPICE ID codes '
     .   //            'for both ''#'' and ''#'' '                )
         CALL ERRCH  ( '#', TARG                                  )
         CALL ERRCH  ( '#', TARG                                  )
         CALL ERRCH  ( '#', OBS                                   )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( 'SPKEZR'                                   )
         RETURN

      END IF

C
C     Now do the same for observer
C
      CALL ZZBODS2C ( SVCTR2, SVOBSN, SVOBSI, SVFND2,
     .                OBS, OBSID, FOUND    )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'The observer, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE toolkit. '
     .   //            'Alternatively you may call SPKEZ '
     .   //            'directly if you know the SPICE ID codes '
     .   //            'for both ''#'' and ''#'' '                )
         CALL ERRCH  ( '#', OBS                                   )
         CALL ERRCH  ( '#', TARG                                  )
         CALL ERRCH  ( '#', OBS                                   )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( 'SPKEZR'                                   )
         RETURN

      END IF

C
C     After all translations are done we can call SPKEZ.
C
      CALL SPKEZ   ( TARGID, ET, REF, ABCORR, OBSID, STARG, LT )


      CALL CHKOUT ( 'SPKEZR' )
      RETURN
      END
