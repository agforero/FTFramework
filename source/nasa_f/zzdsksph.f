C$Procedure ZZDSKSPH ( DSK, bounding spheres for target body )
 
      SUBROUTINE ZZDSKSPH ( BODYID, NSURF, SRFLST, MINRAD, MAXRAD )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Return radii of outer and inner bounding spheres for a given body
C     and surface list. The shape of the body is represented by DSK
C     data. The outputs of this routine are time-independent.
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
C     DSK
C
C$ Keywords
C     
C     DSK
C     GEOMETRY
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'errhnd.inc'
      INCLUDE 'dla.inc'
      INCLUDE 'dsk.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'zzctr.inc'

      INTEGER               BODYID
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      DOUBLE PRECISION      MINRAD
      DOUBLE PRECISION      MAXRAD

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body.
C     NSURF      I   Number of IDs in surface list.
C     SRFLST     I   List of surface IDs.
C     MINRAD     O   Radius of inner bounding sphere for body.
C     MAXRAD     O   Radius of outer bounding sphere for body.
C
C$ Detailed_Input
C
C     BODYID     is the body ID of the target for which radii
C                of bounding spheres are to be generated.
C
C     NSURF,
C     SRFLST     are, respectively, the surface list count and
C                an array containing a list of surface IDs.
C
C                If the count is zero, all surfaces for the body
C                are considered applicable.
C 
C$ Detailed_Output
C
C     MINRAD     is the radius of an inner bounding sphere for the 
C                surface of the body designated by BODYID, NSURF,
C                and SRFLST. The sphere is centered at the target
C                body's center. All points of the body's surface 
C                are outside this sphere.
C              
C                MINRAD is not necessarily a maximum lower bound.
C
C                Units are km.
C
C                
C     MAXRAD     is the radius of an outer bounding sphere for the 
C                surface of the body designated by BODYID, NSURF,
C                and SRFLST. The sphere is centered at the target
C                body's center. All points of the body's surface 
C                are contained in this sphere.
C
C                MINRAD is not necessarily a minimum upper bound.
C
C                Units are km.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If, for a DSK segment's reference frame, frame information
C         cannot be looked up, the error SPICE(NOFRAMEDATA) is
C         signaled.
C
C     2)  If, for a DSK segment's reference frame, the frame's name
C         cannot be looked up, the error SPICE(FRAMENAMENOTFOUND) is
C         signaled.
C
C     3)  If a DSK segment descriptor has an unrecognized coordinate
C         system code, the error SPICE(NOTSUPPORTED) is signaled.
C
C     4)  If no DSK segments are found for the specified target 
C         and surface set, the error SPICE(DSKDATANOTFOUND) is
C         signaled.
C
C     5)  If an error occurs while looking up DSK descriptors,
C         the error will be diagnosed by routines in the call
C         tree of this routine.
C
C     6)  If the surface list size is negative, the error 
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center. The lookup epoch for a segment is the 
C     midpoint of the segment's time coverage interval, so SPK data
C     must be available for that epoch.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup. 
C     
C$ Particulars
C
C     This routine is used as an initialization step by the ZZDSKSBF
C     entry point ZZSUDSKI. 
C
C     The operation of this routine usually will result in physical
C     file reads, and so will be rather slow. 
C
C$ Examples
C
C     See usage in ZZSUDSKI.
C
C$ Restrictions
C
C     1)  This routine does not take time into account when it 
C         selects segments. All segments for the specified body
C         and surface list are selected.
C
C         This functionality might be inappropriate for some 
C         future application. 
C
C     2)  This is a private routine. It is meant to be used only by the
C         DSK subsystem.
C      
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 26-JUL-2016 (NJB) 
C
C        Updated to use reference ellipsoid parameters and altitude
C        bounds to compute candidate values for bounding radii. This
C        applies to the planetodetic coordinate system.
C
C        30-JUN-2016 (NJB) 
C
C        30-JAN-2015 (NJB)
C
C           Updated to provide inner radius in addition to outer
C           radius. Updated to support segment frame centers that don't
C           coincide with the target.
C      
C        14-JAN-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find inner and outer bounding spheres for target body
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      VNORM

      INTEGER               ISRCHI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     EXTERNAL routines
C
      EXTERNAL              ZZDSKBDC

C
C     Local parameters
C
      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

C
C     Local variables
C
      CHARACTER*(LMSGLN)    ERRMSG
      CHARACTER*(FRNMLN)    FRNAME

      DOUBLE PRECISION      BOXCTR ( 3 )
      DOUBLE PRECISION      BOXRAD
      DOUBLE PRECISION      CTRMNR
      DOUBLE PRECISION      DSKDSC ( DSKDSZ )
      DOUBLE PRECISION      F
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      LX
      DOUBLE PRECISION      LY
      DOUBLE PRECISION      LZ
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      MIDTIM
      DOUBLE PRECISION      MINR
      DOUBLE PRECISION      OFFMAG
      DOUBLE PRECISION      OFFSET ( 3 )      
      DOUBLE PRECISION      RE
      DOUBLE PRECISION      RP
      DOUBLE PRECISION      SGMAXR
      DOUBLE PRECISION      SGMINR
      DOUBLE PRECISION      SVMAXR
      DOUBLE PRECISION      SVMINR

      INTEGER               CORSYS
      INTEGER               CTR    ( CTRSIZ )
      INTEGER               DLADSC ( DLADSZ )
      INTEGER               FRAMID
      INTEGER               FRCENT
      INTEGER               FRCLID
      INTEGER               FRCLAS
      INTEGER               HANDLE
      INTEGER               I
      INTEGER               PRVBOD
      INTEGER               PRVFID
      INTEGER               PRVLST ( MAXSRF )
      INTEGER               PRVNLS
      INTEGER               SURFID

      LOGICAL               FIRST
      LOGICAL               FOUND
      LOGICAL               NEWLST
      LOGICAL               SAME
      LOGICAL               SEGFND
      LOGICAL               UPDATE

C
C     Saved variables
C
      SAVE                  CTR
      SAVE                  FIRST
      SAVE                  PRVBOD
      SAVE                  PRVLST
      SAVE                  PRVNLS
      SAVE                  SVMAXR
      SAVE                  SVMINR

C
C     Initial values
C
      DATA                  CTR    / -1, -1      /
      DATA                  FIRST  / .TRUE.      /
      DATA                  PRVFID /  0          /
      DATA                  PRVBOD /  0          /
      DATA                  PRVLST /  MAXSRF * 0 /
      DATA                  PRVNLS / -1          /
      DATA                  SVMAXR / -1.D0       /
      DATA                  SVMINR / -1.D0       /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKSPH' )

      IF ( FIRST ) THEN
         CALL ZZCTRUIN ( CTR )
      END IF

C
C     Check NSURF.
C     
      IF ( NSURF .LT. 0 ) THEN

         CALL SETMSG ( 'NSURF must be non-negative but was #.' )
         CALL ERRINT ( '#', NSURF                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                )
         CALL CHKOUT ( 'ZZDSKSPH'                              )
         RETURN
 
      END IF

C
C     Determine whether the input body surface list matches
C     the previous values. The following code applies whether
C     or not the surface list is non-empty.
C
      NEWLST = .TRUE.

      IF ( .NOT. FIRST ) THEN

         IF ( BODYID .EQ. PRVBOD ) THEN

            IF ( NSURF .EQ. PRVNLS ) THEN

               SAME = .TRUE.
               I    = 1

               DO WHILE ( ( I .LE. NSURF ) .AND.  SAME  )

                  SAME = SRFLST(I) .EQ. PRVLST(I) 
                  I    = I + 1
               END DO
C
C              If SAME is true here, the body and surface list are the
C              same as on the previous call.
C
               NEWLST = .NOT. SAME

            END IF

         END IF

      END IF

C
C     Set PRVNLS to a value that can't match a valid value, so
C     the surface list won't match after an error occurs. We'll
C     reset PRVNLS prior to exit if all goes well.
C     
      PRVNLS = -1 

C
C     Check for DSK update in ZZDSKBSR.
C
      CALL ZZDSKCHK ( CTR, UPDATE )

C
C     Initialize the temporary variables MINR, MAXR. 
C
      MINR = SVMINR
      MAXR = SVMAXR


      IF ( FIRST .OR. UPDATE .OR. NEWLST ) THEN
C
C        Initialize the saved radius data.
C         
         SVMAXR = -1.D0
         SVMINR = DPMAX()

C
C        Prepare to fetch segment data. Initialize the ZZDSKBSR
C        segment list for the body of interest.
C
         CALL ZZDSKBBL ( BODYID )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZDSKSPH' )
            RETURN
         END IF

C
C        Fetch segment DSK descriptors for the indicated body and
C        surface list.
C
         PRVFID = 0
         CALL CLEARD ( 3, OFFSET )

C
C        Re-initialize MINR and MAXR.
C
         MAXR = -1.D0
         MINR = DPMAX()

C
C        Examine all segments for BODYID.
C
         CALL ZZDSKBSS ( BODYID )

         CALL ZZDSKSBD ( BODYID )
         CALL ZZDSKSNS ( ZZDSKBDC, HANDLE, DLADSC, DSKDSC, SEGFND )

         DO WHILE ( SEGFND )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZDSKSPH' )
               RETURN
            END IF

            IF ( NSURF .GT. 0 ) THEN

               SURFID = NINT( DSKDSC(SRFIDX) )

               I      = ISRCHI ( SURFID, NSURF, SRFLST )
            ELSE
               I      = 1
            END IF

            IF ( I .GT. 0 ) THEN
C
C              If we're checking surface IDs, this segment qualifies.
C              Otherwise, we're not checking surface IDs, so the segment
C              qualifies by default.
C
C              Get the frame ID of this segment, and look up the frame's
C              center.
C
               FRAMID = NINT( DSKDSC(FRMIDX) )

               IF ( FRAMID .NE. PRVFID ) THEN
C
C                 Get the frame center for the current segment.
C
                  CALL FRINFO ( FRAMID, FRCENT, FRCLAS, FRCLID, FOUND )
               
                  IF ( .NOT. FOUND ) THEN

                     CALL SETMSG ( 'No frame specification was found '
     .               //            'for frame ID #.'                  )
                     CALL ERRINT ( '#', FRAMID                        )
                     CALL SIGERR ( 'SPICE(NOFRAMEDATA)'               )
                     CALL CHKOUT ( 'ZZDSKSPH'                         )
                     RETURN

                  END IF

                  IF ( FRCENT .EQ. BODYID ) THEN
C
C                    The frame is centered at the target, so
C                    the frame center offset magnitude is zero.
C
                     OFFMAG = 0.D0

                  ELSE

                     CALL FRMNAM ( FRAMID, FRNAME ) 

                     IF ( FAILED() ) THEN
                        CALL CHKOUT ( 'ZZDSKSPH' )
                        RETURN
                     END IF

                     IF ( FRNAME .EQ. ' ' ) THEN

                        CALL SETMSG ( 'No frame name was found '
     .                  //            'for frame ID #.'          )
                        CALL ERRINT ( '#', FRAMID                )
                        CALL SIGERR ( 'SPICE(FRAMENAMENOTFOUND)' )
                        CALL CHKOUT ( 'ZZDSKSPH'                 )
                        RETURN

                     END IF

                     MIDTIM = ( DSKDSC(BTMIDX) + DSKDSC(ETMIDX) ) / 2

                     CALL SPKGPS ( FRCENT, MIDTIM, FRNAME, 
     .                             BODYID, OFFSET, LT     )

                     IF ( FAILED() ) THEN
                        CALL CHKOUT ( 'ZZDSKSPH' )
                        RETURN
                     END IF

                     OFFMAG = VNORM( OFFSET )

                  END IF

               END IF

C
C              Get the segment coordinate system and derive the maximum
C              radius of the segment.
C
               CORSYS = NINT( DSKDSC(SYSIDX) )

C
C              Get bounding radii for the segment relative to the
C              origin of the segment's coordinate system. We'll account
C              for the offset of the origin from the segment's central
C              body as a subsequent step.
C
               IF ( CORSYS .EQ. LATSYS ) THEN

                  SGMINR = DSKDSC(MN3IDX)
                  SGMAXR = DSKDSC(MX3IDX) 


               ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C                 Use the reference spheroid and altitude bounds to
C                 generate initial bounding radii.
C                 
                  RE = DSKDSC(PARIDX  )
                  F  = DSKDSC(PARIDX+1)
                  RP = RE *  (1.D0 - F)

                  IF ( F .GE. 0.D0 ) THEN
C
C                    The spheroid is oblate. The maximum altitude over
C                    the equator is an upper bound for the distance of
C                    any surface point from the origin. The minimum 
C                    altitude over either pole is a lower bound for 
C                    the distance of any surface point from the origin.
C
C                    The DSK descriptor gives us the altitude bounds. 
C
                     SGMAXR = RE + DSKDSC(MX3IDX) 
                     SGMINR = RP + DSKDSC(MN3IDX) 

                  ELSE

C                    The spheroid is prolate. The maximum altitude over
C                    either pole is an upper bound for the distance of
C                    any surface point from the origin.
C
                     SGMAXR = RP + DSKDSC(MX3IDX) 
                     SGMINR = RE + DSKDSC(MN3IDX) 

                  END IF



               ELSE IF ( CORSYS .EQ. RECSYS ) THEN

                  CALL ZZRECBOX ( DSKDSC(MN1IDX), BOXCTR, LX, 
     .                            LY,             LZ,     BOXRAD )
C
C                 SGMINR is a lower bound on the distance of the 
C                 segment from the origin of the coordinate system.
C
                  SGMINR = MAX ( VNORM(BOXCTR) - BOXRAD,  0.0D0 )

                  SGMAXR =       VNORM(BOXCTR) + BOXRAD


               ELSE

                  CALL SETMSG ( 'Coordinate system # is not ' 
     .            //            'currently supported.'        )
                  CALL ERRINT ( '#', CORSYS                   )
                  CALL SIGERR ( 'SPICE(NOTSUPPORTED)'         )
                  CALL CHKOUT ( 'ZZDSKSPH'                    )
                  RETURN

               END IF

C
C              Apply the triangle inequality to derive minimum and
C              maximum values of the distance of the surface from the
C              body center, given the offset between the frame center
C              and the body center, and given bounds on the distance of
C              the surface from the frame's center.
C
               IF ( OFFMAG .LE. SGMINR ) THEN
C
C                 The segment's central body is inside the inner
C                 bounding sphere of the segment.
C
                  CTRMNR = SGMINR - OFFMAG 

               ELSE IF ( OFFMAG .GE. SGMAXR ) THEN
C
C                 The segment's central body is outside the outer
C                 bounding sphere of the segment.
C
                  CTRMNR = OFFMAG - SGMAXR 

               ELSE
C
C                 The segment's central body is between the bounding
C                 spheres. No positive lower radius bound exists.
C
                  CTRMNR = 0.D0

               END IF

C
C              Update the segment's outer bounding radius to 
C              account for the frame center offset (which may
C              be zero).
C
               SGMAXR = SGMAXR + OFFMAG

C
C              Update the global minimum and maximum radii.
C
               MINR = MIN( MINR, CTRMNR )
               MAXR = MAX( MAXR, SGMAXR )

            END IF
C
C           Look at the next segment.
C
            CALL ZZDSKSBD( BODYID )
            CALL ZZDSKSNS( ZZDSKBDC, HANDLE, DLADSC, DSKDSC, SEGFND )

         END DO

         IF (  ( MAXR .GT. 0.D0 ) .AND. (.NOT. FAILED() )  ) THEN
C
C           Update the saved bounds.
C           
            SVMINR = MINR
            SVMAXR = MAXR

         END IF

      END IF


      IF ( MAXR .LT. 0.D0 ) THEN
C
C        We tried to update the radius bounds but didn't find any
C        segments for the specified body.
C
C        We have no radius data for the specified surface list.
C
         IF ( NSURF .EQ. 0 ) THEN

            ERRMSG = 'No segments were found matching the body ID #.'

         ELSE

            ERRMSG = 'No segments were found matching the body ID # '
     .      //       'and the surface list <@>.'

            DO I = 1, NSURF-1

               CALL REPMC  ( ERRMSG, '@', '*, @',    ERRMSG )
               CALL REPMI  ( ERRMSG, '*', SRFLST(I), ERRMSG )

            END DO

            CALL REPMI  ( ERRMSG, '@', SRFLST(NSURF), ERRMSG )

         END IF

         CALL SETMSG ( ERRMSG                   )
         CALL ERRINT ( '#', BODYID              )
         CALL SIGERR ( 'SPICE(DSKDATANOTFOUND)' )
         CALL CHKOUT ( 'ZZDSKSPH'               )
         RETURN

      END IF
      
      IF ( .NOT. FAILED() ) THEN

         FIRST  = .FALSE.

         PRVBOD = BODYID
         PRVNLS = NSURF

         IF ( NEWLST ) THEN

            CALL MOVEI ( SRFLST, NSURF, PRVLST )

         END IF

         MAXRAD = SVMAXR
         MINRAD = SVMINR

      END IF

      CALL CHKOUT ( 'ZZDSKSPH' )
      RETURN
      END
