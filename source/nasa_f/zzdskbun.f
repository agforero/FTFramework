C$Procedure ZZDSKBUN ( DSK, buffered unprioritized normal vector )
 
      SUBROUTINE ZZDSKBUN ( BODYID, NSURF,  SRFLST, ET,     FIXFID,
     .                      NSEG,   HANBUF, DLABUF, DSKBUF, OFFBUF,
     .                      CTRBUF, RADBUF, POINT,  NORMAL         )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find the outward unit normal vector at a specified point on the
C     surface of a body represented by one or more DSK segments. The
C     set of surface IDs to be considered is specified by the caller.
C
C     This routine uses DSK segments in an unprioritized manner.
C     All segments meeting the body, time, and surface constraints
C     are considered.
C
C     Segment descriptor and derived bounding data are passed in.
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
C     GEOMETRY
C     NORMAL
C     RAY
C     SURFACE
C     TOPOGRAPHY
C     VECTOR
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsktol.inc'
      INCLUDE 'zzctr.inc'

      INTEGER               BODYID
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      DOUBLE PRECISION      ET
      INTEGER               FIXFID
      INTEGER               NSEG
      INTEGER               HANBUF ( * )
      INTEGER               DLABUF ( DLADSZ, * )
      DOUBLE PRECISION      DSKBUF ( DSKDSZ, * )
      DOUBLE PRECISION      OFFBUF ( 3,      * )
      DOUBLE PRECISION      CTRBUF ( 3,      * )
      DOUBLE PRECISION      RADBUF ( * )
      DOUBLE PRECISION      POINT  ( 3 )
      DOUBLE PRECISION      NORMAL ( 3 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body.
C     NSURF      I   Number of IDs in surface list.
C     SRFLST     I   List of surface IDs.
C     ET         I   Lookup epoch, expressed as seconds past J2000 TDB.
C     FIXFID     I   Frame ID of body-fixed frame for output vectors.
C     NSEG       I   Number of segments in buffers.
C     HANBUF     I   DSK handle buffer.
C     DLABUF     I   DLA segment descriptor buffer.
C     DSKBUF     I   DSK segment descriptor buffer.
C     OFFBUF     I   Frame center offset buffer.
C     CTRBUF     I   Bounding sphere center buffer.
C     RADBUF     I   Bounding sphere radius buffer.
C     POINT      I   Surface point.
C     NORMAL     O   Outward unit normal vector at surface point.
C
C$ Detailed_Input
C
C     BODYID     is the body ID of the target on which the input
C                point is located.
C
C
C     NSURF,
C     SRFLST     are, respectively, the surface list count and
C                an array containing a list of surface IDs.
C
C                If the count is zero, all surfaces for the body
C                are considered applicable.
C
C
C     ET         is the lookup epoch, specified as seconds past
C                J2000 TDB. Only DSK segments containing ET in
C                their time coverage intervals are considered
C                in the normal vector computation. 
C
C
C     FIXFID     is the frame ID of a body-fixed reference frame
C                centered on the target body.
C
C
C     NSEG       is the number of DSK segments for which the input
C                buffers contain data.
C
C
C     HANBUF     is the DSK handle buffer.
C
C
C     DLABUF     is the DLA segment descriptor buffer.
C
C
C     DSKBUF     is the DSK segment descriptor buffer.
C
C
C     OFFBUF     is the frame offset buffer. For each DSK segment, the
C                entry in this buffer is the position of that segment's
C                reference frame center relative to the target body's
C                center. The vector is expressed in the segment's
C                reference frame. Units are km.
C
C
C     CTRBUF     is the bounding sphere center buffer. For each DSK
C                segment, the entry in this buffer is the position of
C                that segment's outer bounding sphere center relative
C                to the segment's reference frame center. The vector is
C                expressed in the segment's reference frame. Units are
C                km.
C
C
C     RADBUF     is the bounding sphere radius buffer. For each DSK
C                segment, the entry in this buffer is the radius of the
C                segment's outer bounding sphere. Units are km.
C
C
C     POINT      is a point on the surface of the target body. POINT is
C                expressed in the reference frame designated by FIXFID.
C                POINT represents an offset from the target body
C                center.
C
C                POINT must be located on the surface of the target
C                body, within a small tolerance.
C                
C
C$ Detailed_Output
C
C     NORMAL     is the unit length, outward normal vector on the
C                surface of the target body at POINT. NORMAL is
C                expressed in the reference frame designated by FIXFID.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the input segment count is non-positive, the error
C         SPICE(NODSKSEGMENTS) is signaled.
C          
C     2)  If an unrecognized coordinate system code is found
C         in a DSK descriptor, the error SPICE(BADCOORDSYS) is
C         signaled.
C
C     3)  If the input point is contained in more than MAXHIT
C         segments, the error SPICE(TOOMANYHITS) is signaled.
C
C     4)  If the input point is not contained in any DSK segment,
C         the error SPICE(POINTNOTINSEGMENT) is signaled.
C
C     5)  If an unrecognized segment data type is encountered,
C         the error SPICE(NOTSUPPORTED) is signaled.
C
C     6)  If the input point is not sufficiently close to the 
C         surface, the error SPICE(POINTOFFSURFACE) will be 
C         signaled.
C
C     7)  Any other errors that occur while looking up DSK data or
C         mapping the surface point to a normal vector will be signaled
C         by routines in the call tree of this routine.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem.
C
C     Frame kernels may be required in order to perform transformations
C     from DSK segment frames to the output frame. In some cases,
C     additional kernels such as CK kernels and SCLK kernels could be
C     required to support the offset vector lookup.
C     
C$ Particulars
C
C     This routine is meant to be called from the DSK API segment
C     buffering subsystem. All of the input buffers should contain
C     values computed after changes to the set of loaded DSKs.
C
C     This routine takes advantage of the input segment bound
C     information to efficiently determine which segments may contain
C     the input point.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by the DSK
C     subsystem.
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
C-    SPICELIB Version 1.0.0, 22-FEB-2017 (NJB) 
C
C        Added FAILED calls.
C     
C        07-APR-2016 (NJB) 
C
C        Based on first version 17-DEC-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     compute normals on dsk surface
C
C-&



C
C     SPICELIB functions
C
      DOUBLE PRECISION      VDIST

      INTEGER               ISRCHI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C           
      INTEGER               MAXHIT
      PARAMETER           ( MAXHIT = 1000 )


C
C     Local variables
C     
      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      SEGPT  ( 3 )
      DOUBLE PRECISION      SGMARG
      DOUBLE PRECISION      SGXBUF ( 3, 3, MAXHIT )
      DOUBLE PRECISION      VERTS  ( 3, 3 )
      DOUBLE PRECISION      VTEMP  ( 3 )
      DOUBLE PRECISION      XFORM  ( 3, 3 )

      INTEGER               CORSYS
      INTEGER               DTYPE
      INTEGER               I
      INTEGER               J
      INTEGER               NHIT
      INTEGER               PLATE  ( 3 )
      INTEGER               PLID
      INTEGER               PRVFRM
      INTEGER               SEGFID
      INTEGER               SGHIT  ( MAXHIT )
      INTEGER               SURFCE
      INTEGER               WINNER

      LOGICAL               BODYOK
      LOGICAL               DONE
      LOGICAL               INSIDE
      LOGICAL               MULTFR
      LOGICAL               SURFOK 
      LOGICAL               TIMEOK
      LOGICAL               XFND

C
C     Saved variables
C

C
C     Initial values
C

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKBUN' )

C
C     Check the incoming segment count.
C
      IF ( NSEG .LE. 0 ) THEN

         CALL SETMSG ( 'Input segment list was empty. This may '
     .   //            'be due to no DSKs containing data for '
     .   //            'body # having been loaded.'            )
         CALL ERRINT ( '#',  BODYID                            )
         CALL SIGERR ( 'SPICE(NODSKSEGMENTS)'                  )
         CALL CHKOUT ( 'ZZDSKBUN'                              )
         RETURN

      END IF

C
C     Get the segment margin from the tolerance database.
C
      CALL DSKGTL ( KEYSGR, SGMARG )

C
C     Indicate we haven't yet seen a segment frame different
C     from the one designated by FIXFID.
C
      MULTFR = .FALSE.

C
C     Make a local copy of the point. We'll update this copy
C     later if need be.
C
      CALL VEQU ( POINT, SEGPT )

C
C     By default, each segment in the input list must be checked
C     for intersection.
C
C     We start out by trying to eliminate segments from consideration
C     by comparing the point's distance from their centers to the radii
C     of their bounding spheres. Only those segments whose bounding
C     spheres contain the point are examined further.
C
      NHIT   = 0
      PRVFRM = 0

      DO I = 1, NSEG
C
C        BODYOK indicates whether the input body ID matches that 
C        of the current segment.
C
         SURFOK = .FALSE.
         TIMEOK = .FALSE.

         BODYOK = BODYID .EQ. NINT( DSKBUF(CTRIDX,I) )


         IF ( BODYOK ) THEN
C
C           See whether the current segment contains a surface we're
C           supposed to consider. If the surface list is empty, we
C           consider all surfaces. Otherwise, the surface of the
C           segment must be on the surface list in order to qualify.
C
            J = 0

            IF ( NSURF .GT. 0 ) THEN

               SURFCE = NINT( DSKBUF(SRFIDX,I) )
               J      = ISRCHI( SURFCE, NSURF, SRFLST )

            END IF

            SURFOK = ( NSURF .EQ. 0 ) .OR. ( J .GT. 0 )

C
C           See whether the segment covers the input epoch.
C
            TIMEOK =        ( ET .GE. DSKBUF(BTMIDX,I) ) 
     .                .AND. ( ET .LE. DSKBUF(ETMIDX,I) ) 

         END IF


         IF ( BODYOK .AND. SURFOK .AND. TIMEOK  ) THEN
C
C           This segment is to be considered.
C
C           In order to do any geometric comparison, the point must be
C           in the same frame as the segment we're checking. Transform
C           the point if need be.
C
C           Get the segment's frame ID. Get the transformation from the
C           input frame to the output frame if needed.

            SEGFID = NINT( DSKBUF(FRMIDX,I) )
            

            IF ( SEGFID .NE. FIXFID ) THEN
C
C              We have a segment that uses a different frame
C              from that specified by FIXFID.
C              
               MULTFR = .TRUE.

               IF ( SEGFID .NE. PRVFRM ) THEN
C
C                 The frame of the current segment doesn't match
C                 that of the previous segment, so we'll need
C                 to look up the transformation from the input
C                 frame to the segment frame.
C
C                 Otherwise, XFORM already contains the correct
C                 transformation.
C
                  CALL REFCHG ( FIXFID, SEGFID, ET, XFORM ) 
                  
                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'ZZDSKBUN' )
                     RETURN
                  END IF

C
C                 Transform the local copy of the point to the segment's
C                 frame, and shift the local copy of the point so
C                 that it represents an offset relative to center of
C                 the segment's frame.
C
                  CALL MXV  ( XFORM, POINT,  SEGPT  )

C
C                 The direction of the buffered offset is from the body
C                 to the segment frame's center. The offset is
C                 expressed in the segment's frame.
C
                  CALL VSUB ( SEGPT, OFFBUF(1,I), VTEMP )
                  CALL VEQU ( VTEMP,              SEGPT )

               END IF


            ELSE IF ( MULTFR ) THEN
C
C              The input and segment frames are the same for this
C              segment, but the current value of SEGPT needs to be
C              reset.
C
               CALL VEQU ( POINT, SEGPT )

            END IF

C
C           Find the distance of the point from the "center" of the
C           segment's coverage volume.
C
            DIST = VDIST ( CTRBUF(1,I), SEGPT )


            IF ( DIST .LE. RADBUF(I) ) THEN
C
C              The point is inside or on the bounding surface. We'll
C              check the boundary of the segment for an intersection.
C
               CORSYS = NINT( DSKBUF(SYSIDX,I) )

               IF ( CORSYS .EQ. LATSYS ) THEN

                  CALL ZZINLAT ( SEGPT,  DSKBUF(MN1IDX,I),
     .                           SGMARG, 0,                INSIDE )

               ELSE IF ( CORSYS .EQ. RECSYS ) THEN

                  CALL ZZINREC ( SEGPT,  DSKBUF(MN1IDX,I),
     .                           SGMARG, 0,                INSIDE )
 
               ELSE IF ( CORSYS .EQ. PDTSYS ) THEN

                  CALL ZZINPDT ( SEGPT,             DSKBUF(MN1IDX,I), 
     .                           DSKBUF(PARIDX, I), SGMARG,
     .                           0,                 INSIDE           )
               ELSE

                  CALL SETMSG( 'Coordinate system # is not supported.' )
                  CALL ERRINT( '#', CORSYS                             )
                  CALL SIGERR( 'SPICE(BADCOORDSYS)'                    )
                  CALL CHKOUT( 'ZZDSKBUN'                              )
                  RETURN

               END IF              

               IF ( FAILED() ) THEN
                  CALL CHKOUT( 'ZZDSKBUN' )
                  RETURN
               END IF


               IF ( INSIDE ) THEN                   
C
C                 The point in inside the region enclosed by the
C                 boundary of this segment. Save the index of the
C                 segment in the "hit list."
C
                  IF ( NHIT .EQ. MAXHIT ) THEN

                     CALL SETMSG ( 'Too many segments contain the '
     .               //            'input point. Buffer size is #.'  )
                     CALL ERRINT ( '#', MAXHIT                       )
                     CALL SIGERR ( 'SPICE(TOOMANYHITS)'              )
                     CALL CHKOUT ( 'ZZDSKBUN'                        )
                     RETURN

                  END IF

                  NHIT           = NHIT + 1
                  SGHIT ( NHIT ) = I
C
C                 Save the frame transformation for this segment.
C
                  CALL MOVED ( XFORM, 9, SGXBUF(1,1,NHIT) )

               END IF


            END IF
C
C           The current segment matched the input criteria.
C
C           Update the saved segment frame ID to that of the segment
C           we just examined.
C
            PRVFRM = SEGFID

         END IF

      END DO
      
C
C     We have an error if no segments were hit. None of that
C     "not found" nonsense here.
C
      IF ( NHIT .EQ. 0 ) THEN

         CALL SETMSG ( 'Input point (# # #) in frame # does not '
     .   //            'lie inside any segment for the '
     .   //            'specified body (#) and surfaces.'        )
         CALL ERRDP  ( '#', POINT(1)                             )
         CALL ERRDP  ( '#', POINT(2)                             )
         CALL ERRDP  ( '#', POINT(3)                             )
         CALL ERRINT ( '#', FIXFID                               )
         CALL ERRINT ( '#', BODYID                               )
         CALL SIGERR ( 'SPICE(POINTNOTINSEGMENT)'                )
         CALL CHKOUT ( 'ZZDSKBUN' )
         RETURN

      END IF
       
 
C
C     Now process the segments on the hit list. If we find a segment
C     surface on which the input point is located, we compute the
C     normal vector and terminate the search.
C
      I      = 1
      DONE   = .FALSE.
      WINNER = 0
      PRVFRM = 0

      DO WHILE ( .NOT. DONE )
C
C        I is the index in the hit list of the segment
C        we're considering. J is the index of that segment
C        in the parallel input arrays.
C
         J = SGHIT (I)


         SEGFID = NINT( DSKBUF(FRMIDX,J) )

         IF ( SEGFID .NE. FIXFID ) THEN

            IF ( SEGFID .NE. PRVFRM ) THEN
C
C              Transform and shift the input point.
C
C              Here I is an index in the hit list and J is 
C              an index in the input arrays.
C  
               CALL MOVED ( SGXBUF(1,1,I), 9, XFORM )

               CALL MXV  ( XFORM, POINT, SEGPT )
  
               CALL VSUB ( SEGPT, OFFBUF(1,J), VTEMP )
               CALL VEQU ( VTEMP,              SEGPT )

            END IF

         ELSE IF ( MULTFR ) THEN

            CALL VEQU ( POINT, SEGPT )

         END IF

C
C        If the point lies on the surface described by the
C        current segment, find the outward unit normal 
C        vector at the point.
C
         DTYPE = NINT( DSKBUF(TYPIDX,J) )
         XFND  = .FALSE.

         IF ( DTYPE .EQ. 2 ) THEN
C
C           Find the plate on which the point lies, if any.
C
            CALL ZZPTPL02 ( HANBUF(J), DLABUF(1,J), DSKBUF(1,J), SEGPT, 
     .                      PLID,      PLATE,       VERTS,       XFND  )

            IF ( FAILED() ) THEN
               CALL CHKOUT( 'ZZDSKBUN' )
               RETURN
            END IF

            IF ( XFND ) THEN
C
C              Find the unit outward normal at SEGPT. We must
C              convert the output of PLTNRM to unit length.
C               
               CALL PLTNRM ( VERTS(1,1), VERTS(1,2), VERTS(1,3),
     .                       NORMAL                              )
               IF ( FAILED() ) THEN
                  CALL CHKOUT( 'ZZDSKBUN' )
                  RETURN
               END IF

               CALL VHATIP ( NORMAL )
C
C              WINNER is the index in the hit list of the current
C              segment.
C              
               WINNER = I
               DONE   = .TRUE.

            END IF


         ELSE
            
            CALL SETMSG ( 'Segment type is #; this type is not '
     .      //            'currently supported.'                )
            CALL ERRINT ( '#', DTYPE                            )
            CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                 )
            CALL CHKOUT ( 'ZZDSKBUN'                            )
            RETURN
            
         END IF
 

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZDSKBUN' )
            RETURN
         END IF

 

         IF ( .NOT. DONE ) THEN

            IF ( I .EQ. NHIT ) THEN
C
C              We've looked at all of the segments.
C
               DONE = .TRUE.

            ELSE
C
C              Consider the next segment.
C
               I = I + 1

            END IF

         END IF

      END DO

C
C     If we have an result, it may be represented in a frame
C     other than the input point frame. Transform it back to the
C     input frame, if need be.
C
      IF ( XFND ) THEN
C
C        J is the index in the input arrays of the "winning" segment.
C
         J      = SGHIT ( WINNER )         

         SEGFID = NINT( DSKBUF(FRMIDX,J) )

         IF ( SEGFID .NE. FIXFID ) THEN
C
C           The segment frame and input frame differ. The normal
C           vector must be converted back to the input frame.
C
            CALL MOVED( SGXBUF(1,1,WINNER), 9, XFORM )

            CALL MTXV ( XFORM, NORMAL, VTEMP  )
            CALL VEQU ( VTEMP,         NORMAL )

         END IF

      ELSE
C
C        The input point was not legitimate; otherwise we 
C        would have found a solution.
C
         CALL SETMSG ( 'Input point (# # #) in frame # does not '
     .   //            'lie on the surface contained in any '
     .   //            'segment for the specified body (#) '
     .   //            'and surfaces.'                           )
         CALL ERRDP  ( '#', POINT(1)                             )
         CALL ERRDP  ( '#', POINT(2)                             )
         CALL ERRDP  ( '#', POINT(3)                             )
         CALL ERRINT ( '#', FIXFID                               )
         CALL ERRINT ( '#', BODYID                               )
         CALL SIGERR ( 'SPICE(POINTOFFSURFACE)'                  )
         CALL CHKOUT ( 'ZZDSKBUN' )
         RETURN
         
      END IF
 
      CALL CHKOUT ( 'ZZDSKBUN' )
      RETURN
      END
