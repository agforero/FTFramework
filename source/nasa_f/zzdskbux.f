C$Procedure ZZDSKBUX ( DSK, buffered unprioritized ray intercept )
 
      SUBROUTINE ZZDSKBUX ( BODYID, NSURF,  SRFLST, ET,     FIXFID,
     .                      NSEG,   HANBUF, DLABUF, DSKBUF, OFFBUF,
     .                      CTRBUF, RADBUF, VERTEX, RAYDIR, XPT, 
     .                      SEGIDX, DC,     IC,     FOUND          )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find a ray-surface intercept on the surface of a body represented
C     by one or more DSK segments. If multiple intercepts exist, select
C     the one closest to the ray's vertex. The set of surface IDs to be
C     considered is specified by the caller.
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
C     INTERCEPT
C     INTERSECTION
C     RAY
C     SURFACE
C     TOPOGRAPHY
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
      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      XPT    ( 3 )
      INTEGER               SEGIDX
      DOUBLE PRECISION      DC     ( * )
      INTEGER               IC     ( * )
      LOGICAL               FOUND
      
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
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     XPT        O   Intercept point.
C     SEGIDX     O   Index of segment in input buffers.
C     FOUND      O   Found flag. True if and only if intercept exists.
C     SGREED     P   Default segment boundary margin.
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
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray must be represented in the
C                reference frame designated by FIXFID. The vertex is
C                considered to be an offset from the target body.
C                
C
C$ Detailed_Output
C
C
C     XPT        is the surface intercept, if an intercept exists. If
C                multiple intercepts exist, the one closest to the 
C                ray's vertex is selected. XPT is represented in the
C                reference frame designated by FIXFID. XPT is
C                considered to be an offset from the target body.
C
C
C     SEGIDX     is the index, within the input buffer arguments, of
C                the segment providing the surface data on which XPT
C                was found. SEGIDX is valid if and only if FOUND is
C                .TRUE.
C
C
C     FOUND      is a logical flag that is .TRUE. if and only if a 
C                ray-surface intercept was found.
C
C$ Parameters
C
C     SGREED     is the default margin used to determine whether
C                segments should be tested for intersection. Segments
C                are effectively expanded slightly for intersection
C                tests, using this margin. For example, if a segment
C                has latitudinal coordinates, and if the ray hits the
C                eastern longitude boundary of the segment, the
C                latitude of this intercept is within the segment's
C                latitude bounds, RMAX is the segment's maximum radius,
C                and the radius R of the intercept satisfies
C
C                   RMAX < R < RMAX * (1+SGREED)
C
C                then the intercept is considered to be "on" the
C                segment's boundary, and the segment's surface data
C                are tested for intersection.
C
C                See the routines that test points for inclusion in
C                segments for details concerning use of the segment
C                margin. For latitudinal, planetodetic, and 
C                rectangular coordinates respectively, these are
C
C                   ZZINLAT
C                   ZZINPDT
C                   ZZINREC
C
C                See the include file dsktol.inc for the value of
C                this parameter.
C
C                This parameter can be overridden. See the include
C                file dsktol.inc and the routine DSKSTL for details.
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
C         segments, the error SPICE(BUFFERTOOSMALL) is signaled.
C
C     4)  Any other errors that occur while looking up DSK data or
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
C     information to efficiently determine which segments may be
C     intersected by the input ray.
C
C     This routine uses DSK data in an unprioritized fashion: all
C     DSK segments having the required body, surface, and time
C     attributes are considered in the computation.
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
C-    SPICELIB Version 1.0.0, 19-JUL-2016 (NJB) 
C
C        30-JAN-2015 (NJB)
C
C        Updated argument list: now returns segment
C        index and source info arrays. Added header.
C        Added support for rectangular coordinates.
C
C        30-JAN-2015 (NJB)
C
C        Updated to accommodate change in argument list of 
C        ZZDSKSPH.
C
C        28-JAN-2015 (NJB)
C
C        Now uses outer bounding sphere to produce vertex near the
C        target for individual ray-segment intercept computations.
C
C        17-DEC-2014 (NJB)
C
C        Updated to allow DSK segments to use body-fixed frames
C        not centered on the target.
C     
C        10-OCT-2014 (NJB)
C
C        First version.
C
C-&
 
C$ Index_Entries
C
C     buffered unprioritized ray surface intercept 
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VNORM

      INTEGER               ISRCHI

      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local parameters
C     
      DOUBLE PRECISION      RADSCL
      PARAMETER           ( RADSCL = 1.01D0 )

      INTEGER               MAXHIT
      PARAMETER           ( MAXHIT = 1000 )

C
C     Local variables
C     
      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      DMIN
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      MINRAD
      DOUBLE PRECISION      SEGDIR ( 3 )
      DOUBLE PRECISION      SEGVTX ( 3 )
      DOUBLE PRECISION      PNEAR  ( 3 )
      DOUBLE PRECISION      SGDIST ( MAXHIT )
      DOUBLE PRECISION      SGMARG
      DOUBLE PRECISION      SGXBUF ( 3, 3, MAXHIT )
      DOUBLE PRECISION      SPHVTX ( 3 ) 
      DOUBLE PRECISION      SRFX   ( 3 )
      DOUBLE PRECISION      VTEMP  ( 3 )
      DOUBLE PRECISION      XFORM  ( 3, 3 )

      INTEGER               DTYPE
      INTEGER               I
      INTEGER               IORDER ( MAXHIT )
      INTEGER               J
      INTEGER               K
      INTEGER               NHIT
      INTEGER               NXPTS
      INTEGER               PRVFRM
      INTEGER               SEGFID
      INTEGER               SGHIT  ( MAXHIT )
      INTEGER               SURFCE
      INTEGER               WINNER
      
      LOGICAL               BODYOK
      LOGICAL               DONE
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

      CALL CHKIN ( 'ZZDSKBUX' )

C
C     Check the incoming segment count.
C
      IF ( NSEG .LE. 0 ) THEN

         CALL SETMSG ( 'Input segment list was empty. This may '
     .   //            'be due to no DSKs containing data for '
     .   //            'body # having been loaded.'            )
         CALL ERRINT ( '#',  BODYID                            )
         CALL SIGERR ( 'SPICE(NODSKSEGMENTS)'                  )
         CALL CHKOUT ( 'ZZDSKBUX'                              )
         RETURN

      END IF

C
C     No intercept has been found.
C
      SEGIDX = 0
      FOUND  = .FALSE.
      
C
C     Indicate we haven't yet seen a segment frame different
C     from the one designated by FIXFID.
C
      MULTFR = .FALSE.

C
C     Obtain the "greedy" segment margin.
C
      CALL DSKGTL ( KEYSGR, SGMARG )

C
C     Make a local copy of the ray. We'll update this copy
C     later if need be.
C
      CALL VEQU ( VERTEX, SEGVTX )
      CALL VEQU ( RAYDIR, SEGDIR )

C
C     Obtain the radius of an outer bounding sphere for the given body
C     and surface list.
C
      CALL ZZDSKSPH ( BODYID, NSURF, SRFLST, MINRAD, MAXRAD )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZDSKBUX' )
         RETURN
      END IF
C
C     Scale up the bounding sphere to avoid round-off difficulties.
C     We'll use this value in the loop below.
C
      MAXRAD = MAXRAD * RADSCL

C
C     If the ray's vertex is distant from the target, use a vertex on
C     the surface of the outer bounding sphere of the target.
C
C     Note that distant vertices can give rise to large transverse
C     displacements of the ray's intercepts on the segments'
C     boundaries. In cases where the ray intercept is very close to the
C     segment's spatial coverage boundaries, this can cause the ray to
C     miss all plates in type 2 segments, unless a large plate
C     expansion factor is used. Using an intercept on the outer
C     bounding sphere greatly ameliorates this problem.
C        
      IF ( VNORM(SEGVTX) .GT. MAXRAD ) THEN
C
C        Find the intercept of the ray with the outer bounding
C        sphere. We'll use this intercept as the vertex for
C        later computations.
C
         CALL SURFPT ( SEGVTX, SEGDIR, MAXRAD, 
     .                 MAXRAD, MAXRAD, SPHVTX, XFND )

         IF ( FAILED() .OR. ( .NOT. XFND ) ) THEN
C
C           It would be highly unusual for the SURFPT call to
C           fail to produce an intercept. Check anyway.
C
            CALL CHKOUT ( 'ZZDSKBUX' )
            RETURN

         END IF

         CALL VEQU ( SPHVTX, SEGVTX )         
            
      END IF

C
C     By default, each segment in the input list must be checked
C     for intersection.
C
C     We start out by trying to eliminate segments from consideration
C     by comparing the ray's distance from their centers to the radii
C     of their bounding spheres. Only those segments whose bounding
C     spheres are hit are examined further.
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
C           In order to do any geometric comparison, the ray must be in
C           the same frame as the segment we're checking. Transform the
C           input vertex and ray if need be.
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
                     CALL CHKOUT ( 'ZZDSKBUX' )
                     RETURN
                  END IF

C
C                 Transform the local copy of the ray to the segment's
C                 frame, and shift the local copy of the ray's vertex
C                 so that it represents an offset relative to center of
C                 the segment's frame.
C
                  CALL MXV  ( XFORM, RAYDIR, SEGDIR )
                  CALL MXV  ( XFORM, VERTEX, SEGVTX )

C
C                 The direction of the buffered offset is from the body
C                 to the segment frame's center. The offset is
C                 expressed in the segment's frame.
C
                  CALL VSUB ( SEGVTX, OFFBUF(1,I), VTEMP  )
                  CALL VEQU ( VTEMP,               SEGVTX )

               END IF


            ELSE IF ( MULTFR ) THEN
C
C              The input and segment frames are the same for this
C              segment, but the current values of SEGVTX and SEGDIR
C              need to be reset.
C
               CALL VEQU ( VERTEX, SEGVTX )
               CALL VEQU ( RAYDIR, SEGDIR )

            END IF

C
C           Find the distance of the ray from the "center" of the
C           segment's coverage volume.
C
            CALL NPLNPT ( SEGVTX, SEGDIR, CTRBUF(1,I), PNEAR, DIST )


            IF ( DIST .LE. RADBUF(I) ) THEN
C
C              The line containing the ray intersects the bounding
C              surface. We'll check the boundary of the segment for an
C              intersection.
C
               CALL ZZRYTELT ( SEGVTX, SEGDIR, DSKBUF(1,I), 
     .                         SGMARG, NXPTS,  XPT         )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZDSKBUX' )
                  RETURN
               END IF

               IF ( NXPTS .GT. 0 ) THEN
C
C                 The ray hits the boundary of this segment. Save the
C                 index of the segment in the "hit list." Record the
C                 distance from the ray's vertex to the intercept in
C                 the parallel array SGDIST.

                  IF ( NHIT .EQ. MAXHIT ) THEN

                     CALL SETMSG ( 'Too many segments were hit by the '
     .               //            'input ray. Buffer size is #.'      )
                     CALL ERRINT ( '#', MAXHIT                         )
                     CALL SIGERR ( 'SPICE(BUFFERTOOSMALL)'             )
                     CALL CHKOUT ( 'ZZDSKBUX'                          )
                     RETURN

                  END IF

                  NHIT           = NHIT + 1
                  SGHIT ( NHIT ) = I
                  SGDIST( NHIT ) = VDIST( SEGVTX, XPT )    
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
C     Leave now if no segments were hit.
C
      IF ( NHIT .EQ. 0 ) THEN
         CALL CHKOUT ( 'ZZDSKBUX' )
         RETURN
      END IF
 
C
C     Find the order of the segments on the hit list, where
C     the metric is the distance of the ray intercepts from
C     the ray's vertex.
C
      CALL ORDERD ( SGDIST, NHIT, IORDER )
      
C
C     Now process the segments on the hit list in order. If we find a
C     surface intercept (that is, a ray intercept with the surface
C     represented by the segment's data, as opposed to the segment's
C     boundary), compare its distance from the ray's vertex to the
C     vertex-boundary distance of the next segment. If the
C     vertex-surface distance is smaller, we terminate the search,
C     since no other segment can contribute a closer intercept.
C
      I      = 1
      DONE   = .FALSE.
      WINNER = 0
      PRVFRM = 0

      DO WHILE ( .NOT. DONE )
C
C        J is the index in the hit list of the segment
C        we're considering. K is the index of that segment
C        in the parallel input arrays.
C
         J = IORDER(I)
         K = SGHIT (J)


         SEGFID = NINT( DSKBUF(FRMIDX,K) )

         IF ( SEGFID .NE. FIXFID ) THEN

            IF ( SEGFID .NE. PRVFRM ) THEN
C
C              Transform and shift the input ray.
C
C              Here J is an index in the hit list and K is 
C              an index in the input arrays.
C  
               CALL MOVED ( SGXBUF(1,1,J), 9, XFORM )

               CALL MXV  ( XFORM,  VERTEX, SEGVTX )
               CALL MXV  ( XFORM,  RAYDIR, SEGDIR )
  
               CALL VSUB ( SEGVTX, OFFBUF(1,K), VTEMP  )
               CALL VEQU ( VTEMP,               SEGVTX )

            END IF

         ELSE IF ( MULTFR ) THEN

            CALL VEQU ( VERTEX, SEGVTX )
            CALL VEQU ( RAYDIR, SEGDIR )

         END IF

C
C        Find the surface intercept using the segment topography
C        data.
C        
         DTYPE = NINT( DSKBUF(TYPIDX,K) )

         CALL ZZDSKSGX ( HANBUF(K), DLABUF(1,K), DTYPE, ET,
     .                   SEGVTX,    SEGDIR,      SRFX,  DC,
     .                   IC,        XFND                   )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZDSKBUX' )
            RETURN
         END IF


         IF ( XFND ) THEN
C
C           At least one surface intercept exists.
C
            IF ( .NOT. FOUND ) THEN
C
C              The intercept we just found is the first one, and
C              at least for now, it is the winner.
C
C              Save the intercept and vertex-intercept distance
C              for this segment.
C
               FOUND  = .TRUE.
               DMIN   = VDIST( SRFX, SEGVTX )

               CALL VEQU ( SRFX, XPT )

               WINNER = J

            ELSE
C
C              At least one surface intercept was found already.
C
C              Compare the vertex-intercept distance for this segment
C              to the best found so far.
C
               DIST = VDIST( SRFX, SEGVTX )

               IF ( DIST .LT. DMIN ) THEN
C
C                 This intercept is closer to the ray's vertex than
C                 any we've seen yet. We have a new winner.
C
                  DMIN   = DIST

                  CALL VEQU ( SRFX, XPT )

                  WINNER = J
                  
               END IF

            END IF

         END IF

C
C        If there's at least one solution in hand, see whether
C        we can stop looking for better solutions.
C
         IF ( FOUND ) THEN

            IF ( I .LT. NHIT ) THEN
C
C              There are more segments in the hit list. Compare the
C              minimum vertex-intercept distance of the segments we've
C              checked to the vertex-boundary distance of the next
C              segment.
C
               J = IORDER(I+1)

               IF ( DMIN .LE. SGDIST(J) ) THEN
C
C                 The best intercept we've found is closer to the
C                 vertex than any intercept that may exist in the
C                 current segment, or any of the remaining segments in
C                 the hit list.
C
                  DONE  = .TRUE.

               END IF

            END IF

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
C     If we have an intercept, it may be represented in a frame
C     other than the input ray frame. Transform it back to the
C     input frame, and shift it as well if the segment frame and
C     input frame have different centers.
C
      IF ( FOUND ) THEN
C
C        K is the index in the input arrays of the "winning" segment.
C        We'll return this index as SEGIDX.
C
         K      = SGHIT ( WINNER )         
         SEGIDX = K

         SEGFID = NINT( DSKBUF(FRMIDX,K) )

         IF ( SEGFID .NE. FIXFID ) THEN
C
C           The segment frame and input frame differ. The intercept
C           must be converted back to the input frame. It also may
C           need to be shifted so that it is expressed as an offset
C           from the target body. If the shift is done, it must be
C           done before the frame transformation, since the offset
C           is expressed relative to the segment's frame.
C
            IF ( .NOT. VZERO( OFFBUF(1,K) ) ) THEN
C
C              OFFBUF(*,K) contains the offset of the segment frame
C              center from the segment's central body.
C
               CALL VADD ( XPT,   OFFBUF(1,K), VTEMP )
               CALL VEQU ( VTEMP,              XPT   )

            END IF

            CALL MOVED( SGXBUF(1,1,WINNER), 9, XFORM )

            CALL MTXV ( XFORM, XPT, VTEMP )
            CALL VEQU ( VTEMP,      XPT   )

         END IF

      END IF
C
C     FOUND and XPT are set. 
C
      CALL CHKOUT ( 'ZZDSKBUX' )
      RETURN
      END
