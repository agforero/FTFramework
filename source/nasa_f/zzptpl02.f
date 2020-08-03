C$Procedure ZZPTPL02 ( DSK, map point to plate, type 2 )
 
      SUBROUTINE ZZPTPL02 ( HANDLE, DLADSC, DSKDSC, POINT, 
     .                      PLID,   PLATE,  VERTS,  FOUND )
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Map a point to the nearest plate in a specified type 2 DSK
C     segment.
C
C     The point is expressed in the reference frame of the segment and
C     represents an offset from the center of the segment's reference
C     frame.
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
      INCLUDE 'dsk02.inc'
      INCLUDE 'dsktol.inc'

      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      DOUBLE PRECISION      DSKDSC ( * )
      DOUBLE PRECISION      POINT  ( 3 )
      INTEGER               PLID
      INTEGER               PLATE  ( 3 )
      DOUBLE PRECISION      VERTS  ( 3, 3 )
      LOGICAL               FOUND

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DSK file handle.   
C     DLADSC     I   DLA descriptor of segment.
C     DSKDSC     I   DSK descriptor of segment.
C     POINT      I   Input point.
C     PLID       O   Plate ID.
C     PLATE      O   Plate, expressed as an array of three vertex IDs.
C     VERTS      O   Vertices of plate.
C     FOUND      O   Found flag.
C     PTMEMM     P   Default margin for point-plate distance.
C     XFRACT     P   Default plate expansion fraction.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of a DSK file containing a type 2 DSK
C                segment. The input point is to be matched with the
C                closest plate in this segment.
C                
C     DLADSC     is the DLA descriptor of the DSK segment to be used.
C
C     DSKDSC     is the DSK descriptor of the DSK segment to be used.
C
C     POINT      is a point that is on or very near the surface
C                described by the DSK segment. POINT is expressed in
C                the reference frame of the segment and represents an
C                offset from the center of that frame. The frame center
C                may be distinct from the central body of the segment.
C
C$ Detailed_Output
C
C     PLID       is the ID of the plate nearest to the input point, if
C                such a plate is found. A plate ID is the index of a 
C                plate in the set of plates contained in the segment. 
C                PLID is 1-based.
C
C                The distances between POINT and the respective plates
C                are computed using the plate expansion factor XFRACT.
C                See the Parameters section below and the routine
C                PLTEXP for details.
C
C                A plate will be found only if POINT is sufficiently
C                close to at least one plate in the segment, when plate
C                expansion is taken into account. A margin is used to
C                determine whether a plate is close enough to the point
C                to be considered as a solution. See the Parameters
C                section for details.
C
C
C     PLATE      is the plate designated by PLID. It is an array of
C                three integers. The integers are the IDs of the
C                plate's vertices.
C
C     VERTS      is an array containing the plate's vertices. These
C                a 3-dimensional, double precision vectors. The Ith
C                vertex is stored in elements
C
C                   VERTS(1:3,I)
C                
C     FOUND      is a logical flag that is set to .TRUE. if and only
C                if plate nearest to the input point was found
C
C                The outputs PLID, PLATE, and VERTS are valid if and
C                only if FOUND is .TRUE.
C
C$ Parameters
C
C     See the include file 
C
C        dsktol.inc
C
C     for declarations of these parameters. These defaults can be
C     overridden. See DSKSTL for details.
C     
C
C     PTMEMM     is the default value of the point-plate membership
C                parameter, which is used to determine whether a point
C                is close enough to a plate to be considered as a
C                solution. Let S be the current value of this
C                parameter; the bounding radius of the segment is
C                multiplied by (1+S) to produce a distance used for
C                this comparison.
C
C     XFRACT     is the default value of an expansion factor applied to
C                each plate before the distance of POINT from the plate
C                is computed. Let S be the current value of this
C                parameter; each plate is expanded by a factor of (1+S)
C                to produce a plate used for the membership test.
C
C                This expansion is performed to keep results of this
C                routine consistent with those of DSKX02, which also
C                performs plate expansion.
C
C$ Exceptions
C
C     1)  If an invalid voxel edge length is detected, the error  
C         SPICE(VALUEOUTOFRANGE) is signaled.
C          
C     2)  If the coarse voxel scale is zero, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C     3)  If an unrecognized coordinate system is encountered, the
C         error SPICE(NOTSUPPORTED) is signaled.
C         
C     4)  If an error occurs while this routine attempts to read
C         DSK data, the error will be signaled by a routine in
C         the call tree of this routine.
C
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - DSK data: the DSK file designated by HANDLE and containing
C          the segment having the DLA descriptor DLADSC must be loaded
C          at the time this routine is called.
C
C     Kernel data are normally loaded once per program run, NOT every
C     time this routine is called.
C     
C$ Particulars
C
C     This routine supports geometry routines that must associate
C     a computed surface point with the plate to which that point 
C     belongs.
C
C$ Examples
C
C     See usage in ZZDSKNRM.
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
C        23-AUG-2016 (NJB)         
C
C        Bug fix: now saves previous handle and DLA descriptor.
C
C        17-JUN-2016 (NJB) 
C
C           Added support for planetodetic coordinates. Updated
C           bounding radius computation for segments using latitudinal
C           coordinates; this routine now calls ZZSEGBOX. Changed
C           implementation to use ZZINVELT and ZZVOXCVO. Updated
C           parameter descriptions.
C
C        07-OCT-2015 (NJB)
C
C           Now uses plate expansion before testing point-plate 
C           distance.
C
C        15-DEC-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     map point to plate in type 2 dsk segment
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      VDOT

      INTEGER               BRCKTI
      INTEGER               ZZVOX2ID

      LOGICAL               DLASSG
      LOGICAL               FAILED
      LOGICAL               RETURN
      
C
C     Local parameters
C
      INTEGER               BUFSIZ
      PARAMETER           ( BUFSIZ = 1000 )

C
C     Local variables
C     
      DOUBLE PRECISION      BOXCTR ( 3 )
      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      DMIN
      DOUBLE PRECISION      LIMIT
      DOUBLE PRECISION      OFFSET ( 3 )
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      NORMAL ( 3 )
      DOUBLE PRECISION      PNEAR  ( 3 )
      DOUBLE PRECISION      PNTOFF ( 3 )
      DOUBLE PRECISION      PTSRFM
      DOUBLE PRECISION      VOXORI ( 3 )
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VRTTMP ( 3, 3 )
      DOUBLE PRECISION      XPDFRC
      DOUBLE PRECISION      XVERTS ( 3, 3 )

      INTEGER               CGRCOR ( 3 )
      INTEGER               CGREXT ( 3 )
      INTEGER               CGROFF ( 3 )
      INTEGER               CGRPTR
      INTEGER               CGRSCL
      INTEGER               CORSYS
      INTEGER               CGRVID
      INTEGER               I
      INTEGER               J
      INTEGER               K
      INTEGER               N
      INTEGER               NPLATE
      INTEGER               NREAD
      INTEGER               PIDTMP
      INTEGER               PLTBUF ( BUFSIZ )
      INTEGER               PLTPTR
      INTEGER               PLTTMP ( 3 )
      INTEGER               PTRLOC
      INTEGER               PRVDSC ( DLADSZ )
      INTEGER               PRVHAN
      INTEGER               PTROFF
      INTEGER               REMAIN
      INTEGER               START
      INTEGER               VGRCOR ( 3 )
      INTEGER               VGREXT ( 3 )
      INTEGER               VID

      LOGICAL               INSIDE
      LOGICAL               PASS1

C
C     Saved variables
C
      SAVE                  CGRSCL
      SAVE                  CORSYS
      SAVE                  LIMIT
      SAVE                  MAXR
      SAVE                  PASS1
      SAVE                  PRVHAN
      SAVE                  PRVDSC
      SAVE                  VGREXT
      SAVE                  VOXORI
      SAVE                  VOXSIZ

C
C     Initial values
C
      DATA                  LIMIT  / -1.D0      /
      DATA                  PASS1  / .TRUE.     /
      DATA                  PRVDSC / DLADSZ * 0 /
      DATA                  PRVHAN / 0          /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZPTPL02' )

C
C     No plate has been found so far.
C
      FOUND = .FALSE.

C
C     Decide whether we're looking at the segment we saw
C     on the previous call.
C
      IF (  PASS1  .OR. 
     .      .NOT.  DLASSG( HANDLE, PRVHAN, DLADSC, PRVDSC )  ) THEN
C
C        We'll need to look up the voxel grid parameters for this
C        segment.
C
         CALL DSKD02 ( HANDLE, DLADSC, KWVXOR, 1, 3, N, VOXORI )
         CALL DSKD02 ( HANDLE, DLADSC, KWVXSZ, 1, 1, N, VOXSIZ )
         CALL DSKI02 ( HANDLE, DLADSC, KWVGRX, 1, 3, N, VGREXT )
         CALL DSKI02 ( HANDLE, DLADSC, KWCGSC, 1, 1, N, CGRSCL )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZPTPL02' )
            RETURN
         END IF


         IF ( VOXSIZ .EQ. 0.D0 ) THEN

            CALL SETMSG ( 'Voxel edge length is zero; length '
     .      //            'must be positive.'                  )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'             )
            CALL CHKOUT ( 'ZZPTPL02'                           )
            RETURN

         END IF

         IF ( CGRSCL .EQ. 0 ) THEN

            CALL SETMSG ( 'Coarse voxel scale is zero; scale '
     .      //            'must be positive.'                  )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'             )
            CALL CHKOUT ( 'ZZPTPL02'                           )
            RETURN

         END IF

         CORSYS = NINT( DSKDSC(SYSIDX) )

         CALL ZZSEGBOX ( DSKDSC, BOXCTR, MAXR )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZPTPL02' )
            RETURN
         END IF

C
C        We succesfully obtained the desired segment parameters, so we
C        don't need to execute this code again until the segment
C        changes. Save the current handle and DLA descriptor.
C
         PRVHAN = HANDLE
         CALL MOVEI ( DLADSC, DLADSZ, PRVDSC )

         PASS1  = .FALSE.

      END IF

C
C     Look up the point-plate membership margin; compute 
C     the distance limit. This call must be made on every
C     call to ZZPTPL02.
C
      CALL DSKGTL ( KEYPTM, PTSRFM )

      LIMIT = PTSRFM * MAXR

C
C     Look up the plate expansion fraction. This call must be made on
C     every call to ZZPTPL02.
C
      CALL DSKGTL ( KEYXFR, XPDFRC )

C
C     Find out whether the point is within the volume element
C     bounding the segment.
C
      CALL ZZINVELT ( POINT,          CORSYS, DSKDSC(PARIDX), 
     .                DSKDSC(MN1IDX), PTSRFM, 0,              INSIDE )
      
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
      END IF


      IF ( .NOT. INSIDE ) THEN
C
C        The point is too far from the segment to be considered
C        to lie on a plate in that segment.
C
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
         
      END IF

C
C     Map the point to the coordinates of a voxel containing it. If the
C     point is outside the voxel grid, map the point to the closest
C     voxel.
C
      CALL VSUB ( POINT, VOXORI, OFFSET )

      DO I = 1, 3

         J         = INT( OFFSET(I) / VOXSIZ ) + 1

         VGRCOR(I) = BRCKTI( J,  1,  VGREXT(I) )

      END DO

C
C     Compute the coordinates of the coarse voxel containing the fine
C     voxel we just identified. Get the 1-d offset of the fine voxel
C     relative the coarse voxel; this offset gives us the index of the
C     pointer associating the fine voxel with its plate list. The
C     1-d offset PTROFF is 1-based.
C
      CALL ZZVOXCVO ( VGRCOR, VGREXT, CGRSCL, CGRCOR, CGROFF, PTROFF )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
      END IF

C
C     Fetch the pointer from the coarse voxel to the first element of
C     its fine voxel pointer array.
C
C     We'll need the 1-D offset of the coarse voxel from the base of
C     the coarse voxel grid.
C     
      DO I = 1, 3
         CGREXT(I) = VGREXT(I) / CGRSCL
      END DO
 
      CGRVID = ZZVOX2ID ( CGRCOR, CGREXT )

      CALL DSKI02 ( HANDLE, DLADSC, KWCGPT, CGRVID, 1, N, CGRPTR )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
      END IF

      IF ( CGRPTR .LT. 1 ) THEN
C
C        There are no non-empty fine voxels, hence no plates, in the
C        coarse voxel we're looking at.
C        
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
         
      END IF
C
C     Look up the pointer to the plate list for this voxel, and if
C     the pointer is non-null, look up the plate count.
C
      PTRLOC = (CGRPTR - 1) + PTROFF

      CALL DSKI02 ( HANDLE, DLADSC, KWVXPT, PTRLOC, 1, N, PLTPTR )

      IF (  FAILED()  .OR.  ( PLTPTR .LT. 1 )  ) THEN
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
      END IF

      CALL DSKI02 ( HANDLE, DLADSC, KWVXPL, PLTPTR, 1, N, NPLATE )

      IF (  FAILED()  .OR.  ( NPLATE .LT. 1 )  ) THEN
         CALL CHKOUT ( 'ZZPTPL02' )
         RETURN
      END IF

C
C     Loop through the plates, keeping track of the minimum plate-point
C     distance. 
C
      DMIN   = DPMAX()
      REMAIN = NPLATE
      NREAD  = MIN ( REMAIN, BUFSIZ )
      I      = 1
     
      DO WHILE ( REMAIN .GT. 0 ) 
C
C        Look up the current set of plate IDs.
C
         CALL DSKI02 ( HANDLE,   DLADSC, KWVXPL, 
     .                 PLTPTR+I, NREAD,  N,      PLTBUF )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZPTPL02' )
            RETURN
         END IF
C
C        Look up the vertices of each plate in the buffer and find
C        the distance of the point from that plate. Quit if we
C        find a match.
C
         DO J = 1, NREAD

            PIDTMP = PLTBUF(J)
            START  = ( 3*(PIDTMP-1) ) + 1

            CALL DSKI02 ( HANDLE, DLADSC, KWPLAT, START, 3, N, PLTTMP )
            
            DO K = 1, 3
               
               VID   = PLTTMP(K)
               START = ( 3*(VID-1) ) + 1

               CALL DSKD02 ( HANDLE, DLADSC, KWVERT,
     .                       START,  3,      N,      VRTTMP(1,K) )

            END DO
            
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZPTPL02' )
               RETURN
            END IF

C
C           Work with an expanded version of the plate.
C
            CALL PLTEXP ( VRTTMP, XPDFRC, XVERTS )

            CALL PLTNRM ( XVERTS(1,1), XVERTS(1,2), XVERTS(1,3),
     .                    NORMAL                                )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZPTPL02' )
               RETURN
            END IF

            CALL VHATIP ( NORMAL )

            CALL VSUB ( POINT, XVERTS(1,1), PNTOFF )

            IF ( ABS( VDOT(PNTOFF,NORMAL) ) .LE. LIMIT ) THEN
C
C              The input point lies in a narrow region of space
C              bounded by two planes, both of which are parallel
C              to the plate. The plate lies between the planes.
C 
C              This test does not rule out a comparison between POINT
C              and a distant plate, if POINT is close to the plane
C              containing that plate. However, the proportion of such
C              cases will normally be small.
C
               CALL PLTNP ( POINT,       XVERTS(1,1), XVERTS(1,2), 
     .                      XVERTS(1,3), PNEAR,       DIST        )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZPTPL02' )
                  RETURN
               END IF

            ELSE
               DIST = DPMAX()
            END IF

            IF ( DIST .LE. LIMIT ) THEN
C
C              We have a reasonable candidate for the closest plate. 
C
               FOUND = .TRUE.

               IF ( DIST .LT. DMIN ) THEN

                  DMIN  = DIST
                  PLID  = PIDTMP
C
C                 Set the output vertices to the original version.
C              
                  CALL MOVEI ( PLTTMP, 3, PLATE )
                  CALL MOVED ( VRTTMP, 9, VERTS )
C
C                 We'll return the above values if we don't find
C                 a better match.
C                 
               END IF

            END IF

         END DO

C
C        Prepare to read the next set of plate IDs, if any.
C
         REMAIN = REMAIN - NREAD
         I      = I      + NREAD
         NREAD  = MIN ( REMAIN, BUFSIZ )

      END DO

      CALL CHKOUT ( 'ZZPTPL02' )
      RETURN     
      END
