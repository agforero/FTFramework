C$Procedure    ZZVOXSCL ( MKDSK, compute voxel scales )
 
      SUBROUTINE ZZVOXSCL ( EXTENT, AVPLEX, TRGCOR,
     .                      TRGFIN, CORSCL, FINSCL )
  
C$ Abstract
C
C     Compute DSK type 2 coarse and fine voxel scales from plate set
C     extents and target voxel counts.
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
C     FILES
C
C$ Declarations

      IMPLICIT NONE
      INCLUDE 'mkdsk02.inc'
 
      DOUBLE PRECISION      EXTENT ( 2, 3 )
      DOUBLE PRECISION      AVPLEX
      INTEGER               TRGCOR
      INTEGER               TRGFIN
      INTEGER               CORSCL
      DOUBLE PRECISION      FINSCL
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     EXTENT     I   Extents of the vertex set.
C     AVPLEX     I   Average plate extent.
C     TRGCOR     I   Target coarse voxel count.
C     TRGFIN     I   Target fine voxel count.
C     CORSCL     O   Coarse voxel scale.
C     FINSCL     O   Fine voxel scale.
C
C$ Detailed_Input
C
C     EXTENT         is the extent of the input vertex set. Units
C                    are km. See Particulars.
C
C     AVPLEX         is the average extent of the plates in the 
C                    plate set. Units are km. See Particulars. 
C
C     TRGCOR         is the target coarse voxel count. This is an
C                    approximation to the final count; the actual count
C                    will be smaller than, but near, TRGCOR.
C
C     TRGFIN         is the target fine voxel count. This is an
C                    approximation to the final count; the actual count
C                    will be smaller than, but near, TRGFIN.
C
C$ Detailed_Output
C
C     CORSCL         is the coarse voxel scale compatible with the
C                    coarse grid determined by this routine. The grid
C                    is large enough to enclose the input vertex set.
C                    Its sides are aligned with the coordinate axes.
C                    CORSCL is the ratio of the edge length of coarse
C                    voxels to that of fine voxels. CORSCL is an
C                    integer.
C
C                    The edge length of coarse voxels is
C
C                       CORSCL * FINSCL * AVPLEX
C
C
C     FINSCL         is the fine voxel scale compatible with the
C                    coarse grid determined by this routine. 
C
C                    The edge length of fine voxels is 
C
C                       FINSCL * AVPLEX.
C
C
C$ Parameters
C
C     See mkdsk02.inc
C
C$ Exceptions
C
C     This routine is meant to be operated in RETURN SPICE error
C     handling mode. The caller is expected to delete the DSK file if
C     an error occurs during file creation.
C
C
C     1)  If the target coarse voxel count is less than 1 or greater
C         than MAXCGR, the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     2)  If the target fine voxel count is less than TRGCOR or greater
C         than MAXVOX, the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     3)  If the average plate extent is not strictly positive, the
C         error SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The "extent" of a set of vertices is the set of extrema
C     of the coordinates of those vertices.
C
C     The extent of a plate along a coordinate axis is the length of
C     the orthogonal projection of the plate onto that axis.
C
C     The average "extent" of a set of plates is the average of the 
C     extents of the plates along the coordinate axes.
C
C     The coarse voxel grid computed by this routine is large enough to
C     enclose the input vertex set. The grid has sides aligned with the
C     coordinate axes.
C
C     The voxel scales computed by this routine are such that, when the
C     coarse voxel grid is extended by 3 coarse voxels in each
C     dimension, the total number of coarse voxels does not exceed
C     MAXCGR, and the total number of fine voxels does not exceed
C     MAXVOX.
C
C$ Examples
C
C     See usage in MKDSK.
C
C$ Restrictions
C
C     This routine should be called only from within MKDSK.
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
C-    MKDSK Version 1.0.0, 24-FEB-2017 (NJB)
C
C-&
 
C$ Index_Entries
C
C     determine type 2 dsk voxel scales
C
C-&


C
C     SPICELIB functions
C     
      LOGICAL               RETURN

C
C     Local parameters
C

C
C     Local variables
C
      DOUBLE PRECISION      CORSZ
      DOUBLE PRECISION      CORSZ0
      DOUBLE PRECISION      DPCS
      DOUBLE PRECISION      FINSZ
      DOUBLE PRECISION      LX
      DOUBLE PRECISION      LY
      DOUBLE PRECISION      LZ
      DOUBLE PRECISION      Q
      DOUBLE PRECISION      SX
      DOUBLE PRECISION      SY
      DOUBLE PRECISION      SZ
      DOUBLE PRECISION      VOLUME
      
      INTEGER               MAXNF
      INTEGER               NCORSE
      INTEGER               NX
      INTEGER               NY
      INTEGER               NZ


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZVOXSCL' )

C
C     Our target counts must be positive. 
C
      IF (  ( TRGCOR .LT. 1 ) .OR. ( TRGCOR .GT. MAXCGR )  ) THEN

         CALL SETMSG ( 'Target coarse voxel count = #; count must ' 
     .   //            'be in the range 1:#.'                      )
         CALL ERRINT ( '#',  TRGCOR                                )
         CALL ERRINT ( '#',  MAXCGR                                )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                    )
         CALL CHKOUT ( 'ZZVOXSCL'                                  )
         RETURN
         
      END IF

      IF ( TRGFIN .LT. TRGCOR ) THEN

         CALL SETMSG ( 'Target fine voxel count = #; target coarse '
     .   //            'count is #; target fine count must be greater ' 
     .   //            'than or equal to target coarse count.'        )
         CALL ERRINT ( '#',  TRGFIN                                   )
         CALL ERRINT ( '#',  TRGCOR                                   )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                       )
         CALL CHKOUT ( 'ZZVOXSCL'                                     )
         RETURN
         
      END IF

C
C     Check the average plate extent.
C
      IF ( AVPLEX .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Average plate extent = #; extent must ' 
     .   //            'be strictly positive.'                 )
         CALL ERRDP  ( '#',  AVPLEX                            )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                )
         CALL CHKOUT ( 'ZZVOXSCL'                              )
         RETURN
         
      END IF
 
C
C     Derive the coarse voxel size first.
C
C
C     Compute the lengths of the edges of the bounding box of the
C     vertex set. Compute the volume of the box.
C
      LX     = EXTENT(2,1) - EXTENT(1,1)
      LY     = EXTENT(2,2) - EXTENT(1,2)
      LZ     = EXTENT(2,3) - EXTENT(1,3)

C
C     In case we have zero extent in any dimension, assign an arbitrary
C     positive extent.
C
      IF ( LX .EQ. 0.D0 ) THEN
         LX = 1.D0
      END IF

      IF ( LY .EQ. 0.D0 ) THEN
         LY = 1.D0
      END IF

      IF ( LZ .EQ. 0.D0 ) THEN
         LZ = 1.D0
      END IF


      VOLUME = LX * LY * LZ

C
C     Compute the nominal coarse voxel size, given the target
C     coarse voxel count and the box volume. We know
C
C                                3
C        volume = trgcor * corsz0
C
C
      CORSZ0 = ( VOLUME / TRGCOR ) ** ( 1.D0/3.D0 )

C
C     Compute the lengths of the bounding box edges as multiples of
C     CORSZ0.
C
      SX = LX / CORSZ0
      SY = LY / CORSZ0
      SZ = LZ / CORSZ0

C
C     Since 
C
C        SX * SY * SZ = VOLUME / CORSZ0**3 = TRGCOR
C
C     we know that
C
C        INT(SX) * INT(SY) * INT(SZ) <= TRGCOR
C
C     The integers above could serve as valid coarse voxel counts
C     in each dimension. However, we want to allow for as much
C     as 3 coarse voxels of padding in each dimension, so we'll
C     reduce the counts by 3 and compute the coarse voxel edge
C     length based on these reduced counts.
C
C     The portion of MKDSK code where the coarse voxel grid extension
C     takes place is inside the routine ZZMKSPIN. The code fragment
C     that carries out the extension is shown below:
C
C       
C        C
C        C     Extend the coarse voxel grid by at least 1/2
C        C     coarse voxel length along each degree of freedom.
C        C
C              XVMIN = DNINT ( XVMIN - 1.D0 )
C              YVMIN = DNINT ( YVMIN - 1.D0 )
C              ZVMIN = DNINT ( ZVMIN - 1.D0 )
C              XVMAX = DNINT ( XVMAX + 1.D0 )
C              YVMAX = DNINT ( YVMAX + 1.D0 )
C              ZVMAX = DNINT ( ZVMAX + 1.D0 )
C
C     
      NX = MAX(  1,  INT(SX)-3  )
      NY = MAX(  1,  INT(SY)-3  )
      NZ = MAX(  1,  INT(SZ)-3  )
      
C
C     Now we know the coarse voxel count, excluding padding voxels.
C     
      NCORSE = NX * NY * NZ

C
C     Compute the final coarse voxel size based on the bounding
C     box edge lengths and the coarse voxel counts along each 
C     edge.
C
      CORSZ = MAX(  ( LX / NX ), 
     .              ( LY / NY ), 
     .              ( LZ / NZ )  )

C
C     By construction, the coarse voxel grid created by ZZMKSPIN can
C     contain 3 extra coarse voxels in each dimension without exceeding
C     the target coarse voxel count.
C
C     Now compute the fine voxel count. We start with the target fine
C     voxel count and the computed coarse voxel count. From these we
C     can compute an estimated fractional coarse voxel scale.
C
      DPCS = ( DBLE(TRGFIN) / NCORSE )**(1.D0/3)

C
C     The coarse scale is an integer. Round down the fractional scale;
C     this effectively enlarges the fine voxels.
C
C     Never set the coarse scale to a value below 1.
C
      CORSCL = MAX( INT( DPCS ),  1 )

C
C     Compute the maximum number of fine voxels that can occupy
C     the expanded fine grid.
C
      MAXNF  = (NX + 3) * (NY + 3) * (NZ + 3) * (CORSCL**3)


      IF ( MAXNF .GT. MAXVOX ) THEN
C
C        The fine voxel count is too large.
C
C        We don't expect this to happen frequently, but if it does,
C        we'll reduce the coarse voxel scale, which effectively 
C        increases the size of the fine voxels and reduces their 
C        count.
C
         Q      = MAXVOX / DBLE( (NX + 3) * (NY + 3) * (NZ + 3) )

         CORSCL = MIN( CORSCL, INT( Q**(1.D0/3.D0) ) )

      END IF

C
C     The fine voxel edge length is the coarse voxel edge length,
C     scaled down by a factor of CORSCL.
C
      FINSZ  = CORSZ / CORSCL

C
C     MKDSK requires the fine voxel scale, not the edge length.
C
C     The average plate extent times the fine voxel scale is the 
C     fine voxel edge length, so we can derive the scale from 
C     the edge length and extent.
C
      FINSCL = FINSZ / AVPLEX


      CALL CHKOUT ( 'ZZVOXSCL' )
      RETURN
      END
