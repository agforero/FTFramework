C$Procedure  ZZVOXCVO ( Voxel to coarse voxel offset )
 
      SUBROUTINE ZZVOXCVO ( VIXYZ, NVOX, CGRSCL, CGXYZ, CGOFF, CGOF1D )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Map voxel coordinates to coarse voxel coordinates and offset
C     relative to coarse voxel. Offset coordinates are 1-based.
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
C     PRIVATE
C     UTILITY
C
C$ Declarations
 
      IMPLICIT NONE

      INTEGER               VIXYZ  ( 3 )
      INTEGER               NVOX   ( 3 )
      INTEGER               CGRSCL
      INTEGER               CGXYZ  ( 3 )
      INTEGER               CGOFF  ( 3 )
      INTEGER               CGOF1D
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     VIXYX      I   Fine voxel integer coordinates.
C     NVOX       I   Fine voxel grid dimensions.
C     CGRSCL     I   Coarse voxel scale.
C     CGXYZ      O   Coarse voxel coordinates.
C     CGOFF      O   3-D voxel offset relative to coarse voxel.
C     CGOF1D     O   1-D voxel offset relative to coarse voxel.
C
C$ Detailed_Input
C
C     VIXYZ      is an array containing integer Cartesian coordinates
C                of a voxel in the fine voxel grid. Coordinates are
C                1-based.
C
C     NVOX       is an array containing the integer extents of the fine
C                voxel grid. Units are voxels.
C
C     CGRSCL     is the coarse voxel scale. This is the integer ratio
C                of the coarse voxel edge length to the fine voxel
C                edge length.
C          
C$ Detailed_Output
C
C     CGXYZ      is an array containing integer Cartesian coordinates
C                of the coarse voxel that contains the fine voxel
C                indexed by VIXYZ. Coordinates are 1-based.
C
C     CGOFF      is an array containing the integer Cartesian
C                coordinates of the fine voxel relative to the coarse
C                voxel that contains it. Coordinates are 1-based.
C
C     CGOF1D     is the 1-dimensional offset of the input fine voxel
C                relative to the coarse voxel that contains it. The
C                offset is 1-based.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If any input fine voxel grid dimension is non-positive,
C        the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     2) If the coarse voxel scale is non-positive, the error
C        SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine encapsulates a frequently used type 2 DSK address 
C     calculation. The output offset CGOF1D is useful for locating
C     a fine voxel's pointer into the voxel-plate list.
C
C$ Examples
C
C     See usage in the routine ZZMKSPIN.
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
C-    SPICELIB Version 1.0.0, 23-AUG-2016 (NJB)
C
C        $Exceptions section was updated.
C
C     Based on DSKLIB Version 1.0.0, 20-MAR-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     map voxel to coarse voxel offset
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C     
      INTEGER               I


      IF ( RETURN() ) THEN
         RETURN
      END IF

      DO I = 1, 3
         
         IF ( NVOX(I) .LT. 1 ) THEN

            CALL CHKIN  ( 'ZZVOXCVO' )

            CALL SETMSG ( 'Voxel grid dimensions must be positive '
     .      //            'but were # # #.'                        )
            CALL ERRINT ( '#', NVOX(1)                             )
            CALL ERRINT ( '#', NVOX(2)                             )
            CALL ERRINT ( '#', NVOX(3)                             )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )

            CALL CHKOUT ( 'ZZVOXCVO' )
            RETURN

         END IF

      END DO


      DO I = 1, 3

         IF ( ( VIXYZ(I) .LT. 1 ) .OR. ( VIXYZ(I) .GT. NVOX(I) ) ) THEN

            CALL CHKIN  ( 'ZZVOXCVO' )

            CALL SETMSG ( 'Voxel grid coordinates must be inside '
     .      //            'grid having dimensions # x # x # '
     .      //            'but were # # #.'                        )
            CALL ERRINT ( '#', NVOX(1)                             )
            CALL ERRINT ( '#', NVOX(2)                             )
            CALL ERRINT ( '#', NVOX(3)                             )
            CALL ERRINT ( '#', VIXYZ(1)                            )
            CALL ERRINT ( '#', VIXYZ(2)                            )
            CALL ERRINT ( '#', VIXYZ(3)                            )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )

            CALL CHKOUT ( 'ZZVOXCVO' )
            RETURN

         END IF

      END DO


      IF ( CGRSCL .LT. 1 ) THEN

         CALL CHKIN  ( 'ZZVOXCVO' )

         CALL SETMSG ( 'Coarse voxel grid scale must be positive '
     .   //            'but was #.'                               )
         CALL ERRINT ( '#', NVOX(1)                               )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                   )

         CALL CHKOUT ( 'ZZVOXCVO' )
         RETURN

      END IF

      DO I = 1, 3
C
C        Set the Ith coarse grid coordinate. Recall these coordinates
C        are 1-based.
C
         CGXYZ(I) = ( ( VIXYZ(I) - 1 ) /  CGRSCL  )  +  1

C
C        Set the Ith coarse grid coordinate offset. These offsets 
C        are 1-based as well.
C
         CGOFF(I) =  VIXYZ(I)  -  CGRSCL*( CGXYZ(I) - 1 )

      END DO

C
C     Convert the coarse grid-relative offset to a relative
C     ID. The ID is a one-dimensional offset.
C
      CGOF1D =    (CGOFF(3)-1) * CGRSCL * CGRSCL
     .          + (CGOFF(2)-1) * CGRSCL
     .          +  CGOFF(1)
      
      END 
