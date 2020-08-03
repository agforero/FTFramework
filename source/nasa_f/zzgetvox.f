C$Procedure ZZGETVOX ( Coordinates of voxel containing a point )
 
      SUBROUTINE ZZGETVOX ( VOXSIZ, VOXORI, NVOX, XYZ, INBOX, VOXCOR )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     This subroutine returns the 1-based voxel coordinates of the
C     voxel that contains a given point.
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
C     None.
C
C$ Keywords
C
C     PLATE
C
C$ Declarations
 
      IMPLICIT NONE
 
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VOXORI    ( 3 )
      INTEGER               NVOX      ( 3 )
      DOUBLE PRECISION      XYZ       ( 3 )
      LOGICAL               INBOX
      INTEGER               VOXCOR     ( 3 )
 
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     VOXSIZ     I   Voxel size in model units.
C     VOXORI     I   3-Vector location of the voxel grid origin
C                    in model units.
C     NVOX       I   Number of voxel in each dimension.
C     XYZ        I   Model grid coordinates of point.
C     INBOX      O   Logical value is true if point is inside voxel
C                    grid, false if not.
C     VOXCOR     O   3-Index of voxel location with the voxel grid.
C
C$ Detailed_Input
C
C     VOXSIZ     is a double precision scalar value giving the size of
C                each cubical voxel in model units.
C
C     VOXORI     is a double precision 3-vector indicating the location
C                of the voxel grid origin in model units.
C
C     NVOX       is an integer 3-vector giving the length of each
C                side of the voxel grid in voxel units.
C
C     XYZ        is a double precision 3-vector containing the location
C                of the point in question in model units.
C
C$ Detailed_Output
C
C     INBOX      is a logical scalar indicating whether or not the
C                point XYZ lies within the extent of the voxel grid.
C
C     VOXCOR     is an integer 3-vector containing the x, y, and z
C                voxel indices of the voxel containing the point.
C
C                VOXCOR is valid only if INBOX is .TRUE.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If VOXSIZ is less than or equal to zero, the error
C         SPICE(NONPOSITIVEVALUE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This is not a general coordinate conversion routine. This routine
C     returns valid voxel coordinates only for voxels inside the voxel
C     grid.
C
C     Points on the outer surface of the grid, as defined by the voxel
C     origin and the grid dimensions, are considered to be inside the
C     grid.
C
C     Points outside of the grid are not mapped to voxel coordinates.
C
C$ Examples
C
C
C   C
C   C     Some position in model coordinates, i.e. the body fixed
C   C     units and frame coordinates.
C   C
C         DOUBLE PRECISION       POS(3)
C
C               ...
C
C   C
C   C     Retrive the voxel grid geometry description.
C   C
C         CALL VOXDIM ( NVOX, VOXSIZ, VOXORI)
C
C   C
C   C     Convert the coordinate POS to voxel grid indices.
C   C
C         CALL ZZGETVOX ( VOXSIZ, VOXORI, NVOX, POS, INBOX, VOXCOR )
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
C     N.J. Bachcman   (JPL)
C     J.A. Bytof      (JPL)
C     E.D. Wright     (JPL
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 18-MAR-2016 (NJB)
C
C        Added error check for non-positive voxel size.
C        Made small changes to header.
C
C        17-JAN-2016 (NJB)
C
C           Functional change: routine returns valid voxel coordinates
C           only for voxels inside the voxel grid. Points on the outer
C           surface of the grid are considered to be inside the grid.
C
C           Changed check-in to discovery style.
C
C           Corrected rounding for voxels on outer edge of grid.
C           Updated header.
C
C        08-OCT-2009 (NJB)
C
C           Re-ordered header sections.
C
C        14-SEP-2004 (EDW)
C
C           Corrected VOXID(I) functions to allow for
C           negative voxel coordinates. Original functions
C           returned values GE 1. 
C
C           Added Examples section.
C
C        28-JAN-1999 (JAB)
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     returns the voxel grid indices containing a given point
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables.
C 
      DOUBLE PRECISION       TERM
      INTEGER                I
 
C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      IF ( VOXSIZ .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZGETVOX'                            )
         CALL SETMSG ( 'Voxel size was #; must be positive.' )
         CALL ERRDP  ( '#',  VOXSIZ                          )
         CALL SIGERR ( 'SPICE(NONPOSITIVEVALUE)'             )
         CALL CHKOUT ( 'ZZGETVOX'                            )
         RETURN

      END IF

C
C     Initialize 'point in box' flag and voxel coordinates. The
C     coordinates are assigned out-of-range values.
C
      INBOX     = .FALSE.
 
      VOXCOR(1) = 0
      VOXCOR(2) = 0
      VOXCOR(3) = 0

C
C     Scale the point's coordinates to voxel grid space
C     and determine the indices of the voxel that contains it.
C
      DO I = 1, 3
C
C        A Galilean transform. Calculate the voxel coordinate
C        corresponding to the body centered coordinate. This
C        operation performs the same task as TOGRID, but
C        including the operation here improves ZZGETVOX's
C        runtime performance.
C
         TERM = ( XYZ(I) - VOXORI(I) ) / VOXSIZ

C
C        Calculate the voxel index for each degree of freedom
C        corresponding to the voxel coordinate.
C
C        If the point is outside of the grid, return now.
C
         IF ( ( TERM .LT. 0.D0 ) .OR. ( TERM .GT. NVOX(I) ) ) THEN
                        
            RETURN
        
         END IF

C
C        Assign a 1-based value to the Ith component of the voxel's
C        coordinates. The outer surface of the grid is considered part
C        of the grid.
C
C        Note that TERM is non-negative at this point.
C
         IF ( INT(TERM) .LT. NVOX(I) ) THEN

            VOXCOR(I) =  1 + INT(TERM)
            
         ELSE
C
C           TERM is NVOX(I), since the cases 
C
C              TERM > NVOX(I)
C              TERM < NVOX(I)
C
C           have been ruled out.
C
            VOXCOR(I) = NVOX(I)

         END IF

      END DO

      INBOX = .TRUE.
   
      RETURN
      END
