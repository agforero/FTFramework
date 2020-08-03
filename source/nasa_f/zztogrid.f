C$Procedure ZZTOGRID ( Model coordinates to voxel grid coordinates )
 
       SUBROUTINE ZZTOGRID ( MODXYZ, VOXORI, VOXSIZ, GRDXYZ )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Shift and scale plate model coordinates to zero-based, double
C     precision voxel grid coordinates.
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
C     DSK
C     PLATE 
C     TOPOGRAPHY
C     VOXEL
C
C$ Declarations
 
      IMPLICIT NONE
 
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      MODXYZ  ( 3 )
      DOUBLE PRECISION      VOXORI  ( 3 )
      DOUBLE PRECISION      GRDXYZ  ( 3 )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     MODXYZ     I   Coordinates of point in model coordinates.
C     VOXORI     I   Origin of voxel grid in model coordinates.
C     VOXSIZ     I   Voxel size in model coordinates.
C     GRDXYZ     O   Coordinates of point in voxel grid coordinates.
C
C$ Detailed_Input
C
C     MODXYZ     Coordinates of point in model coordinates. The
C                point is expressed as an offset from the center
C                of a body-fixed reference frame associated with
C                a target body. Units are km.
C
C
C     VOXORI     Origin of voxel grid in model coordinates. VOXORI
C                is expressed in the same reference frame as MODXYZ.
C                Units are km.
C
C
C     VOXSIZ     Voxel size in model coordinates. VOXSIZ is the
C                voxel's edge length. Units are km.
C
C
C$ Detailed_Output
C
C     GRDXYZ     Coordinates of the input point, scaled to voxel grid
C                coordinates, and expressed as an offset from the voxel
C                grid origin. These coordinates are zero-based, double
C                precision values. The units are the voxel edge length
C                multiplied by km.
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
C     This routine supports DSKX02.
C
C$ Examples
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
C   C     Retrieve the voxel grid geometry description.
C   C
C         CALL VOXDIM ( NVOX, VOXSIZ, VOXORI)
C
C   C
C   C     Convert the coordinate POS to voxel grid coordinates.
C   C
C         CALL ZZTOGRID ( POS, VOXORI, VOXSIZ, POSVOX )
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
C     N.J. Bachman    (JPL)
C     J.A. Bytof      (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 18-MAY-2016 (NJB) 
C
C        Added error check for non-positive voxel size.
C
C        02-FEB-2016 (NJB) (JAB)
C
C           Renamed routine to ZZTOGRID.
C
C        08-OCT-2009 (NJB)
C
C           Updated header.
C
C        19-OCT-2004 (EDW)
C
C           Added Examples section, edits to comments.
C
C        03-FEB-1999 (JAB)
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     scale plate model coordinates to voxel grid coordinates
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
 
C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      IF ( VOXSIZ .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZTOGRID'                            )
         CALL SETMSG ( 'Voxel size was #; must be positive.' )
         CALL ERRDP  ( '#',  VOXSIZ                          )
         CALL SIGERR ( 'SPICE(NONPOSITIVEVALUE)'             )         
         CALL CHKOUT ( 'ZZTOGRID'                            )
         RETURN

      END IF
      
C
C     Convert model coordinates to voxel grid coordinates
C     via a Galilean transform.
C
      DO I = 1, 3
         GRDXYZ(I) = (MODXYZ(I) - VOXORI(I)) / VOXSIZ
      END DO

      RETURN
      END
