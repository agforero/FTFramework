C$Procedure  ZZINGRD  ( is a voxel inside the grid? )
 
       LOGICAL FUNCTION ZZINGRD ( NVOX, VOXEL )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C 
C     Return true if voxel is inside the voxel grid. This routine
C     operates on integer voxel coordinates.
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
C     VOXEL
C
C$ Declarations
 
      IMPLICIT              NONE
 
      INTEGER               NVOX  ( 3 )
      INTEGER               VOXEL ( 3 )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NVOX       I   Dimensions of voxel grid.
C     VOXEL      I   Coordinates of voxel in voxel grid units.
C
C     The function returns .TRUE. if the voxel is inside the grid.
C
C$ Detailed_Input
C
C     NVOX       Dimensions of voxel grid.
C
C     VOXEL      Coordinates of voxel in voxel grid units. The
C                coordinates are 1-based integers.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if the voxel is inside the grid.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     Error free.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports the SPICELIB routine XDDA.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 02-FEB-2016 (NJB) (JAB)
C
C        Renamed routine.
C
C        Based on DSKLIB Version 1.0.0, 03-FEB-1999 (JAB)
C
C-&
 
C$ Index_Entries
C
C      voxel inside grid
C
C-&
 
 
      INTEGER               I
 
 
      ZZINGRD = .FALSE.
 
C
C     Determine if voxel is outside the voxel grid
C     in any direction.
C
      DO I = 1, 3

         IF (      ( VOXEL(I) .LT. 1       ) 
     .        .OR. ( VOXEL(I) .GT. NVOX(I) )  ) THEN

            RETURN

         END IF

      END DO
 
      ZZINGRD = .TRUE.
 
      RETURN
      END
