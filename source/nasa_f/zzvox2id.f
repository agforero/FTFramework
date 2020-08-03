C$Procedure ZZVOX2ID ( Return voxel index from coords )
 
      INTEGER FUNCTION ZZVOX2ID ( VIXYZ, NVOX )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Given a voxel's grid coordinates and the voxel grid
C     dimensions, return the voxel's 1-dimensional index.
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
C     COORDINATE
C     DSK
C     INDEX
C     VOXEL
C
C$ Declarations
 
      IMPLICIT              NONE
 
      INTEGER               NVOX   ( 3 )
      INTEGER               VIXYZ  ( 3 )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     VIXYZ      I   Voxel's grid coordinates stored in a
C                    three-dimensional vector.
C     NVOX       I   Number of voxels in each coordinate direction.
C
C     The function returns the voxel index.
C
C$ Detailed_Input
C
C     VIXYZ      An integer 3-vector storing the voxel grid coordinates
C                for the point of interest.
C
C     NVOX       An array of three positive integers defining a voxel
C                grid's extents in the X, Y, and Z directions.
C
C$ Detailed_Output
C
C     The function returns the value of the voxel index.
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
C     None.
C
C$ Examples
C
C     See usage in ZZMKSPIN.
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
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C        Updated version info.
C
C        17-JAN-2016 (NJB)
C
C           Re-wrote portions of header.
C
C        08-OCT-2009 (NJB)
C
C           Re-ordered header sections.
C
C        23-SEP-2004 (EDW)
C
C           Recast equation to improve runtime performance.
C           Edited comments for clarity.
C
C           Subroutine renamed from GETXYZ.
C
C        19-FEB-2000 (JAB)
C
C           Removed SAVED variables.
C
C        04-FEB-1999 (JAB)
C
C-&
 
C$ Index_Entries
C
C     given a voxel's coordinates return the voxel index
C
C-&

C
C     Convert from voxel coordinates to voxel index. A more
C     readable form for the following function:
C
C        NX   = NVOX(1)
C        NY   = NVOX(2)
C        NXNY = NX * NY
C
C        ZZVOX2ID = VIXYZ(1) + (VIXYZ(2)-1)*NX + (VIXYZ(3)-1)*NXNY
C
C     Expressing the function in this format improves runtime
C     performance as per Horner's Rule.
C
      ZZVOX2ID = VIXYZ(1) + 
     .           NVOX(1)*(  VIXYZ(2) - 1 + (VIXYZ(3)-1)*NVOX(2)  )

 
      RETURN
      END
