C$Procedure DSKB02 ( DSK, fetch type 2 bookkeeping data )
 
      SUBROUTINE DSKB02 ( HANDLE, DLADSC, NV,     NP,     NVXTOT,
     .                    VTXBDS, VOXSIZ, VOXORI, VGREXT, 
     .                    CGSCAL, VTXNPL, VOXNPT, VOXNPL          )
 
C$ Abstract
C
C     Return bookkeeping data from a DSK type 2 segment.
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
C     DAS
C     DSK
C
C$ Keywords
C
C     DAS
C     DSK
C     FILES
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsk02.inc'
      
      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      INTEGER               NV
      INTEGER               NP
      INTEGER               NVXTOT
      DOUBLE PRECISION      VTXBDS ( 2, 3 )
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VOXORI ( 3 )
      INTEGER               VGREXT ( 3 )
      INTEGER               CGSCAL      
      INTEGER               VTXNPL
      INTEGER               VOXNPT
      INTEGER               VOXNPL
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DSK file handle.
C     DLADSC     I   DLA descriptor.
C     NV         O   Number of vertices in model.
C     NP         O   Number of plates in model.
C     NVXTOT     O   Number of voxels in fine grid.
C     VTXBDS     O   Vertex bounds.
C     VOXSIZ     O   Fine voxel edge length.
C     VOXORI     O   Fine voxel grid origin.
C     VGREXT     O   Fine voxel grid extent.
C     CGSCAL     O   Coarse voxel grid scale.
C     VTXNPL     O   Size of vertex-plate correspondence list.
C     VOXNPT     O   Size of voxel-plate pointer list.
C     VOXNPL     O   Size of voxel-plate correspondence list.
C
C$ Detailed_Input
C
C     HANDLE         is the handle of a DSK file containing a type 2
C                    segment from which data are to be fetched.
C
C     DLADSC         is the DLA descriptor associated with the segment
C                    from which data are to be fetched. 
C
C$ Detailed_Output
C
C     NV         Number of vertices in model.
C
C     NP         Number of plates in model.
C
C     NVXTOT     Total number of voxels in fine grid.
C
C     VTXBDS     Vertex bounds. This is an array of six values giving
C                the minimum and maximum values of each component of
C                the vertex set.  VTXBDS has dimensions ( 2, 3 ).
C                Units are km.
C
C     VOXSIZ     Fine grid voxel size.  DSK voxels are cubes; the edge
C                length of each cube is given by the voxel size.  This
C                size applies to the fine voxel grid.  Units are km.
C
C     VOXORI     Voxel grid origin.  This is the location of the voxel
C                grid origin in the body-fixed frame associated with
C                the target body.  Units are km.
C
C     VGREXT     Voxel grid extent.  This extent is
C                an array of three integers
C                indicating the number of voxels in
C                the X, Y, and Z directions in the
C                fine voxel grid.
C
C     CGSCAL     Coarse voxel grid scale.  The extent
C                of the fine voxel grid is related to
C                the extent of the coarse voxel grid
C                by this scale factor.
C
C     VTXNPL     Vertex-plate correspondence list
C                size.
C
C     VOXNPT     Size of the voxel-to-plate pointer
C                list.
C
C     VOXNPL     Voxel-plate correspondence list size.
C
C
C$ Parameters
C
C     See the INCLUDE files 
C
C         dla.inc
C         dsk02.inc
C         dskdsc.inc
C
C$ Exceptions
C
C     1) If the input handle is invalid, the error will be diagnosed by
C        routines in the call tree of this routine. 
C
C     2) If a file read error occurs, the error will be diagnosed by
C        routines in the call tree of this routine.
C
C     3) If the input DLA descriptor is invalid, the effect of this
C        routine is undefined. The error *may* be diagnosed by routines
C        in the call tree of this routine, but there are no
C        guarantees.
C
C$ Files
C
C     See input argument HANDLE.
C
C$ Particulars
C
C     DSK files are built using the DLA low-level format and
C     the DAS architecture; DLA files are a specialized type of DAS
C     file in which data are organized as a doubly linked list of
C     segments.  Each segment's data belong to contiguous components of
C     character, double precision, and integer type.
C
C$ Examples
C
C     1) Dump several parameters from the first DLA segment of 
C        a DSK file.  The segment is assumed to be of type 2.
C  
C                 PROGRAM EX1 
C                 IMPLICIT NONE
C
C                 INCLUDE 'dla.inc'
C                 INCLUDE 'dskdsc.inc'
C                 INCLUDE 'dsk02.inc'
C
C           C
C           C     Local parameters
C           C 
C                 INTEGER               FILSIZ
C                 PARAMETER           ( FILSIZ = 255 )
C
C                 INTEGER               LNSIZE
C                 PARAMETER           ( LNSIZE = 80 )
C
C           C
C           C     Local variables
C           C 
C                 CHARACTER*(FILSIZ)    DSK
C                 CHARACTER*(LNSIZE)    OUTLIN
C
C                 DOUBLE PRECISION      VOXORI ( 3 )
C                 DOUBLE PRECISION      VOXSIZ
C                 DOUBLE PRECISION      VTXBDS ( 2, 3 )
C
C                 INTEGER               CGSCAL
C                 INTEGER               DLADSC ( DLADSZ )
C                 INTEGER               HANDLE
C                 INTEGER               NP
C                 INTEGER               NV
C                 INTEGER               NVXTOT
C                 INTEGER               VGREXT ( 3 )
C                 INTEGER               VOXNPL
C                 INTEGER               VOXNPT
C                 INTEGER               VTXNPL
C                 
C                 LOGICAL               FOUND
C
C
C           C
C           C     Prompt for the name of the DSK to read.
C           C
C                 CALL PROMPT ( 'Enter DSK name > ', DSK )
C           C
C           C     Open the DSK file for read access.
C           C     We use the DAS-level interface for
C           C     this function.
C           C
C                 CALL DASOPR ( DSK, HANDLE )
C
C           C
C           C     Begin a forward search through the
C           C     kernel, treating the file as a DLA.
C           C     In this example, it's a very short
C           C     search.
C           C
C                 CALL DLABFS ( HANDLE, DLADSC, FOUND )
C
C                 IF ( .NOT. FOUND ) THEN
C           C
C           C        We arrive here only if the kernel
C           C        contains no segments.  This is 
C           C        unexpected, but we're prepared for it.
C           C
C                    CALL SETMSG ( 'No segments found '
C                .   //            'in DSK file #.'    )
C                    CALL ERRCH  ( '#',  DSK           )
C                    CALL SIGERR ( 'SPICE(NODATA)'     )
C
C                 END IF
C
C           C
C           C     If we made it this far, DLADSC is the
C           C     DLA descriptor of the first segment.
C           C
C           C     Read and display type 2 bookkeeping data.
C           C
C                 CALL DSKB02 ( HANDLE, DLADSC, NV,     NP,     NVXTOT,
C                .              VTXBDS, VOXSIZ, VOXORI, VGREXT, CGSCAL,
C                .              VTXNPL, VOXNPT, VOXNPL                 )
C
C           C
C           C     Show vertex and plate counts.
C           C
C                 OUTLIN = 'Number of vertices:                 #'
C                 CALL REPMI  ( OUTLIN, '#', NV, OUTLIN )
C                 CALL TOSTDO ( OUTLIN )
C
C                 OUTLIN = 'Number of plates:                   #'
C                 CALL REPMI  ( OUTLIN, '#', NP, OUTLIN )
C                 CALL TOSTDO ( OUTLIN )
C
C                 OUTLIN = 'Voxel edge length (km):             #'
C                 CALL REPMF  ( OUTLIN, '#', VOXSIZ, 6, 'E', OUTLIN )
C                 CALL TOSTDO ( OUTLIN )
C
C                 OUTLIN = 'Number of voxels:                   #'
C                 CALL REPMI  ( OUTLIN, '#', NVXTOT, OUTLIN )
C                 CALL TOSTDO ( OUTLIN )
C
C                 END 
C
C
C$ Restrictions
C
C     1) The caller must verify that the segment associated with
C        the input DLA descriptor is a DSK type 2 segment.
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
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C        Updated version info.
C
C        23-JAN-2016 (NJB)
C  
C           Removed references to unneeded variables.
C           Updated header comments.
C
C        DSKLIB Version 2.0.0, 05-MAY-2010 (NJB)
C
C           Renamed routine from DSKP02 to DSKB02.
C     
C        DSKLIB Version 1.0.1, 08-OCT-2009 (NJB)
C
C           Updated header.
C
C        Beta Version 1.0.0, 30-OCT-2006 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch parameters from a type 2 dsk segment
C
C-&
 

C
C     SPICELIB functions
C     
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               DBFSIZ
      PARAMETER           ( DBFSIZ = 10 )

      INTEGER               IBFSIZ
      PARAMETER           ( IBFSIZ = 10 )

      INTEGER               VTBSIZ
      PARAMETER           ( VTBSIZ = 6 )

C
C     Local variables
C
      DOUBLE PRECISION      DBUFF  ( DBFSIZ )

      INTEGER               B
      INTEGER               E
      INTEGER               DPBASE
      INTEGER               IBASE
      INTEGER               IBUFF  ( IBFSIZ )


      IF ( RETURN() ) THEN         
         RETURN
      END IF

      CALL CHKIN ( 'DSKB02' )

      DPBASE = DLADSC ( DBSIDX )
      IBASE  = DLADSC ( IBSIDX )
      
C
C     Read the integer parameters first.  These are located at the
C     beginning of the integer component of the segment.
C     
      CALL DASRDI ( HANDLE, IBASE+1, IBASE+IBFSIZ, IBUFF )

      NV     = IBUFF ( IXNV   ) 
      NP     = IBUFF ( IXNP   ) 
      NVXTOT = IBUFF ( IXNVXT ) 
      CGSCAL = IBUFF ( IXCGSC )
      VOXNPT = IBUFF ( IXVXPS )
      VOXNPL = IBUFF ( IXVXLS )
      VTXNPL = IBUFF ( IXVTLS )

      CALL MOVEI ( IBUFF(IXVGRX), 3, VGREXT )

C
C     Read the d.p. parameters.
C
      B = DPBASE + IXVTBD
      E = DPBASE + IXVTBD + DBFSIZ - 1

      CALL DASRDD( HANDLE, B, E, DBUFF )

      CALL MOVED ( DBUFF,   VTBSIZ,  VTXBDS )
      CALL VEQU  ( DBUFF(VTBSIZ+1),  VOXORI )

      VOXSIZ = DBUFF ( IXVXSZ - DSKDSZ )

      CALL CHKOUT ( 'DSKB02' )
      RETURN
      END 

