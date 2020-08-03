C$Procedure DSKZ02 ( DSK, fetch type 2 model size parameters )
 
      SUBROUTINE DSKZ02 ( HANDLE, DLADSC, NV, NP )

C$ Abstract
C
C     Return plate model size parameters---plate count and
C     vertex count---for a type 2 DSK segment.
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
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DSK file handle.
C     DLADSC     I   DLA descriptor.
C     NV         O   Number of vertices.
C     NP         O   Number of plates.
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
C     NV             is the number of vertices belonging to
C                    the specified plate model.
C
C     NP             is the number of plates belonging to the
C                    specified plate model.
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
C        in the call tree of this routine, but there are no guarantees.
C
C$ Files
C
C     See input argument HANDLE.
C
C$ Particulars
C
C     This routine enables SPICE-based user applications to
C     conveniently fetch the plate and vertex counts of a type 2 DSK
C     segment.
C
C     See the routine DSKB02 (DSK, fetch type 2 bookkeeping data)
C     for an interface that returns all type 2 DSK segment
C     bookkeeping data in a single call.
C
C$ Examples
C
C     1) Look up all the vertices associated with each plate
C        of the model contained in a specified type 2 segment. For each
C        plate, display the plate's vertices and normal vector.
C
C        For this example, we'll show the context of this look-up:
C        opening the DSK file for read access, traversing a trivial,
C        one-segment list to obtain the segment of interest.
C
C
C        Example code begins here.
C
C
C           PROGRAM DSKBULK
C           IMPLICIT NONE
C
C           INCLUDE 'dla.inc'
C           INCLUDE 'dsk02.inc'
C
C
C           CHARACTER*(*)         FMT
C           PARAMETER           ( FMT    = '(1X,A,3(1XE16.9))' )
C
C
C           INTEGER               BUFSIZ
C           PARAMETER           ( BUFSIZ = 10000 )
C
C           INTEGER               FILSIZ
C           PARAMETER           ( FILSIZ = 255 )
C
C
C           CHARACTER*(FILSIZ)    DSK
C
C           DOUBLE PRECISION      NORMAL ( 3 )
C           DOUBLE PRECISION      VERTS  ( 3, BUFSIZ )
C
C           INTEGER               DLADSC ( DLADSZ )
C           INTEGER               HANDLE
C           INTEGER               I
C           INTEGER               J
C           INTEGER               N
C           INTEGER               NNORM
C           INTEGER               NP
C           INTEGER               NREAD
C           INTEGER               NV
C           INTEGER               NVTX
C           INTEGER               PLATES  ( 3, BUFSIZ )
C           INTEGER               PLIX
C           INTEGER               REMAIN
C           INTEGER               START
C
C           LOGICAL               FOUND
C
C     C
C     C     Prompt for name of DSK and open file for reading.
C     C
C           CALL PROMPT ( 'Enter DSK name > ', DSK )
C
C           CALL DASOPR ( DSK, HANDLE )
C
C           CALL DLABFS ( HANDLE, DLADSC, FOUND )
C
C           IF ( .NOT. FOUND ) THEN
C
C              CALL SETMSG ( 'No segment found in file #.' )
C              CALL ERRCH  ( '#',  DSK                     )
C              CALL SIGERR ( 'SPICE(NOSEGMENT)'            )
C
C           END IF
C
C     C
C     C     Get segment vertex and plate counts.
C     C
C           CALL DSKZ02 ( HANDLE, DLADSC, NV, NP )
C
C           WRITE (*,*) ' '
C           WRITE (*,*) 'Number of vertices: ', NV
C           WRITE (*,*) 'Number of plates:   ', NP
C     C
C     C     Display the vertices of each plate.
C     C
C           REMAIN = NP
C           START  = 1
C
C           DO WHILE ( REMAIN .GT. 0 )
C     C
C     C        NREAD is the number of plates we'll read on this
C     C        loop pass.
C     C
C              NREAD  = MIN ( BUFSIZ, REMAIN )
C
C              CALL DSKP02 ( HANDLE, DLADSC, START, NREAD, N, PLATES )
C
C              DO I = 1, N
C
C                 PLIX = START + I - 1
C     C
C     C           Read the vertices of the current plate.
C     C
C                 DO J = 1, 3
C                    CALL DSKV02 ( HANDLE, DLADSC, PLATES(J,I),
C          .                       1,      NVTX,   VERTS (1,J)  )
C                 END DO
C     C
C     C           Display the vertices of the current plate:
C     C
C                 WRITE (*,*  ) ' '
C                 WRITE (*,*  ) 'Plate number: ', PLIX
C                 WRITE (*,FMT) '   Vertex 1: ', (VERTS(J,1), J = 1,3)
C                 WRITE (*,FMT) '   Vertex 2: ', (VERTS(J,2), J = 1,3)
C                 WRITE (*,FMT) '   Vertex 3: ', (VERTS(J,3), J = 1,3)
C
C     C
C     C           Display the normal vector of the current plate:
C     C
C                 CALL DSKN02 ( HANDLE, DLADSC, PLIX, NORMAL )
C
C                 WRITE (*,FMT) '   Normal:   ', (NORMAL(J), J = 1,3)
C
C              END DO
C
C              START  = START  + NREAD
C              REMAIN = REMAIN - NREAD
C
C           END DO
C
C     C
C     C     Close the kernel.  This isn't necessary in a stand-
C     C     alone program, but it's good practice in subroutines
C     C     because it frees program and system resources.
C     C
C           CALL DASCLS ( HANDLE )
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran/64bit
C     platform, using a DSK file representing a regular icosahedron,
C     the output was:
C
C
C     Enter DSK name > solid.bds
C
C      Number of vertices:           12
C      Number of plates:             20
C
C      Plate number:            1
C         Vertex 1:   0.000000000E+00  0.000000000E+00  0.117557000E+01
C         Vertex 2:   0.105146000E+01  0.000000000E+00  0.525731000E+00
C         Vertex 3:   0.324920000E+00  0.100000000E+01  0.525731000E+00
C         Normal:     0.491124160E+00  0.356821347E+00  0.794654382E+00
C
C      Plate number:            2
C         Vertex 1:   0.000000000E+00  0.000000000E+00  0.117557000E+01
C         Vertex 2:   0.324920000E+00  0.100000000E+01  0.525731000E+00
C         Vertex 3:  -0.850651000E+00  0.618034000E+00  0.525731000E+00
C         Normal:    -0.187592328E+00  0.577350079E+00  0.794654645E+00
C
C           ... 
C
C      Plate number:           20
C         Vertex 1:   0.850651000E+00 -0.618034000E+00 -0.525731000E+00
C         Vertex 2:   0.000000000E+00  0.000000000E+00 -0.117557000E+01
C         Vertex 3:   0.850651000E+00  0.618034000E+00 -0.525731000E+00
C         Normal:     0.607061680E+00  0.000000000E+00 -0.794654715E+00
C
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 02-JUN-2010 (NJB)
C
C-&
 
C$ Index_Entries
C
C     fetch model size parameters from a type 2 dsk segment
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C     
      INTEGER               N



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKZ02' )

      CALL DSKI02 ( HANDLE, DLADSC, KWNV, 1, 1, N, NV )
      CALL DSKI02 ( HANDLE, DLADSC, KWNP, 1, 1, N, NP )
 
      CALL CHKOUT ( 'DSKZ02' )
      RETURN
      END
