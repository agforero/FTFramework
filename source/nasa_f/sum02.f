C$Procedure   SUM02 ( DSKBRIEF, summarize type 2 segment )
 
      SUBROUTINE SUM02 ( HANDLE, DLADSC, NSIG )
 
C$ Abstract
C
C     Display type 2-specific summary of contents of a SPICE 
C     Digital Shape Kernel (DSK).
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
C     DSKBRIEF
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsk02.inc'

      INTEGER               HANDLE
      INTEGER               DLADSC ( DLADSZ )
      INTEGER               NSIG

 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK file.
C     DLADSC     I   DLA descriptor of segment.
C     NSIG       I   Number of significant digits in floating point
C                    output.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of a DSK file containing a segment
C                to be summarized.
C
C     DLADSC     is the DLA descriptor of a segment to be summarized.
C
C     NSIG       is the number of significant digits in floating point
C                numeric output. The range of NSIG is 6:17.
C
C$ Detailed_Output
C
C     None.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an invalid coarse voxel scale is encountered, this routine
C         signals the error SPICE(VALUEOUTOFRANGE).
C
C$ Files
C
C     See the input HANDLE.
C
C$ Particulars
C
C     This routine displays detailed summary information for a
C     specified type 2 DSK segment. The display is written to 
C     standard output.
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
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    DSKBRIEF Version 1.0.0, 05-OCT-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     summarize type 2 dsk segment
C
C-&



C
C     SPICELIB functions
C     
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 132 )

C
C     Local variables
C
      CHARACTER*(1)         COORD  ( 3 )
      CHARACTER*(LNSIZE)    LABELS ( 3 )
      CHARACTER*(LNSIZE)    OUTLIN
      CHARACTER*(LNSIZE)    TABLE  ( 3 )
      
      DOUBLE PRECISION      VOXORI ( 3 )
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VTXBDS ( 2, 3 )

      INTEGER               CGSCAL
      INTEGER               CGREXT ( 3 )
      INTEGER               I
      INTEGER               NCGR
      INTEGER               NP
      INTEGER               NV
      INTEGER               NVXTOT
      INTEGER               START1
      INTEGER               STARTS ( 2 )
      INTEGER               VOXNPL
      INTEGER               VOXNPT
      INTEGER               VTXNPL
      INTEGER               VGREXT ( 3 )

C
C     Saved variables
C
C     Note:  since this is a main program, SAVE statements are not
C     required.  They are used in case some of this code is later
C     ported to subroutines.
C
      SAVE                  COORD

C
C     Initial values
C     
      DATA COORD  / 'X', 'Y', 'Z' /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'SUM02' )


C
C     Display type 2 parameters.
C
      CALL TOSTDO ( ' ' )
 
      CALL TOSTDO ( '    Type 2 parameters' )
      CALL TOSTDO ( '    -----------------' )
 
      CALL DSKB02 ( HANDLE, DLADSC, NV,     NP,     NVXTOT,
     .              VTXBDS, VOXSIZ, VOXORI, VGREXT, CGSCAL,
     .              VTXNPL, VOXNPT, VOXNPL                  )
 
C
C     Show vertex and plate counts.
C
      OUTLIN = '      Number of vertices:                 #'
      CALL REPMI  ( OUTLIN, '#', NV, OUTLIN )
      CALL TOSTDO ( OUTLIN )
 
      OUTLIN = '      Number of plates:                   #'
      CALL REPMI  ( OUTLIN, '#', NP, OUTLIN )
      CALL TOSTDO ( OUTLIN )
 
C
C     Show voxel size.
C
      OUTLIN = '      Voxel edge length (km):             #'
      CALL REPMF  ( OUTLIN, '#', VOXSIZ, NSIG, 'E', OUTLIN )
      CALL TOSTDO ( OUTLIN )
 
C
C     Show voxel grid count.
C
      OUTLIN = '      Number of voxels:                   #'
      CALL REPMI  ( OUTLIN, '#', NVXTOT, OUTLIN )
      CALL TOSTDO ( OUTLIN )
 
C
C     Show coarse voxel grid count.
C
      IF ( CGSCAL .LT. 1 ) THEN
 
         CALL SETMSG ( 'Coarse voxel grid scale is #; ' //
     .                 'this scale should be an '       //
     .                 'integer > 0'                    )
         CALL ERRINT ( '#',  CGSCAL                     )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
 
      END IF
 
      NCGR = 1
 
      DO I = 1, 3
         CGREXT(I) = VGREXT(I) / CGSCAL
         NCGR      = NCGR * CGREXT(I)
      END DO
 
      OUTLIN = '      Number of coarse voxels:            #'
      CALL REPMI  ( OUTLIN, '#', NCGR, OUTLIN )
      CALL TOSTDO ( OUTLIN )
 
C
C     Show voxel grid extents.
C
      OUTLIN = '      Voxel grid X, Y, Z extents:         # # #'
 
      DO I = 1, 3
         CALL REPMI  ( OUTLIN, '#', VGREXT(I), OUTLIN )
      END DO
 
      CALL TOSTDO ( OUTLIN )
 
 
C
C     Show coarse voxel grid extents.
C
      OUTLIN = '      Coarse voxel grid X, Y, Z extents:  # # #'
 
      DO I = 1, 3
         CALL REPMI  ( OUTLIN, '#', CGREXT(I), OUTLIN )
      END DO
 
      CALL TOSTDO ( OUTLIN )
 
 
C
C     Show vertex bounds.
C
      DO I = 1, 3
 
         LABELS(I) =  '      Min, max vertex # value (km):'
         CALL REPMC( LABELS(I), '#', COORD(I), LABELS(I) )
 
      END DO

      START1 = 43

      CALL CORTAB ( 3, LABELS, START1, NSIG, 2, VTXBDS, STARTS, TABLE )
      
 
      DO I = 1, 3
         CALL TOSTDO( TABLE(I) )
      END DO

      CALL CHKOUT ( 'SUM02' )
      END



