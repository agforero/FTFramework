C$Procedure      ZZINILNK ( Initialize an AB cell linked-list )
 
      SUBROUTINE ZZINILNK ( MAXP, MAXC, NCELL, PNTRS, CELLS )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Initialize an AB cell linked-list structure.
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
C     None
C
C$ Keywords
C
C     DSK
C     Utility
C
C$ Declarations
 
      IMPLICIT NONE
 
      INTEGER               MAXP
      INTEGER               MAXC
      INTEGER               NCELL
      INTEGER               PNTRS  ( * )
      INTEGER               CELLS  ( 2, * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     MAXP       I   Length of A-value array.
C     MAXC       I   Number of cells.
C     NCELL      O   Number of cells in use; initialized to zero.
C     PNTRS      O   Initialized pointer array.
C     CELLS      O   Initialized cell array.
C
C$ Detailed_Input
C
C     MAXP       is the length of the pointer array PNTRS.
C
C     MAXC       is the length of the cell array CELLS. The array has
C                dimensions
C
C                   ( 2, MAXC )
C
C                Elements
C
C                   ( 1, * )
C
C                normally contain "A-values"; elements
C
C                   ( 2, * )
C
C                normally contain pointers.
C
C                
C$ Detailed_Output
C
C     NCELL      is the number of cells in use; this argument is
C                initialized to zero. This is done as a convenience
C                for the calling routine, which usually must keep
C                track of the number of cells in use.
C
C     PNTRS      is the pointer array, initialized so that each
C                element contains -1, which represents a null
C                pointer.
C
C     CELLS      is a cell array, initialized so that each pointer
C                element contains -1, and so that each value 
C                element contains 0.       
C
C$ Parameters
C
C     See the DSK type 2 include file 
C
C        dsk02.inc
C
C$ Exceptions
C
C     1)  If the pointer array size MAXP is non-positive, the 
C         error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     2)  If the cell array size MAXC is non-positive, the 
C         error SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports creation of DSK type 2 segments.
C     It is used for creation of both voxel-plate association
C     data structures and of vertex-plate association data
C     structures.
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C        Updated version info.
C     
C        05-JAN-2016 (NJB)
C
C           Added error checks on input dimension arguments. Re-wrote
C           some of the argument descriptions. Added private routine
C           warning to header abstract.
C
C        08-OCT-2009 (NJB)
C
C           Re-ordered header sections.
C
C        03-FEB-1999 (JAB)
C
C-&
 
C$ Index_Entries
C
C     initialize an AB cell linked-list
C
C-&
 
      INTEGER               I
      LOGICAL               RETURN

      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZINILNK' )

      IF ( MAXP .LT. 1 ) THEN

         CALL SETMSG ( 'Pointer array size MAXP = #; size must ' 
     .   //            'be positive.'                           )
         CALL ERRINT ( '#',  MAXP                               )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZINILNK'                               )
         RETURN  
         
      END IF

      IF ( MAXC .LT. MAXP ) THEN

         CALL SETMSG ( 'Cell array size MAXC = #; size must ' 
     .   //            'be at least as large as pointer array '
     .   //            'size #.'                                )
         CALL ERRINT ( '#',  MAXC                               )
         CALL ERRINT ( '#',  MAXP                               )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZINILNK'                               )
         RETURN  
         
      END IF

C
C     Initialize pointer array and cells.
C 
      DO I = 1, MAXP
         PNTRS(I)   = -1
      END DO
 
      DO I = 1, MAXC
         CELLS(1,I) =  0
         CELLS(2,I) = -1
      END DO

C
C     Set count of cells in use to 0, for the convenience
C     of the calling routine.
C      
      NCELL = 0

      CALL CHKOUT ( 'ZZINILNK' )
      RETURN
      END
