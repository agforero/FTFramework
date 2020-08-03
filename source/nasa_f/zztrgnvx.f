C$Procedure    ZZTRGNVX ( MKDSK, compute target voxel counts )
 
      SUBROUTINE ZZTRGNVX ( NP, TRGCOR, TRGFIN )
  
C$ Abstract
C
C
C     Compute target coarse and fine voxel counts for a DSK type 2
C     segment.
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

      INTEGER               NP
      INTEGER               TRGCOR
      INTEGER               TRGFIN



C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NP         O   Plate count.
C     TRGCOR     O   Target coarse voxel count.
C     TRGFIN     O   Target fine voxel count.
C
C$ Detailed_Input
C
C     NP             is the plate count of a type 2 DSK segment.
C
C$ Detailed_Output
C
C     TRGCOR         is the target coarse voxel count. This is an
C                    approximation to the final count; the actual count
C                    set by MKDSK will be smaller than, but near,
C                    TRGCOR.
C
C     TRGFIN         is the target fine voxel count. This is an
C                    approximation to the final count; the actual count
C                    set by MKDSK will be smaller than, but near,
C                    TRGFIN.
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
C$ Files
C
C     None.
C
C$ Particulars
C
C     The target coarse and fine voxel counts produced by this routine
C     are set by heuristics. The counts are meant to yield reasonable
C     DSK type 2 ray-surface intercept computation speed.
C
C     Actual counts selected by MKDSK will be less than the respective
C     target counts.
C
C     See the source code for details.
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
C-    MKDSK Version 1.0.0, 19-FEB-2017 (NJB)
C
C-&
 
C$ Index_Entries
C
C     compute type 2 dsk target voxel counts
C
C-&


C
C     SPICELIB functions
C
      INTEGER               BRCKTI
      LOGICAL               RETURN

C
C     Local parameters
C

      DOUBLE PRECISION      FS0
      PARAMETER           ( FS0   = 1.2D7 )

      DOUBLE PRECISION      NPLIM
      PARAMETER           ( NPLIM = 1.D6 )

      DOUBLE PRECISION      FAC0
      PARAMETER           ( FAC0  = (MAXPLT-FS0)/(MAXPLT-NPLIM) )

C
C     Local variables
C
      DOUBLE PRECISION      FRAC



      IF ( RETURN() ) THEN
         RETURN   
      END IF
      
      CALL CHKIN ( 'ZZTRGNVX' )

C
C     All of the formulas below for target voxel counts are heuristics.
C     They may be updated if superior formulas are developed.
C
      IF ( NP .GT. NPLIM ) THEN
C
C        For plate counts of NPLIM or more, we go for a full-size
C        coarse voxel grid.
C
         TRGCOR = MAXCGR

C
C        The fine voxel count increases linearly as a function
C        of NP - NPLIM. The count at NP = MAXPLT is 2*MAXPLT.
C
         TRGFIN = NINT(FS0)  +  NINT( ( NP - NPLIM ) * FAC0 * 2 )

      ELSE IF ( NP .GT. 1000 ) THEN
C
C        Scale down the coarse grid proportionally to 
C        the square root of the ratio of NP to NPLIM.
C
         FRAC   = DBLE(NP) / NPLIM

         TRGCOR = NINT( MAXCGR * SQRT(FRAC) )

C
C        Scale down the fine voxel count as well.
C
         TRGFIN = NINT( FS0 * FRAC )


      ELSE IF ( NP .GT. 100 ) THEN

         TRGCOR = 100
         TRGFIN = 1000 

      ELSE

         TRGCOR = 100
         TRGFIN = 100

      END IF

C
C     Bracket the estimates to ensure they're in range.
C
      TRGCOR = BRCKTI ( TRGCOR, 1, MAXCGR )
      TRGFIN = BRCKTI ( TRGFIN, 8, MAXVOX )

            
      CALL CHKOUT ( 'ZZTRGNVX' )
      RETURN
      END
