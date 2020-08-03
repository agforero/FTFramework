C$Procedure ZZRYTELT ( DSK, ray touches coordinate volume element )
 
      SUBROUTINE ZZRYTELT ( VERTEX, RAYDIR, DSKDSC, 
     .                      MARGIN, NXPTS,  XPT    )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find nearest intersection to ray's vertex of ray and
C     a coordinate volume element. If the vertex is inside
C     the element, the vertex is considered to be the solution.
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
C     GEOMETRY
C     INTERCEPT
C     INTERSECTION
C     RAY
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'

      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      DSKDSC ( DSKDSZ )
      DOUBLE PRECISION      MARGIN
      INTEGER               NXPTS
      DOUBLE PRECISION      XPT    ( 3 )
 

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     DSKDSC     I   DSK descriptor of segment.
C     MARGIN     I   Margin used for element expansion.
C     NXPTS      O   Number of intercept points.
C     XPT        O   Intercept.
C
C$ Detailed_Input
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray direction must be represented
C                in the reference frame of the segment to which the
C                volume element boundaries correspond. The vertex is
C                considered to be an offset from the center of the
C                reference frame associated with the segment.
C 
C
C     DSKDSC     is a DSK segment descriptor. The coordinate system
C                and spatial boundaries of a the segment are specified
C                by members of this descriptor.
C
C
C$ Detailed_Output
C
C     XPT        is the intercept of the ray on the boundary of the
C                input volume element, if such an intercept exists. If
C                the ray's vertex is inside the element, XPT is set
C                equal to the vertex. XPT is valid if and only if FOUND
C                is .TRUE.
C
C                XPT is expressed in the reference frame associated
C                with the inputs VERTEX and RAYDIR. XPT represents
C                an offset from the origin of the coordinate system.
C
C                XPT is valid only if NXPTS is set to 1.
C
C
C     NXPTS      is the number of intercept points of the ray and
C                the volume element. 
C
C                Currently there are only two possible values for
C                NXPTS:
C
C                   1 for an intersection
C                   0 for no intersection
C    
C                If the vertex is inside the element, NXPTS is
C                set to 1.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the coordinate system code in the segment descriptor
C        is not recognized, the error SPICE(NOTSUPPORTED) will be
C        signaled.
C              
C$ Files
C
C     None. However, the input segment boundaries normally have
C     been obtained from a loaded DSK file.
C     
C$ Particulars
C
C     This routine sits on top of data DSK type-specific ray-segment
C     intercept routines such as DSKX02.
C
C$ Examples
C
C     See usage in ZZDSKBUX.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by the DSK
C     subsystem.
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
C-    SPICELIB Version 1.0.0, 19-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     find intercept of ray on dsk volume element
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C
      INTEGER               CORSYS


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZRYTELT' )


      CORSYS = NINT( DSKDSC(SYSIDX) )
 
      IF ( CORSYS .EQ. LATSYS ) THEN
 
         CALL ZZRYTLAT ( VERTEX, RAYDIR, DSKDSC(MN1IDX),
     .                   MARGIN, NXPTS,  XPT              )
 

      ELSE IF ( CORSYS .EQ. RECSYS ) THEN
 
         CALL ZZRYTREC ( VERTEX, RAYDIR, DSKDSC(MN1IDX),
     .                   MARGIN, NXPTS,  XPT              )
 

      ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
 
         CALL ZZRYTPDT ( VERTEX,         RAYDIR,
     .                   DSKDSC(MN1IDX), DSKDSC(PARIDX),
     .                   MARGIN,         NXPTS,          XPT )
      ELSE
 
         CALL SETMSG( 'Coordinate system # is not supported.' )
         CALL ERRINT( '#', CORSYS                             )
         CALL SIGERR( 'SPICE(BADCOORDSYS)'                    )
         CALL CHKOUT( 'ZZRYTELT'                              )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZRYTELT' )
      RETURN
      END 
