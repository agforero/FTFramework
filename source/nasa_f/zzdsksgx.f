C$Procedure ZZDSKSGX ( DSK, ray-segment intercept )
 
      SUBROUTINE ZZDSKSGX ( HANDLE, DLADSC, DTYPE, ET, VERTEX,
     .                      RAYDIR, XPT,    DC,    IC, FOUND   )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find the intersection of a ray and the surface described by
C     a single DSK segment.
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
C     GEOMETRY
C     INTERCEPT
C     INTERSECTION
C     RAY
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE
      
      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      INTEGER               DTYPE
      DOUBLE PRECISION      ET
      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      XPT    ( 3 )
      DOUBLE PRECISION      DC     ( * )
      INTEGER               IC     ( * )
      LOGICAL               FOUND

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   DSK file handle.   
C     DLADSC     I   DLA descriptor of segment.
C     DTYPE      I   Data type code.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     XPT        O   Surface intercept, if found.
C     DC         O   D.p. component of source info.
C     IC         O   Integer component of source info.
C     FOUND      O   Found flag. 
C
C$ Detailed_Input
C
C     HANDLE     is the handle of a DSK file containing a segment
C                to be used in a ray-surface intercept computation.
C
C     DLASDC     is the DLA descriptor of the DSK segment to be used.
C
C     DTYPE      is the data type code of the segment. While this
C                information can be retrieved from the DSK descriptor
C                of the segment, the availability of this argument
C                saves the time needed to do so.
C
C     ET         is the epoch of the intersection computation,
C                expressed as seconds past J2000 TDB. This epoch is
C                used for DSK segment selection.
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray's direction vector must be
C                represented in the reference frame of the segment. The
C                vertex is considered to be an offset from the center
C                of the reference frame associated with the segment.
C 
C$ Detailed_Output
C
C     XPT        is the intercept of the ray on the surface described
C                by the segment, if such an intercept exists. If the
C                ray intersects the surface at multiple points, the
C                one closest to the ray's vertex is selected. XPT is
C                valid if and only if FOUND is .TRUE.
C
C                XPT is expressed in the reference frame associated
C                with the specified segment. It represents an offset
C                from the center of this frame. Note that the frame
C                center may differ from the central body of the
C                segment.
C
C
C     DC         is the double precision component of the data
C                source information. Contents are data type-
C                dependent. DC is valid if and only if FOUND 
C                is .TRUE.
C
C
C     IC         is the integer component of the data
C                source information. Contents are data type-
C                dependent. IC is valid if and only if FOUND 
C                is .TRUE.
C
C                For type 2 segments, IC contains just the
C                intercept plate ID in element 1.
C
C
C     FOUND      is a logical flag that is set to .TRUE. if and only
C                if a ray-surface intercept was found.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the DSK segment data type is not recognized, the error
C         SPICE(TYPENOTSUPPORTED) is signaled.
C          
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - DSK data: the DSK file designated by HANDLE and containing
C          the segment having the DLA descriptor DLADSC must be loaded
C          at the time this routine is called.
C
C     Kernel data are normally loaded once per program run, NOT every
C     time this routine is called.
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
C-    SPICELIB Version 1.0.0, 18-FEB-2016 (NJB) 
C
C        Based on first version 20-JAN-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find intercept of ray with surface defined by dsk segment
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      TOUCHD

      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      RETVAL

      INTEGER               PLID


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKSGX' )

C
C     Note: input argument ET is provided to support time-dependent
C     data types.
C     
      RETVAL = TOUCHD( ET )
      DC(1)  = TOUCHD( DC )

      IF ( DTYPE .EQ. 2 ) THEN
C
C        The intercept plate ID is returned in element 1 of
C        IC, if an intercept is found.
C
         CALL DSKX02 ( HANDLE, DLADSC, VERTEX, 
     .                 RAYDIR, PLID,   XPT,    FOUND )

         IF ( FOUND ) THEN
            IC(1) = PLID
         END IF

      ELSE

         CALL SETMSG ( 'DSK ray-surface intercepts are not '
     .   //            'supported for DSK data type #.'     )
         CALL ERRINT ( '#', DTYPE                           )
         CALL SIGERR ( 'SPICE(TYPENOTSUPPORTED)'            )
         CALL CHKOUT ( 'ZZDSKSGX'                           )
         RETURN

      END IF


      CALL CHKOUT ( 'ZZDSKSGX' )
      RETURN
      END
      
