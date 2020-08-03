C$Procedure ZZSINUTL ( DSK, utilities for generalized ray intercept )
 
      SUBROUTINE ZZSINUTL ( TRGCDE, NSURF,  SRFLST, ET,    FIXFID,
     .                      VERTEX, RAYDIR, SPOINT, FOUND, MINRAD,
     .                      MAXRAD, PNEAR,  DIST                  )
  
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     This is the umbrella routine for utilities supporting the 
C     generalized ray-surface intercept capability used by SINCPT.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations

      IMPLICIT NONE
      
      INCLUDE 'dsk.inc'
      
      INTEGER               TRGCDE
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      DOUBLE PRECISION      ET
      INTEGER               FIXFID
      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      SPOINT ( 3 )
      LOGICAL               FOUND
      DOUBLE PRECISION      MINRAD
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      PNEAR  ( 3 )
      DOUBLE PRECISION      DIST
 
C$ Brief_I/O
C
C     Variable  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     TRGCDE     I   ZZSUELIN, ZZSUDSKI
C     NSURF      I   ZZSUELIN, ZZSUDSKI
C     SRFLST     I   ZZSUELIN, ZZSUDSKI
C     ET         I   ZZRAYSFX, ZZRAYNP
C     FIXFID     I   ZZSUDSKI
C     VERTEX     I   ZZRAYSFX, ZZRAYNP
C     RAYDIR     I   ZZRAYSFX, ZZRAYNP
C     SPOINT     O   ZZRAYSFX
C     FOUND      O   ZZRAYSFX
C     MINRAD     O   ZZMINRAD
C     MAXRAD     O   ZZMAXRAD
C     PNEAR      O   ZZRAYNP
C     DIST       O   ZZRAYNP
C
C$ Detailed_Input
C
C     See the entry points.
C
C$ Detailed_Output
C
C     See the entry points. 
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If this routine is called directly, it signals the error 
C         SPICE(BOGUSENTRY).
C        
C     See the entry points for descriptions of errors specific to
C     those routines.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C
C     This routine uses PCK data for target body reference
C     ellipsoids.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine contains the following entry points that support
C     the generalized ray-surface intercept algorithm used by SINCPT:
C
C        ZZSUELIN:   Initialize ray intercept and near point utilities 
C                    for ellipsoidal target surface
C
C
C        ZZSUDSKI:   Initialize ray intercept and near point utilities 
C                    for DSK target surface
C
C        ZZMINRAD:   Return minimum spherical bounding radius for
C                    surface
C
C        ZZMAXRAD:   Return maximum spherical bounding radius for
C                    surface
C
C        ZZRAYSFX:   Compute ray-surface intercept given inputs
C                    set by ZZSUELIN or ZZSUDSKI
C     
C
C        ZZRAYNP:    Compute nearest point on surface to ray given
C                    inputs set by ZZSUELIN or ZZSUDSKI
C
C
C     The geometric computation routines ZZRAYSFX, ZZRAYNP are
C     used as callback routines by the generalized SINCPT algorithm.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1) This is a private routine. 
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
C-    SPICELIB version 1.0.0 01-JUN-2016 (NJB)
C
C        Updated comments. Updated error handling in ZZSUDSKI.
C
C        11-FEB-2016 (NJB) 
C
C           Added header.
C
C        02-NOV-2015 (NJB)
C
C           No longer requires ellipsoidal radii for 
C           initialization for DSK intercept computations.
C
C        10-OCT-2015 (NJB)
C
C           Now includes dsk.inc.
C
C        30-JAN-2015 (NJB)
C
C           Updated to provide an entry point returning a global
C           minimum surface radius. The argument list has been
C           changed to accommodate this new output.
C
C        17-OCT-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     umbrella for generalized sincpt utilities
C
C-&

      
C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
C
C     Target types:
C
      INTEGER               ELLTRG
      PARAMETER           ( ELLTRG = 1 )

      INTEGER               DSKTRG
      PARAMETER           ( DSKTRG = 2 )

C
C     Local variables
C
      DOUBLE PRECISION      SAVMNR
      DOUBLE PRECISION      SAVMXR
      DOUBLE PRECISION      SAVRAD ( 3 )

      INTEGER               N
      INTEGER               SAVFID
      INTEGER               SAVNSF
      INTEGER               SAVSRF ( MAXSRF )
      INTEGER               SAVTRG
      INTEGER               SAVTYP

C
C     Saved variables
C
      SAVE                  SAVFID
      SAVE                  SAVMNR
      SAVE                  SAVMXR
      SAVE                  SAVNSF
      SAVE                  SAVRAD
      SAVE                  SAVSRF
      SAVE                  SAVTRG
      SAVE                  SAVTYP

C
C     Initial values
C

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZSINUTL' )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)' )
      CALL CHKOUT ( 'ZZSINUTL' )
      RETURN




C$Procedure ZZSUELIN ( DSK, initialize SINCPT utilities for ellipsoid )
 
      ENTRY ZZSUELIN ( TRGCDE )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Initialize generalized intercept utilities, which are entry
C     points in this package, by storing the target ID code and
C     ellipsoid radii to be used by the callback entry points of this
C     package.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     INTEGER               TRGCDE
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     TRGCDE     I   Target body ID code.
C
C$ Detailed_Input
C     
C     TRGCDE     is the integer body ID of the target body for which
C                ray-surface intercept or ray-surface near point
C                computations are to be performed by the callback
C                entry points of this package.
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
C     1)  If an error occurs during lookup of the target body's radii,
C         an error will be signaled by a routine in the call tree
C         of this routine.
C        
C     2)  If the body-fixed frame used by the caller of this routine
C         does not have its axes properly aligned with the target
C         body's reference ellipsoid, results from this routine will be
C         invalid.
C       
C         This routine does not have the capability of detecting such
C         an error.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     This routine expects the radii of the target body's reference
C     ellipsoid to be present in the kernel pool.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine supports the generalized ray-surface intercept
C     algorithm used by SINCPT.
C     
C     This routine prepares the callback entry points for computations
C     using an ellipsoidal shape model. It obtains radii for ellipsoid
C     computations: ray-ellipsoid intersection computation using SURFPT
C     and ray-ellipsoid altitude computation using NPEDLN.
C
C     Note the absence of a frame ID from the argument list. This
C     routine assumes the input ray vectors are expressed in a
C     body-fixed frame centered at the target body's center, having its
C     axes aligned with the principal axes of the target body's
C     reference ellipsoid.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1) This is a private routine. 
C
C     2) See exception (2) above.
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
C-    SPICELIB Version 1.0.0, 01-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     initialize intercept utilities for ellipsoidal target
C
C-&


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSUELIN' )

      SAVTYP = ELLTRG

      CALL BODVCD ( TRGCDE, 'RADII', 3, N, SAVRAD )
      
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSUELIN' )
         RETURN
      END IF

      SAVMNR = MIN( SAVRAD(1), SAVRAD(2), SAVRAD(3) )
      SAVMXR = MAX( SAVRAD(1), SAVRAD(2), SAVRAD(3) )
      
      CALL CHKOUT( 'ZZSUELIN' )
      RETURN

     



C$Procedure ZZSUDSKI ( DSK, initialize SINCPT utilities for dsk target )
 
      ENTRY ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FIXFID )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Initialize generalized intercept utilities, which are entry
C     points in this package, by storing the target ID code, DSK
C     surface ID list, and frame ID to be used by the callback entry
C     points of this package. This routine is to be used for 
C     computations in which the target surface is modeled using
C     DSK data.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     INTEGER               TRGCDE
C     INTEGER               NSURF
C     INTEGER               SRFLST ( * )
C     INTEGER               FIXFID
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     TRGCDE     I   Target body ID code.
C     NSURF      I   Number of surface IDs in surface list.
C     SRFLST     I   Surface list.
C     FIXFID     I   Target body-fixed frame ID.
C     MAXSRF     P   Maximum size of surface ID list.
C
C$ Detailed_Input
C     
C     TRGCDE     is the integer body ID of the target body for which
C                ray-surface intercept or ray-surface near point
C                computations are to be performed by the callback
C                entry points of this package.
C
C     NSURF,
C     SRFLST     are, respectively, the count of surface IDs in the
C                surface ID list and the list itself.
C
C                If the list is empty, all surfaces associated with
C                the target body are used.
C
C     FIXFID     is the frame ID code of a body-fixed reference frame
C                associated with the target body. The frame must be
C                centered at the target body's center. 
C
C$ Detailed_Output
C
C     None.
C
C$ Parameters
C
C     MAXSRF     is the maximum supported size of a surface ID list.
C                See the include file dsk.inc for the parameter's
C                value.
C
C$ Exceptions
C
C     1)  If the number of IDs in the surface list is negative or
C         exceeds MAXSRF, the error SPICE(INVALIDCOUNT) will be
C         signaled.
C
C     2)  If an error occurs during calculation of the radii of the
C         inner and outer bounding spheres for the target body,
C         the error will be signaled by a routine in the call tree
C         of this routine.
C
C     3)  The reference frame designated by FIXFID must be centered
C         at the target body's center. This condition must be
C         checked by the caller of this routine.
C        
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine supports the generalized ray-surface intercept
C     algorithm used by SINCPT.
C     
C     This routine prepares the callback entry points for computations
C     using a DSK shape model. It stores the target body ID, DSK
C     surface ID list, and ID of the target body-fixed frame. It 
C     also computes the radii of the inner and output bounding
C     surfaces for the target body; these surfaces are computed using
C     loaded DSK segments for the surfaces indicated by the input
C     body ID and surface ID list.
C
C     The stored quantities are used by the callback entry points of
C     this package.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1) This is a private routine. 
C
C     2) See exception (2) above.
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
C-    SPICELIB Version 1.0.0, 01-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     initialize intercept utilities for dsk target
C
C-&

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSUDSKI' )

      SAVTYP = DSKTRG

      IF (  ( NSURF .LT. 0 ) .OR. ( NSURF .GT. MAXSRF )  ) THEN

         CALL SETMSG ( 'Surface count must be in the range '
     .   //            '0:# but was #.'                     )
         CALL ERRINT ( '#', MAXSRF                          )
         CALL ERRINT ( '#', NSURF                           )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                )
         CALL CHKOUT ( 'ZZSUDSKI'                           )
         RETURN

      END IF

      SAVNSF = NSURF

      CALL MOVEI ( SRFLST, SAVNSF, SAVSRF )

      SAVFID = FIXFID
      SAVTRG = TRGCDE

      CALL CLEARD ( 3, SAVRAD )


      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSUDSKI' )
         RETURN
      END IF

C
C     Fetch minimum and maximum radius of target body surface.
C
      CALL ZZDSKSPH ( TRGCDE, SAVNSF, SAVSRF, SAVMNR, SAVMXR )

      CALL CHKOUT( 'ZZSUDSKI' )
      RETURN





C$Procedure ZZRAYSFX ( DSK, callback for ray-surface intercept )
 
      ENTRY ZZRAYSFX ( VERTEX, RAYDIR, ET, SPOINT, FOUND )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Perform ray-surface intercept using values set via the 
C     initialization entry points of this package. This routine
C     is used as a callback by the generalized ray-surface intercept
C     algorithm.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     DOUBLE PRECISION      VERTEX ( 3 )
C     DOUBLE PRECISION      RAYDIR ( 3 )
C     DOUBLE PRECISION      ET
C     DOUBLE PRECISION      SPOINT ( 3 )
C     LOGICAL               FOUND
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     ET         I   Evaluation epoch, seconds past J2000 TDB.
C     SPOINT     O   Surface intercept point.
C     FOUND      O   Found flag. True if intercept exists.
C
C$ Detailed_Input
C     
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                a ray. When the target's surface is represented by DSK
C                data, Both vectors are expressed in the frame
C                designated by FIXFID, which is set via a call to one
C                of the initialization routines of this package.
C
C                When the target's surface is represented by an
C                ellipsoid, the vectors are presumed to be expressed in
C                a body-fixed frame compatible with that ellipsoid.
C
C     ET         is the epoch for which the computation is to be
C                performed. This epoch is used for DSK segment
C                selection; only segments containing ET in their time
C                coverage interval will be used. ET is expressed as
C                seconds past J2000 TDB.
C
C$ Detailed_Output
C
C     SPOINT     is the surface intercept on the target body 
C                nearest to the ray's vertex, if the intercept
C                exists.
C
C     FOUND      is a logical flag that is set to .TRUE. if and
C                only if an intercept exists.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs while this routine attempts to compute a
C         surface intercept, the error will be signaled by a routine in
C         the call tree of this routine.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C
C     This routine uses PCK data for target body reference ellipsoids.
C
C     If an ellipsoidal target model is selected, this routine expects
C     the radii of the target body's reference ellipsoid to be present
C     in the kernel pool.
C      
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine prepares the local buffers for a ray-surface
C     intercept computation using unprioritized DSK data. It calls
C     ZZSBFXR to perform the computation.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1)  This is a private routine. 
C
C     2)  One of the initialization routines ZZSUDSKI or ZZSUELIN
C         must be called before the first time this routine is called. 
C         Whenever the set of data to be considered changes, an
C         initialization call must be made before this routine may
C         be called.
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
C-    SPICELIB Version 1.0.0, 01-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     callback for generalized ray surface intercept 
C
C-&


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZRAYSFX' )

      IF ( SAVTYP .EQ. ELLTRG ) THEN 

         CALL SURFPT ( VERTEX,    RAYDIR,    SAVRAD(1), 
     .                 SAVRAD(2), SAVRAD(3), SPOINT,    FOUND )


      ELSE IF ( SAVTYP .EQ. DSKTRG ) THEN

         CALL ZZSBFXR ( SAVTRG, SAVNSF, SAVSRF, ET,
     .                  SAVFID, VERTEX, RAYDIR, SPOINT, FOUND )

      ELSE 
         
         CALL SETMSG ( 'Surface type code # is not supported. '
     .   //            'This code branch is not supposed to be '
     .   //            'reached.'                               )
         CALL ERRINT ( '#', SAVTYP                              )
         CALL SIGERR ( 'SPICE(BUG)'                             )
         CALL CHKOUT ( 'ZZRAYSFX'                               )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZRAYSFX' )
      RETURN




C$Procedure ZZMAXRAD ( DSK, maximum radius )
 
      ENTRY ZZMAXRAD ( MAXRAD )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Return the radius of an outer bounding sphere for the target body.
C     The radius is not necessarily a minimum upper bound.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     DOUBLE PRECISION      MAXRAD
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     MAXRAD     O   Radius of outer bounding sphere for target body.
C
C$ Detailed_Input
C     
C     MAXRAD     is the radius of an outer bounding sphere for the
C                target body. The sphere is centered at the center of
C                the target body.
C
C                If the target's surface is modeled as an ellipsoid,
C                this radius is the maximum radius of the ellipsoid. If
C                the target's surface is modeled using DSK data, the
C                radius is determined by the surface list used at
C                initialization time.
C
C                The radius is not necessarily a least upper bound.
C
C                Units are km.
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
C     1)  If a bounding sphere has not been computed due to
C         initialization failure, the value returned by this
C         routine will be invalid.
C
C$ Files
C
C     None. However, see the $Files sections of routines ZZSUELIN
C     and ZZSUDSKI.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine supports the generalized ray-surface intercept
C     algorithm used by SINCPT.
C
C     Note: this routine does not accept an epoch as an input argument,
C     so when the last setup was performed by ZZSUDSKI, this routine
C     computes its result using all DSK segments for the target body
C     (the ID code of which was passed to ZZSUDSKI). This could result
C     in the computed maximum radius being substantially larger than
C     the maximum radius that would be obtained by considering only
C     segments covering the epoch of interest.
C
C     The routine considers all segments for the target in in the
C     interest of efficiency: it simply returns a value that was
C     computed during the last setup call. If the routine were
C     time-dependent, it would need to re-compute the maximum radius
C     each time it was called.
C
C     If a need arises for a time-dependent version of this routine,
C     that routine should be given a new name, and the routine should
C     be added to this package.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1)  This is a private routine. 
C
C     2)  One of the initialization routines ZZSUDSKI or ZZSUELIN
C         must be called before the first time this routine is called. 
C         Whenever the set of data to be considered changes, an
C         initialization call must be made before this routine may
C         be called.
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
C-    SPICELIB Version 1.0.0, 01-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     return radius of maximum bounding sphere for target
C
C-&

      MAXRAD = SAVMXR

      RETURN





C$Procedure ZZMINRAD ( DSK, minimum radius )
 
      ENTRY ZZMINRAD ( MINRAD )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Return the radius of an inner bounding sphere for the target body.
C     The radius is not necessarily a maximum lower bound.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     DOUBLE PRECISION      MINRAD
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     MINRAD     O   Radius of inner bounding sphere for target body.
C
C$ Detailed_Input
C     
C     MINRAD     is the radius of the inner bounding sphere for the
C                target body. The sphere is centered at the center of
C                the target body.
C
C                If the target's surface is modeled as an ellipsoid,
C                this radius is the minimum radius of the ellipsoid. If
C                the target's surface is modeled using DSK data, the
C                radius is determined by the surface list used at
C                initialization time.
C
C                The radius is not necessarily a least upper bound.
C
C                Units are km.
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
C     1)  If a bounding sphere has not been computed due to
C         initialization failure, the value returned by this
C         routine will be invalid.
C
C$ Files
C
C     None. However, see the $Files sections of routines ZZSUELIN
C     and ZZSUDSKI.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine supports the generalized ray-surface intercept
C     algorithm used by SINCPT.
C
C     Note: this routine does not accept an epoch as an input argument,
C     so when the last setup was performed by ZZSUDSKI, this routine
C     computes its result using all DSK segments for the target body
C     (the ID code of which was passed to ZZSUDSKI). This could result
C     in the computed minimum radius being substantially smaller than
C     the minimum radius that would be obtained by considering only
C     segments covering the epoch of interest.
C
C     The routine considers all segments for the target in in the
C     interest of efficiency: it simply returns a value that was
C     computed during the last setup call. If the routine were
C     time-dependent, it would need to re-compute the minimum radius
C     each time it was called.
C
C     If a need arises for a time-dependent version of this routine,
C     that routine should be given a new name, and the routine should
C     be added to this package.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1)  This is a private routine. 
C
C     2)  One of the initialization routines ZZSUDSKI or ZZSUELIN
C         must be called before the first time this routine is called. 
C         Whenever the set of data to be considered changes, an
C         initialization call must be made before this routine may
C         be called.
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
C-    SPICELIB Version 1.0.0, 01-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     return radius of minimum bounding sphere for target
C
C-&


      MINRAD = SAVMNR

      RETURN





C$Procedure ZZRAYNP ( DSK, callback for ray-surface near point )
 
      ENTRY ZZRAYNP ( VERTEX, RAYDIR, ET, PNEAR, DIST )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Compute near point to ray on reference ellipsoid or outer
C     bounding sphere. If the ray intersects the surface, the intercept
C     is returned, and the returned distance is set to zero. This
C     routine is used as a callback by the generalized ray-surface
C     intercept algorithm.
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
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     DOUBLE PRECISION      VERTEX ( 3 )
C     DOUBLE PRECISION      RAYDIR ( 3 )
C     DOUBLE PRECISION      ET
C     DOUBLE PRECISION      PNEAR  ( 3 )
C     DOUBLE PRECISION      DIST
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     ET         I   Evaluation epoch, seconds past J2000 TDB.
C     PNEAR      O   Near point to ray on surface.
C     DIST       O   Distance between near point and ray.
C
C$ Detailed_Input
C     
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                a ray. When the target's surface is represented by DSK
C                data, Both vectors are expressed in the frame
C                designated by FIXFID, which is set via a call to one
C                of the initialization routines of this package.
C
C                When the target's surface is represented by an
C                ellipsoid, the vectors are presumed to be expressed in
C                a body-fixed frame compatible with that ellipsoid.
C
C     ET         is the epoch for which the computation is to be
C                performed. This epoch is used for DSK segment
C                selection; only segments containing ET in their time
C                coverage interval will be used. ET is expressed as
C                seconds past J2000 TDB.
C
C$ Detailed_Output
C
C     PNEAR      is, when the target shape is modeled by a reference
C                ellipsoid, the point on the target body nearest to the
C                ray, if the ray does not intersect the body. If a
C                ray-surface intercept exists, PNEAR is set to the
C                intercept closest to the ray's vertex.
C
C                When the target shape is modeled using DSK data, the
C                computation is performed with an outer bounding sphere
C                for the target used in place of the target's reference
C                ellipsoid. This routine does not attempt to find the
C                closest point on the topographic surface to the ray.
C
C                PNEAR is expressed in the reference frame associated
C                with the input ray's vertex and direction vectors.
C
C     DIST       is the distance between PNEAR and VERTEX. Units are
C                km.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs while this routine attempts to 
C         find the nearest point on the target surface to the
C         input ray, the error will be signaled by a routine
C         in the call tree of this routine.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem.
C
C     If an ellipsoidal target model is selected, this routine expects
C     the radii of the target body's reference ellipsoid to be present
C     in the kernel pool.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine prepares the local buffers for a ray-surface
C     near point computation using unprioritized DSK data. When
C     DSK data are used, an outer bounding sphere for the target
C     is used to model the target shape for the purpose of
C     this computation. The sphere is that computed by the 
C     most recently called initialization entry point.
C
C     This routine calls NPEDLN to perform the computation for
C     both the ellipsoid and DSK cases.
C
C$ Examples
C
C     See usage in SINCPT, ZZSFXCOR.
C
C$ Restrictions
C
C     1)  This is a private routine. 
C
C     2)  One of the initialization routines ZZSUDSKI or ZZSUELIN
C         must be called before the first time this routine is called. 
C         Whenever the set of data to be considered changes, an
C         initialization call must be made before this routine may
C         be called.
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
C-    SPICELIB Version 1.0.0, 01-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     callback for generalized ray surface near point 
C
C-&

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZRAYNP' )


      IF ( SAVTYP .EQ. ELLTRG ) THEN 
C
C        Find the nearest point on the ellipsoid's surface to
C        to the ray.
C
         CALL NPEDLN ( SAVRAD(1), SAVRAD(2), SAVRAD(3), 
     .                 VERTEX,    RAYDIR,    PNEAR,    DIST )


      ELSE IF ( SAVTYP .EQ. DSKTRG ) THEN
C
C        Find the nearest point on the outer bounding sphere to
C        to the ray.
C
         CALL NPEDLN ( SAVMXR, SAVMXR, SAVMXR, 
     .                 VERTEX, RAYDIR, PNEAR,  DIST )

      ELSE 
         
         CALL SETMSG ( 'Surface type code # is not supported. '
     .   //            'This code branch is not supposed to be '
     .   //            'reached.'                               )
         CALL ERRINT ( '#', SAVTYP                              )
         CALL SIGERR ( 'SPICE(BUG)'                             )
         CALL CHKOUT ( 'ZZRAYNP'                                )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZRAYNP' )
      RETURN
      END


      
