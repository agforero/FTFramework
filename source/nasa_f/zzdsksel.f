C$Procedure ZZDSKSEL ( DSK, segment selection callback umbrella )
 
      LOGICAL FUNCTION ZZDSKSEL ( SURFID, NSURF,  SRFLST,
     .                            BODYID, DCLASS, CORSYS, CORPAR,
     .                            COR1,   COR2,   FRAMID, POS,  
     .                            ET,     HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     This is the umbrella routine for DSK segment selection functions
C     that are passed by ZZDSKSNS. Entry points for initializing
C     segment selection functions are included as well.
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
C     SEARCH
C     SEGMENT
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsk.inc'
      INCLUDE 'dskdsc.inc'

      INTEGER               SURFID
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      INTEGER               BODYID
      INTEGER               DCLASS
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      COR1
      DOUBLE PRECISION      COR2
      INTEGER               FRAMID
      DOUBLE PRECISION      POS    ( 3 )
      DOUBLE PRECISION      ET
      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      DOUBLE PRECISION      DSKDSC ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     SURFID     I   ZZDSKMSC, ZZDSKSRC
C     NSURF      I   ZZDSKSIT
C     SRFLST     I   ZZDSKSIT
C     BODYID     I   ZZDSKSBD
C     DCLASS     I   ZZDSKSRC
C     CORSYS     I   ZZDSKMSC
C     CORPAR     I   ZZDSKMSC
C     COR1       I   ZZDSKMSC, ZZDSKUSC
C     COR2       I   ZZDSKMSC, ZZDSKUSC
C     FRAMID     I   ZZDSKMSC, ZZDSKSRC
C     POS        I   ZZDSKSRC
C     ET         I   ZZDSKMSC, ZZDSKSIT, ZZDSKSRC, ZZDSKUSC
C     HANDLE     I   ZZDSKBDC, ZZDSKMMC, ZZDSKUMC, ZZDSKMRC, ZZDSKCIT
C     DLADSC     I   ZZDSKBDC, ZZDSKMMC, ZZDSKUMC, ZZDSKMRC, ZZDSKCIT
C     DSKDSC     I   ZZDSKBDC, ZZDSKMMC, ZZDSKUMC, ZZDSKMRC, ZZDSKCIT
C
C$ Detailed_Input
C
C     See the entry points for descriptions of their inputs.
C                
C$ Detailed_Output
C
C     See the entry points for descriptions of their outputs.
C                 
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If this routine is called directly the error
C         SPICE(BOGUSENTRY) is signaled.
C
C     See the entry points for descriptions of exceptions
C     applicable to those entry points.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem.
C
C     Frame kernels may be required in order to perform transformations
C     from DSK segment frames to the output frame. In some cases,
C     additional kernels such as CK kernels and SCLK kernels could be
C     required to support the offset vector lookup.
C     
C$ Particulars
C
C     The entry points of this routine are used to search for 
C     segments matching specified criteria among the set of DSK
C     segment descriptors stored in the ZZDSKBSR subsystem. 
C           
C     For each kind of search, there is an initialization entry point
C     used to store search parameters, and there is a corresponding
C     segment comparison entry point. The comparison entry points are
C     callback routines that are passed to ZZDSKSNS.
C
C     The supported search types and associated entry points are 
C     listed below.
C
C         Match body ID only
C         ------------------
C
C           Initialization:       ZZDSKSBD
C           Segment comparison:   ZZDSKBDC
C
C
C         Match nothing:
C     
C           Segment comparison:   ZZDSKNOT
C
C           This routine is used to force the ZZDSKSBF package to
C           build a complete segment list for a specified body.
C
C
C     The remaining entry points are included to support the SMAP
C     Alpha DSK Toolkit. They may be useful for DSK type 4 in a
C     future SPICELIB version.
C
C
C         Match body ID, surface list, and time
C         -------------------------------------
C
C           Initialization:       ZZDSKSIT
C           Segment comparison:   ZZDSKCIT
C
C
C         Match body ID, time, and coordinates
C         ------------------------------------
C
C           Initialization:       ZZDSKUSC
C           Segment comparison:   ZZDSKUMC
C
C
C         Match body ID, surface ID, frame ID, coordinate system,
C         coordinate parameters, time, and coordinates
C         ------------------------------------------------------
C
C           Initialization:       ZZDSKMSC
C           Segment comparison:   ZZDSKMMC
C
C
C         Match body ID, surface ID, frame ID, data class, time,
C         and whether coverage includes a specified position vector
C         ---------------------------------------------------------
C
C           Initialization:       ZZDSKSRC
C           Segment comparison:   ZZDSKMRC
C
C
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C        Added entry point ZZDSKNOT.
C
C        05-FEB-2016 (NJB) 
C
C
C        29-APR-2015 (NJB)
C
C           Updated description of ZZDSKUMC to specify that the body,
C           not the surface, is matched.
C
C        10-OCT-2014 (NJB)
C
C           Added time-independent body search routines
C
C              ZZDSKSBD, ZZDSKCBD
C
C        29-SEP-2014 (NJB)
C
C           Updated ZZDSKUSC, ZZDSKUMC to use body IDs.
C
C        18-SEP-2014 (NJB)
C
C           Updated ZZDSKCIT, ZZDSKSIT to use list of surfaces as well
C           as a body ID. Updated umbrella ZZDSKSEL to support new
C           arguments.
C
C        15-SEP-2014 (NJB)
C
C           Utilities for selecting segments via 
C           calls to ZZDSKSNS
C
C        13-MAY-2010 (NJB)
C  
C           Developed to support SMAP Alpha DSK Toolkit.
C
C-&
 
C$ Index_Entries
C
C     umbrella for dsk segment selection callback 
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      TWOPI

      INTEGER               BSRCHI
      INTEGER               TOUCHI

      LOGICAL               ZZDSKBDC
      LOGICAL               ZZDSKCIT
      LOGICAL               ZZDSKMMC
      LOGICAL               ZZDSKMRC
      LOGICAL               ZZDSKMSC
      LOGICAL               ZZDSKNOT
      LOGICAL               ZZDSKSBD
      LOGICAL               ZZDSKSIT
      LOGICAL               ZZDSKSRC
      LOGICAL               ZZDSKUMC
      LOGICAL               ZZDSKUSC

C
C     Local parameters
C
      DOUBLE PRECISION      MARGIN
      PARAMETER           ( MARGIN = 1.D-12 )

C
C     Local variables
C     
      DOUBLE PRECISION      ALT
      DOUBLE PRECISION      CO1MAX
      DOUBLE PRECISION      CO1MIN
      DOUBLE PRECISION      CO2MAX
      DOUBLE PRECISION      CO2MIN
      DOUBLE PRECISION      F
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      LOCCOR ( 1 )
      DOUBLE PRECISION      LOCPOS ( 3 )
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      PI2
      DOUBLE PRECISION      R
      DOUBLE PRECISION      RE
      DOUBLE PRECISION      RMAT   ( 3, 3 )
      DOUBLE PRECISION      SAVCO1
      DOUBLE PRECISION      SAVCO2
      DOUBLE PRECISION      SAVET
      DOUBLE PRECISION      SAVPAR ( NSYPAR )
      DOUBLE PRECISION      SAVPOS ( 3 )
      DOUBLE PRECISION      SCALE

      INTEGER               I
      INTEGER               SAVBID
      INTEGER               SAVCLS
      INTEGER               SAVFID
      INTEGER               SAVNSF
      INTEGER               SAVSRF ( MAXSRF )
      INTEGER               SAVSYS
      INTEGER               SAVTRG
      INTEGER               SEGFID
      INTEGER               SEGSYS
      INTEGER               SURF

      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  PI2
      SAVE                  SAVBID
      SAVE                  SAVCLS
      SAVE                  SAVCO1
      SAVE                  SAVCO2
      SAVE                  SAVET
      SAVE                  SAVFID
      SAVE                  SAVNSF
      SAVE                  SAVPAR
      SAVE                  SAVPOS
      SAVE                  SAVSRF
      SAVE                  SAVSYS
      SAVE                  SAVTRG

C
C     Initial values
C
      DATA                  FIRST / .TRUE. /

C
C     Set return value.
C
      ZZDSKSEL = .FALSE.

C
C     The following no-op calls are used to suppress compiler 
C     warnings.
C         
      I = TOUCHI ( HANDLE )
      I = TOUCHI ( CORSYS )
      I = TOUCHI ( DLADSC )
      I = TOUCHI ( FRAMID )

      CALL CHKIN  ( 'ZZDSKSEL' )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)' )
      CALL CHKOUT ( 'ZZDSKSEL' )
      RETURN



C$Procedure ZZDSKSBD ( DSK, set body ID )
 
      ENTRY ZZDSKSBD ( BODYID )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Set body ID only for segment matching.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               BODYID
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body.
C
C     The function always returns the value .FALSE.
C
C$ Detailed_Input
C
C     BODYID     is the body ID against which DSK segments' body IDs
C                are to be compared.
C
C$ Detailed_Output
C
C     The function always returns the value .FALSE. An output
C     is necessary because this routine is an entry point of
C     a function. The returned value is meaningless.
C
C     This function operates by side effects.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is meant to prepare for DSK segment searches
C     using the DSK segment comparison callback 
C
C        ZZDSKBDC
C
C     That routine will indicate a match for segments having
C     body ID equal to the input BODYID.
C
C$ Examples
C
C     See usage in ZZDSKSBA.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     set body id for body id selection callback 
C
C-&

C
C     Set the function's return value. This value has no
C     particular meaning and is not meant to be used by
C     the caller. This setting suppresses compiler warnings.
C
      ZZDSKSBD = .FALSE.
C
C     Save input value.
C
      SAVBID = BODYID

      RETURN
      



C$Procedure ZZDSKBDC ( DSK, check segment's body ID )
 
      ENTRY ZZDSKBDC ( HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Check DSK segment for body ID match.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK containing segment to be checked.
C     DLADSC     I   DLA descriptor of segment to be checked.
C     DSKDSC     I   DLA descriptor of segment to be checked.
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKSBD.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of the DSK containing a segment to
C                be checked.
C
C     DLADSC, 
C     DSKDSC     are, respectively, the DLA descriptor and DSK
C                descriptor of a segment to be checked.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKSBD.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is used to search for DSK segments that
C     have body ID codes matching that set via the latest
C     call to 
C
C        ZZDSKSBD
C
C     These searches are conducted by ZZDSKSNS. This function
C     is passed to that routine as a callback argument.
C
C$ Examples
C
C     See usage in ZZDSKSBA.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk body id segment selection callback 
C
C-&

C
C     The body ID is the only DSK segment attribute that must match.
C
      ZZDSKBDC = NINT( DSKDSC(CTRIDX) ) .EQ. SAVBID

      RETURN
      



C$Procedure ZZDSKNOT ( DSK, match nothing )
 
      ENTRY ZZDSKNOT ( HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Return .FALSE. for any segment.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK containing segment to be checked.
C     DLADSC     I   DLA descriptor of segment to be checked.
C     DSKDSC     I   DLA descriptor of segment to be checked.
C
C     The function returns .FALSE. for all segments.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of the DSK containing a segment to
C                be checked.
C
C     DLADSC, 
C     DSKDSC     are, respectively, the DLA descriptor and DSK
C                descriptor of a segment to be checked.
C
C$ Detailed_Output
C
C     The function returns .FALSE. for all segments.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This function is passed to ZZDSKSNS as a callback argument.
C
C     This routine makes it possible to force ZZDSKSNS to build
C     a complete segment list for a specified body in linear time.
C
C     By indicating "no match" for any segment, this routine 
C     prevents a search for matching segments from terminating 
C     until all available files have been searched.
C
C$ Examples
C
C     See usage in ZZDSKBBL.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk no-match callback 
C
C-&

C
C     Whatever's in this segment, it doesn't match.
C
      ZZDSKNOT = .FALSE.

      RETURN
      




C$Procedure ZZDSKSIT ( DSK, set body ID, surfaces, and time )
 
      ENTRY ZZDSKSIT ( BODYID, NSURF, SRFLST, ET )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Set body ID code, surface list, and time for segment matching.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               BODYID
C     INTEGER               NSURF
C     INTEGER               SRFLST ( * )
C     DOUBLE PRECISION      ET
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body.
C     NSURF      I   Surface ID count.
C     SRFLST     I   Surface ID list.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C
C     The function always returns the value .FALSE.
C
C$ Detailed_Input
C
C     BODYID     is the body ID against which DSK segments' body IDs
C                are to be compared.
C
C     NSURF,
C     SRFLST     are, respectively, a surface ID count and a list of
C                surface IDs. Matching segments must have a surface
C                belonging to this list. 
C
C                If the list is empty, all surfaces are considered
C                to be matches.
C
C     ET         is an epoch, expressed as seconds past J2000 TDB.
C                All matching segments must contain ET in their
C                time coverage intervals.
C
C$ Detailed_Output
C
C     The function always returns the value .FALSE. An output
C     is necessary because this routine is an entry point of
C     a function. The returned value is meaningless.
C
C     This function operates by side effects.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is meant to prepare for DSK segment searches
C     using the DSK segment comparison callback 
C
C        ZZDSKCIT
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     set body id, surfaces, and time for corresponding callback 
C
C-&


C
C     Set return value.
C
      ZZDSKSIT = .FALSE.
C
C     Save input values.
C
      SAVBID = BODYID
      SAVET  = ET

      IF ( NSURF .GT. MAXSRF ) THEN

         CALL CHKIN ( 'ZZDSKSIT' )

         CALL SETMSG ( 'Maximum allowed surface ID count is #; '
     .   //            'input count was #.'                     )
         CALL ERRINT ( '#', MAXSRF                              )
         CALL ERRINT ( '#', NSURF                               )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZDSKSIT' )
         RETURN

      END IF

      SAVNSF = NSURF

      DO I = 1, NSURF
         SAVSRF(I) = SRFLST(I)
      END DO

      CALL SHELLI ( NSURF, SAVSRF )

      RETURN






C$Procedure ZZDSKCIT ( DSK, check body ID, surface, and time )
 
      ENTRY ZZDSKCIT ( HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Check DSK segment for body ID, surface, and time coverage match.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK containing segment to be checked.
C     DLADSC     I   DLA descriptor of segment to be checked.
C     DSKDSC     I   DLA descriptor of segment to be checked.
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKSIT.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of the DSK containing a segment to
C                be checked.
C
C     DLADSC, 
C     DSKDSC     are, respectively, the DLA descriptor and DSK
C                descriptor of a segment to be checked.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKSIT.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is used to search for DSK segments that
C     have body ID codes and time coverage intervals matching the
C     values set via the latest call to 
C
C        ZZDSKSIT
C
C     These searches are conducted by ZZDSKSNS. This function
C     is passed to that routine as a callback argument.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk body id and time segment selection callback 
C
C-&


C
C     Currently, we don't need to access the DSK file in order
C     to make this test, but the required inputs are available.
C
      ZZDSKCIT = .FALSE.

      IF ( SAVBID  .EQ. NINT( DSKDSC(CTRIDX) )  ) THEN

         IF (      ( SAVET .GE. DSKDSC(BTMIDX) )
     .       .AND. ( SAVET .LE. DSKDSC(ETMIDX) )  ) THEN

            IF ( SAVNSF .LT. 1 ) THEN
C
C              There are no surface ID constraints; we have
C              a match.
C
               ZZDSKCIT = .TRUE.

            ELSE
C
C              We have a match if and only if the surface ID of this
C              segment is on the list of allowed surface IDs.
C
               SURF     = NINT( DSKDSC(SRFIDX) )

               ZZDSKCIT = BSRCHI( SURF, SAVNSF, SAVSRF ) .GT. 0

            END IF

         END IF

      END IF

      RETURN






C$Procedure ZZDSKUSC ( DSK, set body ID, coordinates, and time )
 
      ENTRY ZZDSKUSC ( BODYID, ET, COR1, COR2 )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Set body ID, time, and coordinates for "raw" unchecked search.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               BODYID
C     DOUBLE PRECISION      ET
C     DOUBLE PRECISION      COR1
C     DOUBLE PRECISION      COR2
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     COR1       I   First coordinate.
C     COR2       I   Second coordinate.
C
C     The function always returns the value .FALSE.
C
C$ Detailed_Input
C
C     BODYID     is the body ID against which DSK segments' body IDs
C                are to be compared.
C
C     ET         is an epoch, expressed as seconds past J2000 TDB.
C                All matching segments must contain ET in their
C                time coverage intervals.
C
C     COR1,
C     COR2       are two coordinates that must be included in the
C                spatial region covered by a matching segment.
C
C                These are "domain" coordinates for class 1
C                segments.
C
C$ Detailed_Output
C
C     The function always returns the value .FALSE. An output
C     is necessary because this routine is an entry point of
C     a function. The returned value is meaningless.
C
C     This function operates by side effects.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is meant to prepare for DSK segment searches
C     using the DSK segment comparison callback 
C
C        ZZDSKUMC
C
C     This routine is applicable only for cases where the loaded
C     segments have a common reference frame and coordinate system,
C     and where the segment class is 1 ("single-valued function").
C     
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     set body id, time, and coordinates for callback 
C
C-&


      ZZDSKUSC = .FALSE.

      SAVBID = BODYID
      SAVET  = ET
      SAVCO1 = COR1
      SAVCO2 = COR2

      RETURN





C$Procedure ZZDSKUMC ( DSK, check body ID, coordinates, and time )
 
      ENTRY ZZDSKUMC ( HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Perform an unchecked segment match using body ID, time,
C     and coordinates. 
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK containing segment to be checked.
C     DLADSC     I   DLA descriptor of segment to be checked.
C     DSKDSC     I   DLA descriptor of segment to be checked.
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKUSC.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of the DSK containing a segment to
C                be checked.
C
C     DLADSC, 
C     DSKDSC     are, respectively, the DLA descriptor and DSK
C                descriptor of a segment to be checked.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID, time, and coordinates set via ZZDSKUSC.
C
C     A time and coordinate match means that the input values are
C     contained in the respective time interval and spatial coverage
C     region of the segment.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is used to search for DSK segments that
C     have body ID codes and time coverage intervals matching the
C     values set via the latest call to 
C
C        ZZDSKUSC.
C
C     These searches are conducted by ZZDSKSNS. This function
C     is passed to that routine as a callback argument.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk body id, time, and coordinate selection callback 
C
C-&



C
C     We don't have a match to begin with.
C
      ZZDSKUMC = .FALSE.

      IF ( FIRST ) THEN
         PI2   = TWOPI()
         FIRST = .FALSE.
      END IF

C
C     Check the body ID first.
C
      IF ( SAVBID  .EQ. NINT( DSKDSC(CTRIDX) )  ) THEN
C
C        The body ID matches; check the time coverage.
C
         IF (      ( SAVET .GE. DSKDSC(BTMIDX) )
     .       .AND. ( SAVET .LE. DSKDSC(ETMIDX) )  ) THEN
C
C           Check the coordinates. Note that we don't 
C           know whether the frame and coordinates are
C           reasonable. It's up to the user to ensure this.
C
C           We do need to know whether the first coordinate
C           is a longitude, since we may need to adjust it to
C           get it into the range of the segment.
C   
            SEGSYS = NINT( DSKDSC(SYSIDX) )
C
C           Make a local copy of the first coordinate.
C
            LOCCOR(1) = SAVCO1

            IF (      ( SEGSYS .EQ. LATSYS )
     .           .OR. ( SEGSYS .EQ. PDTSYS ) ) THEN
C
C              Adjust segment bounds using a small margin.
C
               CO1MIN = DSKDSC( MN1IDX ) - MARGIN
               CO1MAX = DSKDSC( MX1IDX ) + MARGIN
               CO2MIN = DSKDSC( MN2IDX ) - MARGIN
               CO2MAX = DSKDSC( MX2IDX ) + MARGIN
C
C              Move longitude into range.
C
               IF ( LOCCOR(1) .LT. CO1MIN ) THEN

                  LOCCOR(1) = LOCCOR(1) + PI2

               ELSE IF ( LOCCOR(1) .GT. CO1MAX ) THEN

                  LOCCOR(1) = LOCCOR(1) - PI2

               END IF

            ELSE
               
               SCALE  = ( 1.D0 + MARGIN )

               CO1MIN = DSKDSC( MN1IDX ) - SCALE * ABS( DSKDSC(MN1IDX) )
               CO1MAX = DSKDSC( MX1IDX ) + SCALE * ABS( DSKDSC(MX1IDX) )
               CO2MIN = DSKDSC( MN2IDX ) - SCALE * ABS( DSKDSC(MN2IDX) )
               CO2MAX = DSKDSC( MX2IDX ) + SCALE * ABS( DSKDSC(MX2IDX) )

            END IF

C
C           Check the first coordinate against the segment's
C           corresponding coverage range.
C
            IF (     ( LOCCOR(1) .LT. CO1MIN )
     .          .OR. ( LOCCOR(1) .GT. CO1MAX ) ) THEN
C
C              The first input coordinate is not covered by this
C              segment. 
C
               RETURN

            END IF

C
C           Compare the second coordinate against the segment's
C           corresponding coverage range.
C
            IF (     ( SAVCO2 .LT. CO2MIN )
     .          .OR. ( SAVCO2 .GT. CO2MAX ) ) THEN
C
C              The second input coordinate is not covered by this
C              segment. 
C
               RETURN

            END IF

C
C           At this point we have a match.
C
            ZZDSKUMC = .TRUE.

         END IF
C
C        This is the end of the time check block.
C 
      END IF
C
C     This is the end of the surface ID check block.
C 
      RETURN
      




C$Procedure ZZDSKMSC ( DSK, setup matched attribute search )
 
      ENTRY ZZDSKMSC ( BODYID, SURFID, FRAMID, CORSYS, 
     .                 CORPAR, ET,     COR1,   COR2   )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Set inputs for a matched attribute search, also known as a
C     "checked" search.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               BODYID
C     INTEGER               SURFID
C     INTEGER               FRAMID
C     INTEGER               CORSYS
C     DOUBLE PRECISION      CORPAR ( * )
C     DOUBLE PRECISION      ET
C     DOUBLE PRECISION      COR1
C     DOUBLE PRECISION      COR2
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body.
C     SURFID     I   ID code of target surface.
C     FRAMID     I   ID code of body-fixed reference frame.
C     CORSYS     I   Coordinate system code.
C     CORPAR     I   Coordinate system parameters.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     COR1       I   First coordinate.
C     COR2       I   Second coordinate.
C
C     The function always returns the value .FALSE.
C
C$ Detailed_Input
C
C     BODYID     is the body ID against which DSK segments' body IDs
C                are to be compared.
C
C     SURFID     is the surface ID against which DSK segments' surface
C                IDs are to be compared.
C
C     FRAMID     is the reference frame ID against which DSK segments'
C                frame IDs are to be compared.
C
C     CORSYS     is the coordinate system code against which DSK
C                segments' coordinate systems are to be compared.
C
C     CORPAR     is an array containing the coordinate system
C                parameters against which DSK segments' coordinate
C                system parameters are to be compared.
C
C                A small margin is used for parameter comparison.
C
C
C     ET         is an epoch, expressed as seconds past J2000 TDB.
C                All matching segments must contain ET in their
C                time coverage intervals.
C
C     COR1,
C     COR2       are two coordinates that must be included in the
C                spatial region covered by a matching segment.
C
C                These are "domain" coordinates for class 1
C                segments.
C
C$ Detailed_Output
C
C     The function always returns the value .FALSE. An output
C     is necessary because this routine is an entry point of
C     a function. The returned value is meaningless.
C
C     This function operates by side effects.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is meant to prepare for DSK segment searches
C     using the DSK segment comparison callback 
C
C        ZZDSKMMC
C
C     This routine supports "matched attribute" searches. 
C     Segment attributes are checked explicitly.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     set inputs for matched attribute callback 
C
C-&


      ZZDSKMSC = .FALSE.

      SAVTRG = BODYID
      SAVBID = SURFID
      SAVFID = FRAMID
      SAVSYS = CORSYS
      
      CALL MOVED ( CORPAR, NSYPAR, SAVPAR )

      SAVET  = ET
      SAVCO1 = COR1
      SAVCO2 = COR2

      RETURN






C$Procedure ZZDSKMMC ( DSK, matched segment attribute comparison )
 
      ENTRY ZZDSKMMC ( HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C
C     Perform a matched attribute comparison using target ID, surface
C     ID, frame ID, coordinate system, coordinate parameters, time, and
C     coordinates.
C
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK containing segment to be checked.
C     DLADSC     I   DLA descriptor of segment to be checked.
C     DSKDSC     I   DLA descriptor of segment to be checked.
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKMSC.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of the DSK containing a segment to
C                be checked.
C
C     DLADSC, 
C     DSKDSC     are, respectively, the DLA descriptor and DSK
C                descriptor of a segment to be checked.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the 
C
C        body ID
C        surface ID
C        frame ID
C        coordinate system
C        coordinate parameters
C        time
C        coordinates
C
C     set via ZZDSKMSC.
C
C     A time and coordinate match means that the input values are
C     contained in the respective time interval and spatial coverage
C     region of the segment.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is used to search for DSK segments that have
C     attributes matching the values set via the latest call to
C
C        ZZDSKMSC
C
C     This is a "matched attribute" search, also known as a
C     "checked" search. See Detailed_Output above.
C
C     These searches are conducted by ZZDSKSNS. This function
C     is passed to that routine as a callback argument.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk matched segment attribute callback 
C
C-&

C
C     We don't have a match to begin with.
C
      ZZDSKMMC = .FALSE.

      IF ( FIRST ) THEN
         PI2   = TWOPI()
         FIRST = .FALSE.
      END IF

C
C     Check the target ID.
C
      IF ( SAVTRG  .NE. NINT( DSKDSC(CTRIDX) )  ) THEN
         RETURN
      END IF

C
C     Check the surface ID.
C
      IF ( SAVBID  .NE. NINT( DSKDSC(SRFIDX) )  ) THEN
         RETURN
      END IF

C
C     Check the frame ID.
C
      IF ( SAVFID  .NE. NINT( DSKDSC(FRMIDX) )  ) THEN
         RETURN
      END IF

C
C     Check the coordinate system.
C
      IF ( SAVSYS  .NE. NINT( DSKDSC(SYSIDX) )  ) THEN
         RETURN
      END IF

C
C     If the system is planetodetic, check the reference
C     ellipsoid parameters.
C
      IF ( SAVSYS .EQ. PDTSYS ) THEN


         IF (  ABS( SAVPAR(1) - DSKDSC(PARIDX)   ) .GT. MARGIN ) THEN
            RETURN
         END IF

         IF (  ABS( SAVPAR(2) - DSKDSC(PARIDX+1) ) .GT. MARGIN ) THEN
            RETURN
         END IF

      END IF

C
C     The segment attributes match; check the time coverage.
C
      IF (      ( SAVET .GE. DSKDSC(BTMIDX) )
     .    .AND. ( SAVET .LE. DSKDSC(ETMIDX) )  ) THEN
C
C        Check the coordinates.  
C
         SEGSYS = NINT( DSKDSC(SYSIDX) )
C
C        Make a local copy of the first coordinate.
C
         LOCCOR(1) = SAVCO1

         IF (      ( SEGSYS .EQ. LATSYS )
     .        .OR. ( SEGSYS .EQ. PDTSYS ) ) THEN
C
C           Adjust segment bounds using a small margin.
C
            CO1MIN = DSKDSC( MN1IDX ) - MARGIN
            CO1MAX = DSKDSC( MX1IDX ) + MARGIN
            CO2MIN = DSKDSC( MN2IDX ) - MARGIN
            CO2MAX = DSKDSC( MX2IDX ) + MARGIN

C
C           Move longitude into range.
C
            IF ( LOCCOR(1) .LT. CO1MIN ) THEN

               LOCCOR(1) = LOCCOR(1) + PI2

            ELSE IF ( LOCCOR(1) .GT. CO1MAX ) THEN

               LOCCOR(1) = LOCCOR(1) - PI2

            END IF

         ELSE
               
            SCALE  = ( 1.D0 + MARGIN )

            CO1MIN = DSKDSC( MN1IDX ) - SCALE * ABS( DSKDSC(MN1IDX) )
            CO1MAX = DSKDSC( MX1IDX ) + SCALE * ABS( DSKDSC(MX1IDX) )
            CO2MIN = DSKDSC( MN2IDX ) - SCALE * ABS( DSKDSC(MN2IDX) )
            CO2MAX = DSKDSC( MX2IDX ) + SCALE * ABS( DSKDSC(MX2IDX) )

         END IF

C
C        Check the first coordinate against the segment's
C        corresponding coverage range.
C
         IF (     ( LOCCOR(1) .LT. CO1MIN )
     .       .OR. ( LOCCOR(1) .GT. CO1MAX ) ) THEN
C
C           The first input coordinate is not covered by this
C           segment. 
C
            RETURN

         END IF

C
C        Compare the second coordinate against the segment's
C        corresponding coverage range.
C
         IF (     ( SAVCO2 .LT. CO2MIN )
     .       .OR. ( SAVCO2 .GT. CO2MAX ) ) THEN
C
C           The second input coordinate is not covered by this
C           segment. 
C
            RETURN

         END IF

C
C        At this point we have a match.
C
         ZZDSKMMC = .TRUE.

      END IF

      RETURN
      






C$Procedure ZZDSKSRC ( DSK, setup rectangular coordinate search )
 
      ENTRY ZZDSKSRC ( SURFID, BODYID, DCLASS, ET, FRAMID, POS )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Set surface and target ID, data class, ET, and rectangular
C     coordinates for a coverage match search. Set the frame ID
C     as well; this is needed to define the position vector.
C     Matching segments need not use the same frame.
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               SURFID
C     INTEGER               BODYID
C     INTEGER               DCLASS
C     DOUBLE PRECISION      ET
C     INTEGER               FRAMID
C     DOUBLE PRECISION      POS    ( 3 )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     SURFID     I   ID code of target surface.
C     BODYID     I   ID code of target body.
C     DCLASS     I   Data class of segment.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     FRAMID     I   ID code of body-fixed reference frame.
C     POS        I   Cartesian position vector.
C
C     The function always returns the value .FALSE.
C
C$ Detailed_Input
C
C     SURFID     is the surface ID against which DSK segments' surface
C                IDs are to be compared.
C
C     BODYID     is the body ID against which DSK segments' body IDs
C                are to be compared.
C
C     DCLASS     is the data class against which DSK segments' data
C                classes are to be compared.
C
C     ET         is an epoch, expressed as seconds past J2000 TDB.
C                All matching segments must contain ET in their
C                time coverage intervals.
C
C    
C     FRAMID     is the reference frame ID of the frame in which
C                the input vector POS is expressed. FRAMID is
C                not used for comparison with segments' frame IDs.
C
C     POS        is a Cartesian position vector, expressed in
C                the frame designated by FRAMID. POS represents
C                an offset from the target body's center.
C
C                Segments are considered to "match" POS if the
C                latitude and longitude of POS are within the
C                spatial boundaries of those segments.
C
C$ Detailed_Output
C
C     The function always returns the value .FALSE. An output
C     is necessary because this routine is an entry point of
C     a function. The returned value is meaningless.
C
C     This function operates by side effects.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is meant to prepare for DSK segment searches
C     using the DSK segment comparison callback 
C
C        ZZDSKMRC
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     1)  This is a private routine. It is meant to be used only by the
C         DSK subsystem.
C
C     2)  This routine works only with segments that use latitudinal or
C         planetodetic coordinates.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     set inputs for dsk rectangular coordinate callback 
C
C-&

      ZZDSKSRC = .FALSE.

      SAVBID = SURFID
      SAVTRG = BODYID
      SAVCLS = DCLASS
      SAVET  = ET
      SAVFID = FRAMID
      CALL VEQU ( POS, SAVPOS )
      
      RETURN




C$Procedure ZZDSKMRC ( DSK, rectangular coordinate comparison )
 
      ENTRY ZZDSKMRC ( HANDLE, DLADSC, DSKDSC )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Perform a segment match using surface ID, target, data
C     class, time, and rectangular coordinates. 
C
C     A search using this function must be initialized
C     by a call to ZZDSKSRC.
C
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
C     BODY
C     DSK
C     SEARCH
C     SEGMENT
C
C$ Declarations
C
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK containing segment to be checked.
C     DLADSC     I   DLA descriptor of segment to be checked.
C     DSKDSC     I   DLA descriptor of segment to be checked.
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the body ID set via ZZDSKSRC.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of the DSK containing a segment to
C                be checked.
C
C     DLADSC, 
C     DSKDSC     are, respectively, the DLA descriptor and DSK
C                descriptor of a segment to be checked.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if and only if the specified segment
C     matches the 
C
C        surface ID
C        body ID
C        class
C        time
C        frame ID
C        coordinates
C
C     set via ZZDSKSRC.
C
C     A time and coordinate match means that the input values are
C     contained in the respective time interval and spatial coverage
C     region of the segment.
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
C     This routine does not use SPICE kernels directly, but it 
C     supports searches for DSK segments. 
C     
C$ Particulars
C
C     This routine is used to search for DSK segments that have
C     attributes matching the values set via the latest call to
C
C        ZZDSKSRC
C
C     These searches are conducted by ZZDSKSNS. This function
C     is passed to that routine as a callback argument.
C
C$ Examples
C
C     None.
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
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk match rectangular coordinate callback 
C
C-&



C
C     We don't have a match to begin with.
C
      ZZDSKMRC = .FALSE.

      IF ( FIRST ) THEN
         PI2   = TWOPI()
         FIRST = .FALSE.
      END IF


C     Check the surface ID.
C
      IF ( SAVBID  .NE. NINT( DSKDSC(SRFIDX) )  ) THEN
         RETURN
      END IF

C
C     Reject any segment whose target (center) ID doesn't match
C     the request ID.
C
      IF ( SAVTRG  .NE. NINT( DSKDSC(CTRIDX) )  ) THEN
         RETURN
      END IF

C
C     Check the time coverage.
C
      IF (      ( SAVET .LT. DSKDSC(BTMIDX) )
     .     .OR. ( SAVET .GT. DSKDSC(ETMIDX) )  ) THEN
         RETURN
      END IF

C
C     Check the data class.
C
      IF ( SAVCLS  .NE. NINT( DSKDSC(CLSIDX) )  ) THEN
         RETURN
      END IF


C
C     Check whether the position vector is covered by the segment.
C     In order to determine this, we need to transform the vector
C     into the frame of the segment, if the frames differ.
C
      SEGFID = NINT ( DSKDSC(FRMIDX) )

      IF ( SAVFID .EQ. SEGFID ) THEN
C
C        The request frame and segment frame match. Just copy
C        the saved vector.
C
         CALL VEQU ( SAVPOS, LOCPOS )

      ELSE
C
C        Transform the saved vector to the frame of the
C        segment. The transformation epoch is the saved
C        value of ET.
C
         CALL REFCHG ( SAVFID, SEGFID, SAVET, RMAT   )
         CALL MXV    ( RMAT,   SAVPOS,        LOCPOS ) 
 
      END IF


C     We do need to know whether the first coordinate is a longitude,
C     since we may need to adjust it to get it into the range of the
C     segment.
C   
      SEGSYS = NINT( DSKDSC(SYSIDX) )

      IF (      ( SEGSYS .EQ. LATSYS )
     .     .OR. ( SEGSYS .EQ. PDTSYS ) ) THEN
C
C        Find the latitude and longitude of the input point, expressed
C        in the coordinate system of this segment.
C
         IF ( SEGSYS .EQ. LATSYS ) THEN

            CALL RECLAT ( LOCPOS, R, LON, LAT )

         ELSE IF ( SEGSYS .EQ. PDTSYS ) THEN

            RE = DSKDSC(PARIDX  )
            F  = DSKDSC(PARIDX+1)

            CALL RECGEO ( LOCPOS, RE, F, LON, LAT, ALT )

         ELSE

            CALL CHKIN  ( 'ZZDSKMRC'                      )
            CALL SETMSG ( 'Backstop error (0): this code '
     .      //            'should be unreachable.'        )
            CALL SIGERR ( 'SPICE(BUG)'                    )
            CALL CHKOUT ( 'ZZDSKMRC'                      )
            RETURN
               
         END IF

C
C        Adjust segment bounds using a small margin.
C
         CO1MIN = DSKDSC( MN1IDX ) - MARGIN
         CO1MAX = DSKDSC( MX1IDX ) + MARGIN
         CO2MIN = DSKDSC( MN2IDX ) - MARGIN
         CO2MAX = DSKDSC( MX2IDX ) + MARGIN

C        Move longitude into range.
C
         IF ( LON .LT. CO1MIN ) THEN

            LON = LON + PI2

         ELSE IF ( LON .GT. CO1MAX ) THEN

            LON = LON - PI2

         END IF

      ELSE
               
         CALL CHKIN  ( 'ZZDSKMRC'                        )
         CALL SETMSG ( 'Only planetocentric and planeto' 
     .   //            'detic coordinates are supported ' 
     .   //            'by this entry point. Segment '
     .   //            'coordinate system was #.'        )
         CALL ERRINT ( '#',  SEGSYS                      )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED'              )
         CALL CHKOUT ( 'ZZDSKMRC'                        )
         RETURN
 
      END IF

C
C     Check the first coordinate against the segment's
C     corresponding coverage range.
C
      IF (     ( LON .LT. CO1MIN )
     .    .OR. ( LON .GT. CO1MAX ) ) THEN
C
C        The first input coordinate is not covered by this
C        segment. 
C
         RETURN

      END IF

C
C     Compare the second coordinate against the segment's
C     corresponding coverage range.
C
      IF (     ( LAT .LT. CO2MIN )
     .    .OR. ( LAT .GT. CO2MAX ) ) THEN
C
C        The second input coordinate is not covered by this
C        segment. 
C
         RETURN

      END IF

C
C     At this point we have a match.
C
      ZZDSKMRC = .TRUE.

      RETURN     
      END
