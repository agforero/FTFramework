C$Procedure DSKTOL ( DSK, tolerance umbrella )
 
      SUBROUTINE DSKTOL ( KEYWRD, DPVAL )
 
C$ Abstract
C
C     Umbrella routine for DSK tolerance and margin parameter access
C     routines.
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
C     MARGIN
C     NUMERIC
C     TOLERANCE
C     
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dsktol.inc'
      
      INTEGER               KEYWRD
      DOUBLE PRECISION      DPVAL
      
C$ Brief_I/O
C
C     Variable  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     KEYWRD     I   DSKGTL, DSKSTL
C     DPVAL     I-O  DSKGTL, DSKSTL
C
C$ Detailed_Input
C
C     See the entry points for descriptions of their input arguments.
C
C$ Detailed_Output
C
C     See the entry points for descriptions of their output arguments.
C
C$ Parameters
C
C     See the include file 
C
C        dsktol.inc 
C
C     for descriptions and values of the tolerance or margin parameters
C     accessed by the entry points of this routine, and of the keyword
C     parameters used to refer to them.
C
C$ Exceptions
C
C     1)  If this routine is called directly, the error
C         SPICE(BOGUSENTRY) is signaled.
C
C     See the entry points for descriptions of exceptions specific to
C     those entry points.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The routines in this package serve to centralize numeric
C     tolerance and margin values used by the DSK subsystem. The
C     subsystem retrieves values from this package to use at run time.
C
C     The entry points of this routine are
C
C        DSKGTL {DSK, get tolerance value}
C        DSKSTL {DSK, set tolerance value}
C 
C     To minimize run time overhead, the "keywords" used by these
C     routines to identify parameters are actually integer codes.
C
C     SPICE users may override certain values maintained by this
C     package; others values are fixed.
C     
C     One use of this system would be to disable the "greedy"
C     algorithms used to prevent "false miss" ray-surface intercepts
C     results. Setting to zero the margins used for this purpose would
C     accomplish this.
C
C     It is recommended that any change to the tolerance values made at
C     run time be programmed only by SPICE experts.
C
C$ Examples
C
C     See the entry points.
C
C$ Restrictions
C
C     1) The default settings used by the DSK subsystem should
C        be overridden only by expert SPICE users.
C
C     2) The entry points of this routine do not check the 
C        validity of new parameter values supplied by the
C        calling application.
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
C-    SPICELIB Version 1.0.0, 01-AUG-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     dsk tolerance and margin umbrella
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               NDSPAR
      PARAMETER           ( NDSPAR = 6 )

      INTEGER               NAMLEN
      PARAMETER           ( NAMLEN = 6 )

C
C     Local variables
C      
      CHARACTER*(NAMLEN)    NAMES  ( NDSPAR )

      DOUBLE PRECISION      DPPARS ( NDSPAR )

      LOGICAL               ISFIXD ( NDSPAR )

C
C     Saved values
C
      SAVE                  DPPARS
      SAVE                  ISFIXD
      SAVE                  NAMES

C
C     Initial values
C
C     Caution: this array must be declared in 
C     the order defined by the keywords declared
C     in dsktol.inc.
C
      DATA                  DPPARS /
     .                      XFRACT, SGREED, SGPADM,
     .                      PTMEMM, ANGMRG, LONALI       
     .                             /

      DATA                  ISFIXD /
     .                      .FALSE., .FALSE., .FALSE.,
     .                      .FALSE., .TRUE.,  .TRUE.
     .                             /

      DATA                  NAMES  /
     .                      'XFRACT', 'SGREED', 'SGPADM',
     .                      'PTMEMM', 'ANGMRG', 'LONALI'
     .                             /

      CALL CHKIN  ( 'DSKTOL' )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)' )
      CALL CHKOUT ( 'DSKTOL' )
      RETURN



C$Procedure DSKGTL ( DSK, get tolerance )
 
      ENTRY DSKGTL ( KEYWRD, DPVAL )
 
C$ Abstract
C
C     Retrieve the value of a specified DSK tolerance or margin
C     parameter.
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
C     MARGIN
C     NUMERIC
C     TOLERANCE
C     
C$ Declarations
C      
C     INTEGER               KEYWRD
C     DOUBLE PRECISION      DPVAL
C     
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     KEYWRD     I   Code specifying parameter to retrieve.
C     DPVAL      O   Value of parameter.
C
C$ Detailed_Input
C
C     KEYWRD     is an integer code specifying the parameter to
C                retrieve. See the include file dsktol.inc for
C                a description of the possible keywords.
C
C$ Detailed_Output
C
C     DPVAL      is the value of the parameter specified by KEYWRD.
C
C$ Parameters
C
C     See the include file 
C
C        dsktol.inc 
C
C     for descriptions and values of the tolerance or margin parameters
C     accessed by this routine, and of the keyword parameters used to
C     refer to them.
C
C$ Exceptions
C
C     1)  If the input keyword is not recognized, the error
C         SPICE(INDEXOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The routines in this package serve to centralize numeric
C     tolerance and margin values used by the DSK subsystem.
C     The subsystem retrieves values from this package to use
C     at run time. 
C
C     The entry points of this routine are
C
C        DSKGTL {DSK, get tolerance value}
C        DSKSTL {DSK, set tolerance value}
C 
C     To minimize run time overhead, the "keywords" used by
C     these routines to identify parameters are actually 
C     integer codes.
C
C     SPICE users may override certain values maintained by
C     this package; others values are fixed. It is recommended
C     that any change to the tolerance values made at run 
C     time be performed only by expert SPICE users. 
C
C$ Examples
C
C     1) Obtain the DSK type 2 plate expansion fraction:
C
C           DOUBLE PRECISION      F
C
C             ...
C
C           CALL DSKGTL ( KEYXFR, F )
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
C-    SPICELIB Version 1.0.0, 01-AUG-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     retrieve dsk tolerance or margin parameters
C
C-&

C
C     Use discovery check-in.
C     
      IF (  ( KEYWRD .LT. 1 ) .OR. ( KEYWRD .GT. NDSPAR )  ) THEN

         CALL CHKIN  ( 'DSKGTL'                                     )
         CALL SETMSG ( 'Valid keyword range is 1:#; keyword was #.' )
         CALL ERRINT ( '#', NDSPAR                                  )
         CALL ERRINT ( '#', KEYWRD                                  )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'                     )
         CALL CHKOUT ( 'DSKGTL'                                     )
         RETURN

      END IF

      DPVAL = DPPARS( KEYWRD )
      
      RETURN





C$Procedure DSKSTL ( DSK, set tolerance )
 
      ENTRY DSKSTL ( KEYWRD, DPVAL )
 
C$ Abstract
C
C     Set the value of a specified DSK tolerance or margin parameter.
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
C     MARGIN
C     NUMERIC
C     TOLERANCE
C     
C$ Declarations
C      
C     INTEGER               KEYWRD
C     DOUBLE PRECISION      DPVAL
C     
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     KEYWRD     I   Code specifying parameter to modify.
C     DPVAL      I   New value of parameter.
C
C$ Detailed_Input
C
C     KEYWRD     is an integer code specifying the parameter to
C                set. See the include file dsktol.inc for
C                a description of the possible keywords.
C
C     DPVAL      is the new value of the parameter specified by KEYWRD.
C                This value will be retrieved by future calls to
C                to DSKGTL made to retrieve the specified parameter.
C
C                <<< Use extreme caution. This routine performs no
C                checks on DPVAL. >>>
C
C$ Detailed_Output
C
C     None. This routine operates by side effects.
C
C$ Parameters
C
C     See the include file 
C
C        dsktol.inc 
C
C     for descriptions and values of the tolerance or margin parameters
C     accessed by this routine, and of the keyword parameters used to
C     refer to them.
C
C$ Exceptions
C
C     1)  If the input keyword is not recognized, the error
C         SPICE(INDEXOUTOFRANGE) is signaled.
C
C     2)  If an attempt is made to modify a fixed parameter,
C         the error SPICE(IMMUTABLEVALUE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The routines in this package serve to centralize numeric
C     tolerance and margin values used by the DSK subsystem.
C     The subsystem retrieves values from this package to use
C     at run time. 
C
C     The entry points of this routine are
C
C        DSKGTL {DSK, get tolerance value}
C        DSKSTL {DSK, set tolerance value}
C 
C     To minimize run time overhead, the "keywords" used by
C     these routines to identify parameters are actually 
C     integer codes.
C
C     SPICE users may override certain values maintained by
C     this package; others values are fixed. It is recommended
C     that any change to the tolerance values made at run 
C     time be performed only by expert SPICE users. 
C
C$ Examples
C
C     1) Change the DSK type 2 plate expansion fraction to
C        1.e-8:
C
C            CALL DSKSTL ( KEYXFR, 1.D-8 )
C
C$ Restrictions
C
C     1) The default settings used by the DSK subsystem should
C        be overridden only by expert SPICE users.
C
C     2) The entry points of this routine do not check the 
C        validity of new parameter values supplied by the
C        calling application.
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
C-    SPICELIB Version 1.0.0, 01-AUG-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     set dsk tolerance or margin parameters
C
C-&
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'DSKSTL' )


      IF (  ( KEYWRD .LT. 1 ) .OR. ( KEYWRD .GT. NDSPAR )  ) THEN

         CALL SETMSG ( 'Valid keyword range is 1:#; keyword was #.' )
         CALL ERRINT ( '#', NDSPAR                                  )
         CALL ERRINT ( '#', KEYWRD                                  )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'                     )
         CALL CHKOUT ( 'DSKSTL'                                     )
         RETURN

      END IF

      IF ( ISFIXD(KEYWRD) ) THEN
         
         CALL SETMSG ( 'The parameter # cannot be modified.' )
         CALL ERRCH  ( '#', NAMES(KEYWRD)                    )
         CALL SIGERR ( 'SPICE(IMMUTABLEVALUE)'               )
         CALL CHKOUT ( 'DSKSTL'                              )
         RETURN
         
      END IF
      
C
C     We have a valid parameter index. We don't check
C     the new parameter value; the user presumably knows 
C     the reason for change.
C
      DPPARS(KEYWRD) = DPVAL

      CALL CHKOUT ( 'DSKSTL' )
      RETURN
      END
