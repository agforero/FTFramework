C$Procedure      ZZCKE06 ( C-Kernel, evaluate, type 6 )
 
      SUBROUTINE ZZCKE06 ( RECORD, QSTATE, CLKOUT )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Evaluate a single data record from a type 6 CK segment. The
C     output is expressed as an interpolated unit quaternion and
C     quaternion derivative rather than as a C-matrix and angular
C     velocity vector.
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
C     CK
C
C$ Keywords
C
C     POINTING
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'ck06.inc'
      INCLUDE 'ckparam.inc'
      
 
      DOUBLE PRECISION      RECORD   ( * )
      DOUBLE PRECISION      QSTATE   ( * )
      DOUBLE PRECISION      CLKOUT

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     RECORD    I-O  Data type 6 record.
C     QSTATE     O   Interpolated output record.
C     CLKOUT     O   SCLK associated with input record.
C
C$ Detailed_Input
C
C     RECORD      is a record from a type 6 CK segment which, when
C                 evaluated at the epoch contained in its first
C                 element, will give the attitude and angular velocity
C                 of a spacecraft structure or instrument relative to a
C                 base reference frame.
C
C                 The structure of the record is as follows:
C
C                    +----------------------+
C                    | evaluation epoch     |
C                    +----------------------+
C                    | subtype code         |
C                    +----------------------+
C                    | number of packets (n)|
C                    +----------------------+
C                    | nominal SCLK rate    |
C                    +----------------------+
C                    | packet 1             |
C                    +----------------------+
C                    | packet 2             |
C                    +----------------------+
C                             .
C                             .
C                             .
C                    +----------------------+
C                    | packet n             |
C                    +----------------------+
C                    | epochs 1--n          |
C                    +----------------------+
C
C                See the CK Required Reading or the include file
C                ck06.inc for details on CK type 6 packet contents.
C
C
C$ Detailed_Output
C
C     RECORD     has been modified due to its use as a workspace array.
C                The contents are undefined.
C
C     QSTATE     is an interpolated output record, represented as a unit
C                quaternion and its derivative with respect to time.
C
C     CLKOUT     is the encoded SCLK associated with the returned
C                C-matrix and angular velocity vector.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the input record contains an unrecognized subtype code,
C         the error SPICE(NOTSUPPORTED) is signaled.
C
C     2)  If the record subtype is one for which quaternion derivatives
C         are stored (subtypes 0 and 2), and if the Ith quaternion in
C         the input record is farther than its negative from the (I-1)st
C         quaternion in the record, the error SPICE(BADQUATSIGN)
C         is signaled.
C
C         For subtypes 1 and 3, this condition is not considered an
C         error: the closer to the preceding quaternion of the two
C         quaternion representations is used for interpolation.
C
C     3)  If a zero-magnitude quaternion is produced as a result
C         of interpolating quaternions, the error SPICE(DIVIDEBYZERO)
C         is signaled.
C
C     4)  If the input record contains a non-positive SCLK rate value,
C         the error SPICE(INVALIDSCLKRATE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine returns quaternion states that, presuming validity
C     of the input record, are suitable for interpolation, unlike
C     quaternion states obtained by calling CKE06 and then M2Q. The
C     latter method of obtaining quaternions is subject to branch
C     singularities.
C   
C     The exact format and structure of CK type 6 (MEX/Rosetta Attitude
C     file interpolation) CK segments is described in the CK Required
C     Reading.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     1)  This routine performs minimal error checking. The input data
C         are assumed to have been checked when the source CK file was
C         created.
C
C     2)  With the exception of the check described in item 2 of 
C         the Exceptions section above, the input data are assumed to
C         be suitable for the interpolation method specified by the
C         input record's subtype and packet count (which implies an
C         interpolating polynomial degree).
C              
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C     B.V. Semenov   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 10-AUG-2015 (NJB) (BVS)
C
C-&


C$ Index_Entries
C
C     evaluate type_6 ck_segment
C
C-&
 
C$ Revisions
C
C     None.
C
C-&

 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      LGRINT
      DOUBLE PRECISION      VDISTG
      DOUBLE PRECISION      VDOTG
      DOUBLE PRECISION      VNORMG

      LOGICAL               RETURN
 
C
C     Local parameters
C

C
C     Index of evaluation epoch in record:
C
      INTEGER               EPCIDX
      PARAMETER           ( EPCIDX = 1 )

C
C     Index of subtype code in record:
C
      INTEGER               SBTIDX
      PARAMETER           ( SBTIDX = 2 )

C
C     Index of packet count in record:
C
      INTEGER               CNTIDX
      PARAMETER           ( CNTIDX = 3 )

C
C     Index at which packets start; packet base:
C
      INTEGER               PKTIDX
      PARAMETER           ( PKTIDX = 5 )

      INTEGER               PKTBAS
      PARAMETER           ( PKTBAS = PKTIDX - 1 )

C
C     Local variables
C
      DOUBLE PRECISION      AV     ( 3 )
      DOUBLE PRECISION      DQ     ( 0 : 3 )
      DOUBLE PRECISION      DS     ( 0 : 3 )
      DOUBLE PRECISION      LOCREC ( CKMRSZ )
      DOUBLE PRECISION      MAGS
      DOUBLE PRECISION      Q      ( 0 : 3 )
      DOUBLE PRECISION      QAV    ( 0 : 3 )
      DOUBLE PRECISION      QNEG   ( 0 : 3 )
      DOUBLE PRECISION      RADTRM ( 0 : 3 )
      DOUBLE PRECISION      RATE
      DOUBLE PRECISION      SCLDDQ ( 0 : 3 )
      DOUBLE PRECISION      SCLKDP
      DOUBLE PRECISION      STATE  ( 8 )
      DOUBLE PRECISION      VBUFF  ( 6 )
      DOUBLE PRECISION      WORK   ( CKMRSZ * 2,  2 )

      INTEGER               FROM
      INTEGER               I
      INTEGER               J
      INTEGER               N
      INTEGER               NEWPTR
      INTEGER               PACKSZ
      INTEGER               PKTSZS ( 0 : C06NST-1 )
      INTEGER               PRVPTR
      INTEGER               SUBTYP
      INTEGER               TO
      INTEGER               XSTART
      INTEGER               YSTART

C
C     Saved variables
C
      SAVE                  PKTSZS

C
C     Initial values
C
      DATA                  PKTSZS / C06PS0, C06PS1, C06PS2, C06PS3 /

C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'ZZCKE06' )

C
C     Transfer the input record's epoch to the output epoch.
C     
      CLKOUT = RECORD(1)
 
C
C     Capture the subtype from the record and set the packet size
C     accordingly.
C
      SUBTYP =  NINT( RECORD(SBTIDX) )

      
      IF (  ( SUBTYP .LT. 0 ) .OR. ( SUBTYP .GE. C06NST)  )  THEN
         
         CALL SETMSG ( 'Unexpected CK type 6 subtype # found in ' 
     .   //            'type 6 segment.'                          )
         CALL ERRINT ( '#',  SUBTYP                               )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                      )
         CALL CHKOUT ( 'ZZCKE06'                                  )
         RETURN
      
      ELSE

         PACKSZ = PKTSZS( SUBTYP )

      END IF

C
C     Get the packet count and epoch.
C
      N       =  NINT( RECORD(CNTIDX) )
      SCLKDP  =        RECORD(EPCIDX)

C
C     Get the nominal clock rate.
C
      RATE = RECORD(4)

      IF ( RATE .LE. 0.D0 ) THEN

         CALL SETMSG ( 'SCLK rate is #; rate must be positive.' )
         CALL ERRDP  ( '#', RATE                                )
         CALL SIGERR ( 'SPICE(INVALIDSCLKRATE)'                 )
         CALL CHKOUT ( 'ZZCKE06'                                )
         RETURN

      END IF

C
C     Adjust quaternion "signs" as necessary to minimize distance
C     between successive quaternions. This adjustment is performed
C     only for subtypes that don't store quaternion derivatives
C     (these are the Lagrange subtypes).
C
      IF (  ( SUBTYP .EQ. C06TP1 ) .OR. ( SUBTYP .EQ. C06TP3 )  ) THEN
C
C        For these subtypes, only the quaternions themselves need be
C        adjusted.  
C
C        PRVPTR is the index of the "previous" quaternion---the one to
C        which the successor and its negative will be compared.
C
         PRVPTR = PKTIDX

         DO I = 2, N
C
C           NEWPTR points to the quaternion ahead of the one
C           pointed to by PRVPTR.
C
            NEWPTR = PKTIDX + PACKSZ*(I-1)

            CALL VMINUG ( RECORD(NEWPTR), 4,  QNEG )

C
C           Replace the Ith quaternion with QNEG if QNEG is closer
C           than the current quaternion to the previous quaternion.
C
            IF (     VDISTG( RECORD(PRVPTR), QNEG,           4 ) 
     .          .LT. VDISTG( RECORD(PRVPTR), RECORD(NEWPTR), 4 ) ) THEN

               CALL MOVED ( QNEG, 4, RECORD( NEWPTR ) )

            END IF

            PRVPTR = NEWPTR
            
         END DO


      ELSE
C
C        For the Hermite types, if the quaternions need to be adjusted,
C        we have an error condition.
C
C        PRVPTR is the index of the "previous" quaternion---the one to
C        which the successor and its negative will be compared.
C
         PRVPTR = PKTIDX

         DO I = 2, N
C
C           NEWPTR points to the quaternion ahead of the one
C           pointed to by PRVPTR. 
C
            NEWPTR = PKTIDX + PACKSZ*(I-1)

            CALL VMINUG ( RECORD(NEWPTR), 4,  QNEG )
C
C           For the Hermite subtypes, it's an error for the current
C           quaternion to be closer to QNEG than to the previous
C           quaternion.
C
            IF (     VDISTG( RECORD(PRVPTR), QNEG,           4 ) 
     .          .LT. VDISTG( RECORD(PRVPTR), RECORD(NEWPTR), 4 ) ) THEN

               CALL SETMSG ( 'Quaternion sign error: quaternion at '
     .         //            'index # in the input record is farther '
     .         //            'than its negative from the preceding '
     .         //            'quaternion in the record. Quaternion '
     .         //            'is (#, #, #, #); predecessor is ' 
     .         //            '(#, #, #, #). This makes the quaternion '
     .         //            'sequence unsuitable for Hermite '
     .         //            'interpolation. The quaternions, and '
     .         //            'if applicable, their derivatives, '
     .         //            'must be adjusted before they are '
     .         //            'passed to this routine.'                 )
               CALL ERRINT ( '#', I                                    )
               CALL ERRDP  ( '#', RECORD(NEWPTR  )                     )
               CALL ERRDP  ( '#', RECORD(NEWPTR+1)                     )
               CALL ERRDP  ( '#', RECORD(NEWPTR+2)                     )
               CALL ERRDP  ( '#', RECORD(NEWPTR+3)                     )
               CALL ERRDP  ( '#', RECORD(PRVPTR  )                     )
               CALL ERRDP  ( '#', RECORD(PRVPTR+1)                     )
               CALL ERRDP  ( '#', RECORD(PRVPTR+2)                     )
               CALL ERRDP  ( '#', RECORD(PRVPTR+3)                     )
               CALL SIGERR ( 'SPICE(BADQUATSIGN)'                      )
               CALL CHKOUT ( 'ZZCKE06'                                 )
               RETURN

            END IF
            
            PRVPTR = NEWPTR

         END DO

      END IF



      IF ( SUBTYP .EQ. C06TP1 ) THEN
C
C        We perform Lagrange interpolation on each quaternion 
C        component, and obtain quaternion derivatives from the
C        interpolating polynomials. 
C
C        We'll transpose the pointing information in the input record so
C        that contiguous pieces of it can be shoved directly into the
C        interpolation routine LGRIND. 
C
         N   =  NINT( RECORD(CNTIDX) )

         CALL XPSGIP ( PACKSZ, N, RECORD(PKTIDX) ) 

C
C        We interpolate each state component in turn.
C
         XSTART   =   PKTIDX   +  N * PACKSZ
 
         DO I = 1, PACKSZ
 
            YSTART    =   PKTIDX  +  N * (I-1)
 

            CALL LGRIND ( N,
     .                    RECORD(XSTART),
     .                    RECORD(YSTART),
     .                    WORK,
     .                    SCLKDP,
     .                    STATE(I),
     .                    STATE(I+4)  )
     
         END DO  

C
C        The output quaternion is a unitized version of the 
C        interpolated state.
C
         MAGS = VNORMG( STATE, 4 )

         IF ( MAGS .EQ. 0.D0 ) THEN

            CALL SETMSG ( 'Quaternion magnitude at SCLK # was zero.' )
            CALL ERRDP  ( '#',  SCLKDP                               )
            CALL SIGERR ( 'SPICE(DIVIDEBYZERO)'                      )
            CALL CHKOUT ( 'ZZCKE06'                                  )
            RETURN

         END IF

         CALL VSCLG ( 1.D0/MAGS, STATE, 4, Q )

C
C        Find the time derivative of the unit quaternion:
C        Letting S represent the quaternion portion of STATE, we
C        have
C
C           Q = S/||S||
C
C
C        Then letting < , > denote the 4-dimensional inner product
C        operator, we have
C
C
C                      d(S)/dt      < Q, d(S)/dt >         
C           d(Q)/dt =  -------  -   -------------- * Q
C                       ||S||            ||S||        
C
C
         CALL MOVED ( STATE(5), 4, DS )

         CALL VSCLG ( 1.D0          / MAGS,  DS, 4,  SCLDDQ )
         CALL VSCLG ( VDOTG(Q,DS,4) / MAGS,  Q,  4,  RADTRM )

         CALL VSUBG ( SCLDDQ, RADTRM, 4, DQ )
 
C
C        Scale the derivative from 1/tick to 1/second.
C
         CALL VSCLG ( 1.D0/RATE,  DQ,  4,  SCLDDQ )

         CALL MOVED ( Q,      4, QSTATE(1) )
         CALL MOVED ( SCLDDQ, 4, QSTATE(5) )


      
      ELSE IF ( SUBTYP .EQ. C06TP3 ) THEN
C
C        This is the easiest case:  we perform Lagrange interpolation
C        on each quaternion or angular velocity component.
C
C        We'll transpose the pointing information in the input record so
C        that contiguous pieces of it can be shoved directly into the
C        interpolation routine LGRINT.  We allow LGRINT to overwrite
C        the state values in the input record, since this saves local
C        storage and does no harm.  (See the header of LGRINT for a
C        description of its work space usage.)
C
         N   =  NINT( RECORD(CNTIDX) )
 
         CALL XPSGIP ( PACKSZ, N, RECORD(PKTIDX) ) 

C
C        We interpolate each state component in turn.
C
         XSTART  =  PKTIDX   +  N * PACKSZ
 
         DO I = 1, PACKSZ
 
            YSTART    =   PKTIDX  +  N * (I-1)
 
            STATE(I)  =   LGRINT ( N,
     .                             RECORD(XSTART),
     .                             RECORD(YSTART),
     .                             LOCREC,
     .                             SCLKDP          )
         END DO  


         MAGS = VNORMG( STATE, 4 )

         IF ( MAGS .EQ. 0.D0 ) THEN

            CALL SETMSG ( 'Quaternion magnitude at SCLK # was zero.' )
            CALL ERRDP  ( '#',  SCLKDP                               )
            CALL SIGERR ( 'SPICE(DIVIDEBYZERO)'                      )
            CALL CHKOUT ( 'ZZCKE06'                                  )
            RETURN

         END IF

         CALL VSCLG ( 1.D0/MAGS, STATE, 4, Q )

C
C        The angular velocity already is in units of radians/second.
C
         CALL VEQU ( STATE(5), AV )

C
C        Convert AV to a quaternion derivative. We have from
C        the header of QDQ2AV
C
C                       *    
C           AV =  -2 * Q  * DQ  
C
C        so
C
C           DQ =  -1/2 * Q * AV 
C
C        
         CALL VSCLIP ( -0.5D0, AV )

         QAV(0) = 0.D0
         CALL VEQU ( AV, QAV(1) )

         CALL QXQ ( Q, QAV, DQ )


         CALL MOVED ( Q,  4, QSTATE    )
         CALL MOVED ( DQ, 4, QSTATE(5) )

      ELSE
C
C        We have a Hermite-style subtype.  Whether it's subtype 0
C        or 2, we perform Hermite interpolation on the quaternions.
C
C        We interpolate each quaternion component in turn.  Attitude and
C        angular velocity are interpolated separately.
C
         XSTART   =   PKTIDX  +  PACKSZ * N 
   
         DO I = 1, 4
 
            DO J = 1, N
C
C              For the Jth input packet, copy the Ith position and
C              velocity components into the local record buffer RECORD.
C
C              In order to perform Hermite interpolation, the
C              quaternions and quaternion derivatives must have a
C              common time scale. So prior to interpolation, we scale
C              the units of the quaternion derivatives from radians/sec
C              to radians/tick.
C
               FROM         = PKTBAS + PACKSZ*(J-1) + I
               TO           =              2 * J    - 1
            
               LOCREC(TO  ) = RECORD ( FROM     )
               LOCREC(TO+1) = RECORD ( FROM + 4 ) * RATE

            END DO

C
C           Interpolate the Ith quaternion and quaternion derivative
C           components.
C        
            CALL HRMINT ( N, 
     .                    RECORD(XSTART),
     .                    LOCREC,
     .                    SCLKDP,           
     .                    WORK,
     .                    STATE(I  ),
     .                    STATE(I+4)      ) 

         END DO

C
C        The output quaternion is a unitized version of the 
C        interpolated state.
C
         MAGS = VNORMG( STATE, 4 )

         IF ( MAGS .EQ. 0.D0 ) THEN

            CALL SETMSG ( 'Quaternion magnitude at SCLK # was zero.' )
            CALL ERRDP  ( '#',  SCLKDP                               )
            CALL SIGERR ( 'SPICE(DIVIDEBYZERO)'                      )
            CALL CHKOUT ( 'ZZCKE06'                                  )
            RETURN

         END IF

         CALL VSCLG ( 1.D0/MAGS, STATE, 4, Q )


 
         IF ( SUBTYP .EQ. C06TP0 ) THEN
C
C           Find the time derivative of the unit quaternion:
C           Letting S represent the quaternion portion of STATE, we
C           have
C
C              Q = S/||S||
C
C
C           Then letting < , > denote the 4-dimensional inner product
C           operator, we have
C
C
C                         d(S)/dt      < Q, d(S)/dt >         
C              d(Q)/dt =  -------  -   -------------- * Q
C                          ||S||            ||S||        
C
C
            CALL MOVED ( STATE(5), 4, DS )

            CALL VSCLG ( 1.D0          / MAGS,  DS, 4,  SCLDDQ )
            CALL VSCLG ( VDOTG(Q,DS,4) / MAGS,  Q,  4,  RADTRM )

            CALL VSUBG ( SCLDDQ, RADTRM, 4, DQ )
C
C           Scale the derivative from radians/tick to
C           radians/second.
C
            CALL VSCLG ( 1.D0/RATE, DQ, 4, SCLDDQ )

C
C           Store Q and DQ in QSTATE. In the process, 
C
            CALL MOVED ( Q,      4, QSTATE    )
            CALL MOVED ( SCLDDQ, 4, QSTATE(5) )


         ELSE  
C
C           This is subtype 2; we perform Hermite interpolation on
C           the angular velocity and its derivative.
C
C           Now interpolate angular velocity, using separate angular
C           velocity data and angular acceleration.
C
            DO I = 1, 3
  
               DO J = 1, N
C
C                 For the Jth input packet, copy the Ith angular
C                 velocity and angular acceleration components into the
C                 local record buffer LOCREC.  Note that, as with
C                 quaternion derivatives, we must scale angular
C                 acceleration from radians/sec**2 to
C                 radians/(sec*tick) before interpolating. We would
C                 need to scale the angular acceleration to
C                 radians/sec**2 for output, if we were returning this
C                 quantity. However, we're returning only angular
C                 velocity, which is already in the correct units of
C                 radians/second.
C
                  FROM         = PKTBAS + PACKSZ*(J-1) +  8  +  I
                  TO           =              2 * J    -  1
         
                  LOCREC(TO  ) = RECORD ( FROM     )
                  LOCREC(TO+1) = RECORD ( FROM + 3 ) * RATE
            
               END DO
 
C
C              Interpolate the Ith angular velocity and angular
C              acceleration components of the attitude. We'll
C              capture the result in a temporary buffer, then
C              transfer the velocity to the output argument AV.
C
               CALL HRMINT ( N, 
     .                       RECORD(XSTART),
     .                       LOCREC,
     .                       SCLKDP,           
     .                       WORK,
     .                       VBUFF(I  ),
     .                       VBUFF(I+3)     ) 
 
            END DO
C
C           Fill in the angular velocity in the output angular
C           velocity vector using the results of interpolating
C           velocity and acceleration.
C
C           The angular velocity is already in units of
C           radians/second.
C
            CALL VEQU ( VBUFF, AV )
C
C           Convert AV to a quaternion derivative. We have from
C           the header of QDQ2AV
C
C                          *    
C              AV =  -2 * Q  * DQ  
C
C           so
C
C              DQ =  -1/2 * Q * AV 
C
C        
            CALL VSCLIP ( -0.5D0, AV )

            QAV(0) = 0.D0
            CALL VEQU ( AV, QAV(1) )

            CALL QXQ ( Q, QAV, DQ )

            CALL MOVED ( Q,  4, QSTATE    )
            CALL MOVED ( DQ, 4, QSTATE(5) )

         END IF
C
C        We've handled the type 0 and type 2 cases.
C              
C
C        We've computed the angular velocity AV for the Hermite
C        subtypes, if a.v. was requested.
C
      END IF
C
C     We've handled all four subtypes.
C

      CALL CHKOUT ( 'ZZCKE06' )
      RETURN
      END
 
