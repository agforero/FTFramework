program main

!*****************************************************************************80
!
!! MAIN is the main program for MIXTURE.
!
!  Discussion:
!
!    MIXTURE carries out the mixture simulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: acid_num = 20
  integer ( kind = 4 ), parameter :: comp_max = 15
  integer ( kind = 4 ), parameter :: iunit = 1

  character acid_sym(acid_num)
  real ( kind = 8 ) alpha(acid_num,comp_max)
  real ( kind = 8 ) alpha_sum(comp_max)
  integer ( kind = 4 ) comp_label(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_weight(comp_max)
  real ( kind = 8 ) comp_weight_post(comp_max)
  integer ( kind = 4 ) ierror
  character ( len = 30 ) mixture_file_name
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) site_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MIXTURE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Last modified on 09 May 2002.' 
!
!  Read mixture parameters from a file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Read mixture data from file.'
  write ( *, '(a)' ) ' '

  mixture_file_name = 'mixture.dat'

  open ( unit = iunit, file = mixture_file_name, form = 'formatted' )

  call mixture_read ( acid_num, acid_sym, alpha, alpha_sum, &
    comp_label, comp_max, comp_num, comp_weight, ierror, iunit )

  close ( unit = iunit )
!
!  Print the amino acid parameters.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Numeric key for amino acid abbreviations.'
  write ( *, '(a)' ) ' '

  call amino_print ( acid_num, acid_sym )
!
!  Print the component parameters.
!
  call comp_param_print ( acid_num, acid_sym, comp_max, comp_num, &
    alpha, alpha_sum, comp_weight )
!
!  Initialize the Dirichlet parameter estimates.
!
  call weight_init ( comp_weight_post, comp_num )
!
!  Repeatedly observe the process and update the parameter estimates.
!
  sample_num = 10000
  site_num = 10

  call observe ( acid_num, alpha, alpha_sum, comp_max, comp_num, comp_weight, &
    comp_weight_post, sample_num, site_num )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MIXTURE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine amino_print ( acid_num, acid_sym )

!*****************************************************************************80
!
!! AMINO_PRINT prints the amino acid parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
  implicit none

  integer ( kind = 4 ) acid_num

  integer ( kind = 4 ) acid_i
  character ( len = 27 ) acid_name
  character acid_sym(acid_num)
  character c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  Amino Acid Symbol'
  write ( *, '(a)' ) ' '
  do acid_i = 1, acid_num
    c = acid_sym(acid_i)
    call ch_to_amino_name ( c, acid_name )
    write ( *, '(i3,2x,a,2x,a)' ) acid_i, acid_sym(acid_i), acid_name
  end do

  return
end
subroutine binomial_sample ( a, b, x )

!*****************************************************************************80
!
!! BINOMIAL_SAMPLE samples the Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Algorithm BU,
!    William Kennedy and James Gentle,
!    Statistical Computing,
!    Dekker, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) u
  integer ( kind = 4 ) x

  x = 0

  do i = 1, a

    call r8_random ( 0.0D+00, 1.0D+00, u )

    if ( u <= b ) then
      x = x + 1
    end if

  end do

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_next ( line, cval, done )

!*****************************************************************************80
!
!! CH_NEXT "reads" space-separated characters from a string, one at a time.
!
!  Example:
!
!    Input:
!
!      LINE = ' A  B, C    DE  F'
!
!    Output:
!
!      'A', 'B', 'C', 'D', 'E', 'F', and then blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    characters, possibly separated by spaces or commas.
!
!    Output, character CVAL.  If DONE is FALSE, then CVAL contains the
!    "next" character read from LINE.  If DONE is TRUE, then
!    CVAL is blank.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another character
!    was read, or TRUE if no more characters could be read.
!
  implicit none

  character cval
  logical done
  integer ( kind = 4 ) i
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1

  if ( done ) then
    next = 1
    done = .false.
  end if

  do i = next, len(line)

    if ( line(i:i) /= ' ' .and. line(i:i) /= ',' ) then
      cval = line(i:i)
      next = i + 1
      return
    end if

  end do

  done = .true.
  next = 1
  cval = ' '

  return
end
subroutine ch_to_amino_name ( c, amino_name )

!*****************************************************************************80
!
!! CH_TO_AMINO_NAME converts a character to an amino acid name.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl Branden, John Tooze,
!    Introduction to Protein Structure,
!    Garland Publishing, 1991.
!
!  Parameters:
!
!    Input, character C, the one letter code for an amino acid.
!    Lower and upper case letters are treated the same.
!
!    Output, character ( len = * ) AMINO_NAME, the full name of the
!    corresponding amino acid.  The longest name is 27 characters.  If 
!    the input code is not recognized, then AMINO_NAME will be set to '???'.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 23

  character ( len = * ) amino_name
  character ( len = 27 ), dimension ( n ) :: amino_table = (/ &
    'Alanine                    ', &
    'Aspartic acid or Asparagine', &
    'Cysteine                   ', &
    'Aspartic acid              ', &
    'Glutamic acid              ', &
    'Phenylalanine              ', &
    'Glycine                    ', &
    'Histidine                  ', &
    'Isoleucine                 ', &
    'Lysine                     ', &
    'Leucine                    ', &
    'Methionine                 ', &
    'Asparagine                 ', &
    'Proline                    ', &
    'Glutamine                  ', &
    'Arginine                   ', &
    'Serine                     ', &
    'Threonine                  ', &
    'Valine                     ', &
    'Tryptophan                 ', &
    'Undetermined amino acid    ', &
    'Tyrosine                   ', &
    'Glutamic acid or Glutamine ' /)
  character c
  logical ch_eqi
  character, dimension ( n ) :: c_table = (/ &
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
    'X', 'Y', 'Z' /)
  integer ( kind = 4 ) i

  do i = 1, n
    if ( ch_eqi ( c, c_table(i) ) ) then
      amino_name = amino_table(i)
      return
    end if
  end do

  amino_name = '???'

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine comp_param_print ( acid_num, acid_sym, comp_max, comp_num, &
  beta, beta_sum, comp_weight )

!*****************************************************************************80
!
!! COMP_PARAM_PRINT prints the parameters for the mixture components.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    Input, real ( kind = 8 ) BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the 
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    Input, real ( kind = 8 ) BETA_SUM(COMP_MAX), the sum of the values of
!    BETA(ACID_I,COMP_I) for a given component COMP_I.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weight of each
!    component.  These values should be nonnegative, and sum to 1.  They 
!    represent the relative proportion of each component in the mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) beta(acid_num,comp_max)
  real ( kind = 8 ) beta_sum(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_weight(comp_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of components = ', comp_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(''Compon:'',20i8)' ) ( comp_i, comp_i = 1, comp_num )
  write ( *, '(''Weight:'',20f8.4)' ) comp_weight(1:comp_num)
  write ( *, '(a)' ) ' '

  do acid_i = 1, acid_num
    write ( *, '(i2,2x,a1,2x,20f8.4)' )acid_i, acid_sym(acid_i), &
      beta(acid_i,1:comp_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a3,4x,20f8.4)' ) 'Sum', beta_sum(1:comp_num)

  return
end
subroutine comp_stats_print ( acid_num, comp_max, comp_num, comp_mean, &
  comp_variance )

!*****************************************************************************80
!
!! COMP_STATS_PRINT prints the mean and variance for the mixture components.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_mean(comp_max,acid_num)
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_variance(comp_max,acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected means for each component PDF:'
  write ( *, '(a)' ) ' '
  do acid_i = 1, acid_num
    write ( *, '(9f8.4)' ) comp_mean(1:comp_num,acid_i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected variances for each component PDF:'
  write ( *, '(a)' ) ' '
  do acid_i = 1, acid_num
    write ( *, '(9f8.4)' ) comp_variance(1:comp_num,acid_i)
  end do

  return
end
subroutine discrete_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! DISCRETE_CDF_INV inverts the Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes 
!    1 through A.  Each entry must be nonnegative.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument for which
!    CDF(X-1) < CDF <= CDF(X)
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  b_sum = sum ( b(1:a) )

  cum = 0.0D+00

  do j = 1, a

    cum = cum + b(j) / b_sum

    if ( cdf <= cum ) then
      x = j
      return
    end if

  end do

  x = a

  return
end
subroutine discrete_sample ( a, b, x )

!*****************************************************************************80
!
!! DISCRETE_SAMPLE samples the Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes 
!    1 through A.  Each entry must be nonnegative.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  b_sum = sum ( b(1:a) )

  call r8_random ( 0.0D+00, 1.0D+00, cdf )

  call discrete_cdf_inv ( cdf, a, b, x )

  return
end
subroutine favor_ratio_compute ( acid_num, alpha, alpha_sum, comp_max, &
  comp_num, comp_weight, favor )

!*****************************************************************************80
!
!! FAVOR_RATIO_COMPUTE computes the amino acid favor ratio.
!
!  Discussion:
!
!    This routine computes the ratio by which a component density 
!    favors an amino acid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, real ( kind = 8 ) ALPHA(COMP_MAX,ACID_NUM); ALPHA(I,J) is the
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    Input, real ( kind = 8 ) ALPHA_SUM(COMP_MAX), the sum of the values of 
!    ALPHA(COMP_I,ACID_I) for a given component COMP_I.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components 
!    in the Dirichlet mixture.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weight of each
!    component.  These values should be nonnegative, and sum to 1.  They 
!    represent the relative proportion of each component in the mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  real ( kind = 8 ) alpha(acid_num,comp_max)
  real ( kind = 8 ) alpha_sum(comp_max)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) factor(acid_num)
  real ( kind = 8 ) favor(comp_max,acid_num)
  real ( kind = 8 ) comp_weight(comp_num)

  do acid_i = 1, acid_num

    factor(acid_i) = 0.0D+00
    do comp_i = 1, comp_num
      factor(acid_i) = factor(acid_i) + comp_weight(comp_i) * &
        alpha(acid_i,comp_i) / alpha_sum(comp_i)
    end do
  end do

  do acid_i = 1, acid_num
    do comp_i = 1, comp_num
      favor(comp_i,acid_i) = alpha(acid_i,comp_i) / alpha_sum(comp_i) &
        / factor(acid_i)
    end do
  end do

  return
end
subroutine favor_ratio_print ( acid_num, acid_sym, comp_max, comp_num, favor )

!*****************************************************************************80
!
!! FAVOR_RATIO_PRINT prints the favor ratios.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in 
!    the Dirichlet mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) favor(comp_max,acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Favor ratios:'
  write ( *, '(a)' ) ' '
  write ( *, '(''Component: '',9i6)' ) ( comp_i, comp_i = 1, comp_num )
  write ( *, '(a)' ) ' '

  do acid_i = 1, acid_num
    write ( *, '(i2,2x,a1,2x,9f6.2)' ) acid_i, acid_sym(acid_i), &
      ( favor(comp_i,acid_i), comp_i = 1, comp_num )
  end do

  return
end
subroutine i4_next ( line, ival, done )

!*****************************************************************************80
!
!! I4_NEXT "reads" integers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    integer ( kind = 4 )s.  These may be separated by spaces or commas.
!
!    Output, integer ( kind = 4 ) IVAL.  If DONE is FALSE, then IVAL contains
!    the "next" integer read from LINE.  If DONE is TRUE, then
!    IVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another integer
!    was read, or TRUE if no more integers could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1

  ival = 0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( next > len(line) ) then
    done = .true.
    return
  end if

  call s_to_i4 ( line(next:), ival, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine mixture_print ( acid_num, acid_sym, alpha, alpha_sum, comp_label, &
  comp_max, comp_num, comp_weight )

!*****************************************************************************80
!
!! MIXTURE_PRINT prints the Dirichlet mixture parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Input, real ( kind = 8 ) ALPHA(COMP_MAX,ACID_NUM); ALPHA(I,J) is the 
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    Input, real ( kind = 8 ) ALPHA_SUM(COMP_MAX), the sum of the values of 
!    ALPHA(COMP_I,ACID_I) for a given component COMP_I.
!
!    Input, integer ( kind = 4 ) COMP_LABEL(COMP_NUM), the label of 
!    each component.  Normally, component I has label I.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weight of each
!    component.  These values should be nonnegative, and sum to 1.  They 
!    represent the relative proportion of each component in the mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  real ( kind = 8 ) alpha(acid_num,comp_max)
  real ( kind = 8 ) alpha_sum(comp_max)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_label(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_weight(comp_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'Number of components in the Dirichlet mixture:', &
    comp_num

  do comp_i = 1, comp_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Component ', comp_i
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Label 	 = ', comp_label(comp_i)
    write ( *, '(a,g14.6)' ) '  Mixture weight = ', comp_weight(comp_i)
    write ( *, '(a,g14.6)' ) '  Parameter sum  = ', alpha_sum(comp_i)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Parameters:'
    write ( *, '(a)' ) ' '
    do acid_i = 1, acid_num
      write ( *, '(i2,2x,a1,2x,g14.6)' ) acid_i, acid_sym(acid_i), &
        alpha(acid_i,comp_i)
    end do

  end do

  return
end
subroutine mixture_read ( acid_num, acid_sym, beta, beta_sum, comp_label, &
  comp_max, comp_num, comp_weight, ierror, iunit )

!*****************************************************************************80
!
!! MIXTURE_READ reads the Dirichlet mixture parameters from a file.
!
!  Discussion:
!
!    The data in the file is delimited by keywords.
!
!    The first lines (not necessarily in order!) may include
!
!      ClassName = string
!      NumDistr = N           the number of components in the mixture.
!      Alphabet = string
!      Order = A C D E ...    the order of the amino acids.
!      AlphaChar = 20
!      NumDistr = 9           the number of distributions
!      EndClassName = string
!
!    For each component, there are four lines:
!
!      Number= N              the component number, starting with 0
!      Mixture= N             the mixture weight, out of a total of 1.0
!      Alpha=  |A| A1 A2 ...  the parameter sum, and individual parameters
!      Comment=               a comment, which describes the frequencies.
!
!    In the comment, the symbol "><" indicates the mean background frequency;
!    residues to the left of that symbol occur more frequently
!    than background, residues to the right less frequently.  Commas separate
!    residues differing in frequency by a factor of 2.
!
!    For example, the comment
!      S A T , C G P >< N V M , Q H R I K F L D W , E Y
!    indicates that for this component, the frequency of
!    proline is just above the mean, and serine, alanine and
!    threonine are twice as frequent in this component than they
!    are on average.  By contrast, tyrosine and glutamic acid are
!    between 4 and 8 times less likely in this component than on
!    average.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Output, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Output, real ( kind = 8 ) BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the 
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    Output, real ( kind = 8 ) BETA_SUM(COMP_MAX), the sum of the values of
!    BETA(ACID_I,COMP_I) for a given component COMP_I.
!
!    Output, integer ( kind = 4 ) COMP_LABEL(COMP_NUM), the label of each 
!    component.  Normally, component I has label I.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet
!    mixture components.
!
!    Output, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    Output, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weight of each
!    component.  These values should be nonnegative, and sum to 1.  They 
!    represent the relative proportion of each component in the mixture.
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0: no error occurred; nonzero: an error occurred.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which the data
!    is to be read.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  real ( kind = 8 ) beta(acid_num,comp_max)
  real ( kind = 8 ) beta_sum(comp_max)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_label(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_weight(comp_max)
  logical done
  integer ( kind = 4 ) iequal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ngoofy
  integer ( kind = 4 ) nrec
  logical s_begin
  character ( len = 500 ) string

  ierror = 0
  comp_i = 0
  comp_num = 0
  nrec = 0
  ngoofy = 0

  do

    read ( iunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    nrec = nrec + 1
!
!  Ignore blank lines.
!
    if ( string == ' ' ) then
!
!  Ignore the CLASSNAME field.
!
    else if ( s_begin ( string, 'CLASSNAME' ) ) then
!
!  Ignore the ENDCLASSNAME field.
!
    else if ( s_begin ( string, 'ENDCLASSNAME' ) ) then
!
!  Ignore the NAME field.
!
    else if ( s_begin ( string, 'NAME' ) ) then
!
!  Ignore the ALPHABET field.
!
    else if ( s_begin ( string, 'ALPHABET' ) ) then
!
!  Read the ORDER field, since it tells us how to interpret the ALPHA's.
!
    else if ( s_begin ( string, 'ORDER' ) ) then

      iequal = index ( string, '=' )
      done = .true.
      do acid_i = 1, acid_num
        call ch_next ( string(iequal+1:), acid_sym(acid_i), done )
      end do
!
!  Ignore the ALPHACHAR field.
!
    else if ( s_begin ( string, 'ALPHACHAR' ) ) then
!
!  Read the NUMDISTR field.
!
    else if ( s_begin ( string, 'NUMDISTR' ) ) then

      iequal = index ( string, '=' )
      done = .true.
      call i4_next ( string(iequal+1:), comp_num, done )

      if ( comp_num < 1 ) then
        ierror = 1
        return
      else if ( comp_num > comp_max ) then
        ierror = 2
        return
      end if
!
!  Read the NUMBER field.
!
    else if ( s_begin ( string, 'NUMBER' ) ) then

      comp_i = comp_i + 1

      if ( comp_i > comp_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MIXTURE_READ - Fatal error!'
        write ( *, '(a,i6)' ) '  Number of components = ', comp_i
        write ( *, '(a,i6)' ) '  exceeding reported value of ', comp_num
        stop
      end if

      iequal = index ( string, '=' )
      done = .true.
      call i4_next ( string(iequal+1:), comp_label(comp_i), done )
!
!  Read the MIXTURE field.
!
    else if ( s_begin ( string, 'MIXTURE' ) ) then

      iequal = index ( string, '=' )
      done = .true.
      call r8_next ( string(iequal+1:), comp_weight(comp_i), done )
!
!  Read the ALPHA field.
!
    else if ( s_begin ( string, 'ALPHA' ) ) then

      iequal = index ( string, '=' )
      done = .true.
      call r8_next ( string(iequal+1:), beta_sum(comp_i), done )
      do acid_i = 1, acid_num
        call r8_next ( string(iequal+1:), beta(acid_i,comp_i), done )
      end do
!
!  Ignore the COMMENT field.
!
    else if ( s_begin ( string, 'COMMENT' ) ) then
!
!  Unexpected field:
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIXTURE_READ - Warning!'
      write ( *, '(a)' ) '  Goofy record: '
      write ( *, '(a)' ) trim ( string )

      ngoofy = ngoofy + 1
      stop

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MIXTURE_READ - Note:'
  write ( *, '(a,i6)' ) '  Number of records read was ', nrec
  write ( *, '(a,i6)' ) '  Number of goofy records was ', ngoofy

  return
end
subroutine multinomial_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! MULTINOMIAL_PDF computes a Multinomial PDF.
!
!  Discussion:
!
!    PDF(X)(A,B,C) = Comb(A,B,X) * Product ( 1 <= I <= B ) C(I)^X(I)
!
!    where Comb(A,B,X) is the multinomial coefficient
!      C( A; X(1), X(2), ..., X(B) )
!
!    PDF(X)(A,B,C) is the probability that in A trials there
!    will be exactly X(I) occurrences of event I, whose probability
!    on one trial is C(I), for I from 1 to B.
!
!    As soon as A or B gets large, the number of possible X's explodes,
!    and the probability of any particular X can become extremely small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    Input, integer ( kind = 4 ) A, the total number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of different possible outcomes on
!    one trial.
!
!    Input, integer ( kind = 4 ) C(B); C(I) is the probability of outcome I on
!    any one trial.
!
!    Output, real ( kind = 8 ) PDF, the value of the multinomial PDF.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) c(b)
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  real ( kind = 8 ) pdf_log
  real ( kind = 8 ) r8_gamma_log
  integer ( kind = 4 ) x(b)
!
!  To try to avoid overflow, do the calculation in terms of logarithms.
!  Note that Gamma(A+1) = A factorial.
!
  arg = real ( a + 1, kind = 8 )
  pdf_log = r8_gamma_log ( arg )

  do i = 1, b
    arg = real ( x(i) + 1, kind = 8 )
    pdf_log = pdf_log + x(i) * log ( c(i) ) - r8_gamma_log ( arg )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine multinomial_sample ( a, b, c, x )

!*****************************************************************************80
!
!! MULTINOMIAL_SAMPLE samples the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer-Verlag, New York, 1986, page 559.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the total number of trials.
!    0 <= A.
!
!    Input, integer ( kind = 4 ) B, the number of outcomes possible
!    on one trial.  1 <= B.
!
!    Input, real ( kind = 8 ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    SUM ( 1 <= I <= B) C(I) = 1.0.
!
!    Output, integer ( kind = 4 ) X(B); X(I) is the number of
!    occurrences of event I during the N trials.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifactor
  integer ( kind = 4 ) ntot
  real ( kind = 8 ) prob
  real ( kind = 8 ) sum2
  integer ( kind = 4 ) x(b)

  ntot = a

  sum2 = 1.0D+00

  x(1:b) = 0

  do ifactor = 1, b - 1

    prob = c(ifactor) / sum2
!
!  Generate a binomial random deviate for NTOT trials with
!  single trial success probability PROB.
!
    call binomial_sample ( ntot, prob, x(ifactor) )

    ntot = ntot - x(ifactor)
    if ( ntot <= 0 ) then
      return
    end if

    sum2 = sum2 - c(ifactor)

  end do
!
!  The last factor gets what's left.
!
  x(b) = ntot

  return
end
subroutine observe ( acid_num, alpha, alpha_sum, comp_max, comp_num, &
  comp_weight, comp_weight_post, sample_num, site_num )

!*****************************************************************************80
!
!! OBSERVE repeatedly observes the process and updates the parameter estimates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, real ( kind = 8 ) ALPHA(COMP_MAX,ACID_NUM); ALPHA(I,J) is the 
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    Input, real ( kind = 8 ) ALPHA_SUM(COMP_MAX), the sum of the values of 
!    ALPHA(COMP_I,ACID_I) for a given component COMP_I.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the
!    Dirichlet mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  real ( kind = 8 ) alpha(acid_num,comp_max)
  real ( kind = 8 ) alpha_sum(comp_max)
  integer ( kind = 4 ) comp
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_weight(comp_max)
  real ( kind = 8 ) comp_weight_post(comp_max)
  real ( kind = 8 ) comp_weight_prior(comp_max)
  integer ( kind = 4 ) count(comp_max)
  real ( kind = 8 ) prob(acid_num,comp_max)
  integer ( kind = 4 ) sample_i
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) site_num
  real ( kind = 8 ) sum2
  integer ( kind = 4 ) x(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'Site numbers = ', site_num
  write ( *, '(a)' ) ' '

  count(1:comp_num) = 0

  do sample_i = 1, sample_num

    comp_weight_prior(1:comp_num) = comp_weight_post(1:comp_num)

    sum2 = sum ( comp_weight_post(1:comp_num) )
!
!  Choose a particular density component COMP.
!
    call discrete_sample ( comp_num, comp_weight, comp )

    count(comp) = count(comp) + 1

    if ( sample_i <= 10 .or. mod ( 10 * sample_i, sample_num ) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Sample number ', sample_i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  WEIGHTS:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  True, Estimated, Est norm, Count:'
      write ( *, '(a)' ) ' '
      do comp_i = 1, comp_num
        write ( *, '(i6,5g14.6)' ) comp_i, comp_weight(comp_i), &
          comp_weight_post(comp_i), comp_weight_post(comp_i) / sum2, &
          real ( count(comp_i), kind = 8 ) / real ( sample_i, kind = 8 )
      end do
    end if
!
!  Sample the density number COMP.
!
!  Well, I know the damn ALPHA's don't have to add up to 1, but
!  probabilities do, so here's yet another shot in the dark.
!
    do comp_i = 1, comp_num
      do acid_i = 1, acid_num
        prob(acid_i,comp_i) = alpha(acid_i,comp_i) / alpha_sum(comp_i)
      end do
    end do

    call multinomial_sample ( site_num, acid_num, prob(1,comp), x )
!
!  Update the Dirichlet parameter estimates.
!
    call weight_update ( acid_num, comp_num, comp_weight_post, &
      comp_weight_prior, prob, site_num, x )

  end do

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2.  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0 are similar to those in
!    reference 1, but are unpublished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ) :: d1 = -5.772156649015328605195174D-01
  real ( kind = 8 ) :: d2 = 4.227843350984671393993777D-01
  real ( kind = 8 ) :: d4 = 1.791759469228055000094023D+00
  real ( kind = 8 ), parameter :: frtbig = 2.25D+76
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 2.55D+305
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.79D+308
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq

  y = x

  if ( 0.0D+00 < y .and. y <= xbig ) then

    if ( y <= epsilon ( y ) ) then

      res = - log ( y )
!
!  EPS < X <= 1.5.
!
    else if ( y <= 1.5D+00 ) then

      if ( y < 0.6796875D+00 ) then
        corr = -log ( y )
        xm1 = y
      else
        corr = 0.0D+00
        xm1 = ( y - 0.5D+00 ) - 0.5D+00
      end if

      if ( y <= 0.5D+00 .or. 0.6796875D+00 <= y ) then

        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm1 + p1(i)
          xden = xden * xm1 + q1(i)
        end do

        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

      else

        xm2 = ( y - 0.5D+00 ) - 0.5D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

      end if
!
!  1.5 < X <= 4.0.
!
    else if ( y <= 4.0D+00 ) then

      xm2 = y - 2.0D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
!
!  4.0 < X <= 12.0.
!
    else if ( y <= 12.0D+00 ) then

      xm4 = y - 4.0D+00
      xden = -1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm4 + p4(i)
        xden = xden * xm4 + q4(i)
      end do

      res = d4 + xm4 * ( xnum / xden )
!
!  Evaluate for 12 <= argument.
!
    else

      res = 0.0D+00

      if ( y <= frtbig ) then

        res = c(7)
        ysq = y * y

        do i = 1, 6
          res = res / ysq + c(i)
        end do

      end if

      res = res / y
      corr = log ( y )
      res = res + sqrtpi - 0.5D+00 * corr
      res = res + y * ( corr - 1.0D+00 )

    end if
!
!  Return for bad arguments.
!
  else

    res = xinf

  end if
!
!  Final adjustments and return.
!
  r8_gamma_log = res

  return
end
subroutine r8_next ( line, rval, done )

!*****************************************************************************80
!
!! R8_NEXT "reads" real numbers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    Output, real ( kind = 8 ) RVAL.  If DONE is FALSE, then RVAL contains the
!    "next" real value read from LINE.  If DONE is TRUE, then
!    RVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  real ( kind = 8 ) rval

  rval = 0.0D+00

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( len ( line ) < next ) then
    done = .true.
    return
  end if

  call s_to_r8 ( line(next:), rval, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine r8_random ( rlo, rhi, r )

!*****************************************************************************80
!
!! R8_RANDOM returns a random real in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RLO, RHI, the minimum and maximum values.
!
!    Output, real ( kind = 8 ) R, the randomly chosen value.
!
  implicit none

  logical, save :: seed = .false.
  real ( kind = 8 ) r
  real ( kind = 8 ) rhi
  real ( kind = 8 ) rlo
  real ( kind = 8 ) t

  if ( .not. seed ) then
    call random_seed
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0D+00 - t ) * rlo + t * rhi

  return
end
subroutine r8row_mean ( lda, m, n, a, mean )

!*****************************************************************************80
!
!! R8ROW_MEAN returns the means of rows of a real array.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      2
!      5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input, real ( kind = 8 ) A(LDA,N), the array to be examined.
!
!    Output, real ( kind = 8 ) MEAN(M), the means, or averages, of the rows.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mean(m)

  do i = 1, m

    mean(i) = sum ( a(i,1:n) ) / real ( n, kind = 8 )

  end do

  return
end
subroutine r8row_variance ( lda, m, n, a, variance )

!*****************************************************************************80
!
!! R8ROW_VARIANCE returns the variances of the rows of a real array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(LDA,N), the array whose variances are desired.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variances of the rows.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance(m)

  do i = 1, m

    mean = sum ( a(i,1:n) ) / real ( n, kind = 8 )

    variance(i) = 0.0D+00
    do j = 1, n
      variance(i) = variance(i) + ( a(i,j) - mean )**2
    end do

    if ( 1 < n ) then
      variance(i) = variance(i) / real ( n - 1, kind = 8 )
    else
      variance(i) = 0.0D+00
    end if

  end do

  return
end
function s_begin ( s1, s2 )

!*****************************************************************************80
!
!! S_BEGIN is TRUE if one string matches the beginning of the other.
!
!  Example:
!
!     S1              S2      S_BEGIN
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     FALSE    (Do not allow a blank to match
!                                       anything but another blank string.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be compared.
!
!    Output, logical S_BEGIN, .TRUE. if the strings match up to
!    the end of the shorter string, ignoring case,
!    FALSE otherwise.
!
  implicit none

  logical ch_eqi
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  logical s_begin
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len_trim ( s1 )
  len2 = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not
!  what most people want.
!
  if ( len1 == 0 .or. len2 == 0 ) then

    if ( len1 == 0 .and. len2 == 0 ) then
      s_begin = .true.
    else
      s_begin = .false.
    end if

    return

  end if

  i1 = 0
  i2 = 0
!
!  Find the next nonblank in S1.
!
  do

    do

      i1 = i1 + 1

      if ( i1 > len1 ) then
        s_begin = .true.
        return
      end if

      if ( s1(i1:i1) /= ' ' ) then
        exit
      end if

    end do
!
!  Find the next nonblank in S2.
!
    do

      i2 = i2 + 1

      if ( i2 > len2 ) then
        s_begin = .true.
        return
      end if

      if ( s2(i2:i2) /= ' ' ) then
        exit
      end if

    end do
!
!  If the characters match, get the next pair.
!
    if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
      exit
    end if

  end do

  s_begin = .false.

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads a real number from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine sample_analyze ( acid_num, acid_sym, mixture_mean, sample_num, &
  sample_x )

!*****************************************************************************80
!
!! SAMPLE_ANALYZE analyzes the samples from the Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of samples to take.
!
!    Input, real ( kind = 8 ) SAMPLE_X(ACID_NUM,SAMPLE_NUM), the sampled data.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) sample_num

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  real ( kind = 8 ) mixture_mean(acid_num)
  real ( kind = 8 ) sample_mean(acid_num)
  real ( kind = 8 ) sample_variance(acid_num)
  real ( kind = 8 ) sample_x(acid_num,sample_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE_ANALYZE analyzes the samples of'
  write ( *, '(a)' ) '    the Dirichlet mixture distribution;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of samples taken = ', sample_num

  call r8row_mean ( acid_num, acid_num, sample_num, sample_x, sample_mean )

  call r8row_variance ( acid_num, acid_num, sample_num, sample_x, &
    sample_variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Acid  Expected mean, Observed mean'
  write ( *, '(a)' ) ' '

  do acid_i = 1, acid_num
    write ( *, '(i2,2x,a1,2x,2g14.6)' ) acid_i, acid_sym(acid_i), &
      mixture_mean(acid_i), sample_mean(acid_i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Acid  Observed variance'
  write ( *, '(a)' ) ' '

  do acid_i = 1, acid_num
    write ( *, '(i2,2x,a1,2x,2g14.6)' ) acid_i, acid_sym(acid_i), &
      sample_variance(acid_i)
  end do

  return
end
subroutine sample_project1 ( acid_num, comp_max, comp_num, comp_mean, &
  comp_variance, sample_comp, sample_num, sample_x )

!*****************************************************************************80
!
!! SAMPLE_PROJECT1 computes the projection of a sample onto mixture components.
!
!  Discussion:
!
!    This method used the distance from the sample to each of the means,
!    normalized by the variance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the
!    Dirichlet mixture.
!
!    Input, integer ( kind = 4 ) SAMPLE_COMP(SAMPLE_NUM), the index of the
!    mixture component that generated each sample.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of samples to take.
!
!    Output, real ( kind = 8 ) SAMPLE_X(ACID_NUM,SAMPLE_NUM), the sampled data.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max
  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) sample_num

  integer ( kind = 4 ) acid_i
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_mean(comp_max,acid_num)
  real ( kind = 8 ) comp_variance(comp_max,acid_num)
  real ( kind = 8 ) cv
  real ( kind = 8 ) cx
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_tot
  integer ( kind = 4 ) sample_comp(sample_num)
  real ( kind = 8 ) sample_comp_dist(comp_num)
  integer ( kind = 4 ) sample_i
  real ( kind = 8 ) sample_x(acid_num,sample_num)
  real ( kind = 8 ) sx
!
!  For each sample instantiation:
!
  do sample_i = 1, min ( sample_num, 10 )
!
!  For each component PDF:
!
    dist_tot = 0.0D+00

    do comp_i = 1, comp_num

      dist = 0.0D+00
      do acid_i = 1, acid_num
        sx = sample_x(acid_i,sample_i)
        cx = comp_mean(comp_i,acid_i)
        cv = comp_variance(comp_i,acid_i)
        dist = dist + ( sx - cx )**2 / cv
      end do

      sample_comp_dist(comp_i) = sqrt ( dist )

      dist_tot = dist_tot + sqrt ( dist )

    end do

    do comp_i = 1, comp_num
      sample_comp_dist(comp_i) = sample_comp_dist(comp_i) / dist_tot
    end do

    write ( *, '(i8)' ) sample_comp(sample_i)
    write ( *, '(9f8.4)' ) sample_comp_dist(1:comp_num)

  end do

  return
end
subroutine sample_project ( acid_num, alpha, comp_max, comp_num, sample_comp, &
  sample_num, sample_x )

!*****************************************************************************80
!
!! SAMPLE_PROJECT computes the projection of a sample onto mixture components.
!
!  Discussion:
!
!    This method uses the relative probabilities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, integer ( kind = 4 ) COMP_MAX, the maximum number of Dirichlet
!    mixture components.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    Input, integer ( kind = 4 ) SAMPLE_COMP(SAMPLE_NUM), the index of the 
!    mixture component that generated each sample.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of samples to take.
!
!    Output, real ( kind = 8 ) SAMPLE_X(ACID_NUM,SAMPLE_NUM), the sampled data.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max
  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) sample_num

  real ( kind = 8 ) a_prod
  real ( kind = 8 ) a_sum
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  real ( kind = 8 ) a(acid_num)
  integer ( kind = 4 ) acid_i
  real ( kind = 8 ) alpha(acid_num,comp_max)
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) dist_tot
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) sample_comp(sample_num)
  real ( kind = 8 ) sample_comp_dist(comp_num)
  integer ( kind = 4 ) sample_i
  real ( kind = 8 ) sample_x(acid_num,sample_num)
  real ( kind = 8 ) x(acid_num)
!
!  For each sample instantiation:
!
  do sample_i = 1, min ( sample_num, 1 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'X:'
    write ( *, '(a)' ) ' '
    do acid_i = 1, acid_num
      x(acid_i) = sample_x(acid_i,sample_i)
      write ( *, '(g14.6)' ) x(acid_i)
    end do
!
!  For each component PDF:
!
    dist_tot = 0.0D+00

    do comp_i = 1, comp_num

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'ALPHAs for Component ', comp_i
      write ( *, '(a)' ) ' '
      do acid_i = 1, acid_num
        a(acid_i) = alpha(acid_i,comp_i)
        write ( *, '(g14.6)' ) a(acid_i)
      end do
  
      a_sum = sum ( a(1:acid_num) )

      a_prod = 1.0D+00
      do i = 1, acid_num
        a_prod = a_prod * gamma ( a(i) )
      end do

      pdf = gamma ( a_sum ) / a_prod
      do i = 1, acid_num
        pdf = pdf * x(i)**a(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 'PDF = ', pdf

      sample_comp_dist(comp_i) = pdf

      dist_tot = dist_tot + pdf

    end do

    write ( *, '(a,g14.6)' ) 'DIST_TOT = ', dist_tot

    do comp_i = 1, comp_num
      sample_comp_dist(comp_i) = sample_comp_dist(comp_i) / dist_tot
    end do

    write ( *, '(i8)' ) sample_comp(sample_i)
    write ( *, '(9f8.4)' ) sample_comp_dist(1:comp_num)

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine weight_init ( comp_weight_post, comp_num )

!*****************************************************************************80
!
!! WEIGHT_INIT initializes the estimated weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COMP_WEIGHT_POST(COMP_NUM), the initial weights.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the
!    Dirichlet mixture.
!
  implicit none

  integer ( kind = 4 ) comp_num

  real ( kind = 8 ) comp_weight_post(comp_num)
  integer ( kind = 4 ) comp_i

  comp_weight_post(1:comp_num) = 1.0D+00

  return
end
subroutine weight_update ( acid_num, comp_num, comp_weight_post, &
  comp_weight_prior, prob, site_num, x )

!*****************************************************************************80
!
!! WEIGHT_UPDATE updates the estimated weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ACID_NUM, the number of amino acids.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    Output, real ( kind = 8 ) COMP_WEIGHT_POST(COMP_NUM), the updated weights.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT_PRIOR(COMP_NUM), the old weights.
!
!    Input, real ( kind = 8 ) PROB(ACID_NUM,COMP_NUM), probabilities for each amino
!    acid at each component.
!
!    Input, integer ( kind = 4 ) SITE_NUM, ?
!
!    Input, integer ( kind = 4 ) X(ACID_NUM), the observed event.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_num

  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_pdf
  real ( kind = 8 ) comp_weight_post(comp_num)
  real ( kind = 8 ) comp_weight_prior(comp_num)
  real ( kind = 8 ) prob(acid_num,comp_num)
  integer ( kind = 4 ) site_num
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp(comp_num)
  integer ( kind = 4 ) x(acid_num)

  do comp_i = 1, comp_num
!
!  Compute the probability of the observed event X, 
!  supposing PDF component I was used
!
    call multinomial_pdf ( x, site_num, acid_num, prob(1,comp_i), comp_pdf )

    temp(comp_i) = comp_pdf

  end do
!
!  Normalize the relative probabilities.
!
  sum2 = sum ( temp(1:comp_num) )

  do comp_i = 1, comp_num
    comp_weight_post(comp_i) = comp_weight_prior(comp_i) + temp(comp_i) / sum2
  end do

  return
end
