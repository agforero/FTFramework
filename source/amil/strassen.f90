! From: "James Van Buskirk" <torsop@ix.netcom.com>
! Newsgroups: comp.lang.fortran
! Subject: Re: ? Strassen algorithm for multiplying matrix
! Date: Thu, 14 Sep 2000 06:27:20 -0600

! Latest revision - 16 September 2000

module mykinds
   implicit none
   integer, parameter :: wp = selected_real_kind(12, 60)
   integer            :: threshold
end module mykinds



module matrix_multiply
   use mykinds
   implicit none

contains

recursive subroutine strassen(a, b, c, N)
   integer, intent(in)   :: N
   real(wp), intent(in)  :: a(N,N), b(N,N)
   real(wp), intent(out) :: c(N,N)

   real(wp), dimension(N/2,N/2) :: a11,a21,a12,a22,b11,b21,b12,b22
   real(wp), dimension(N/2,N/2) :: q1,q2,q3,q4,q5,q6,q7

   if(iand(N,1) /= 0 .OR. N < threshold) then
      c = matmul(a,b)
   else
      a11 = a(:N/2,:N/2)
      a21 = a(N/2+1:,:N/2)
      a12 = a(:N/2,N/2+1:)
      a22 = a(N/2+1:,N/2+1:)
      b11 = b(:N/2,:N/2)
      b21 = b(N/2+1:,:N/2)
      b12 = b(:N/2,N/2+1:)
      b22 = b(N/2+1:,N/2+1:)
      call strassen(a11+a22, b11+b22, q1, N/2)
      call strassen(a21+a22, b11, q2, N/2)
      call strassen(a11, b12-b22, q3, N/2)
      call strassen(a22, -b11+b21, q4, N/2)
      call strassen(a11+a12, b22, q5, N/2)
      call strassen(-a11+a21, b11+b12, q6, N/2)
      call strassen(a12-a22, b21+b22, q7, N/2)
      c(:N/2,:N/2) = q1+q4-q5+q7
      c(N/2+1:,:N/2) = q2+q4
      c(:N/2,N/2+1:) = q3+q5
      c(N/2+1:,N/2+1:) = q1+q3-q2+q6
   end if

   return
end subroutine strassen


recursive subroutine simple(a, b, c, N)
   use mykinds
   implicit none
   integer, intent(in)   :: N
   real(wp), intent(in)  :: a(N,N), b(N,N)
   real(wp), intent(out) :: c(N,N)
   real(wp), dimension(N/2,N/2) :: a11,a21,a12,a22,b11,b21,b12,b22,q1,q2

   if(iand(N,1) /= 0 .OR. N < threshold) then
      c = matmul(a,b)
   else
      a11 = a(:N/2,:N/2)
      a21 = a(N/2+1:,:N/2)
      a12 = a(:N/2,N/2+1:)
      a22 = a(N/2+1:,N/2+1:)
      b11 = b(:N/2,:N/2)
      b21 = b(N/2+1:,:N/2)
      b12 = b(:N/2,N/2+1:)
      b22 = b(N/2+1:,N/2+1:)
      call simple(a11, b11, q1, N/2)
      call simple(a12, b21, q2, N/2)
      c(:N/2,:N/2) = q1+q2
      call simple(a21, b11, q1, N/2)
      call simple(a22, b21, q2, N/2)
      c(N/2+1:,:N/2) = q1+q2
      call simple(a11, b12, q1, N/2)
      call simple(a12, b22, q2, N/2)
      c(:N/2,N/2+1:) = q1+q2
      call simple(a21, b12, q1, N/2)
      call simple(a22, b22, q2, N/2)
      c(N/2+1:,N/2+1:) = q1+q2
   end if

   return
end subroutine simple

end module matrix_multiply


program strassen_test
   use mykinds
   use matrix_multiply
   implicit none
   real(wp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), e(:,:)
   integer  :: N
   real     :: t0, t1

   write(*, '(a)', advance='no') ' Enter the order of the matrices:> '
   read(*,*) N
   write(*, '(a)', advance='no') ' Enter the threshold:> '
   read(*,*) threshold
   allocate(a(N,N), b(N,N), c(N,N), d(N,N), e(N,N))
   call random_seed()
   call random_number(a)
   call random_number(b)
   call cpu_time(t0)
   c = matmul(a,b)
   call cpu_time(t1)
   write(*,'(a, f10.3)') ' Time for matmul        = ', t1-t0
   call cpu_time(t0)
   call strassen(a, b, d, N)
   call cpu_time(t1)
   write(*,'(a, f10.3)') ' Time for strassen      = ', t1-t0
   write(*,'(a, g12.4)') ' RMS error for strassen = ', sqrt(sum((c-d)**2))/N
   call cpu_time(t0)
   call simple(a,b,e,N)
   call cpu_time(t1)
   write(*,'(a, f10.3)') ' Time for simple        = ', t1-t0
   write(*,'(a, g12.4)') ' RMS error for simple   = ', sqrt(sum((c-e)**2))/N
   stop
end program strassen_test

