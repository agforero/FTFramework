
 implicit none

 complex :: a

 a=(1.,2.)
 print*,a,real(a),aimag(a)
 a=a*(0.,1.)
 print*,a

 end
