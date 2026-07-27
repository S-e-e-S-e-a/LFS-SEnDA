module mod_states
! Modelstate definition for LICOM
   use mod_dimensions
   type states
      real h0(nx,ny)                           ! Sea surface Height unit:meter
      real u(nx,ny,nz)                         ! 3-D u-velocity
      real v(nx,ny,nz)                         ! 3-D v-velocity
      real t(nx,ny,nz)                         ! 3-D Temperature
      real s(nx,ny,nz)                         ! 3-D Temperature 
   end type states
   integer, parameter ::  global_ndim = 4*nx*ny*nz+nx*ny  ! Dimension of states fuww

   type states4
      real*4 h0(nx,ny)
      real*4 u(nx,ny,nz)
      real*4 v(nx,ny,nz)
      real*4 t(nx,ny,nz)
      real*4 s(nx,ny,nz)
   end type states4


   type sub_states
      real h0
      real u(nz)                               ! 1-D u-velocity
      real v(nz)                               ! 1-D v-velocity
      real t(nz)                               ! 1-D Temperature
      real s(nz)                               ! 1-D Temperature
   end type sub_states

   type sub_states4
      real*4 h0
      real*4 u(nz)                               ! 1-D u-velocity
      real*4 v(nz)                               ! 1-D v-velocity
      real*4 t(nz)                               ! 1-D Temperature
      real*4 s(nz)                               ! 1-D Temperature
   end type sub_states4

   integer, parameter ::  local_ndim=4*nz+1    ! Dimension of sub_states fuww


! Overloaded and generic operators
!   interface operator(+)
!!      module procedure add_states
!   end interface

 !  interface operator(-)
!      module procedure subtract_states
!   end interface

!   interface operator(*)
!      module procedure states_real_mult,&
!                       real_states_mult,&
!                       states_states_mult
!   end interface

!   interface operator(/)
!      module procedure divide_states
!   end interface

!   interface assignment(=)
!      module procedure assign_states
!      module procedure states4to8
!      module procedure states8to4
!   end interface

contains

   subroutine states4to8(A,B)
      type(states), intent(out) :: A
      type(states4), intent(in)  :: B
      A%h0=DBLE(B%h0)
      A%u=DBLE(B%u)
      A%v=DBLE(B%v)
      A%t=DBLE(B%t)
      A%s=DBLE(B%s)
   end subroutine states4to8

   subroutine states8to4(A,B)
      type(states), intent(in)  :: B
      type(states4),  intent(out) :: A
      A%h0=real(B%h0)
      A%u=real(B%u)
      A%v=real(B%v)
      A%t=real(B%t)
      A%s=real(B%s)
   end subroutine states8to4

   subroutine substates4to8(A,B)
      type(sub_states), intent(out) :: A
      type(sub_states4), intent(in)  :: B
      A%h0=DBLE(B%h0)
      A%u=DBLE(B%u)
      A%v=DBLE(B%v)
      A%t=DBLE(B%t)
      A%s=DBLE(B%s)
   end subroutine substates4to8

   subroutine substates8to4(A,B)
      type(sub_states4), intent(out) :: A
      type(sub_states), intent(in)  :: B
      A%h0=real(B%h0)
      A%u=real(B%u)
      A%v=real(B%v)
      A%t=real(B%t)
      A%s=real(B%s)
   end subroutine substates8to4



   type (sub_states) function getA(A,i,j)
      implicit none
      type(states), intent(in)     :: A
      integer, intent(in) :: i,j
      getA%h0 = A%h0(i,j)
      getA%u(:)=A%u(i,j,:)
      getA%v(:)=A%v(i,j,:)
      getA%t(:)=A%t(i,j,:)
      getA%s(:)=A%s(i,j,:)
   end function getA

   type (sub_states4) function getA4(A,i,j)
      implicit none
      type(states4), intent(in)     :: A
      integer, intent(in) :: i,j
      getA4%h0 = A%h0(i,j)
      getA4%u(:)=A%u(i,j,:)
      getA4%v(:)=A%v(i,j,:)
      getA4%t(:)=A%t(i,j,:)
      getA4%s(:)=A%s(i,j,:)
   end function getA4

   subroutine putA(subA,A,i,j)
      implicit none
      type(sub_states), intent(in) :: subA
      type(states), intent(inout)   :: A
      integer, intent(in) :: i,j
      A%h0(i,j) =subA%h0
      A%u(i,j,:)=subA%u(:)
      A%v(i,j,:)=subA%v(:)
      A%t(i,j,:)=subA%t(:)
      A%s(i,j,:)=subA%s(:)
   end subroutine putA



!   function add_states(A,B)
!      type(states) add_states
!      type(states), intent(in) :: A
!      type(states), intent(in) :: B
!       add_states%h0 = A%h0 + B%h0
!       add_states%u = A%u + B%u
!       add_states%v = A%v + B%v
!       add_states%t = A%t + B%t
!       add_states%s = A%s + B%s
!   end function add_states


   function subtract_states8(A,B) result(C)
      type(sub_states) C
      type(sub_states), intent(in) :: A
      type(sub_states), intent(in) :: B
       C%h0 = A%h0 - B%h0
       C%u = A%u - B%u
       C%v = A%v - B%v
       C%t = A%t - B%t
       C%s = A%s - B%s
   end function subtract_states8

   function add_states8(A,B) result(C)
      type(sub_states) C
      type(sub_states), intent(in) :: A
      type(sub_states), intent(in) :: B
       C%h0 = A%h0 + B%h0
       C%u = A%u + B%u
       C%v = A%v + B%v
       C%t = A%t + B%t
       C%s = A%s + B%s
   end function add_states8

   function multi_states8(A, scalar) result(C)
      type(sub_states) :: C
      type(sub_states), intent(in) :: A
      real*8, intent(in) :: scalar
       C%h0 = A%h0 * scalar
       C%u  = A%u  * scalar
       C%v  = A%v  * scalar
       C%t  = A%t  * scalar
       C%s  = A%s  * scalar
   end function multi_states8

   function assign_states8(r) result(A)
    type(sub_states) :: A
    real*8, intent(in) :: r

    A%h0 = r
    A%u = r
    A%v = r
    A%t = r
    A%s = r
   end function assign_states8

   function subtract_states(A,B) result(C)
      type(sub_states4) C
      type(sub_states4), intent(in) :: A
      type(sub_states4), intent(in) :: B
       C%h0 = A%h0 - B%h0
       C%u = A%u - B%u
       C%v = A%v - B%v
       C%t = A%t - B%t
       C%s = A%s - B%s
   end function subtract_states

   function add_states(A,B) result(C)
      type(sub_states4) C
      type(sub_states4), intent(in) :: A
      type(sub_states4), intent(in) :: B
       C%h0 = A%h0 + B%h0
       C%u = A%u + B%u
       C%v = A%v + B%v
       C%t = A%t + B%t
       C%s = A%s + B%s
   end function add_states
   
   function multi_states(A, scalar) result(C)
      type(sub_states4) :: C
      type(sub_states4), intent(in) :: A
      real, intent(in) :: scalar
       C%h0 = A%h0 * scalar
       C%u  = A%u  * scalar
       C%v  = A%v  * scalar
       C%t  = A%t  * scalar
       C%s  = A%s  * scalar
   end function multi_states

   function assign_states(r) result(A)
    type(sub_states4) :: A
    real, intent(in) :: r
    
    A%h0 = r
    A%u = r
    A%v = r
    A%t = r
    A%s = r
   end function assign_states

  function nan_states(A, missing_val) result(B)
    type(sub_states4), intent(in) :: A
    real, intent(in) :: missing_val
    type(sub_states4) :: B
    
    B = A
    if (B%h0 > missing_val) B%h0 = 0.0
    where (B%u > missing_val)
        B%u = 0.0
    end where
    where (B%v > missing_val)
        B%v = 0.0
    end where
    where (B%t > missing_val)
        B%t = 0.0
    end where
    where (B%s > missing_val)
        B%s = 0.0
    end where
   end function nan_states
!   function states_real_mult(A,B)
!      type(states) states_real_mult
!      type(states), intent(in) :: A
!      real, intent(in) :: B
!       states_real_mult%h0 = B*A%h0
!       states_real_mult%u = B*A%u
!       states_real_mult%v = B*A%v
!       states_real_mult%t = B*A%t
!       states_real_mult%s = B*A%s
!   end function states_real_mult
!
!   function real_states_mult(B,A)
!      type(states) real_states_mult
!      type(states), intent(in) :: A
!      real, intent(in) :: B
!       real_states_mult%h0 = B*A%h0
!       real_states_mult%u = B*A%u
!       real_states_mult%v = B*A%v
!       real_states_mult%t = B*A%t
!       real_states_mult%s = B*A%s
!   end function real_states_mult
!
!   function states_states_mult(A,B)
!      type(states) states_states_mult
!      type(states), intent(in) :: A
!      type(states), intent(in) :: B
!       states_states_mult%h0 = A%h0 * B%h0
!       states_states_mult%u = A%u * B%u
!       states_states_mult%v = A%v * B%v
!       states_states_mult%t = A%t * B%t
!       states_states_mult%s = A%s * B%s
!   end function states_states_mult
!
!
!   subroutine assign_states(A,r)
!      type(states), intent(out) :: A
!      real, intent(in) :: r
!       A%h0= r
!       A%u = r
!       A%v = r
!       A%t = r
!       A%s = r
!   end subroutine assign_states

end module mod_states

