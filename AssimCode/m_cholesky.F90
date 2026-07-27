module m_chol
contains


subroutine chol_inv(A, Ainv, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in)  :: A(n,n)
    real(8), intent(out) :: Ainv(n,n)

    real(8) :: tmp(n,n)
    integer :: info, i, j

    tmp = A

    call dpotrf('U', n, tmp, n, info)
    if (info /= 0) then
        print *, 'DPOTRF failed, info=', info
        stop
    end if

    call dpotri('U', n, tmp, n, info)
    if (info /= 0) then
        print *, 'DPOTRI failed, info=', info
        stop
    end if

    do i = 1, n
        do j = 1, i-1
            tmp(i,j) = tmp(j,i)
        end do
    end do

    Ainv = tmp

end subroutine chol_inv

subroutine chol_inv_pivoted(A, Ainv, n, cutoff_ratio)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in)  :: A(n,n)
    real(8), intent(in)  :: cutoff_ratio 
    real(8), intent(out) :: Ainv(n,n)

    real(8) :: tmp(n,n), work(3*n), max_diag, tol
    integer :: piv(n), info, i, j, r

    tmp = A
    Ainv = 0.0d0
    piv = 0 

    max_diag = 0.0d0
    do i = 1, n
        if (tmp(i,i) > max_diag) max_diag = tmp(i,i)
    end do
    tol = max_diag * cutoff_ratio

    call dpstrf('U', n, tmp, n, piv, r, tol, work, info)

    if (r < n) then
        tmp(r+1:n, r+1:n) = 0.0d0 
    end if

    call dtrtri('U', 'N', r, tmp, n, info)
    if (info /= 0) then
        print *, 'DTRTRI failed, info=', info
        return
    end if

    call dlauum('U', r, tmp, n, info)
    do j = 1, r
        do i = 1, r
            Ainv(piv(i), piv(j)) = tmp(min(i,j), max(i,j))
        end do
    end do

end subroutine chol_inv_pivoted

!      <><><><><><><><><> old <><><><><><><><><><>
  subroutine cholesky(Aa,nn) 
    integer, intent(in) :: nn 
    real, intent(inout), dimension(1:nn,1:nn) :: Aa 
    integer :: i, j, k 
    real :: sum=0 
!    print*, '<><><><><><><> Cholesky is working for you<><><><><><><>'
    do i=1, nn 
       do j=1, i-1 
          sum=0 
          do k=1, j-1 
             sum=sum + Aa(i,k)*Aa(j,k) 
          end do
          Aa(i,j)=(1/Aa(j,j))*(Aa(i,j)-sum)  
       end do
       sum=0 
       do k=1, i-1 
          sum=sum+(Aa(i,k))**2 
       end do
       Aa(i,i)=sqrt(Aa(i,i)-sum)
       do j=i+1, nn ! Lower triangle sqrt matrix
          Aa(i,j)=0 
       end do
    end do
!    print*, '<><><><><><><> Cholesky made good job!! <><><><><><><>'
    return  
  end subroutine cholesky

!         Exponential covariance function
  function exponential(h)
    real, intent(in) :: h
    real exponential

    if(h<0) stop 'm_cholesky exponential covariance: negative distance'
    exponential=exp(-h)
  end function exponential

!         Gaussian covariance function
  function gaussian(h)
    real, intent(in) :: h
    real gaussian

    if(h<0) stop 'm_cholesky gaussian covariance: negative distance'
    gaussian=0.999*exp(-h*h)
    if (h.eq.0.0) gaussian=1.0
  end function gaussian

!         spherical covariance function
  function spherical(h)
    real, intent(in) :: h
    real spherical

    if(h<0) stop 'm_cholesky spherical covariance: negative distance'
    if(h>1.) then
       spherical=0.
    else
       spherical=1.-1.5*h+0.5*h*h*h
    endif
  end function spherical

end module m_chol
