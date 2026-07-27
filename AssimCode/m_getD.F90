module m_getD
contains
subroutine getDvec(D,subD,nrobs,lobs,nobs)
! Returns the subD vector of active measurements
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nobs
   real*4,    intent(in)  :: D(nrobs)
   logical, intent(in)  :: lobs(nrobs)
   real*4,    intent(out) :: subD(nobs)

   integer j,m

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         j=j+1
         subD(j)=D(m)
      endif
   enddo
end subroutine getDvec

subroutine getDvecnew(D,subD,nrobs,lobs,wlobs,nobs)
! Returns the subD vector of active measurements
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nobs
   real*4,    intent(in)  :: D(nrobs)
   logical, intent(in)  :: lobs(nrobs)
   real*4,    intent(in)  :: wlobs(nrobs)
   real*4,    intent(out) :: subD(nobs)

   integer j,m

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         j=j+1
         subD(j)=D(m)*wlobs(m)
!lzt
!         print*,'D(m),wlobs(m),subD(j)=',D(m),wlobs(m),subD(j)
      endif
   enddo
end subroutine getDvecnew

subroutine getS(S_local,subS,S_k,nrobs,nrens,lobs,nobs,h,k) ! zyz Incorrect!
! Returns the subD matrix corresponding to active measurements
   implicit none
   integer, intent(in)  :: nrobs,h,k,S_k
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nobs
   real*4,    intent(in)  :: S_local(nrens,nrobs)
   logical, intent(in)  :: lobs(nrobs)
   real*4,    intent(out) :: subS(nobs,nrens)

   integer j,m,n,a

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         a=0
         do n=1,S_k
            if(S_local(nrens+1,n)==m) then
                a=1
                j=j+1
                subS(j,:)=S_local(1:nrens,n)
            endif
         enddo
         if(a==0) print*,m,'is not found,i=,j=',h,k
      endif
   enddo
end subroutine getS

subroutine getD(D,subD,nrobs,nrens,lobs,nobs)
! Returns the subD matrix corresponding to active measurements
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nobs
   real*4,    intent(in)  :: D(nrens,nrobs)
   logical, intent(in)  :: lobs(nrobs)
   real*4,    intent(out) :: subD(nobs,nrens)

   integer j,m

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         j=j+1
         subD(j,:)=D(:,m)
      endif
   enddo

end subroutine getD

subroutine getD4S(D,subD,nrobs,nrens,lobs,wlobs,nobs)
! Returns the subD matrix corresponding to active measurements
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nobs
   real*4,    intent(in)  :: D(nrobs,nrens)
   logical, intent(in)  :: lobs(nrobs)
   real*4,    intent(in)  :: wlobs(nrobs)
   real*4,    intent(out) :: subD(nobs,nrens)

   integer j,m

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         j=j+1
         subD(j,:)=D(m,:)*wlobs(m)
      endif
   enddo

end subroutine getD4S
end module m_getD
