module m_modsubstate_point
contains

      real function modsubstate_point(cvar,mem,meanssh,ix,iy,iz,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
      use mod_dimensions
      use mod_states
!      use mod_meanssh
      implicit none

      integer         , intent(in) :: nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2 
      character(len=*), intent(in) :: cvar
      type(sub_states)    , intent(in) :: mem(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2)
      integer         , intent(in) :: ix,iy,iz
      real            ,  intent(in):: meanssh(nx,ny)
      integer ixx,iyy

      if(ix-nx_begin+prep_rx>nx) then
        ixx=ix-nx_begin+prep_rx-nx
      else if(ix-nx_begin+prep_rx<1) then
         ixx=ix-nx_begin+prep_rx+nx
      else
        ixx=ix-nx_begin+prep_rx
      endif
      
      iyy=iy-ny_begin+ry1 
    
      select case (trim(cvar))
!#if defined (ICE) 
!      case ('ficem')
!         modstate_point=mstate%fice
!      case ('hicem')
!         modstate_point=mstate%hice
!#endif
      case ('ut')
         modsubstate_point=mem(ixx,iyy)%u(iz)
      case ('vt')
         modsubstate_point=mem(ixx,iyy)%v(iz)
      case ('TEM')
         modsubstate_point=mem(ixx,iyy)%t(iz)
      case ('SAL')
         modsubstate_point=mem(ixx,iyy)%s(iz)
      case ('SST')
        modsubstate_point=mem(ixx,iyy)%t(1)
      case ('SSH')
         modsubstate_point=mem(ixx,iyy)%h0 - meanssh(ix,iy)
      case ('SLA')
         modsubstate_point=mem(ixx,iyy)%h0 - meanssh(ix,iy)
      case default
         print *,'No match in modstate_point'
         print *,cvar
         stop
      end select
   end function modsubstate_point

end module m_modsubstate_point
