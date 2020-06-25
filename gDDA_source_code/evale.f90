 subroutine evale(cxE00,akd,dx,x0,ixyz0,iocc,flgwav,xyzc0,w0,mxnat,mxn3,nat,nat0,nx,ny,nz,cxE,CFLAVG)
   use ddprecision,only : WP
   use ddcommon_0,only: CFLPAR
   implicit none

   complex(WP),parameter :: cxi=(0.0_WP,1.0_WP)
!*** arguments:
   integer,intent(in) :: mxn3,mxnat,nat,nat0,nx,ny,nz
   integer,intent(in),dimension(nat0,3) :: ixyz0
   integer(kind=2),intent(in),dimension(mxnat) :: iocc
   real(WP),intent(in),dimension(3) :: akd,dx,x0
   complex(WP),intent(in),dimension(3) ::cxE00
   complex(WP),intent(out),dimension(nat,3) :: cxE
   ! new input arguments defining the gaussian beam
   integer,intent(in)::flgwav              ! Option for wavefront: 0--planewave; 1--gaussian
   real(WP),intent(in),dimension(3)::xyzc0 ! Gaussian beam center; in unit of d (dipole spacing)
   real(WP),intent(in)::w0                 ! Gaussian beam waist; in unit of d (dipole spacing)
   character(len=12)::CFLAVG
! note: cxE should be dimensioned to cxE(nat,3) in this routine
!       so that first 3*nat elements of cxE are employed.

!***  local variables:
   complex(WP) :: cxfac
   real(WP) :: x,x1,x2
   integer :: ia,ix,iy,iz,m
   integer,dimension(3) :: shift
!***********************************************************************
! subroutine evale

! given:   cxE00(1-3)=incident E field at origin (complex) at t=0
!          akd(1-3)=(kx,ky,kz)*d for incident wave (d=effective
!                    lattice spacing)
!          dx(1-3)=(dx/d,dy/d,dz/d) for lattice (dx,dy,dz=lattice
!                   spacings in x,y,z directions, d=(dx*dy*dz)**(1/3)
!          x0(1-3)=(x,y,z)location/(d*dx(1-3)) in TF of lattice site
!                  with ix=0,iy=0,iz=0
!          ixyz0(1-nat0,3)=[x-x0(1)]/dx, [y-x0(2)]/dy, [z-x0(3)]/dz
!                  for each of nat0 physical dipoles
!          mxnat,mxn3=dimensioning information
!          nat0=number of dipoles in physical target
!          nat=number of locations at which to calculate cxE

! returns: cxE(1-nat,3)=incident E field at nat locations at t=0

! B.T.Draine, Princeton Univ. Obs., 88.05.09
! History:
! 90.11.06 (btd): modified to pass array dimension.
! 90.11.29 (btd): modified to allow option for either
!                   physical locations only (nat=nat0), or
!                   extended dipole array (nat>nat0)
! 90.11.30 (btd): corrected error for case nat>nat0
! 90.11.30 (btd): corrected another error for case nat>nat0
! 90.12.03 (btd): change ordering of xyz0 and cxE
! 90.12.05 (btd): corrected error in dimensioning of cxE
! 90.12.10 (btd): remove xyz0, replace with ixyz0
! 97.11.02 (btd): add dx to argument list to allow use with
!                 noncubic lattices.
! 07.06.20 (btd): add x0 to the argument list to specify location
!                 in tf corresponding to ix=0,iy=0,iz=0
! 07.09.11 (btd): changed ixyz0 from integer*2 to integer
! 08.03.14 (btd): v7.05
!                 corrected dimensioning
!                 ixyz0(mxnat,3) -> ixyzo(nat0,3)
! copyright (c) 1993,1997,2007 b.t. draine and p.j. flatau
! this code is covered by the gnu general public license.
!***********************************************************************

! evaluate electric field vector at each dipole location.

! if nat=nat0, then evaluate e only at occupied sites.
! if nat>nat0, then evaluate e at all sites.

    open(unit=2332,file='temp-shift.txt',status='old',action='read')
        read(2332,*) shift
    close(2332)
    open(321,file='Einc_'//CFLAVG(1:4)//'_'//CFLPAR,status='unknown',action='write') ! for cxE
       if(nat==nat0)then
          do ia=1,nat0
             x=0._wp
             do m=1,3
                x=x+akd(m)*dx(m)*(real(ixyz0(ia,m),kind=WP)+x0(m))
             end do
             if (flgwav==0) then
                 cxfac=exp(cxi*x)
             else if (flgwav == 1) then
                 call gaussian_beam(akd,w0,xyzc0,real(ixyz0(ia,1),kind=WP), &
                  real(ixyz0(ia,2),kind=WP),real(ixyz0(ia,3),kind=WP),cxfac)
             end if
             do m=1,3
                cxE(ia,m)=cxE00(m)*cxfac
             end do
             if (iocc(ia)==1) then
                 write(321,'(3i7,6es20.6e3)') ixyz0(ia,1:3)-shift,cxE(ia,1:3) ! write cxE to file
             end if
          end do
       else
          ia=0
          do iz=1,nz
             if (flgwav==0) then
                 x1=akd(3)*dx(3)*(real(iz,kind=wp)+x0(3))
             end if
             do iy=1,ny
                if (flgwav==0) then
                    x2=x1+akd(2)*dx(2)*(real(iy,kind=wp)+x0(2))
                end if
                do ix=1,nx
                   ia=ia+1
                   if (flgwav==0) then
                       x=x2+akd(1)*dx(1)*(real(ix,kind=wp)+x0(1))
                       cxfac=exp(cxi*x)
                   else if (flgwav==1) then
                       call gaussian_beam(akd,w0,xyzc0,real(ix,kind=wp)+x0(1),&
                          real(iy,kind=wp)+x0(2),real(iz,kind=wp)+x0(3),cxfac)
                   end if
                   do m=1,3
                      cxE(ia,m)=cxE00(m)*cxfac
                   end do
                   if (iocc(ia)==1) then
                      write(321,'(3i7,6es20.6e3)') ix-shift(1),iy-shift(2),iz-shift(3),cxE(ia,1:3) ! write cxE to file
                   end if
                end do
             end do
          end do
       end if
    close(321)
 end subroutine evale

 subroutine gaussian_beam(k0,w0,xyzc0,x0,y0,z0,E)
   use ddprecision,only:WP
   implicit none
   
   complex(WP),parameter::cxI=(0.0_WP,1.0_WP)
   
   real(WP),intent(in),dimension(3)::k0 ! wavevector = 2*PI/lambda
   real(WP),intent(in)::w0 ! waist
   real(WP),intent(in),dimension(3)::xyzc0 ! center of Gaussian beam
   real(WP),intent(in)::x0,y0,z0 ! point coordinates at target
   complex(WP),intent(out)::E  ! output complex electric field [1]
   
   real(WP)::r,zR,wz,invRz,phi
   real(WP)::xc0,yc0,zc0
   real(WP)::theta_k0,phi_k0,k,kx,ky,kz,xc,yc,zc,x,y,z
 
   theta_k0=atan(sqrt(k0(1)**2+k0(2)**2)/k0(3))
   phi_k0=atan2(k0(2),k0(1))
   
   xc0=xyzc0(1); yc0=xyzc0(2); zc0=xyzc0(3)
   xc=xc0*cos(theta_k0)*cos(phi_k0)+yc0*cos(theta_k0)*sin(phi_k0)-zc0*sin(theta_k0)
   yc=-xc0*sin(phi_k0)+yc0*cos(phi_k0)
   zc=xc0*sin(theta_k0)*cos(phi_k0)+yc0*sin(theta_k0)*sin(phi_k0)+zc0*cos(theta_k0)
 
   x=x0*cos(theta_k0)*cos(phi_k0)+y0*cos(theta_k0)*sin(phi_k0)-z0*sin(theta_k0)
   y=-x0*sin(phi_k0)+y0*cos(phi_k0)
   z=x0*sin(theta_k0)*cos(phi_k0)+y0*sin(theta_k0)*sin(phi_k0)+z0*cos(theta_k0)
 
   kx=k0(1)*cos(theta_k0)*cos(phi_k0)+k0(2)*cos(theta_k0)*sin(phi_k0)-k0(3)*sin(theta_k0)
   ky=-k0(1)*sin(phi_k0)+k0(2)*cos(phi_k0)
   kz=k0(1)*sin(theta_k0)*cos(phi_k0)+k0(2)*sin(theta_k0)*sin(phi_k0)+k0(3)*cos(theta_k0)
 
   k=kz
   r=sqrt((x-xc)**2+(y-yc)**2)
   zR=k*w0**2/2
   wz=w0*sqrt(1+((z-zc)/zR)**2)
   invRz=(z-zc)/((z-zc)**2+zR**2)
   phi=atan((z-zc)/zR)
   
   E=w0/wz*exp(-r**2/wz**2)*exp(cxI*(k*(z-zc)+k*r**2/2*invRz-phi))
 end subroutine gaussian_beam
