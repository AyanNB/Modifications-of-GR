module tovmod1
  implicit none

  integer :: nzones = 0
  integer :: isurf = 0
  real*8  :: rho_c = 1.0d0,rho_e=1.0d0
  real*8  :: rmax  = 1.0d0
  real*8  :: dr = 0.0d0
  real*8  :: rho_crit=1.0d0

  real*8,allocatable  :: rs(:)
  real*8,allocatable  :: ri(:)
  real*8,allocatable  :: rho(:),rho_cit(:),rho_eff(:)
  real*8,allocatable  :: press(:), press_e(:)
  real*8,allocatable  :: mgrav(:),mgrav_e(:),mrest(:)

  

  ! other thermodynamics vars



  ! some constants
  real*8,parameter    :: ggrav = 6.673d-8
  real*8,parameter    :: clite = 2.99792458d10
  real*8,parameter    :: pi = 3.14159265358979d0
  real*8,parameter    :: g_c2 = ggrav/clite**2
  real*8,parameter    :: msun = 1.98892d33

  real*8    :: a(1:24),ksi(1:7),ksi_e(1:7) 




contains 
  subroutine allocate
    implicit none

    allocate(rs(nzones))
    allocate(ri(nzones))
    allocate(rho(nzones))
    allocate(press(nzones))
    allocate(press_e(nzones))
    allocate(mgrav(nzones))
    allocate(mgrav_e(nzones))
    allocate(mrest(nzones))

  end subroutine allocate

  subroutine initialize
    implicit none

    rs     = 0.0d0
    ri     = 0.0d0
    rho    = 0.0d0
    press  = 0.0d0
    mgrav  = 0.0d0
    press_e= 0.0d0
    mgrav_e= 0.0d0

  end subroutine initialize

end module tovmod1



program tov
  use tovmod1
  implicit none

  ! number of zones to use
  nzones = 4000
  call allocate

  ! set up outer radius and zone size
  rmax  = 1.5d6
  dr = rmax/dble(nzones-1)

  call make_a_tov_star

end program tov

subroutine make_a_tov_star

  use tovmod1
  implicit none
  integer :: g,n,k,l
  real*8,allocatable :: rb(:),M(:),M_e(:)

  n=50
  allocate(rb(n))
  allocate(M(n))
  allocate(M_e(n))
  allocate(rho_cit(n))
  allocate(rho_eff(n))

  a(1)=4.078d0
  a(2)=7.587d0
  a(3)=0.00839d0
  a(4)=0.21695d0
  a(5)=3.614d0
  a(6)=11.942d0
  a(7)=13.751d0
  a(8)=1.3373d0
  a(9)=3.606d0
  a(10)=a(6)
  a(11)=-22.996d0
  a(12)=1.6229d0
  a(13)=4.88d0
  a(14)=14.274d0
  a(15)=23.560d0
  a(16)=-1.5564d0
  a(17)=2.095d0
  a(18)=15.294d0
  a(19)=0.084d0
  a(20)=6.36d0
  a(21)=11.67d0
  a(22)=-0.042d0
  a(23)=14.8d0
  a(24)=14.18d0

 
  do k=1,n
     rho_cit(k) =2.5d14 + k*0.5d14
     rho_crit=1.0d17
     rho_eff(k) = rho_cit(k) - rho_cit(k)**2/rho_crit
  ! initialize everything to 0
     call initialize

     g=1
     

     call tov_integrate(k)
     do while (press_e(g)>1.0d25)
          g=g+1
          isurf=g-1
     end do

     

     rb(k)=rs(isurf)
     M(k)=mgrav(isurf)  

  enddo

  open(unit=444,file="Mrb3",status='unknown')
  do l=1,n
    write(444,"(i8,1P10E15.6)") l,rb(l)/1.0d5,M(l)/msun,rho_eff(l)
  enddo
  close(444) 
  
        


end subroutine make_a_tov_star




subroutine tov_integrate(k)
  use tovmod1
  implicit none

  ! local variables
  integer :: k 
  integer :: i,l,p,t
  real*8  :: logrho,logpress,logrho_e,logpress_e
  real*8  :: old_data(5)
  real*8  :: new_data(5)
  real*8 :: min_press = 10.0d0
  rho_c=rho_cit(k)
  rho_e=rho_eff(k)
    
  logrho=log10(rho_c)
  logrho_e=log10(rho_e)
  ksi(1)=(a(1)+a(2)*logrho+a(3)*logrho**3)/(1+a(4)*logrho)/(exp(a(5)*(logrho-a(6)))+1)
  do i=7,15,4
    ksi((i+1)/4)=(a(i)+a(i+1)*logrho)/(exp(a(i+2)*(a(i+3)-logrho))+1d0)
  enddo
  do i=19,22,3
     ksi((i-4)/3)=a(i)/(1+(a(i+1)*(logrho-a(i+2)))**2)
  enddo
  ksi(7)=ksi(1)+ksi(2)+ksi(3)+ksi(4)+ksi(5)+ksi(6) 
  
  ksi_e(1)=(a(1)+a(2)*logrho_e+a(3)*logrho_e**3)/(1+a(4)*logrho_e)/(exp(a(5)*(logrho_e-a(6)))+1)
  do i=7,15,4
    ksi_e((i+1)/4)=(a(i)+a(i+1)*logrho_e)/(exp(a(i+2)*(a(i+3)-logrho_e))+1d0)
  enddo
  do i=19,22,3
     ksi_e((i-4)/3)=a(i)/(1+(a(i+1)*(logrho_e-a(i+2)))**2)
  enddo
  ksi_e(7)=ksi_e(1)+ksi_e(2)+ksi_e(3)+ksi_e(4)+ksi_e(5)+ksi_e(6)
  
  press(1) = 10**ksi(7) 
  press_e(1) = 10**ksi_e(7)

  rho(1) = rho_e
  rs(1) = 0.01
  mgrav(1) = (4*pi*rho_c*rs(1)**3)/3.0d0
  mgrav_e(1) = (4*pi*rho_e*rs(1)**3)/3.0d0
  mrest(1) = (4*pi*rho_c*rs(1)**3)/3.0d0/sqrt(1-2*g_c2*mgrav(1)/rs(1))
  do i=2,nzones
     rs(i) = rs(i-1) + dr
  enddo


  
  do i=2,nzones
     old_data(1) = max(press_e(i-1),min_press)
     old_data(2) = mgrav_e(i-1)
     old_data(3) = mgrav(i-1)
     old_data(4) = mrest(i-1)
     call tov_RK3(old_data,new_data,rs(i-1))     
     
     press_e(i) = new_data(1)
     mgrav_e(i) = new_data(2)
     mgrav(i) = new_data(3)
     mrest(i) = new_data(4)

     

     
  enddo
   
end subroutine tov_integrate


subroutine tov_RK3(old_data,new_data,rad)

  implicit none

  real*8, intent(IN) :: old_data(5)
  real*8, intent(OUT):: new_data(5)

  real*8, intent(IN) :: rad
  real*8  :: source_data(5)
  real*8  :: tmp_data(5)
  real*8  :: ktemp(5,3)
  real*8  :: dr    
  real*8  :: rmax
  integer :: nzones
  
  nzones=4000
  rmax  = 1.5d6

  dr = rmax/dble(nzones-1)  
  ktemp(1:5,1:3) = 0.0d0

  tmp_data(1:5) = old_data(1:5)
  call tov_RHS(tmp_data,source_data,rad)
  ktemp(1:5,1) = dr * source_data(1:5)

  tmp_data(1:5) = old_data(1:5) + 0.5d0*ktemp(1:5,1)
  call tov_RHS(tmp_data,source_data,rad + 0.5d0*dr)
  ktemp(1:5,2) = dr * source_data(1:5)

  tmp_data(1:5) = old_data(1:5) - ktemp(1:5,1) &
       + 3.0d0*ktemp(1:5,2)
  call tov_RHS(tmp_data,source_data,rad + dr)
  ktemp(1:5,3) = dr * source_data(1:5)

  new_data(1:5) = old_data(1:5) + 1.0d0/6.0d0 * &
       (ktemp(1:5,1) + 4.0d0*ktemp(1:5,2) + ktemp(1:5,3))

end subroutine tov_RK3

subroutine tov_RHS(old_data,source_data,r)
  use tovmod1, only: rho_crit, pi, g_c2, ggrav, &
       clite
  implicit none
  
  real*8, intent(IN)  :: r
  real*8, intent(IN)  :: old_data(5)
  real*8, intent(OUT) :: source_data(5)

  ! local variables
  real*8 :: mgrav,mgrav_e,press
  real*8 :: rho,rho_e,press_e
  real*8 :: min_press = 10.0d0

! define the RHS for the integration from n to n+1
! Store results in indices 1 and 2 of array source_data
! note that you can use the constants defined in tovmod!

  ! input data, assign for convenience:

  ! this is to make sure, we don't have neg. press.
  press_e = max(old_data(1),min_press)
  mgrav_e = old_data(2)
  mgrav = old_data(3)
  
 

  call newraph(press_e,rho_e) 
  rho = rho_crit/2*(1-(1-4*rho_e/rho_crit)**0.5)

  
  source_data(1) = ggrav*(rho_e+press_e/clite**2)*(mgrav_e+4*pi*(r**3)*press_e/clite**2)/(r*(2*g_c2*mgrav_e-r))
  source_data(2) = 4*pi*(r**2)*rho_e
  source_data(3) = 4*pi*(r**2)*rho
  source_data(4) = 4*pi*(r**2)*rho/sqrt(1-2*g_c2*mgrav/r)
end subroutine tov_RHS

subroutine newraph(press,rho)
  use tovmod1, only: a,ksi
  implicit none
  real*8, intent(IN)  :: press
  real*8, intent(OUT) :: rho
  real*8 :: b,c,t,m,l,d,f
  real*8 :: psi(1:6),tsi(1:6)
  real*8, allocatable :: x(:)
  real*8, allocatable :: eqn(:),deqn(:)
  integer   :: i,n,k  

  n=50
  allocate(x(n))
  allocate(eqn(n))
  allocate(deqn(n))


  x = 0.0d0 
  eqn = 0.0d0
  deqn = 0.0d0
 
  x(1) =log10(1.0d15)
  
  do i=1,n

      m = (a(1)+a(2)*x(i)+a(3)*x(i)**3)  
      t = exp(a(5)*(x(i)-a(6)))+1.0d0
      l = (1+a(4)*x(i))
      psi(1) = m/l/t
      tsi(1)=(m*(-a(4)*t-a(5)*l*(t-1.0d0))+t*l*(a(2)+3.0d0*a(3)*x(i)**2))/(l*t)**2
      do k=7,15,4
         b=(a(k)+a(k+1)*x(i))
         c=(exp(a(k+2)*(a(k+3)-x(i)))+1d0)
         psi((k+1)/4)=b/c
         tsi((k+1)/4)=a(k+2)*b*(c-1)/c**2 + a(k+1)/c
      enddo
      do k=19,22,3
         f=(x(i)-a(k+2))*a(k+1)
         d=1+(f)**2
         psi((k-4)/3)=a(k)/d
         tsi((k-4)/3)=(2.0d0*a(k)*a(k+1)*f)/d**2
      enddo
      do k=1,6
         eqn(i) = eqn(i)+psi(k)
      enddo
      eqn(i)=eqn(i)-log10(press)
      do k=1,6
         deqn(i) = deqn(i)+tsi(k)
      enddo  
       
      
    
      x(i+1)=x(i) - eqn(i)/deqn(i)

  enddo

  rho=10**x(n)  

end subroutine newraph   


