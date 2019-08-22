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
  real*8,allocatable  :: eps(:)
  real*8,allocatable  :: mgrav(:),mgrav_e(:),mrest(:)


  ! other thermodynamics vars



  ! some constants
  real*8,parameter    :: ggrav = 6.673d-8
  real*8,parameter    :: clite = 2.99792458d10
  real*8,parameter    :: pi = 3.14159265358979d0
  real*8,parameter    :: g_c2 = ggrav/clite**2
  real*8,parameter    :: msun = 1.98892d33

  ! EOS parameters
  real*8    :: gamma = 2.5d0
  real*8    :: eosk = 0.007452876

contains 
  subroutine allocate
    implicit none

    allocate(rs(nzones))
    allocate(ri(nzones))
    allocate(rho(nzones))
    allocate(press(nzones))
    allocate(press_e(nzones))
    allocate(eps(nzones))
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
    eps    = 0.0d0
    mgrav  = 0.0d0


  end subroutine initialize

end module tovmod1



program tov
  use tovmod1
  implicit none

  ! number of zones to use
  nzones = 10000
  call allocate

  ! set up outer radius and zone size
  rmax  = 1.5d6
  dr = rmax/dble(nzones-1)

  call make_a_tov_star

end program tov

subroutine make_a_tov_star

  use tovmod1
  implicit none
  integer :: i,n,k,l,g,un=444
  real*8,allocatable :: rb(:),M(:),Mr(:),Mg(:)
  real*8 :: rho_cr(3)
  character(12)::filename
  n=50
  allocate(rb(n))
  allocate(M(n))
  allocate(Mr(n))
  allocate(Mg(n))
  allocate(rho_cit(n))
  allocate(rho_eff(n))
  do i=1,3
     rho_cr(i)=(i**2-i+2)*1.0d17
     do k=1,n
     	rho_cit(k) =0.5d14 + k*0.5d14
     	rho_crit=rho_cr(i)
     	rho_eff(k) = rho_cit(k) - rho_cit(k)**2/rho_crit
  ! initialize everything to 0
     	call initialize
     	g=1
  ! perform the integration
     	call tov_integrate(k)
     	do while (press_e(g)>10)
           g=g+1
           isurf=g-1
     	end do



     	rb(k)=rs(isurf)
     	M(k)=mgrav(isurf)  
     	Mr(k)=mrest(isurf)
        Mg(k)=mgrav_e(isurf)
     enddo


 
     write(filename,'("Mrb",i0)'),i
     un=un+1;
     open(unit=un,file=filename,status='unknown')
     do l=1,n
     write(un,"(i8,1P10E15.6)") l,rb(l)/1.0d5,M(l)/msun,Mr(l)/msun,rho_cit(l),Mg(l)/msun
     enddo
     close(unit=un)  
   enddo      


end subroutine make_a_tov_star




subroutine tov_integrate(k)
  use tovmod1
  implicit none

  ! local variables
  integer :: k 
  integer :: i
  real*8  :: rsrifac
  real*8  :: old_data(5)
  real*8  :: new_data(5)
  real*8 :: min_press = 1.0d-30

  rho_e=rho_eff(k)
  rho_c=rho_cit(k)    


  press(1) = eosk*(rho_c)**gamma  ! EOS needed!
  press_e(1) = press(1)-clite**2*rho_c**2/rho_crit 

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
  
  nzones=10000
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
  use tovmod1, only: rho_crit,gamma, eosk, pi, g_c2, ggrav, &
       clite
  implicit none
  
  real*8, intent(IN)  :: r
  real*8, intent(IN)  :: old_data(5)
  real*8, intent(OUT) :: source_data(5)

  ! local variables
  real*8 :: mgrav,mgrav_e,press
  real*8 :: rho,mu,eps,rho_e,press_e,mu_e,eps_e
  real*8 :: min_press = 1.0d-30

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
  press = eosk*(rho)**gamma
!  eps_e = press_e /(gamma - 1.0d0)/rho_e
  eps = press/(gamma - 1.0d0)/rho
  ! set up energy density for convenience
  mu_e  = rho_e + rho_e*eps/clite**2
  mu  = rho + rho*eps/clite**2
  ! implement dp/dr and dm/dr here:
  
  source_data(1) = ggrav*(mu_e+press_e/clite**2)*(mgrav_e+4*pi*(r**3)*press_e/clite**2)/(r*(2*g_c2*mgrav_e-r))
  source_data(2) = 4*pi*(r**2)*mu_e
  source_data(3) = 4*pi*(r**2)*mu
  source_data(4) = 4*pi*(r**2)*rho/sqrt(1-2*g_c2*mgrav/r)
end subroutine tov_RHS

subroutine newraph(press_e,rho_e)
  use tovmod1, only: rho_crit,gamma, eosk, clite
  implicit none
  real*8, intent(IN)  :: press_e
  real*8, intent(OUT) :: rho_e
  real*8 :: h
  real*8, allocatable :: x(:),t(:),m(:),l(:)
  real*8, allocatable :: eqn(:),deqn(:)
  integer   :: i,n  

  n=10
  allocate(x(n))
  allocate(eqn(n))
  allocate(deqn(n))
  allocate(t(n))
  allocate(m(n))
  allocate(l(n))

  x = 0.0d0 
  eqn = 0.0d0
  deqn = 0.0d0
 
  x(1) = 1.0d15





  do i=1,n
     t(i) =  rho_crit/2*(1-(1-4*x(i)/rho_crit)**0.5)
     m(i) =  rho_crit*(1-(1-4*x(i)/rho_crit)**0.5)
     l(i) = (1-4*x(i)/rho_crit)**0.5
     h = 0.2973017*gamma*eosk
     eqn(i) = eosk*t(i)**gamma - (clite*t(i))**2/rho_crit-press_e
     deqn(i) = ((h*m(i)**(gamma-1)) - m(i)*clite**2/rho_crit)/l(i)

     x(i+1)=x(i) - eqn(i)/deqn(i)
 
  enddo

  rho_e=x(n)  

end subroutine newraph  

