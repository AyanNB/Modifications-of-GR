module tovmod1
  implicit none

  integer :: nzones = 0
  integer :: isurf = 0
  real*8  :: rho_c = 1.0d0
  real*8  :: rmax  = 1.0d0
  real*8  :: dr = 0.0d0,rho_cit=1.0d0


  real*8,allocatable  :: rs(:)
  real*8,allocatable  :: ri(:)
  real*8,allocatable  :: rho(:)
  real*8,allocatable  :: press(:)

  real*8,allocatable  :: mgrav(:),mrest(:)


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

    allocate(mgrav(nzones))
    allocate(mrest(nzones))


  end subroutine allocate

  subroutine initialize
    implicit none

    rs     = 0.0d0
    ri     = 0.0d0
    rho    = 0.0d0
    press  = 0.0d0
   
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
  integer :: l,g
  real*8 :: pr,rr

 
     
     
     rho_cit=1.65d15-(1.65d15)**2/(4d16)
     pr=eosk*(1.65d15)**gamma
     rr=1.65d15
  ! initialize everything to 0
     call initialize
     g=1
  ! perform the integration
     call tov_integrate
     do while (press(g)>10)
        g=g+1
        isurf=g-1
     end do






 

     open(unit=444,file="Mrb4d16",status='unknown')
     do l=1,isurf
     write(444,"(i8,1P10E15.6)") l,rs(l)/1.0d5,mgrav(l)/msun,mrest(l)/msun,press(l)/pr,rho(l)/rr
     enddo
     close(unit=444)  
     


end subroutine make_a_tov_star




subroutine tov_integrate
  use tovmod1
  implicit none

  ! local variables
  integer :: k 
  integer :: i
  real*8  :: rsrifac
  real*8  :: old_data(5)
  real*8  :: new_data(5)
  real*8 :: min_press = 1.0d-30


  rho_c=rho_cit    


  press(1) = eosk*(rho_c)**gamma  ! EOS needed!


  rho(1) = rho_c
  rs(1) = 0.01
  mgrav(1) = (4*pi*rho_c*rs(1)**3)/3.0d0

  mrest(1) = (4*pi*rho(1)*rs(1)**3)/3.0d0/sqrt(1-2*g_c2*mgrav(1)/rs(1))
  do i=2,nzones
     rs(i) = rs(i-1) + dr
  enddo


  
  do i=2,nzones
     old_data(1) = max(press(i-1),min_press)
     old_data(2) = mgrav(i-1)
     old_data(3) = mrest(i-1)
     call tov_RK3(old_data,new_data,rs(i-1))     
     
     press(i) = new_data(1)
     mgrav(i) = new_data(2)
     mrest(i) = new_data(3)
     rho(i) = (press(i)/eosk)**(1.0d0/gamma)
     
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
  use tovmod1, only: gamma, eosk, pi, g_c2, ggrav, &
       clite
  implicit none
  
  real*8, intent(IN)  :: r
  real*8, intent(IN)  :: old_data(5)
  real*8, intent(OUT) :: source_data(5)

  ! local variables
  real*8 :: mgrav,press
  real*8 :: rho,mu,eps
  real*8 :: min_press = 1.0d-30

! define the RHS for the integration from n to n+1
! Store results in indices 1 and 2 of array source_data
! note that you can use the constants defined in tovmod!

  ! input data, assign for convenience:

  ! this is to make sure, we don't have neg. press.
  press = max(old_data(1),min_press)
  mgrav = old_data(2)


  rho = (press/eosk)**(1.0d0/gamma)


  eps = press/(gamma - 1.0d0)/rho
  ! set up energy density for convenience

  mu  = rho + rho*eps/clite**2
  ! implement dp/dr and dm/dr here:
  
  source_data(1) = ggrav*(mu+press/clite**2)*(mgrav+4*pi*(r**3)*press/clite**2)/(r*(2*g_c2*mgrav-r))
  source_data(2) = 4*pi*(r**2)*mu
  source_data(3) = 4*pi*(r**2)*rho/sqrt(1-2*g_c2*mgrav/r)
end subroutine tov_RHS


