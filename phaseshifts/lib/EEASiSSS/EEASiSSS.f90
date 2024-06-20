! program EEASiSSS_2015_03_01.
! Elastic Electron-Atom Scattering in Solids and Solid Surfaces,
! author: John Rundgren, jru@kth.se .
! The author appreciates acknowledgement in publications by citation:
! J. Rundgren, Phys.Rev.B68,125405(2003),
!              Phys.Rev.B76,195441(2007).
!
!==============================================================================
module param
logical :: prnt
!basic constants.
!alfa=fine-structure constant, c=speed of light, cm2=1/c**2.
real(8) :: pi
real(8),parameter :: rydb=13.60569172d0,bohr=0.5291772083d0,&
  alfa=7.297352533d-03,c_light=2.d0/alfa,cm2=alfa*alfa*0.25d0
real(8) :: fpi3,facmy
character(len=2),parameter :: Elem(92)=(/&
'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si',&
'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',&
'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo',&
'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba',&
'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',&
'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',&
'At','Rn','Fr','Ra','Ac','Th','Pa','U '/)

 !select program operations.
 character(len=1) :: SelVCoul,SelVxc,SelVMadl,SelMTrad,SelOutp,SpinPS
 character(len=1) :: Omega,Theta,Sherman,CrossSection,PSvsE,UvsR
 character(len=2) :: RhoPot
 character(len=80) :: compound,W2kData
 character(len=1),allocatable :: BmSel(:)
 integer,allocatable :: SelBm(:)
 integer :: Nsel

!Mattheiss superposition prescription:
!nshell>=2 is no. of neighbor atomic shells used in Mattheiss scheme.
integer,parameter :: nshell=7

!idW   =atom identifier for WIEN2k,
!id_at =atom identifier,
!radfac=factor in Mattheiss prescription (in superposition),
 integer,allocatable :: idW(:)
 character,allocatable :: id_at(:)*7
 real(8),parameter :: radfac=0.98d0

!parameters of the integration routine ode.f90.
!node  =max. no. of integration steps.
!wave_start=start of u1(kappa,r) integration,
!eta_cutoff=phase shift cut-off,
!relerr0   =temporary relative error for investigating the magnitude
!           of wave amplitude and phase,
!relerr1   =relative error of WaveODE
!relerr2   =relative error of phase_integrate_ode
integer,parameter :: node=4000
real(8),parameter :: wave_start=1.d-12
real(8),parameter :: eta_cutoff=1.d-10
real(8),parameter :: relerr0=1.d-02,relerr1=1.d-08,relerr2=1.d-08
end module param

!==============================================================================
module maths

! subroutine DE(errsub,XCmin,XCmax,bestmem_XC)
! Differential Evolution for Optimal Control Problems.
!
! Written by Dr. Feng-Sheng Wang, Department of Chemical Engineering,
! National Chung Cheng University, Chia-Yi 621, Taiwan.
! Ref.:
! R. STORN and K.V. PRICE, J. Global Optimization 11, 341-59 (1997).
! http://www1.icsi.berkeley.edu/~storn/code.html
!
! intent(in) from external module:
!   Dim_XC : Dimension of decision parameters.
! intent(in):
!   NP : Population size.
!   itermax : Maximum number of iteration.
!   strategy=1,2,...,6 : Strategy of mutation operations, 2 is default.
!   method(1) = 0, Fixed mutation scaling factors F_XC,
!                  = 1, Random mutation scaling factors F_XC=[0, 1],
!                  = 2, Random mutation scaling factors F_XC=[-1, 1].
!   method(2) = 1, Random combination factor F_CR used for
!                         strategy = 6 in mutation operation,
!                   = 0, fixed combined factor provided by the user .
!   method(3) = 1, Saving results in external data file.
!                   = other, displaying results only.
!   errf : User provided subroutine errf(xc,fitness), where "xc" is
!       decision parameter vector and "fitness" is fitness value.
!   XCmin(Dim_XC) : Lower bound of decision parameters.
!   XCmax(Dim_XC) : Upper bound of decision parameters.
!   VTR : Desired fitness value.
!   CR_XC : Crossover factor for decision parameters.
!   refresh : Intermediate output produced after "refresh" iterations,
!       with "refresh < 1" no intermediate output.
!
! intent(in) and modified (depending on method):
!   F_XC : Mutation scaling factor for decision parameters.
!   F_CR : Mutation scaling factor for crossover factor, with strategy=6.
!
! intent(out):
!   bestval : Best value of fitness in error function.
!   bestmen_XC(Dim_XC) : Best decision parameters.
!   nfeval : Number of function calls.
!   itval : Iteration for which the bestval occurred.
implicit none
integer :: Dim_XC, NP, NPfac, itermax, strategy, refresh, method(3), &
  nfeval, itval
real(8) :: VTR,CR_XC,F_XC,F_CR,bestval
contains

  subroutine DE(errsub,XCmin,XCmax,bestmem_XC)
  use maths
  implicit none
  real(8),intent(in) :: XCmin(Dim_XC),XCmax(Dim_XC)
  real(8), dimension(Dim_XC), intent(out) :: bestmem_XC
  !
  integer :: i, ibest, iter
  integer, dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
  integer, dimension(4) :: ind
  real(8),parameter :: bohr=0.5291772083d0 !,rydb=13.60569172d0
  real(8) :: tempval
  real(8), dimension(:,:) :: pop_XC,bm_XC,mui_XC,mpo_XC,popold_XC,rand_XC,ui_XC !LD size(NP,Dim_XC)
  real(8), dimension(NP) :: val
  real(8), dimension(:) :: bestmemit_XC, rand_C1, tmp                           !LD size(Dim_XC)
!f2py intent(callback) errsub
  external errsub
!f2py integer np
!f2py intent(in,out) val
!f2py errsub(np,val)

  ! Initialize a population.
  pop_XC=0.0d0
  do i=1,NP
    call random_number(rand_C1)
    pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
  end do

  ! Evaluate fitness functions and find the best member.
  val=0.0d0
  nfeval=0
  ibest=1
  tmp(:)=pop_XC(1,:)
  call errsub(tmp, val(1))
  bestval=val(1)                                      !bestval=fitness
  nfeval=nfeval+1
  do i=2,NP
    tmp(:)=pop_XC(i,:)
    call errsub(tmp,val(i))
    nfeval=nfeval+1
    if (val(i) < bestval) then
      ibest=i
      bestval=val(i)                                  !bestval=fitness
    end if
  end do
  bestmemit_XC=pop_XC(ibest,:)
  bestmem_XC=bestmemit_XC
  itval=1

  ! Perform evolutionary computation.
  bm_XC=0.0d0
  rot=(/(i,i=0,NP-1)/)
  iter=1
  do while (iter <= itermax)
    popold_XC=pop_XC
    ! Mutation operation.
    ind=randperm(4)
    a1=randperm(NP)
    rt=mod(rot+ind(1),NP)
    a2=a1(rt+1)
    rt=mod(rot+ind(2),NP)
    a3=a2(rt+1)
    rt=mod(rot+ind(3),NP)
    a4=a3(rt+1)
    rt=mod(rot+ind(4),NP)
    a5=a4(rt+1)
    bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)

    ! Generating a random sacling factor.
    select case (method(1))
      case (1)
        call random_number(F_XC)
      case(2)
        call random_number(F_XC)
        F_XC=2.0d0*F_XC-1.0d0
    end select

    ! select a mutation strategy.
    select case (strategy)
    case (1)
      ui_XC=bm_XC + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:))
    case default
      ui_XC=popold_XC(a3,:) + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:))
    case (3)
      ui_XC=popold_XC + &
          F_XC*(bm_XC-popold_XC + popold_XC(a1,:)-popold_XC(a2,:))
    case (4)
      ui_XC=bm_XC + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:) + &
          popold_XC(a3,:)-popold_XC(a4,:))
    case (5)
      ui_XC=popold_XC(a5,:) + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:) + &
          popold_XC(a3,:)-popold_XC(a4,:))
    case (6) ! A linear crossover combination of bm_XC and popold_XC
      if (method(2) == 1) call random_number(F_CR)
      ui_XC=popold_XC + &
          F_CR*(bm_XC - popold_XC) + &
          F_XC*(popold_XC(a1,:) - popold_XC(a2,:))
    end select

    ! Crossover operation.
    call random_number(rand_XC)
    mui_XC=0.0d0
    mpo_XC=0.0d0
    where (rand_XC < CR_XC)
      mui_XC=1.0d0
      ! mpo_XC=0.0d0
    elsewhere
      ! mui_XC=0.0d0
      mpo_XC=1.0d0
    end where
    ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC

    ! Evaluate fitness functions and find the best member.
    do i=1,NP
      ! Confine each of feasible individuals in the lower-upper bound.
      ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
      tmp(:)=ui_XC(i,:)
      call errsub(tmp,tempval)
      nfeval=nfeval+1
      if (tempval < val(i)) then
        pop_XC(i,:)=ui_XC(i,:)
        val(i)=tempval
        if (tempval < bestval) then
          bestval=tempval                          !bestval=fitness
          itval=iter
          bestmem_XC=ui_XC(i,:)
        end if
      end if
    end do
    bestmemit_XC=bestmem_XC
    if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
      if (method(3)==1) write(61,203) iter
      write(*,203) iter
      do i=1,Dim_XC
        if (method(3)==1) write(61,202) i, bestmem_XC(i)*bohr
        write(*,202) i,bestmem_XC(i)*bohr
      end do
      if (method(3)==1) write(61,201) bestval
        write(*,201) bestval
    end if
    iter=iter+1
    if(bestval <= VTR )then    ! .and. refresh > 0) then
      exit
    endif
  end do
  return
  201 format(2x, 'bestval =', es9.1 /)
  202 format(5x, 'bestmem_XC(', I3, ') =', es10.3)
  203 format(2x, 'No. of iteration  =', I8)
  end subroutine DE
  !-------------------------------------------------------------------
  function randperm(num)
  implicit none
  integer, intent(in) :: num
  integer :: number, i, j, k
  integer, dimension(num) :: randperm
  real(8), dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
          number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
          number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
  end function randperm
  !--------------------------------------------------------------------
  subroutine sph_bessel_j(lmax,x,bj)
  !SPH_BESSEL_J is spherical bessel function j(l,x).
  implicit none
  integer :: l,ll,lx,lmax,llp1,llp3
  real(8)  :: bj(0:lmax),x,xi,xl,x2h,x2hh,cl,f,s,almax
  real(8),allocatable :: aux(:)
  almax=lmax
  !small arguments: limit determined by mashine precision 2.22E-16 .
  if(x<1.0188d-02)then
    x2h=x*x*0.5d0; x2hh=x2h*0.5d0
    bj(0)=1.d0-(x2h/3.d0)*(1.d0-x2hh/5.d0)
    xl=1.d0; cl=1.d0
    do l=1,lmax
      xl=xl*x; ll=l+l; cl=cl/(ll+1)
      bj(l)=cl*xl*(1.d0-(x2h/(ll+3))*(1.d0-x2hh/(ll+5)))
    enddo
    return
  endif
  !large arguments.
  xi=1.d0/x
  if(x>almax)then
    !upward recurrence.
    bj(0)=sin(x)*xi
    if(lmax==0) return
    bj(1)=(bj(0)-cos(x))*xi
    if(lmax==1) return
    ll=3
    do l=2,lmax
      bj(l)=(ll*xi)*bj(l-1)-bj(l-2)
      ll=ll+2
    enddo
  else
    !downward recurrence.
    if(lmax<100)then
      lx=lmax+10+int(x*(1.d0-almax*0.0075d0))
    else
      lx=lmax+10+int(x*0.25d0)
    endif
    allocate(aux(0:lx)) !lx variable.
    aux(lx)=0.d0
    aux(lx-1)=1.d-30
    llp1=lx+lx-1; llp3=lx+lx+1
    s=llp1*aux(lx-1)**2
    do l=lx-2,0,-1
      llp1=llp1-2; llp3=llp3-2
      aux(l)=(llp3*xi)*aux(l+1)-aux(l+2)
      s=s+llp1*aux(l)*aux(l)
    enddo
    !normalisation by sum (2*l+1)*j(l,x)**2=1, l=0:lx.
    f=1.d0/sqrt(s)
    bj(0:lmax)=aux(0:lmax)*f
  endif
  return
  end subroutine sph_bessel_j
  !--------------------------------------------------------------------
  subroutine sph_bessel_y(lmax,x,by)
  !SPH_BESSEL_Y is spherical bessel function y(l,x).
  implicit none
  integer :: l,ll,lmax
  real(8)  :: by(0:lmax),x,xi
  xi=1.d0/x
  by(0)=-cos(x)*xi
  if(lmax==0) return
  by(1)=(by(0)-sin(x))*xi
  if(lmax==1) return
  ll=3
  do l=2,lmax
    by(l)=(ll*xi)*by(l-1)-by(l-2)
    ll=ll+2
  enddo
  return
  end subroutine sph_bessel_y
  !--------------------------------------------------------------------
  function integral_egt(x1,x2,m,n,r,d,f)
  !INTEGRAL_EGT is exponential-grid trapezoidal integration;
  !f(r)*(r**m) is integrated.
  !Ref.: Abramowitz and Stegun, chap. 25.2.
  implicit none
  real(8) :: integral_egt
  integer,intent(in) :: n,m
  real(8),intent(in) :: x1,x2,r(n),d,f(n)
  integer :: i1,i2,i
  real(8) :: f1,f2,p1,p2,s,w1,w2
  !p is a subscript value on a continuous exponential-grid.
  p1=log(x1/r(1))/d+1.d0
  p2=log(x2/r(1))/d+1.d0
  i1=max(int(p1),2); i1=min(int(p1),n-1)
  i2=max(int(p2),2); i2=min(int(p2),n-1)
  !trapezoidal integration on the interval [r(i1),r(i2)].
  s=0.d0
  w1=f(i1)*r(i1)**m
  do i=i1+1,i2
    w2=f(i)*r(i)**m
    s=s+w1+w2
    w1=w2
  enddo
  !end-point corrections:
  !lagrangian 3-point interpolation at x1 and x2, followed by
  !trapezoidal integration on [r(i1),x1] (subtraction) and on
  ![r(i2),x2] (addition).
  if(i1==1) i1=2
  p1=p1-i1
  p2=p2-i2
  f1=0.5d0*p1*(p1-1.d0)*f(i1-1)+(1.d0-p1*p1)*f(i1)+0.5d0*p1*(p1+1.d0)*&
     f(i1+1)
  f2=0.5d0*p2*(p2-1.d0)*f(i2-1)+(1.d0-p2*p2)*f(i2)+0.5d0*p2*(p2+1.d0)*&
     f(i2+1)
  s=s-(f(i1)*r(i1)**m+f1*x1**m)*p1+(f(i2)*r(i2)**m+ f2*x2**m)*p2
  integral_egt=0.5d0*d*s
  return
  end function integral_egt
  !--------------------------------------------------------------------
  function interpol_egl3(x,n,r,d,f,messg)
  !INTERPOL_EGL3 is 'exponential-grid lagrangian 3-point interpolation.
  !Ref.: Abramowitz and Stegun, Sec. 25.2.11.
  implicit none
  real(8) :: interpol_egl3
  integer,intent(in) :: n
  real(8),intent(in) :: x,r(n),d,f(n)
  character(len=*),intent(in) :: messg
  integer :: i1
  real(8) :: p1,p
  if(x-r(n)>0.d0)then
    write(61,900)'interpol_egl3 stop: x-r(n),r(n)=',x-r(n),r(n)
    write(61,902)'error location: ',messg
    stop
  endif
  900 format(a,es9.2,f7.3)
  902 format(2a)
  p1=(log(x/r(1)))/d+1.d0
  i1=nint(p1)
  if(i1==n) i1=i1-1
  p=p1-i1
  interpol_egl3=&
    0.5d0*p*(p-1.d0)*f(i1-1)+(1.d0-p*p)*f(i1)+0.5d0*p*(p+1.d0)*f(i1+1)
  return
  end function interpol_egl3
  !--------------------------------------------------------------------
  function dfdx_egl3(x,n,r,d,f,messg)
  !DFDX_EGL3 is 'exponential-grid lagrangian 3-point differentiation
  !with respect to radius.
  !Ref.: Abramowitz and Stegun, Sec. 25.2.11.
  implicit none
  real(8) :: dfdx_egl3
  integer,intent(in) :: n
  real(8),intent(in) :: x,r(n),d,f(n)
  character(len=*),intent(in) :: messg
  integer :: i1
  real(8) :: p1,p
  if(x-r(n)>0.d0)then
    write(61,'(a,es9.2)')'dfdx_egl3 stop: x-r(n)=',x-r(n)
    write(61,'(2a)')'error location: ',messg
    stop
  endif
  p1=(log(x/r(1)))/d+1.d0
  i1=nint(p1)
  if(i1==n) i1=i1-1
  p=p1-i1
  dfdx_egl3=((p-0.5d0)*f(i1-1)-2.d0*p*f(i1)+(p+0.5d0)*f(i1+1))/(x*d)
  return
  end function dfdx_egl3
  !--------------------------------------------------------------------
  function simpson_s(rx,dx,nx,FUNCT,JSTART,JEND)
  !simpson's 3/8 rule, even point at start.
  !from http://www.WIEN2k.at with permission.
  implicit none
  integer :: nx,jstart,jend,js,jf,j,npoint
  real(8) :: simpson_s,rx(nx),dx,FUNCT(nx)
  NPOINT=JEND-JSTART+1
  JS=JSTART+3-3*MOD(NPOINT,2)
  JF=JEND-2
  simpson_s=0.D0
  do J=JS,JF,2
    simpson_s=simpson_s&
      +FUNCT(J)*RX(J)+4.*FUNCT(J+1)*RX(J+1)+FUNCT(J+2)*RX(J+2)
  enddo
  simpson_s=simpson_s*dx/3.0D0
  if(npoint.eq.2)then
    simpson_s=simpson_s&
      +(1-MOD(NPOINT,2))*(FUNCT(JSTART)+FUNCT(JStart+1))/2.D0&
      *(RX(JStart+1)-RX(JSTART))
  elseif(mod(npoint,2).eq.0)then
    simpson_s=simpson_s&
          +(FUNCT(jstart)*RX(jstart)&
      +3.d0*FUNCT(jstart+1)*RX(jstart+1)&
      +3.d0*FUNCT(jstart+2)*RX(jstart+2)&
           +FUNCT(jstart+3)*RX(jstart+3))*dx*0.375d0
  endif
  return
  end function simpson_s
  !-------------------------------------------------------------------
  subroutine poisson(rx,dx,nx,z,rho4pir2,QvsR,VvsR)
  implicit none
  integer :: nx,j
  real(8) :: rx(nx),rho4pir2(nx),QvsR(nx),VvsR(nx),f(nx),s(nx)
  real(8) :: z,dx,a
  do j=1,nx
    QvsR(j)=simpson_s(rx,dx,nx,rho4pir2,1,j)
  enddo
  f=rho4pir2/rx
  do j=1,nx
    s(j)=simpson_s(rx,dx,nx,f,1,j)
  enddo
  a=s(nx)
  VvsR=2.d0*(QvsR/rx+a-s)-2.d0*z/rx
  return
  end subroutine poisson
!----------------------------------------------------------------------
subroutine spcoef(n,x,f,prim1,primn,b,c,d,name)
!SPCOEF calculates coefficients defining a smooth cubic interpolatory
!spline; written by L.F. Shampine.
!Ref.: L.F. Shampine, R.C. Allen, and S. Pruess, Fundamentals of
!      numerical computing (Wiley, 1997);
!Ref.: anoymous ftp, ftp.wiley.com, public/college/math/sapcodes.

!input arguments:
!   n    = number of data points.
!   x    = vector of n values of the independent variable ordered
!          so that  x(1) < x(2) < ... < x(n).
!   f    = vector of values of the dependent variable.
!   prim1= df/dx at point 1,
!   primn= df/dx at point n,
!          external df/dx used if |df/dx|<1.d+30, else internal
!          calculation,
!    name= name of calling routine.
!output arguments:
!    b   = vector of s'(x(i)) values.
!    c   = vector of s''(x(i))/2 values.
!    d   = vector of s'''(x(i)+)/6 values (i < n).
!   flag =  0  normal return.
!        = -1  one of x, f, b, c, or d has fewer than n entries.
!        = -2  x vector is incorrectly ordered.

      implicit none
      character(len=*)   :: name
      integer :: n,flag,i,k
      real(8),dimension(n) :: x,f,b,c,d
      real(8) :: fp1,fpn,p,prim1,primn
      real(8),dimension(n-1) :: diff_f,h

      flag = 0
      if ((size(x,dim=1) < n) .or. (size(f,dim=1) < n) .or. &
          (size(b,dim=1) < n) .or. (size(c,dim=1) < n) .or. &
          (size(d,dim=1) < n)) then
        flag = -1
        write(61,*)'routine ',name,': spcoef stop, flag=',flag
        stop
      end if
      h = diff(x)
      if (any(h <= 0)) then
        flag = -2
        write(61,*)'routine ',name,': spcoef stop, flag=',flag
        stop
      end if

!calculate coefficients for the tridiagonal system: store
!sub-diagonal in b, diagonal in d, difference quotient in c.
      b(1:n-1) = h
      diff_f = diff(f)
      c(1:n-1) = diff_f/h
      if (n == 2) then
        b(1) = c(1)
        c(1) = 0.0d0
        d(1) = 0.0d0
        b(2) = b(1)
        c(2) = 0.0d0
        return
      end if
      d(1) = 2.0d0*b(1)
      do i = 2,n-1
        d(i) = 2.0d0*(b(i) + b(i-1))
      end do
      d(n) = 2.0d0*b(n-1)

!calculate estimates for the end slopes using polynomials
!that interpolate the data nearest the end.
      if(abs(prim1).lt.1.d+30)then
        fp1=prim1
      else
        fp1 = c(1) - b(1)*(c(2) - c(1))/(b(1) + b(2))
        if (n > 3) then
          fp1 = fp1 + b(1)*((b(1) + b(2))*(c(3) - c(2))/ &
                (b(2) + b(3)) - c(2) + c(1))/(x(4) - x(1))
        end if
      endif

      if(abs(primn).lt.1.d+30)then
        fpn=primn
      else
        fpn = c(n-1) + b(n-1)*(c(n-1) - c(n-2))/(b(n-2) + b(n-1))
        if (n > 3) then
          fpn = fpn + b(n-1)*(c(n-1) - c(n-2) - (b(n-2) + b(n-1))* &
                (c(n-2) - c(n-3))/(b(n-2) + b(n-3)))/(x(n) - x(n-3))
        end if
      endif

!calculate the right hand side and store it in c.
      c(n) = 3.0d0*(fpn - c(n-1))
      do i = n-1,2,-1
        c(i) = 3.0d0*(c(i) - c(i-1))
      end do
      c(1) = 3.0d0*(c(1) - fp1)

!solve the tridiagonal system.
      do k = 2,n
        p = b(k-1)/d(k-1)
        d(k) = d(k) - p*b(k-1)
        c(k) = c(k) - p*c(k-1)
      end do
      c(n) = c(n)/d(n)
      do k = n-1,1,-1
        c(k) = (c(k) - b(k)*c(k+1))/d(k)
      end do

!calculate the coefficients defining the spline.
      d(1:n-1) = diff(c)/(3.0d0 * h)
      b(1:n-1) = diff_f/h - h*(c(1:n-1) + h*d(1:n-1))
      b(n) = b(n-1) + h(n-1)*(2.0d0*c(n-1) + h(n-1)*3.0d0*d(n-1))
      return
      contains

        function diff(v)
        !auxiliary function to compute the forward difference of data
        !stored in a vector v.
        real(8),dimension(:),intent(in) :: v
        real(8),dimension(size(v)-1)    :: diff
        integer :: n
        n=size(v)
        diff=v(2:n)-v(1:n-1)
        return
        end function diff
end subroutine spcoef
!-----------------------------------------------------------------------
subroutine svalue(n,x,f,b,c,d,t,s,s1,last_interval,name)
!SVALUE evaluates the spline s and s1 at t using coefficients from
!SPCOEF written by L.F. Shampine.
!there are n nodes ordered so that x(1) < x(2) < ... < x(n).
!input arguments:
!    x, f, b, c, d are defined as in spcoef.
!    t             point where the spline s is to be evaluated,
!    name is name of calling routine.
!output arguments:
!    s        =  value of spline at t
!    s1       =  ds/dt at t
!    flag     =  0  normal return
!             =  1  t < x(1)
!             =  2  t > x(n)
      implicit none
      character(len=*)    :: name
      integer,intent(in)  :: n
      integer :: flag
      integer             :: interval,last_interval,j
      real(8),intent(in)   :: t
      real(8),intent(out)  :: s,s1
      real(8)              :: dt
      real(8),dimension(n),intent(in) :: x,f,b,c,d

!search for correct interval for t.
      if(t < x(1))then
        flag = 1
        interval = 1
      elseif(t < x(n))then       ! x(1) <= t < x(n).
        flag = 0
        if (t >= x(last_interval)) then
          do j = last_interval,n-1
            if (t < x(j+1)) then
              interval = j
              exit
            endif
          enddo
        else
          do j = last_interval-1,1,-1
            if (t >= x(j)) then
              interval = j
              exit
            endif
          enddo
        endif
      elseif (t > x(n)) then
        flag = 2
        interval = n - 1
      else                          ! t = x(n)
        flag = 0
        interval = n - 1
      end if
      last_interval = interval
      if(name.ne.'rho_center'.and.flag.ne.0)then
        write(61,*)'routine ',name,': svalue warning, flag=',flag
      endif

!evaluate cubic polynomial.
      dt = t - x(interval)
      s = f(interval)+dt*(b(interval)+dt*(c(interval)+dt*d(interval)))
      s1= b(interval)+dt*(c(interval)*2.d0+dt*d(interval)*3.d0)
      return
end subroutine svalue
!-----------------------------------------------------------------------
!subroutine ODE
!integrates a system of neqn first order ordinary differential eqs.:
!  dy(i)/dt = f(t,y(1:neqn))
!  y(i) given at an initial point t.
!
!Ref.: L.F. Shampine and M.K. Gordon, Computer solution of ordinary
!      differential equations: the initial value problem (Freeman,1975).
!Ref.: http//www.netlib.org, search for ode/ode.f.
!
!the subroutine integrates from  t  to  tout .  on return the
!parameters in the call list are set for continuing the integration.
!the user has only to define a new value  tout  and call  ode  again.
!
!the differential equations are actually solved by a suite of codes
!de ,  step , and  intrp .  ode  allocates virtual storage in the
!arrays  work  and  iwork  and calls  de .  de  is a supervisor which
!directs the solution.  it calls on the routines  step  and  intrp
!to advance the integration and to interpolate at output points.
!step  uses a modified divided difference form of the adams pece
!formulas and local extrapolation.  it adjusts the order and step
!size to control the local error per unit step in a generalized
!sense.  normally each call to  step  advances the solution one step
!in the direction of  tout .  for reasons of efficiency  de
!integrates beyond  tout  internally, though never beyond
!t+10*(tout-t), and calls  intrp  to interpolate the solution at
!tout .  an option is provided to stop the integration at  tout  but
!it should be used only if it is impossible to continue the
!integration beyond  tout .
!
!the parameters represent:
!   f -- double precision subroutine f(t,y,yp) to evaluate
!             derivatives yp(i)=dy(i)/dt
!   neqn (integer*4)-- number of equations to be integrated,
!   y(*) (real(8))-- solution vector at t,
!   t (real(8))-- independent variable,
!   tout (real(8))-- point at which solution is desired,
!   relerr,abserr (real(8))-- relative and absolute error tolerances for
!        local error test .  at each step the code requires
!          dabs(local error) .le. dabs(y)*relerr + abserr
!        for each component of the local error and solution vectors,
!   iflag (integer*4)-- indicates status of integration,
!   work(*) (real(8))  -- arrays to hold information internal to
!   iwork(*) (integer*4)    which is necessary for subsequent calls,
!   OR(*) (real(8)) -- integration points chosen by the program,
!   OU(*) (real(8)) -- the corresponding values of y(1),
!   NOX (integer*4) -- the number of integration points,
!   NODE (integer*4) -- dimension of OR and OU.
!
!first call to ode:
!the user must provide storage in his calling program for the arrays
!in the call list,
!   y(neqn), work(100+21*neqn), iwork(5),
!declare  f  in an external statement, supply the double precision
!subroutine f(t,y,yp)  to evaluate
!   dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
!and initialize the parameters:
!   neqn  -- number of equations to be integrated
!   y(*)  -- vector of initial conditions
!   t     -- starting point of integration
!   tout  -- point at which solution is desired
!   relerr,abserr -- relative and absolute local error tolerances
!   iflag -- +1,-1.  indicator to initialize the code.  normal input
!            is +1.  the user should set iflag=-1 only if it is
!            impossible to continue the integration beyond  tout .
!all parameters except  f ,  neqn  and  tout  may be altered by the
!code on output so must be variables in the calling program.
!
!output from ode:
!   neqn -- unchanged
!   y(*) -- solution at  t
!   t    -- last point reached in integration. normal return has t=tout.
!   tout -- unchanged
!   relerr,abserr -- normal return has tolerances unchanged.  iflag=3
!        signals tolerances increased
!   iflag = 2 -- normal return.  integration reached  tout
!         = 3 -- integration did not reach  tout  because error
!                tolerances too small.  relerr ,  abserr  increased
!                appropriately for continuing
!         = 4 -- integration did not reach  tout  because more than
!                NODE steps needed
!         = 5 -- integration did not reach  tout  because equations
!                appear to be stiff
!         = 6 -- invalid input parameters (fatal error)
!        the value of  iflag  is returned negative when the input
!        value is negative and the integration does not reach  tout ,
!        i.e., -3, -4, -5.
!   work(*),iwork(*) -- information generally of no interest to the
!                       user but necessary for subsequent calls.
!
!subsequent calls to ode:
!subroutine  ode  returns with all information needed to continue
!the integration.  if the integration reached  tout , the user need
!only define a new  tout  and call again.  if the integration did not
!reach  tout  and the user wants to continue, he just calls again.
!the output value of  iflag  is the appropriate input value for
!subsequent calls.  the only situation in which it should be altered
!is to stop the integration internally at the new  tout , i.e.,
!change output  iflag=2  to input  iflag=-2 .  error tolerances may
!be changed by the user before continuing.  all other parameters must
!remain unchanged.
subroutine ode(f,neqn,y,t,tout,relerr,abserr,iflag,nwork,work,&
               iwork,OH,OR,OU,NOX,NODE)
      implicit none
      logical :: start,phase1,nornd
      integer :: neqn,nwork,NODE,iyy,iwt,ip,iyp,iypout,iphi,NOX,&
        istart,iphase,iwork(5),ialpha,ibeta,isig,iv,iw,ig,ipsi,ix,ih,&
        ihold,itold,idelsn,iflag
      real(8) :: y(neqn),t,tout,relerr,abserr,work(nwork),&
        OH,OR(NODE),OU(2,NODE),twou,fouru
!f2py intent(callback) f
      external f
!     The following lines define the signature of function f for F2PY:
!f2py real*8 y
!f2py f(y)
!
      data ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart,&
        itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/
      twou=2.d0*epsilon(1.d0); fouru=2.d0*twou
      iyy = 100
      iwt = iyy + neqn
      ip = iwt + neqn
      iyp = ip + neqn
      iypout = iyp + neqn
      iphi = iypout + neqn
      NOX=0
      if(iabs(iflag) .eq. 1) go to 1
      start = work(istart) .gt. 0.0d0
      phase1 = work(iphase) .gt. 0.0d0
      nornd = iwork(2) .ne. -1
    1 call de(f,neqn,y,t,tout,relerr,abserr,&
        iflag,work(iyy),work(iwt),work(ip),work(iyp),work(iypout),&
        work(iphi),work(ialpha),work(ibeta),work(isig),work(iv),&
        work(iw),work(ig),phase1,work(ipsi),work(ix),work(ih),&
        work(ihold),start,work(itold),work(idelsn),iwork(1),nornd,&
        iwork(3),iwork(4),iwork(5),OH,OR,OU,NOX,NODE,twou,fouru)
      work(istart) = -1.0d0
      if(start) work(istart) = 1.0d0
      work(iphase) = -1.0d0
      if(phase1) work(iphase) = 1.0d0
      iwork(2) = -1
      if(nornd) iwork(2) = 1
      return
contains

!-----------------------------------------------------------------------
subroutine de(f,neqn,y,t,tout,relerr,abserr,iflag,yy,wt,p,yp,&
              ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold,start,&
              told,delsgn,ns,nornd,k,kold,isnold,OH,OR,OU,NOX,NODE,&
              twou,fouru)

!ode  merely allocates storage for  de  to relieve the user of the
!inconvenience of a long call list.
!NODE  is the maximum number of steps allowed in one call to ode .
      implicit none
      logical :: stiff,crash,start,phase1,nornd
      integer :: neqn,NODE,iflag,isn,kle4,isnold,NOX,kold,k,ns
      real(8) :: t,tout,relerr,abserr,eps,told,del,absdel,tend,&
        releps,abseps,delsgn,x,h,twou,fouru,hold,OH
      real(8) :: y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),&
        yp(neqn),ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),&
        w(12),g(13),OR(NODE),OU(2,NODE)
!f2py intent(callback) f
      external f
!test for improper parameters.
      if(neqn .lt. 1) go to 10
      if(t .eq. tout) go to 10
      if(relerr .lt. 0.0d0  .or.  abserr .lt. 0.0d0) go to 10
      eps = dmax1(relerr,abserr)
      if(eps .le. 0.0d0) go to 10
      if(iflag .eq. 0) go to 10
      isn = isign(1,iflag)
      iflag = iabs(iflag)
      if(iflag .eq. 1) go to 20
      if(t .ne. told) go to 10
      if(iflag .ge. 2  .and.  iflag .le. 5) go to 20
   10 iflag = 6
      return
!on each call set interval of integration and counter for number of
!steps. adjust input error tolerances to define weight vector for
!subroutine  step.
   20 continue
      del = tout - t
      absdel = dabs(del)
      tend = t + 10.0d0*del
      if(isn .lt. 0) tend = tout
      kle4 = 0
      stiff = .false.
      releps = relerr/eps
      abseps = abserr/eps
      if(iflag .eq. 1) go to 30
      if(isnold .lt. 0) go to 30
      if(delsgn*del .gt. 0.0d0) goto 50
!on start and restart also set work variables x and yy(*), store the
!direction of integration and initialize the step size.
   30 continue
      start = .true.
      x = t
      yy(1:neqn)=y(1:neqn)
      delsgn = dsign(1.0d0,del)
      h = dsign(dmax1(dabs(tout-x),fouru*dabs(x)),tout-x)
      if(OH.gt.0.d0) h=min(h,OH)
      if(NOX.eq.0)then
        NOX=NOX+1; OR(NOX)=x; OU(1:neqn,NOX)=y(1:neqn)
      endif
!if already past output point, interpolate and return.
   50 continue
      if(dabs(x-t) .lt. absdel) goto 60
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      !NOX is not incremented when OR and OU replace overshoot values.
      OR(NOX)=tout; OU(1:neqn,NOX)=y(1:neqn)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
!if cannot go past output point and sufficiently close,
!extrapolate and return.
   60 continue
      if(isn .gt. 0  .or.  dabs(tout-x) .ge. fouru*dabs(x)) goto 80
      h = tout - x
      call f(x,yy,yp)
      y(1:neqn)=yy(1:neqn)+h*yp(1:neqn)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
!test for too many steps.
   80 continue
      if(NOX .lt. NODE-2) goto 100
      iflag = isn*4
      if(stiff) iflag = isn*5
      write(61,'(a,i4,a,i2)')&
           'ode: more than ',NODE,' steps, iflag=',iflag
      y(1:neqn)=yy(1:neqn)
      t = x
      told = t
      isnold = 1
      return
!limit step size, set weight vector and take a step.
  100 continue
      h = dsign(dmin1(dabs(h),dabs(tend-x)),h)
      wt(1:neqn) = releps*dabs(yy(1:neqn)) + abseps
      call step(x,yy,f,neqn,h,eps,wt,start,hold,k,kold,crash,phi,p,yp,&
                psi,alpha,beta,sig,v,w,g,phase1,ns,nornd,twou,fouru)
      NOX=NOX+1; OR(NOX)=x; OU(1:neqn,NOX)=yy(1:neqn)
!test for tolerances too small.
      if(.not.crash) go to 130
      iflag = isn*3
      relerr = eps*releps
      abserr = eps*abseps
      y(1:neqn)=yy(1:neqn)
      t = x
      told = t
      isnold = 1
      return
!augment counter on number of steps and test for stiffness.
  130 continue
      kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      goto 50
end subroutine de
!-----------------------------------------------------------------------
!SUBROUTINE STEP integrates a system of first order ordinary
!differential equations one step, normally from x to x+h, using a
!modified divided difference form of the adams pece formulas.  local
!extrapolation is used to improve absolute stability and accuracy.
!the code adjusts its order and step size to control the local error
!per unit step in a generalized sense.  special devices are included
!to control roundoff error and to detect when the user is requesting
!too much accuracy.
!
!the parameters represent:
!   x     -- independent variable (real(8))
!   y(*)  -- solution vector at x (real(8))
!   yp(*) -- derivative of solution vector at  x  after successful
!            step (real(8))
!   neqn  -- number of equations to be integrated (integer*4)
!   h     -- appropriate step size for next step. normally determined by
!            code (real(8))
!   eps   -- local error tolerance.  must be variable (real(8))
!   wt(*) -- vector of weights for error criterion (real(8))
!   start -- logical variable set .true. for first step,  .false.
!            otherwise (logical*4)
!   hold  -- step size used for last successful step (real(8))
!   k     -- appropriate order for next step (determined by code)
!   kold  -- order used for last successful step
!   crash -- logical variable set .true. when no step can be taken,
!            .false. otherwise.
!the arrays  phi, psi  are required for the interpolation subroutine
!intrp.  the array p is internal to the code.  all are real(8)
!
!input to  step:
!
!   first call:
!the user must provide storage in his driver program for all arrays
!in the call list, namely,
!   y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
!
!the user must also declare  start  and  crash  logical variables
!and  f  an external subroutine, supply the subroutine  f(x,y,yp)
!to evaluate
!   dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
!and initialize only the following parameters:
!   x     -- initial value of the independent variable
!   y(*)  -- vector of initial values of dependent variables
!   neqn  -- number of equations to be integrated
!   h     -- nominal step size indicating direction of integration
!            and maximum size of step.  must be variable
!   eps   -- local error tolerance per step.  must be variable
!   wt(*) -- vector of non-zero weights for error criterion
!   start -- .true.
!
!step  requires the l2 norm of the vector with components
!local error(l)/wt(l)  be less than  eps  for a successful step.  the
!array  wt  allows the user to specify an error test appropriate
!for his problem.  for example,
!   wt(l) = 1.0  specifies absolute error,
!         = dabs(y(l))  error relative to the most recent value of
!              the l-th component of the solution,
!         = dabs(yp(l))  error relative to the most recent value of
!              the l-th component of the derivative,
!         = dmax1(wt(l),dabs(y(l)))  error relative to the largest
!              magnitude of l-th component obtained so far,
!         = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed
!              relative-absolute test where  relerr  is relative
!              error,  abserr  is absolute error and  eps =
!              dmax1(relerr,abserr) .
!
!   subsequent calls:
!subroutine  step  is designed so that all information needed to
!continue the integration, including the step size  h  and the order
!k , is returned with each step.  with the exception of the step
!size, the error tolerance, and the weights, none of the parameters
!should be altered.  the array  wt  must be updated after each step
!to maintain relative error tests like those above.  normally the
!integration is continued just beyond the desired endpoint and the
!solution interpolated there with subroutine  intrp .  if it is
!impossible to integrate beyond the endpoint, the step size may be
!reduced to hit the endpoint since the code will not take a step
!larger than the  h  input.  changing the direction of integration,
!i.e., the sign of  h , requires the user set  start = .true. before
!calling  step  again.  this is the only situation in which  start
!should be altered.
!
!output from  step:
!
!   successful step:
!the subroutine returns after each successful step with  start  and
!crash  set .false. .  x  represents the independent variable
!advanced one step of length  hold  from its value on input and  y
!the solution vector at the new value of  x .  all other parameters
!represent information corresponding to the new  x  needed to
!continue the integration.
!
!   unsuccessful step:
!when the error tolerance is too small for the machine precision,
!the subroutine returns without taking a step and  crash = .true. .
!an appropriate step size and error tolerance for continuing are
!estimated and all other information is restored as upon input
!before returning.  to continue with the larger tolerance, the user
!just calls the code again.  a restart is neither required nor
!desirable.
subroutine step(x,y,f,neqn,h,eps,wt,start,hold,k,kold,crash,phi,p,yp,&
                psi,alpha,beta,sig,v,w,g,phase1,ns,nornd,twou,fouru)
      implicit none
      logical :: start,crash,phase1,nornd
      integer :: neqn,l,k,kold,ifail,kp1,kp2,km1,km2,ns,nsp1,i,im1,&
        iq,nsm2,j,limit1,nsp2,limit2,ip1,knew
      real(8) :: y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),&
        psi(12),alpha(12),beta(12),sig(13),w(12),v(12),g(13),gstr(13),&
        two(13)
      real(8) :: x,h,eps,hold,twou,fouru,p5eps,round,sum,absh,&
        realns,temp1,temp2,temp3,temp4,temp5,temp6,reali,tau,xold,erk,&
        erkm1,erkm2,erkp1,err,rho,hnew,r
      external f
      data two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,&
        512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
      data gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0,&
        0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,&
        0.00468d0/

!begin BLOCK 0
!check if step size/error tolerance is too small for machine precision.
!if first step, initialize phi array and estimate a starting step size.
!
!if step size is too small, determine an acceptable one.
      crash = .true.
      if(dabs(h) .ge. fouru*dabs(x)) go to 5
      h = dsign(fouru*dabs(x),h)
      return
    5 p5eps = 0.5d0*eps
!if error tolerance is too small, increase it to an acceptable value.
      round = 0.0d0
      do 10 l = 1,neqn
        round = round + (y(l)/wt(l))**2
   10 continue
      round = twou*dsqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      return
   15 crash = .false.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      if(.not.start) go to 99
!initialize. compute appropriate step size for first step.
      call f(x,y,yp)
      sum = 0.0d0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
        sum = sum + (yp(l)/wt(l))**2
   20 continue
      sum = dsqrt(sum)
      absh = dabs(h)
      if(eps .lt. 16.0d0*sum*h*h) absh = 0.25d0*dsqrt(eps/sum)
      h = dsign(dmax1(absh,fouru*dabs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0d0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
        phi(l,15) = 0.0d0
   25 continue
   99 ifail = 0
!end BLOCK 0
!
!begin BLOCK 1
!compute coefficients of formulas for this step. avoid computing
!those quantities not changed when step size is not changed.
!                  ***
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
!ns is the number of steps taken with size h, including the current one.
!when k.lt.ns, no coefficients change.
      if(h .ne. hold) ns = 0
      if(ns.le.kold)   ns=ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
!compute those components of alpha(*),beta(*),psi(*),sig(*) which are
!changed.
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
        sig(i+1) = reali*alpha(i)*sig(i)
  105 continue
  110 psi(k) = temp1
!compute coefficients g(*). initialize v(*) and set w(*).
!g(2) is set in data statement.
      if(ns .gt. 1) go to 120
      do iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
        w(iq) = v(iq)
      enddo
      go to 140
!if order was raised, update diagonal part of v(*).
  120 if(k .le. kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      if(nsm2 .lt. 1) go to 130
      do 125 j = 1,nsm2
        i = k-j
        v(i) = v(i) - alpha(j+1)*v(i+1)
  125 continue
!update v(*) and set w(*).
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
        w(iq) = v(iq)
  135 continue
      g(nsp1) = w(1)
!compute the g(*) in the work vector w(*).
  140 nsp2 = ns + 2
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
          w(iq) = w(iq) - temp6*w(iq+1)
  145   continue
        g(i) = w(1)
  150 continue
  199 continue
!end BLOCK 1
!
!begin BLOCK 2
!predict a solution p(*), evaluate derivatives using predicted solution,
!estimate local error at order k and errors at orders k, k-1, k-2
!as if constant step size were used.
!
!change phi to phi star.
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
          phi(l,i) = temp1*phi(l,i)
  205   continue
  210 continue
!predict solution and differences.
  215 do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
        p(l) = 0.0d0
  220 continue
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
          phi(l,i) = phi(l,i) + phi(l,ip1)
  225   continue
  230 continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
        phi(l,16) = (p(l) - y(l)) - tau
  235 continue
      go to 250
  240 do 245 l = 1,neqn
        p(l) = y(l) + h*p(l)
  245 continue
  250 xold = x
      x = x + h
      absh = dabs(h)
      call f(x,p,yp)
!estimate errors at orders k,k-1,k-2.
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 264 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
!!      if(km2)265,260,255
!!255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
!!260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
!!265   erk = erk + (temp4*temp3)**2
        if(km2>0)then
          erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
          erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
        elseif(km2==0)then
          erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
        endif
        erk = erk + (temp4*temp3)**2
  264 continue
!!    if(km2)280,275,270
!!270 erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
!!275 erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
!!280 temp5 = absh*dsqrt(erk)
      if(km2>0)then
        erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
        erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
      elseif(km2==0)then
        erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
      endif
      temp5 = absh*dsqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
!test if order should be lowered.
!!    if(km2)299,290,285
!!285 if(dmax1(erkm1,erkm2) .le. erk) knew = km1
!!    go to 299
!!290 if(erkm1 .le. 0.5d0*erk) knew = km1
!!299 if(err .le. eps) go to 400
      if(km2>0)then
        if(dmax1(erkm1,erkm2) .le. erk) knew = km1
      elseif(km2==0)then
        if(erkm1 .le. 0.5d0*erk) knew = km1
      endif
      if(err .le. eps) go to 400
!test if step successful.
!end BLOCK 2
!
!begin BLOCK 3
!the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
!if third consecutive failure, set order to one. if step fails more than
!three times, consider an optimal step size.  double error tolerance and
!return, if estimated step size is too small for machine precision.
!
!restore x, phi(*,*) and psi(*).
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
          phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  305   continue
  310 continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
        psi(i-1) = psi(i) - h
  315 continue
!on third failure, set order to one.  thereafter, use optimal step size.
  320 ifail = ifail + 1
      temp2 = 0.5d0
!!    if(ifail - 3) 335,330,325
!!325 if(p5eps .lt. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
!!330 knew = 1
!!335 h = temp2*h
      if(ifail-3>0)then
        if(p5eps .lt. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
        knew = 1
      elseif(ifail-3==0)then
        knew=1
      endif
      h = temp2*h
      k = knew
      if(dabs(h) .ge. fouru*dabs(x)) go to 340
      crash = .true.
      h = dsign(fouru*dabs(x),h)
      eps = eps + eps
      return
  340 go to 100
!end BLOCK 3
!
!begin BLOCK 4
!the step is successful.  correct the predicted solution, evaluate
!the derivatives using the corrected solution and update the
!differences. determine best order and step size for next step.
  400 kold = k
      hold = h
!correct and evaluate.
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
        phi(l,15) = (y(l) - p(l)) - rho
      enddo
      go to 420
  410 do l=1,neqn
        y(l)=p(l)+temp1*(yp(l)-phi(l,1))
      enddo
  420 continue
      call f(x,y,yp)
!update differences for next step.
      do l=1,neqn
        phi(l,kp1)=yp(l)-phi(l,1)
        phi(l,kp2)=phi(l,kp1)-phi(l,kp2)
      enddo
      do i=1,k
        do l=1,neqn
          phi(l,i)=phi(l,i)+phi(l,kp1)
        enddo
      enddo
!estimate error at order k+1 unless:
!in first phase when always raise order,
!already decided to lower order,
!step size not constant so estimate unreliable.
      erkp1 = 0.0d0
      if(knew .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(knew .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
        erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
  440 continue
      erkp1 = absh*gstr(kp1)*dsqrt(erkp1)
!using estimated error at order k+1, determine appropriate order for
!next step.
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5d0*erk) go to 460
      go to 450
  445 if(erkm1 .le. dmin1(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
!here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have been
!lowered in block 2. thus order is to be raised. raise order.
  450 k = kp1
      erk = erkp1
      go to 460
!lower order.
  455 k = km1
      erk = erkm1
!with new order determine appropriate step size for next step.
  460 hnew = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hnew = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*dmax1(0.5d0,dmin1(0.9d0,r))
      hnew = dsign(dmax1(hnew,fouru*dabs(x)),h)
  465 h = hnew
      return
!end BLOCK 4
end subroutine step
!-----------------------------------------------------------------------
subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
!
!the methods in subroutine  step  approximate the solution near  x
!by a polynomial.  subroutine  intrp  approximates the solution at
!xout  by evaluating the polynomial there.
!
!input to intrp:
!   the user provides storage in the calling program for the arrays
!       y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
!   xout     -- point at which solution is desired.
!
!output from  intrp:
!   yout(*)  -- solution at  xout
!   ypout(*) -- derivative of solution at  xout
!the remaining parameters are returned unaltered from their input
!values.  integration with  step  may be continued.
      implicit none
      integer :: neqn,ki,kold,kip1,i,j,jm1,limit1,l
      real(8) :: y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),&
        psi(12),g(13),w(13),rho(13),x,xout,hi,term,temp1,temp2,temp3,&
        psijm1,gamma,eta
      data g(1)/1.0d0/,rho(1)/1.0d0/
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
!initialize w(*) for computing g(*)
      do 5 i = 1,ki
        temp1 = i
        w(i) = 1.0d0/temp1
    5 continue
      term = 0.0d0
!compute g(*)
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
          w(i) = gamma*w(i) - eta*w(i+1)
   10   continue
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
        term = psijm1
   15 continue
!interpolate
      do 20 l = 1,neqn
        ypout(l) = 0.0d0
        yout(l) = 0.0d0
   20 continue
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
          ypout(l) = ypout(l) + temp3*phi(l,i)
   25   continue
   30 continue
      do 35 l = 1,neqn
        yout(l) = y(l) + hi*yout(l)
   35 continue
      return
end subroutine intrp
end subroutine ode

end module maths

!==============================================================================
module space
implicit none
!nieq   = # inequivalent atoms in the structure,
!neq(i) = # equivalent atoms for i=1:nieq,
!nlatp  = # lattice points in the unit cell,
integer :: nieq,nlatp
integer,allocatable :: neq(:),ieq(:)
!rx(i,ir)  )= radial grid rx(i,ir)=r(1,ir)*exp(dx(ir)*(i-1)),
!             where ir=1:nieq and i=1:nx(ir),
!dx(:)   = logarithmic increment in radial grid,
!rmt(:)  = MT radius,
!rmtn(:) = min MT radius,
!rmtx(:) = max MT radius,
!nx(:)   = no. of grid points in atomic calculation,
!nxz(:)  = outermost radial grid point set in ChargeAndPoisson,
real(8),allocatable :: rx(:,:),dx(:),rmt(:),rmt0(:),&
  rmtn(:),rmtx(:),rmt1(:,:)
integer,allocatable :: nx(:),nxz(:)
integer :: nxx
!rc(i,j) =the i'th coordinate of the j'th axis of the unit cell,
!rk(i,j) =the i'th coordinate of the j'th atom in the unit cell,
!volUC   =volume of unit cell (inclusive vacuum slab),
!volXC   =volume of bulk unit cell (filled with crystal electrons).
integer,allocatable :: na(:,:),ia(:,:),ncon(:)
real(8),allocatable :: ad(:,:),rk(:,:)
real(8) :: rcX(3,3),volUC,volXC
end module space

!==============================================================================
module enrgy
use space,only : rx,dx,nxz,ad,ia
real(8) :: V0,mVstp,rsis
real(8),allocatable,dimension(:) :: z,cls,Vstp,Vrmt,V01
real(8),allocatable,dimension(:,:) :: Vtot,VC,VCit,Vxc,QvsR,&
  atrho4pir2,myxc,rs
real(8),allocatable :: Vtot1(:,:,:)
real(8) :: E_DE
!for self-energy data in separate file: next two lines.
!integer :: nsr,nsp
!real(8),allocatable :: sp(:),sr(:),sdat(:,:)
!
!Madelung items.
!ewii = Ewald matrix for inequivalent atoms,
!dpii = dipole moment matrix for inequivalent atoms.
logical :: wrMadl
real(8),allocatable :: ewii(:,:),dpii(:,:)
!Energy grid items.
!eev(1:ne) = energy grid (in eV),
![ee1,ee2] = energy interval for phase shift calculation,
!de1       = first energy step in energy grid,
!dee       = uniform energy step in interpolated phase shift table,
!rMTweight = weight of rerr relative to verr in subroutine errMT,
!dVeV = accuracy of Vrmt in subroutine MTopt.
integer :: ne,ie
real(8)  :: ee1,ee2,dee
real(8),allocatable :: eev(:)
real(8),parameter :: de1=1.d0 !(in eV)
real(8)  :: rMTweight,dVeV,rerr,verr
!Interstitial potential.
!V0vsE1,2 = input files for ATA purpose,
!v0ev     = energy grid of V0.
real(8),allocatable :: v0ev(:),V0ata(:)
!Phase shift quantum numbers etc.
!iatom =atom subscript,
!kappa =spin-orbit coupling quantum number (dirac eq.):
!      =-l-1, l=0,1,2,... and =l, l=1,2,...(dirac eq.),
!      =l, l=0,1,2,... (schroedinger eq.),
!emv0  =Eprimary-V0(inner potential),
!q     =wave number corresponding to emv0,
integer :: iatom,kappa,lmax
real(8) :: emv0,q
!
!Self-energy data by the courtesy of B.E.Sernelius.
!Ref.: K.W.Shung, B.E.Sernelius, and G.D.Mahan, PRB 36,4499 (1987).
integer,parameter :: nsr=17,nsp=151
real(8),dimension(nsr),parameter :: sr=(/&
0.00d0,0.01d0,0.02d0,0.05d0,0.10d0,0.20d0,0.30d0,0.40d0,0.50d0,0.70d0,1.00d0,1.50d0,2.00d0,3.00d0,&
4.00d0,5.00d0,6.00d0/)
real(8),dimension(nsp),parameter :: sp=(/&
0.00d0,0.02d0,0.04d0,0.06d0,0.08d0,0.10d0,0.12d0,0.14d0,0.16d0,0.18d0,0.20d0,0.22d0,0.24d0,0.26d0,&
0.28d0,0.30d0,0.32d0,0.34d0,0.36d0,0.38d0,0.40d0,0.42d0,0.44d0,0.46d0,0.48d0,0.50d0,0.52d0,0.54d0,&
0.56d0,0.58d0,0.60d0,0.62d0,0.64d0,0.66d0,0.68d0,0.70d0,0.72d0,0.74d0,0.76d0,0.78d0,0.80d0,0.82d0,&
0.84d0,0.86d0,0.88d0,0.90d0,0.92d0,0.94d0,0.96d0,0.98d0,1.00d0,1.02d0,1.04d0,1.06d0,1.08d0,1.10d0,&
1.12d0,1.14d0,1.16d0,1.18d0,1.20d0,1.22d0,1.24d0,1.26d0,1.28d0,1.30d0,1.32d0,1.34d0,1.36d0,1.38d0,&
1.40d0,1.42d0,1.44d0,1.46d0,1.48d0,1.50d0,1.52d0,1.54d0,1.56d0,1.58d0,1.60d0,1.62d0,1.64d0,1.66d0,&
1.68d0,1.70d0,1.72d0,1.74d0,1.76d0,1.78d0,1.80d0,1.82d0,1.84d0,1.86d0,1.88d0,1.90d0,1.92d0,1.94d0,&
1.96d0,1.98d0,2.00d0,2.02d0,2.04d0,2.06d0,2.08d0,2.10d0,2.12d0,2.14d0,2.16d0,2.18d0,2.20d0,2.22d0,&
2.24d0,2.26d0,2.28d0,2.30d0,2.32d0,2.34d0,2.36d0,2.38d0,2.40d0,2.42d0,2.44d0,2.46d0,2.48d0,2.50d0,&
2.52d0,2.54d0,2.56d0,2.58d0,2.60d0,2.62d0,2.64d0,2.66d0,2.68d0,2.70d0,2.72d0,2.74d0,2.76d0,2.78d0,&
2.80d0,2.82d0,2.84d0,2.86d0,2.88d0,2.90d0,2.92d0,2.94d0,2.96d0,2.98d0,3.00d0/)
!sdat=reshape(sdat1,shape=(/nsp,nsr/)) in subroutine PartialWaveMethod.
real(8),allocatable,dimension(:,:) :: sdat
real(8),dimension(nsp*nsr),parameter :: sdat1=(/&
-0.66340d0,-0.66330d0,-0.66310d0,-0.66260d0,-0.66200d0,-0.66120d0,-0.66020d0,-0.65910d0,-0.65770d0,-0.65620d0,-0.65450d0,&
-0.65260d0,-0.65050d0,-0.64830d0,-0.64580d0,-0.64320d0,-0.64030d0,-0.63720d0,-0.63400d0,-0.63050d0,-0.62680d0,-0.62290d0,&
-0.61880d0,-0.61450d0,-0.60990d0,-0.60500d0,-0.60000d0,-0.59460d0,-0.58900d0,-0.58310d0,-0.57700d0,-0.57050d0,-0.56370d0,&
-0.55660d0,-0.54920d0,-0.54130d0,-0.53310d0,-0.52450d0,-0.51540d0,-0.50580d0,-0.49570d0,-0.48500d0,-0.47370d0,-0.46160d0,&
-0.44870d0,-0.43480d0,-0.41970d0,-0.40310d0,-0.38440d0,-0.36250d0,-0.33170d0,-0.30140d0,-0.28060d0,-0.26330d0,-0.24850d0,&
-0.23530d0,-0.22350d0,-0.21290d0,-0.20310d0,-0.19420d0,-0.18590d0,-0.17820d0,-0.17110d0,-0.16450d0,-0.15820d0,-0.15240d0,&
-0.14690d0,-0.14180d0,-0.13690d0,-0.13230d0,-0.12790d0,-0.12380d0,-0.11990d0,-0.11620d0,-0.11260d0,-0.10930d0,-0.10610d0,&
-0.10300d0,-0.10010d0,-0.09727d0,-0.09459d0,-0.09203d0,-0.08958d0,-0.08722d0,-0.08496d0,-0.08280d0,-0.08071d0,-0.07870d0,&
-0.07678d0,-0.07493d0,-0.07314d0,-0.07142d0,-0.06977d0,-0.06817d0,-0.06662d0,-0.06513d0,-0.06369d0,-0.06230d0,-0.06096d0,&
-0.05965d0,-0.05840d0,-0.05718d0,-0.05600d0,-0.05486d0,-0.05375d0,-0.05267d0,-0.05163d0,-0.05062d0,-0.04964d0,-0.04869d0,&
-0.04777d0,-0.04687d0,-0.04600d0,-0.04515d0,-0.04433d0,-0.04353d0,-0.04275d0,-0.04199d0,-0.04125d0,-0.04054d0,-0.03984d0,&
-0.03916d0,-0.03849d0,-0.03785d0,-0.03722d0,-0.03660d0,-0.03600d0,-0.03542d0,-0.03485d0,-0.03429d0,-0.03375d0,-0.03322d0,&
-0.03270d0,-0.03220d0,-0.03170d0,-0.03122d0,-0.03075d0,-0.03029d0,-0.02984d0,-0.02940d0,-0.02897d0,-0.02855d0,-0.02814d0,&
-0.02773d0,-0.02734d0,-0.02696d0,-0.02658d0,-0.02621d0,-0.02585d0,-0.02549d0,-0.02515d0,-0.61980d0,-0.61970d0,-0.61950d0,&
-0.61920d0,-0.61870d0,-0.61800d0,-0.61720d0,-0.61620d0,-0.61500d0,-0.61360d0,-0.61210d0,-0.61040d0,-0.60840d0,-0.60630d0,&
-0.60400d0,-0.60150d0,-0.59880d0,-0.59590d0,-0.59280d0,-0.58950d0,-0.58600d0,-0.58220d0,-0.57830d0,-0.57410d0,-0.56970d0,&
-0.56500d0,-0.56010d0,-0.55500d0,-0.54960d0,-0.54390d0,-0.53790d0,-0.53170d0,-0.52510d0,-0.51820d0,-0.51100d0,-0.50340d0,&
-0.49550d0,-0.48710d0,-0.47830d0,-0.46910d0,-0.45940d0,-0.44910d0,-0.43830d0,-0.42690d0,-0.41490d0,-0.40230d0,-0.38920d0,&
-0.37560d0,-0.36180d0,-0.34750d0,-0.33290d0,-0.31990d0,-0.30740d0,-0.29630d0,-0.28700d0,-0.27060d0,-0.25900d0,-0.24830d0,&
-0.23830d0,-0.22900d0,-0.22050d0,-0.21250d0,-0.20500d0,-0.19800d0,-0.19150d0,-0.18530d0,-0.17950d0,-0.17400d0,-0.16880d0,&
-0.16390d0,-0.15930d0,-0.15490d0,-0.15060d0,-0.14660d0,-0.14280d0,-0.13920d0,-0.13570d0,-0.13240d0,-0.12920d0,-0.12620d0,&
-0.12330d0,-0.12050d0,-0.11780d0,-0.11520d0,-0.11270d0,-0.11030d0,-0.10800d0,-0.10580d0,-0.10370d0,-0.10170d0,-0.09967d0,&
-0.09776d0,-0.09591d0,-0.09412d0,-0.09240d0,-0.09073d0,-0.08911d0,-0.08754d0,-0.08603d0,-0.08456d0,-0.08314d0,-0.08176d0,&
-0.08042d0,-0.07912d0,-0.07786d0,-0.07663d0,-0.07544d0,-0.07428d0,-0.07316d0,-0.07206d0,-0.07100d0,-0.06996d0,-0.06895d0,&
-0.06797d0,-0.06701d0,-0.06608d0,-0.06517d0,-0.06429d0,-0.06343d0,-0.06259d0,-0.06177d0,-0.06097d0,-0.06018d0,-0.05942d0,&
-0.05860d0,-0.05794d0,-0.05723d0,-0.05653d0,-0.05585d0,-0.05518d0,-0.05453d0,-0.05389d0,-0.05327d0,-0.05266d0,-0.05206d0,&
-0.05147d0,-0.05090d0,-0.05034d0,-0.04979d0,-0.04925d0,-0.04873d0,-0.04821d0,-0.04773d0,-0.04723d0,-0.04675d0,-0.04628d0,&
-0.04581d0,-0.04536d0,-0.04493d0,-0.04450d0,-0.04407d0,-0.60270d0,-0.60270d0,-0.60250d0,-0.60210d0,-0.60170d0,-0.60110d0,&
-0.60030d0,-0.59930d0,-0.59820d0,-0.59690d0,-0.59540d0,-0.59380d0,-0.59190d0,-0.58990d0,-0.58770d0,-0.58530d0,-0.58270d0,&
-0.57980d0,-0.57680d0,-0.57360d0,-0.57020d0,-0.56650d0,-0.56270d0,-0.55860d0,-0.55430d0,-0.54970d0,-0.54490d0,-0.53990d0,&
-0.53450d0,-0.52900d0,-0.52310d0,-0.51700d0,-0.51050d0,-0.50380d0,-0.49670d0,-0.48930d0,-0.48150d0,-0.47330d0,-0.46480d0,&
-0.45580d0,-0.44640d0,-0.43660d0,-0.42630d0,-0.41570d0,-0.40470d0,-0.39340d0,-0.38190d0,-0.37030d0,-0.35840d0,-0.34630d0,&
-0.33390d0,-0.32240d0,-0.31130d0,-0.30100d0,-0.29130d0,-0.28190d0,-0.27220d0,-0.26220d0,-0.25250d0,-0.24330d0,-0.23470d0,&
-0.22660d0,-0.21900d0,-0.21190d0,-0.20520d0,-0.19880d0,-0.19290d0,-0.18730d0,-0.18190d0,-0.17690d0,-0.17210d0,-0.16750d0,&
-0.16320d0,-0.15900d0,-0.15510d0,-0.15130d0,-0.14770d0,-0.14430d0,-0.14100d0,-0.13790d0,-0.13480d0,-0.13190d0,-0.12910d0,&
-0.12650d0,-0.12390d0,-0.12140d0,-0.11900d0,-0.11670d0,-0.11450d0,-0.11230d0,-0.11030d0,-0.10830d0,-0.10630d0,-0.10450d0,&
-0.10270d0,-0.10090d0,-0.09924d0,-0.09760d0,-0.09601d0,-0.09440d0,-0.09297d0,-0.09150d0,-0.09012d0,-0.08875d0,-0.08742d0,&
-0.08613d0,-0.08488d0,-0.08366d0,-0.08247d0,-0.08132d0,-0.08019d0,-0.07910d0,-0.07803d0,-0.07699d0,-0.07598d0,-0.07499d0,&
-0.07403d0,-0.07309d0,-0.07218d0,-0.07128d0,-0.07041d0,-0.06956d0,-0.06872d0,-0.06791d0,-0.06712d0,-0.06634d0,-0.06558d0,&
-0.06484d0,-0.06411d0,-0.06340d0,-0.06271d0,-0.06203d0,-0.06136d0,-0.06071d0,-0.06008d0,-0.05945d0,-0.05884d0,-0.05824d0,&
-0.05765d0,-0.05707d0,-0.05651d0,-0.05596d0,-0.05541d0,-0.05488d0,-0.05435d0,-0.05384d0,-0.05333d0,-0.05284d0,-0.05235d0,&
-0.05187d0,-0.05141d0,-0.57140d0,-0.57130d0,-0.57110d0,-0.57080d0,-0.57040d0,-0.56980d0,-0.56910d0,-0.56820d0,-0.56720d0,&
-0.56600d0,-0.56470d0,-0.56310d0,-0.56140d0,-0.55950d0,-0.55750d0,-0.55520d0,-0.55280d0,-0.55010d0,-0.54730d0,-0.54430d0,&
-0.54100d0,-0.53760d0,-0.53390d0,-0.53000d0,-0.52590d0,-0.52160d0,-0.51700d0,-0.51220d0,-0.50720d0,-0.50190d0,-0.49630d0,&
-0.49050d0,-0.48440d0,-0.47800d0,-0.47130d0,-0.46440d0,-0.45710d0,-0.44960d0,-0.44180d0,-0.43380d0,-0.42560d0,-0.41710d0,&
-0.40860d0,-0.39990d0,-0.39120d0,-0.38240d0,-0.37350d0,-0.36440d0,-0.35530d0,-0.34590d0,-0.33650d0,-0.32740d0,-0.31860d0,&
-0.31020d0,-0.30230d0,-0.29460d0,-0.28780d0,-0.28130d0,-0.27630d0,-0.26990d0,-0.26200d0,-0.25410d0,-0.24650d0,-0.23920d0,&
-0.23230d0,-0.22580d0,-0.21950d0,-0.21360d0,-0.20800d0,-0.20260d0,-0.19750d0,-0.19270d0,-0.18800d0,-0.18360d0,-0.17940d0,&
-0.17540d0,-0.17150d0,-0.16780d0,-0.16430d0,-0.16090d0,-0.15760d0,-0.15450d0,-0.15150d0,-0.14860d0,-0.14580d0,-0.14320d0,&
-0.14060d0,-0.13810d0,-0.13570d0,-0.13330d0,-0.13110d0,-0.12890d0,-0.12680d0,-0.12480d0,-0.12280d0,-0.12090d0,-0.11910d0,&
-0.11730d0,-0.11550d0,-0.11390d0,-0.11220d0,-0.11060d0,-0.10910d0,-0.10760d0,-0.10610d0,-0.10470d0,-0.10330d0,-0.10200d0,&
-0.10070d0,-0.09939d0,-0.09815d0,-0.09694d0,-0.09575d0,-0.09460d0,-0.09348d0,-0.09238d0,-0.09131d0,-0.09026d0,-0.08924d0,&
-0.08825d0,-0.08727d0,-0.08632d0,-0.08539d0,-0.08448d0,-0.08359d0,-0.08272d0,-0.08186d0,-0.08103d0,-0.08021d0,-0.07941d0,&
-0.07863d0,-0.07786d0,-0.07711d0,-0.07638d0,-0.07566d0,-0.07495d0,-0.07425d0,-0.07357d0,-0.07290d0,-0.07225d0,-0.07161d0,&
-0.07098d0,-0.07036d0,-0.06975d0,-0.06915d0,-0.06857d0,-0.06799d0,-0.06743d0,-0.06687d0,-0.06633d0,-0.06579d0,-0.54000d0,&
-0.54000d0,-0.53980d0,-0.53950d0,-0.53910d0,-0.53850d0,-0.53780d0,-0.53700d0,-0.53600d0,-0.53490d0,-0.53370d0,-0.53230d0,&
-0.53070d0,-0.52890d0,-0.52700d0,-0.52490d0,-0.52270d0,-0.52020d0,-0.51760d0,-0.51480d0,-0.51170d0,-0.50850d0,-0.50510d0,&
-0.50150d0,-0.49770d0,-0.49370d0,-0.48940d0,-0.48500d0,-0.48030d0,-0.47540d0,-0.47030d0,-0.46500d0,-0.45950d0,-0.45380d0,&
-0.44790d0,-0.44180d0,-0.43560d0,-0.42920d0,-0.42270d0,-0.41620d0,-0.40950d0,-0.40290d0,-0.39620d0,-0.38940d0,-0.38260d0,&
-0.37580d0,-0.36880d0,-0.36180d0,-0.35470d0,-0.34740d0,-0.34010d0,-0.33300d0,-0.32610d0,-0.31940d0,-0.31310d0,-0.30700d0,&
-0.30140d0,-0.29600d0,-0.29120d0,-0.28650d0,-0.28240d0,-0.28030d0,-0.27580d0,-0.26920d0,-0.26250d0,-0.25590d0,-0.24950d0,&
-0.24330d0,-0.23740d0,-0.23170d0,-0.22630d0,-0.22110d0,-0.21610d0,-0.21140d0,-0.20680d0,-0.20250d0,-0.19830d0,-0.19430d0,&
-0.19050d0,-0.18680d0,-0.18330d0,-0.17990d0,-0.17660d0,-0.17350d0,-0.17050d0,-0.16750d0,-0.16470d0,-0.16200d0,-0.15940d0,&
-0.15680d0,-0.15440d0,-0.15200d0,-0.14970d0,-0.14750d0,-0.14530d0,-0.14320d0,-0.14120d0,-0.13920d0,-0.13730d0,-0.13550d0,&
-0.13370d0,-0.13190d0,-0.13020d0,-0.12850d0,-0.12690d0,-0.12530d0,-0.12380d0,-0.12230d0,-0.12090d0,-0.11950d0,-0.11810d0,&
-0.11670d0,-0.11540d0,-0.11410d0,-0.11290d0,-0.11170d0,-0.11050d0,-0.10930d0,-0.10820d0,-0.10700d0,-0.10590d0,-0.10490d0,&
-0.10380d0,-0.10280d0,-0.10180d0,-0.10080d0,-0.09988d0,-0.09894d0,-0.09802d0,-0.09712d0,-0.09623d0,-0.09537d0,-0.09452d0,&
-0.09369d0,-0.09287d0,-0.09207d0,-0.09129d0,-0.09051d0,-0.08976d0,-0.08902d0,-0.08829d0,-0.08757d0,-0.08687d0,-0.08618d0,&
-0.08550d0,-0.08483d0,-0.08417d0,-0.08353d0,-0.08289d0,-0.08227d0,-0.08166d0,-0.50230d0,-0.50220d0,-0.50210d0,-0.50180d0,&
-0.50140d0,-0.50090d0,-0.50020d0,-0.49940d0,-0.49860d0,-0.49750d0,-0.49640d0,-0.49510d0,-0.49370d0,-0.49210d0,-0.49040d0,&
-0.48860d0,-0.48650d0,-0.48440d0,-0.48200d0,-0.47960d0,-0.47690d0,-0.47410d0,-0.47110d0,-0.46800d0,-0.46470d0,-0.46120d0,&
-0.45760d0,-0.45390d0,-0.44990d0,-0.44590d0,-0.44170d0,-0.43750d0,-0.43310d0,-0.42860d0,-0.42410d0,-0.41950d0,-0.41490d0,&
-0.41020d0,-0.40560d0,-0.40090d0,-0.39610d0,-0.39140d0,-0.38660d0,-0.38180d0,-0.37690d0,-0.37200d0,-0.36700d0,-0.36190d0,&
-0.35670d0,-0.35150d0,-0.34620d0,-0.34100d0,-0.33600d0,-0.33110d0,-0.32640d0,-0.32190d0,-0.31770d0,-0.31370d0,-0.30990d0,&
-0.30640d0,-0.30310d0,-0.30040d0,-0.29800d0,-0.29590d0,-0.29480d0,-0.29560d0,-0.28990d0,-0.28420d0,-0.27840d0,-0.27260d0,&
-0.26680d0,-0.26130d0,-0.25590d0,-0.25070d0,-0.24570d0,-0.24090d0,-0.23630d0,-0.23190d0,-0.22760d0,-0.22350d0,-0.21960d0,&
-0.21580d0,-0.21210d0,-0.20860d0,-0.20520d0,-0.20190d0,-0.19880d0,-0.19570d0,-0.19270d0,-0.18990d0,-0.18710d0,-0.18440d0,&
-0.18180d0,-0.17930d0,-0.17690d0,-0.17450d0,-0.17220d0,-0.17000d0,-0.16780d0,-0.16570d0,-0.16370d0,-0.16170d0,-0.15970d0,&
-0.15780d0,-0.15600d0,-0.15420d0,-0.15250d0,-0.15080d0,-0.14910d0,-0.14750d0,-0.14590d0,-0.14440d0,-0.14280d0,-0.14140d0,&
-0.13990d0,-0.13850d0,-0.13720d0,-0.13580d0,-0.13450d0,-0.13320d0,-0.13200d0,-0.13070d0,-0.12950d0,-0.12830d0,-0.12720d0,&
-0.12600d0,-0.12490d0,-0.12380d0,-0.12280d0,-0.12170d0,-0.12070d0,-0.11970d0,-0.11870d0,-0.11770d0,-0.11680d0,-0.11590d0,&
-0.11490d0,-0.11400d0,-0.11320d0,-0.11230d0,-0.11140d0,-0.11060d0,-0.10980d0,-0.10900d0,-0.10820d0,-0.10740d0,-0.10660d0,&
-0.10590d0,-0.10510d0,-0.10440d0,-0.10370d0,-0.47820d0,-0.47810d0,-0.47800d0,-0.47770d0,-0.47730d0,-0.47680d0,-0.47620d0,&
-0.47550d0,-0.47470d0,-0.47370d0,-0.47270d0,-0.47150d0,-0.47020d0,-0.46880d0,-0.46720d0,-0.46560d0,-0.46380d0,-0.46180d0,&
-0.45980d0,-0.45760d0,-0.45530d0,-0.45280d0,-0.45030d0,-0.44760d0,-0.44480d0,-0.44190d0,-0.43880d0,-0.43570d0,-0.43250d0,&
-0.42920d0,-0.42590d0,-0.42250d0,-0.41900d0,-0.41550d0,-0.41200d0,-0.40850d0,-0.40490d0,-0.40130d0,-0.39780d0,-0.39410d0,&
-0.39050d0,-0.38690d0,-0.38320d0,-0.37940d0,-0.37560d0,-0.37180d0,-0.36780d0,-0.36380d0,-0.35980d0,-0.35560d0,-0.35150d0,&
-0.34740d0,-0.34340d0,-0.33950d0,-0.33580d0,-0.33220d0,-0.32880d0,-0.32560d0,-0.32260d0,-0.31970d0,-0.31700d0,-0.31480d0,&
-0.31270d0,-0.31080d0,-0.30940d0,-0.30850d0,-0.30780d0,-0.30900d0,-0.30780d0,-0.30260d0,-0.29710d0,-0.29160d0,-0.28610d0,&
-0.28070d0,-0.27540d0,-0.27030d0,-0.26540d0,-0.26060d0,-0.25610d0,-0.25160d0,-0.24740d0,-0.24330d0,-0.23930d0,-0.23550d0,&
-0.23180d0,-0.22830d0,-0.22480d0,-0.22150d0,-0.21830d0,-0.21520d0,-0.21210d0,-0.20920d0,-0.20640d0,-0.20370d0,-0.20100d0,&
-0.19840d0,-0.19590d0,-0.19350d0,-0.19110d0,-0.18880d0,-0.18650d0,-0.18440d0,-0.18220d0,-0.18020d0,-0.17810d0,-0.17620d0,&
-0.17430d0,-0.17240d0,-0.17060d0,-0.16880d0,-0.16700d0,-0.16530d0,-0.16370d0,-0.16210d0,-0.16050d0,-0.15890d0,-0.15740d0,&
-0.15590d0,-0.15450d0,-0.15310d0,-0.15170d0,-0.15030d0,-0.14900d0,-0.14770d0,-0.14640d0,-0.14520d0,-0.14390d0,-0.14270d0,&
-0.14150d0,-0.14040d0,-0.13920d0,-0.13810d0,-0.13700d0,-0.13600d0,-0.13490d0,-0.13390d0,-0.13280d0,-0.13180d0,-0.13090d0,&
-0.12990d0,-0.12890d0,-0.12800d0,-0.12710d0,-0.12620d0,-0.12530d0,-0.12440d0,-0.12360d0,-0.12270d0,-0.12190d0,-0.12110d0,&
-0.12030d0,-0.46080d0,-0.46080d0,-0.46060d0,-0.46040d0,-0.46000d0,-0.45960d0,-0.45900d0,-0.45830d0,-0.45760d0,-0.45670d0,&
-0.45580d0,-0.45470d0,-0.45350d0,-0.45220d0,-0.45090d0,-0.44940d0,-0.44780d0,-0.44610d0,-0.44430d0,-0.44240d0,-0.44040d0,&
-0.43830d0,-0.43610d0,-0.43380d0,-0.43140d0,-0.42900d0,-0.42650d0,-0.42390d0,-0.42130d0,-0.41870d0,-0.41600d0,-0.41320d0,&
-0.41050d0,-0.40770d0,-0.40490d0,-0.40210d0,-0.39930d0,-0.39650d0,-0.39370d0,-0.39080d0,-0.38790d0,-0.38500d0,-0.38200d0,&
-0.37900d0,-0.37590d0,-0.37280d0,-0.36960d0,-0.36640d0,-0.36300d0,-0.35970d0,-0.35620d0,-0.35290d0,-0.34960d0,-0.34640d0,&
-0.34330d0,-0.34040d0,-0.33760d0,-0.33490d0,-0.33240d0,-0.33000d0,-0.32780d0,-0.32590d0,-0.32420d0,-0.32260d0,-0.32130d0,&
-0.32040d0,-0.31960d0,-0.31940d0,-0.31990d0,-0.32140d0,-0.32150d0,-0.31610d0,-0.31090d0,-0.30550d0,-0.30020d0,-0.29490d0,&
-0.28980d0,-0.28480d0,-0.28000d0,-0.27530d0,-0.27080d0,-0.26640d0,-0.26220d0,-0.25820d0,-0.25420d0,-0.25040d0,-0.24680d0,&
-0.24320d0,-0.23980d0,-0.23640d0,-0.23320d0,-0.23010d0,-0.22700d0,-0.22410d0,-0.22130d0,-0.21850d0,-0.21580d0,-0.21320d0,&
-0.21060d0,-0.20810d0,-0.20570d0,-0.20340d0,-0.20110d0,-0.19890d0,-0.19670d0,-0.19460d0,-0.19250d0,-0.19050d0,-0.18860d0,&
-0.18660d0,-0.18480d0,-0.18290d0,-0.18110d0,-0.17940d0,-0.17770d0,-0.17600d0,-0.17440d0,-0.17280d0,-0.17120d0,-0.16970d0,&
-0.16820d0,-0.16670d0,-0.16530d0,-0.16390d0,-0.16250d0,-0.16110d0,-0.15980d0,-0.15850d0,-0.15720d0,-0.15600d0,-0.15470d0,&
-0.15350d0,-0.15230d0,-0.15120d0,-0.15000d0,-0.14890d0,-0.14780d0,-0.14670d0,-0.14560d0,-0.14460d0,-0.14360d0,-0.14260d0,&
-0.14160d0,-0.14060d0,-0.13960d0,-0.13870d0,-0.13770d0,-0.13680d0,-0.13590d0,-0.13500d0,-0.13410d0,-0.44770d0,-0.44760d0,&
-0.44750d0,-0.44730d0,-0.44690d0,-0.44650d0,-0.44600d0,-0.44540d0,-0.44470d0,-0.44390d0,-0.44310d0,-0.44210d0,-0.44100d0,&
-0.43990d0,-0.43870d0,-0.43740d0,-0.43600d0,-0.43450d0,-0.43290d0,-0.43130d0,-0.42960d0,-0.42780d0,-0.42590d0,-0.42400d0,&
-0.42210d0,-0.42000d0,-0.41800d0,-0.41590d0,-0.41380d0,-0.41160d0,-0.40940d0,-0.40720d0,-0.40500d0,-0.40280d0,-0.40060d0,&
-0.39830d0,-0.39610d0,-0.39380d0,-0.39150d0,-0.38920d0,-0.38690d0,-0.38450d0,-0.38210d0,-0.37960d0,-0.37700d0,-0.37450d0,&
-0.37180d0,-0.36910d0,-0.36630d0,-0.36350d0,-0.36060d0,-0.35780d0,-0.35510d0,-0.35240d0,-0.34980d0,-0.34730d0,-0.34500d0,&
-0.34270d0,-0.34060d0,-0.33860d0,-0.33680d0,-0.33510d0,-0.33370d0,-0.33230d0,-0.33130d0,-0.33050d0,-0.32970d0,-0.32950d0,&
-0.32960d0,-0.32980d0,-0.33100d0,-0.33440d0,-0.33210d0,-0.32670d0,-0.32150d0,-0.31630d0,-0.31110d0,-0.30600d0,-0.30100d0,&
-0.29610d0,-0.29140d0,-0.28680d0,-0.28240d0,-0.27810d0,-0.27400d0,-0.27000d0,-0.26610d0,-0.26230d0,-0.25870d0,-0.25520d0,&
-0.25170d0,-0.24840d0,-0.24520d0,-0.24210d0,-0.23910d0,-0.23610d0,-0.23330d0,-0.23050d0,-0.22780d0,-0.22520d0,-0.22260d0,&
-0.22010d0,-0.21770d0,-0.21530d0,-0.21300d0,-0.21080d0,-0.20860d0,-0.20640d0,-0.20430d0,-0.20230d0,-0.20030d0,-0.19840d0,&
-0.19650d0,-0.19460d0,-0.19280d0,-0.19100d0,-0.18930d0,-0.18760d0,-0.18590d0,-0.18430d0,-0.18270d0,-0.18110d0,-0.17950d0,&
-0.17800d0,-0.17660d0,-0.17510d0,-0.17370d0,-0.17230d0,-0.17090d0,-0.16960d0,-0.16830d0,-0.16700d0,-0.16570d0,-0.16450d0,&
-0.16330d0,-0.16210d0,-0.16090d0,-0.15970d0,-0.15860d0,-0.15750d0,-0.15640d0,-0.15530d0,-0.15420d0,-0.15320d0,-0.15210d0,&
-0.15110d0,-0.15010d0,-0.14910d0,-0.14810d0,-0.14720d0,-0.14630d0,-0.42930d0,-0.42920d0,-0.42910d0,-0.42890d0,-0.42870d0,&
-0.42830d0,-0.42790d0,-0.42740d0,-0.42690d0,-0.42620d0,-0.42560d0,-0.42480d0,-0.42400d0,-0.42310d0,-0.42220d0,-0.42120d0,&
-0.42010d0,-0.41910d0,-0.41790d0,-0.41670d0,-0.41550d0,-0.41430d0,-0.41300d0,-0.41170d0,-0.41030d0,-0.40900d0,-0.40760d0,&
-0.40620d0,-0.40480d0,-0.40340d0,-0.40200d0,-0.40060d0,-0.39920d0,-0.39770d0,-0.39630d0,-0.39480d0,-0.39330d0,-0.39180d0,&
-0.39030d0,-0.38880d0,-0.38720d0,-0.38550d0,-0.38390d0,-0.38220d0,-0.38040d0,-0.37860d0,-0.37670d0,-0.37470d0,-0.37270d0,&
-0.37070d0,-0.36860d0,-0.36660d0,-0.36460d0,-0.36260d0,-0.36070d0,-0.35890d0,-0.35720d0,-0.35560d0,-0.35400d0,-0.35260d0,&
-0.35130d0,-0.35010d0,-0.34910d0,-0.34810d0,-0.34740d0,-0.34680d0,-0.34630d0,-0.34620d0,-0.34620d0,-0.34630d0,-0.34670d0,&
-0.34740d0,-0.34850d0,-0.35060d0,-0.35660d0,-0.35290d0,-0.34720d0,-0.34200d0,-0.33690d0,-0.33180d0,-0.32680d0,-0.32200d0,&
-0.31720d0,-0.31260d0,-0.30810d0,-0.30370d0,-0.29950d0,-0.29540d0,-0.29140d0,-0.28760d0,-0.28380d0,-0.28020d0,-0.27670d0,&
-0.27320d0,-0.26990d0,-0.26670d0,-0.26350d0,-0.26050d0,-0.25750d0,-0.25460d0,-0.25180d0,-0.24900d0,-0.24630d0,-0.24370d0,&
-0.24120d0,-0.23870d0,-0.23630d0,-0.23390d0,-0.23160d0,-0.22930d0,-0.22710d0,-0.22500d0,-0.22290d0,-0.22080d0,-0.21880d0,&
-0.21680d0,-0.21490d0,-0.21300d0,-0.21120d0,-0.20940d0,-0.20760d0,-0.20590d0,-0.20410d0,-0.20250d0,-0.20080d0,-0.19920d0,&
-0.19770d0,-0.19610d0,-0.19460d0,-0.19310d0,-0.19160d0,-0.19020d0,-0.18880d0,-0.18740d0,-0.18600d0,-0.18470d0,-0.18340d0,&
-0.18210d0,-0.18080d0,-0.17960d0,-0.17840d0,-0.17710d0,-0.17600d0,-0.17480d0,-0.17360d0,-0.17250d0,-0.17140d0,-0.17030d0,&
-0.16920d0,-0.16810d0,-0.16710d0,-0.41320d0,-0.41320d0,-0.41310d0,-0.41300d0,-0.41280d0,-0.41260d0,-0.41230d0,-0.41200d0,&
-0.41160d0,-0.41120d0,-0.41080d0,-0.41030d0,-0.40980d0,-0.40920d0,-0.40870d0,-0.40810d0,-0.40750d0,-0.40680d0,-0.40620d0,&
-0.40550d0,-0.40490d0,-0.40420d0,-0.40350d0,-0.40280d0,-0.40210d0,-0.40150d0,-0.40080d0,-0.40010d0,-0.39940d0,-0.39870d0,&
-0.39800d0,-0.39730d0,-0.39660d0,-0.39590d0,-0.39520d0,-0.39450d0,-0.39370d0,-0.39290d0,-0.39210d0,-0.39130d0,-0.39040d0,&
-0.38950d0,-0.38850d0,-0.38750d0,-0.38650d0,-0.38540d0,-0.38420d0,-0.38300d0,-0.38180d0,-0.38050d0,-0.37920d0,-0.37790d0,&
-0.37660d0,-0.37530d0,-0.37410d0,-0.37290d0,-0.37190d0,-0.37080d0,-0.36990d0,-0.36900d0,-0.36820d0,-0.36750d0,-0.36690d0,&
-0.36640d0,-0.36600d0,-0.36580d0,-0.36560d0,-0.36560d0,-0.36580d0,-0.36600d0,-0.36650d0,-0.36700d0,-0.36770d0,-0.36860d0,&
-0.36980d0,-0.37140d0,-0.37390d0,-0.37810d0,-0.38520d0,-0.37700d0,-0.37110d0,-0.36580d0,-0.36070d0,-0.35570d0,-0.35090d0,&
-0.34610d0,-0.34150d0,-0.33700d0,-0.33260d0,-0.32830d0,-0.32410d0,-0.32010d0,-0.31620d0,-0.31240d0,-0.30860d0,-0.30500d0,&
-0.30150d0,-0.29810d0,-0.29480d0,-0.29150d0,-0.28840d0,-0.28530d0,-0.28230d0,-0.27940d0,-0.27650d0,-0.27370d0,-0.27100d0,&
-0.26840d0,-0.26580d0,-0.26320d0,-0.26080d0,-0.25830d0,-0.25600d0,-0.25370d0,-0.25140d0,-0.24920d0,-0.24700d0,-0.24490d0,&
-0.24280d0,-0.24080d0,-0.23880d0,-0.23680d0,-0.23490d0,-0.23300d0,-0.23120d0,-0.22940d0,-0.22760d0,-0.22580d0,-0.22410d0,&
-0.22240d0,-0.22080d0,-0.21920d0,-0.21760d0,-0.21600d0,-0.21450d0,-0.21300d0,-0.21150d0,-0.21000d0,-0.20860d0,-0.20720d0,&
-0.20580d0,-0.20440d0,-0.20310d0,-0.20180d0,-0.20040d0,-0.19920d0,-0.19790d0,-0.19670d0,-0.19540d0,-0.19420d0,-0.19300d0,&
-0.40190d0,-0.40190d0,-0.40190d0,-0.40180d0,-0.40180d0,-0.40170d0,-0.40160d0,-0.40150d0,-0.40130d0,-0.40120d0,-0.40110d0,&
-0.40090d0,-0.40080d0,-0.40060d0,-0.40050d0,-0.40040d0,-0.40020d0,-0.40010d0,-0.40000d0,-0.39990d0,-0.39980d0,-0.39970d0,&
-0.39960d0,-0.39960d0,-0.39950d0,-0.39950d0,-0.39940d0,-0.39940d0,-0.39940d0,-0.39940d0,-0.39930d0,-0.39930d0,-0.39930d0,&
-0.39920d0,-0.39920d0,-0.39910d0,-0.39900d0,-0.39890d0,-0.39880d0,-0.39860d0,-0.39840d0,-0.39820d0,-0.39790d0,-0.39760d0,&
-0.39720d0,-0.39680d0,-0.39640d0,-0.39590d0,-0.39540d0,-0.39480d0,-0.39420d0,-0.39360d0,-0.39310d0,-0.39250d0,-0.39200d0,&
-0.39150d0,-0.39100d0,-0.39060d0,-0.39020d0,-0.38990d0,-0.38970d0,-0.38950d0,-0.38940d0,-0.38930d0,-0.38930d0,-0.38950d0,&
-0.38960d0,-0.39000d0,-0.39040d0,-0.39090d0,-0.39150d0,-0.39210d0,-0.39280d0,-0.39360d0,-0.39460d0,-0.39560d0,-0.39680d0,&
-0.39820d0,-0.40000d0,-0.40220d0,-0.40520d0,-0.40980d0,-0.42080d0,-0.41930d0,-0.41160d0,-0.40520d0,-0.39960d0,-0.39430d0,&
-0.38920d0,-0.38430d0,-0.37950d0,-0.37490d0,-0.37040d0,-0.36610d0,-0.36180d0,-0.35770d0,-0.35360d0,-0.34970d0,-0.34590d0,&
-0.34210d0,-0.33850d0,-0.33500d0,-0.33150d0,-0.32810d0,-0.32480d0,-0.32160d0,-0.31850d0,-0.31540d0,-0.31240d0,-0.30950d0,&
-0.30660d0,-0.30390d0,-0.30110d0,-0.29840d0,-0.29580d0,-0.29330d0,-0.29080d0,-0.28830d0,-0.28590d0,-0.28350d0,-0.28120d0,&
-0.27900d0,-0.27670d0,-0.27460d0,-0.27240d0,-0.27030d0,-0.26830d0,-0.26630d0,-0.26430d0,-0.26230d0,-0.26040d0,-0.25850d0,&
-0.25670d0,-0.25490d0,-0.25310d0,-0.25130d0,-0.24960d0,-0.24790d0,-0.24630d0,-0.24460d0,-0.24300d0,-0.24140d0,-0.23990d0,&
-0.23830d0,-0.23680d0,-0.23530d0,-0.23380d0,-0.23240d0,-0.23100d0,-0.22960d0,-0.22820d0,-0.39990d0,-0.39990d0,-0.39990d0,&
-0.39990d0,-0.39990d0,-0.39990d0,-0.39990d0,-0.40000d0,-0.40000d0,-0.40000d0,-0.40010d0,-0.40020d0,-0.40020d0,-0.40030d0,&
-0.40040d0,-0.40060d0,-0.40070d0,-0.40090d0,-0.40110d0,-0.40130d0,-0.40150d0,-0.40170d0,-0.40200d0,-0.40230d0,-0.40260d0,&
-0.40290d0,-0.40320d0,-0.40350d0,-0.40380d0,-0.40420d0,-0.40450d0,-0.40480d0,-0.40520d0,-0.40550d0,-0.40580d0,-0.40610d0,&
-0.40640d0,-0.40660d0,-0.40690d0,-0.40710d0,-0.40720d0,-0.40740d0,-0.40750d0,-0.40760d0,-0.40770d0,-0.40770d0,-0.40770d0,&
-0.40760d0,-0.40750d0,-0.40740d0,-0.40720d0,-0.40710d0,-0.40690d0,-0.40680d0,-0.40670d0,-0.40660d0,-0.40650d0,-0.40650d0,&
-0.40650d0,-0.40650d0,-0.40660d0,-0.40670d0,-0.40690d0,-0.40720d0,-0.40750d0,-0.40780d0,-0.40820d0,-0.40880d0,-0.40940d0,&
-0.41000d0,-0.41080d0,-0.41160d0,-0.41240d0,-0.41330d0,-0.41420d0,-0.41520d0,-0.41630d0,-0.41750d0,-0.41890d0,-0.42050d0,&
-0.42220d0,-0.42430d0,-0.42680d0,-0.43000d0,-0.43430d0,-0.44120d0,-0.45510d0,-0.44820d0,-0.44070d0,-0.43400d0,-0.42790d0,&
-0.42230d0,-0.41700d0,-0.41190d0,-0.40700d0,-0.40230d0,-0.39780d0,-0.39330d0,-0.38900d0,-0.38480d0,-0.38080d0,-0.37680d0,&
-0.37290d0,-0.36910d0,-0.36550d0,-0.36190d0,-0.35840d0,-0.35490d0,-0.35160d0,-0.34830d0,-0.34510d0,-0.34200d0,-0.33890d0,&
-0.33590d0,-0.33300d0,-0.33010d0,-0.32730d0,-0.32460d0,-0.32190d0,-0.31930d0,-0.31670d0,-0.31410d0,-0.31170d0,-0.30920d0,&
-0.30680d0,-0.30450d0,-0.30220d0,-0.29990d0,-0.29770d0,-0.29550d0,-0.29340d0,-0.29130d0,-0.28920d0,-0.28720d0,-0.28520d0,&
-0.28320d0,-0.28130d0,-0.27940d0,-0.27750d0,-0.27570d0,-0.27390d0,-0.27210d0,-0.27040d0,-0.26870d0,-0.26700d0,-0.26530d0,&
-0.26370d0,-0.26200d0,-0.26040d0,-0.25890d0,-0.25730d0,-0.40670d0,-0.40670d0,-0.40680d0,-0.40680d0,-0.40690d0,-0.40700d0,&
-0.40720d0,-0.40730d0,-0.40750d0,-0.40780d0,-0.40800d0,-0.40830d0,-0.40860d0,-0.40890d0,-0.40930d0,-0.40970d0,-0.41010d0,&
-0.41060d0,-0.41100d0,-0.41150d0,-0.41210d0,-0.41260d0,-0.41320d0,-0.41380d0,-0.41440d0,-0.41500d0,-0.41570d0,-0.41630d0,&
-0.41700d0,-0.41760d0,-0.41830d0,-0.41900d0,-0.41970d0,-0.42030d0,-0.42100d0,-0.42160d0,-0.42230d0,-0.42290d0,-0.42360d0,&
-0.42420d0,-0.42470d0,-0.42530d0,-0.42590d0,-0.42640d0,-0.42690d0,-0.42730d0,-0.42780d0,-0.42820d0,-0.42860d0,-0.42890d0,&
-0.42930d0,-0.42960d0,-0.43000d0,-0.43030d0,-0.43070d0,-0.43100d0,-0.43140d0,-0.43180d0,-0.43230d0,-0.43270d0,-0.43320d0,&
-0.43370d0,-0.43430d0,-0.43490d0,-0.43550d0,-0.43620d0,-0.43690d0,-0.43770d0,-0.43850d0,-0.43940d0,-0.44030d0,-0.44130d0,&
-0.44230d0,-0.44340d0,-0.44440d0,-0.44550d0,-0.44670d0,-0.44790d0,-0.44920d0,-0.45050d0,-0.45200d0,-0.45350d0,-0.45520d0,&
-0.45700d0,-0.45900d0,-0.46120d0,-0.46380d0,-0.46690d0,-0.47050d0,-0.47510d0,-0.48150d0,-0.49340d0,-0.50360d0,-0.49810d0,&
-0.49100d0,-0.48400d0,-0.47730d0,-0.47090d0,-0.46500d0,-0.45940d0,-0.45400d0,-0.44890d0,-0.44400d0,-0.43920d0,-0.43460d0,&
-0.43020d0,-0.42590d0,-0.42170d0,-0.41760d0,-0.41360d0,-0.40980d0,-0.40600d0,-0.40230d0,-0.39870d0,-0.39520d0,-0.39170d0,&
-0.38840d0,-0.38510d0,-0.38190d0,-0.37870d0,-0.37560d0,-0.37260d0,-0.36970d0,-0.36670d0,-0.36390d0,-0.36110d0,-0.35840d0,&
-0.35570d0,-0.35300d0,-0.35050d0,-0.34790d0,-0.34540d0,-0.34300d0,-0.34050d0,-0.33820d0,-0.33590d0,-0.33360d0,-0.33130d0,&
-0.32910d0,-0.32690d0,-0.32480d0,-0.32270d0,-0.32060d0,-0.31860d0,-0.31660d0,-0.31460d0,-0.31260d0,-0.31070d0,-0.30880d0,&
-0.30700d0,-0.30510d0,-0.41860d0,-0.41860d0,-0.41870d0,-0.41880d0,-0.41890d0,-0.41910d0,-0.41920d0,-0.41950d0,-0.41970d0,&
-0.42000d0,-0.42040d0,-0.42070d0,-0.42110d0,-0.42150d0,-0.42200d0,-0.42250d0,-0.42300d0,-0.42360d0,-0.42410d0,-0.42480d0,&
-0.42540d0,-0.42610d0,-0.42680d0,-0.42750d0,-0.42820d0,-0.42900d0,-0.42970d0,-0.43050d0,-0.43130d0,-0.43210d0,-0.43290d0,&
-0.43370d0,-0.43450d0,-0.43540d0,-0.43620d0,-0.43700d0,-0.43780d0,-0.43860d0,-0.43940d0,-0.44020d0,-0.44100d0,-0.44180d0,&
-0.44250d0,-0.44320d0,-0.44400d0,-0.44470d0,-0.44530d0,-0.44600d0,-0.44660d0,-0.44730d0,-0.44790d0,-0.44850d0,-0.44910d0,&
-0.44970d0,-0.45030d0,-0.45090d0,-0.45160d0,-0.45220d0,-0.45290d0,-0.45360d0,-0.45430d0,-0.45500d0,-0.45580d0,-0.45660d0,&
-0.45740d0,-0.45830d0,-0.45920d0,-0.46010d0,-0.46110d0,-0.46210d0,-0.46320d0,-0.46430d0,-0.46540d0,-0.46650d0,-0.46770d0,&
-0.46890d0,-0.47010d0,-0.47130d0,-0.47260d0,-0.47400d0,-0.47540d0,-0.47690d0,-0.47840d0,-0.48010d0,-0.48180d0,-0.48360d0,&
-0.48550d0,-0.48770d0,-0.49010d0,-0.49270d0,-0.49560d0,-0.49900d0,-0.50300d0,-0.50780d0,-0.51420d0,-0.52380d0,-0.54550d0,&
-0.54160d0,-0.53660d0,-0.53010d0,-0.52330d0,-0.51660d0,-0.51000d0,-0.50370d0,-0.49760d0,-0.49190d0,-0.48630d0,-0.48100d0,&
-0.47590d0,-0.47100d0,-0.46620d0,-0.46160d0,-0.45720d0,-0.45280d0,-0.44860d0,-0.44450d0,-0.44060d0,-0.43670d0,-0.43290d0,&
-0.42920d0,-0.42560d0,-0.42200d0,-0.41860d0,-0.41520d0,-0.41190d0,-0.40860d0,-0.40550d0,-0.40240d0,-0.39930d0,-0.39630d0,&
-0.39340d0,-0.39050d0,-0.38770d0,-0.38490d0,-0.38220d0,-0.37950d0,-0.37690d0,-0.37430d0,-0.37180d0,-0.36930d0,-0.36680d0,&
-0.36440d0,-0.36200d0,-0.35970d0,-0.35740d0,-0.35520d0,-0.35290d0,-0.35070d0,-0.34860d0,-0.34650d0,-0.34440d0,-0.43170d0,&
-0.43170d0,-0.43180d0,-0.43190d0,-0.43200d0,-0.43220d0,-0.43240d0,-0.43270d0,-0.43290d0,-0.43330d0,-0.43360d0,-0.43400d0,&
-0.43440d0,-0.43490d0,-0.43540d0,-0.43590d0,-0.43650d0,-0.43710d0,-0.43770d0,-0.43840d0,-0.43900d0,-0.43970d0,-0.44050d0,&
-0.44120d0,-0.44200d0,-0.44280d0,-0.44360d0,-0.44450d0,-0.44530d0,-0.44620d0,-0.44700d0,-0.44790d0,-0.44880d0,-0.44970d0,&
-0.45060d0,-0.45150d0,-0.45240d0,-0.45330d0,-0.45420d0,-0.45500d0,-0.45590d0,-0.45680d0,-0.45770d0,-0.45850d0,-0.45930d0,&
-0.46020d0,-0.46100d0,-0.46180d0,-0.46260d0,-0.46330d0,-0.46410d0,-0.46490d0,-0.46560d0,-0.46640d0,-0.46720d0,-0.46790d0,&
-0.46870d0,-0.46950d0,-0.47040d0,-0.47120d0,-0.47210d0,-0.47290d0,-0.47380d0,-0.47470d0,-0.47570d0,-0.47670d0,-0.47770d0,&
-0.47870d0,-0.47980d0,-0.48090d0,-0.48210d0,-0.48320d0,-0.48440d0,-0.48560d0,-0.48680d0,-0.48810d0,-0.48940d0,-0.49060d0,&
-0.49200d0,-0.49340d0,-0.49480d0,-0.49620d0,-0.49770d0,-0.49930d0,-0.50090d0,-0.50260d0,-0.50440d0,-0.50640d0,-0.50840d0,&
-0.51060d0,-0.51290d0,-0.51540d0,-0.51820d0,-0.52120d0,-0.52470d0,-0.52870d0,-0.53330d0,-0.53910d0,-0.54690d0,-0.55950d0,&
-0.57810d0,-0.57410d0,-0.57020d0,-0.56420d0,-0.55790d0,-0.55130d0,-0.54480d0,-0.53830d0,-0.53200d0,-0.52590d0,-0.52000d0,&
-0.51430d0,-0.50880d0,-0.50350d0,-0.49840d0,-0.49350d0,-0.48870d0,-0.48400d0,-0.47950d0,-0.47510d0,-0.47090d0,-0.46680d0,&
-0.46270d0,-0.45880d0,-0.45490d0,-0.45120d0,-0.44750d0,-0.44390d0,-0.44040d0,-0.43700d0,-0.43370d0,-0.43040d0,-0.42720d0,&
-0.42400d0,-0.42090d0,-0.41790d0,-0.41490d0,-0.41200d0,-0.40910d0,-0.40630d0,-0.40350d0,-0.40080d0,-0.39810d0,-0.39550d0,&
-0.39290d0,-0.39040d0,-0.38790d0,-0.38540d0,-0.38300d0,-0.38060d0,-0.37820d0,-0.44480d0,-0.44480d0,-0.44490d0,-0.44500d0,&
-0.44510d0,-0.44530d0,-0.44550d0,-0.44570d0,-0.44600d0,-0.44640d0,-0.44670d0,-0.44710d0,-0.44760d0,-0.44800d0,-0.44860d0,&
-0.44910d0,-0.44970d0,-0.45030d0,-0.45090d0,-0.45160d0,-0.45220d0,-0.45300d0,-0.45370d0,-0.45450d0,-0.45530d0,-0.45610d0,&
-0.45690d0,-0.45780d0,-0.45860d0,-0.45950d0,-0.46040d0,-0.46130d0,-0.46220d0,-0.46310d0,-0.46410d0,-0.46500d0,-0.46590d0,&
-0.46690d0,-0.46780d0,-0.46870d0,-0.46970d0,-0.47060d0,-0.47150d0,-0.47240d0,-0.47330d0,-0.47420d0,-0.47510d0,-0.47600d0,&
-0.47690d0,-0.47780d0,-0.47860d0,-0.47950d0,-0.48030d0,-0.48120d0,-0.48210d0,-0.48300d0,-0.48390d0,-0.48480d0,-0.48570d0,&
-0.48660d0,-0.48760d0,-0.48850d0,-0.48950d0,-0.49050d0,-0.49150d0,-0.49260d0,-0.49370d0,-0.49480d0,-0.49590d0,-0.49710d0,&
-0.49830d0,-0.49950d0,-0.50080d0,-0.50200d0,-0.50330d0,-0.50460d0,-0.50590d0,-0.50720d0,-0.50860d0,-0.51000d0,-0.51140d0,&
-0.51290d0,-0.51440d0,-0.51590d0,-0.51750d0,-0.51910d0,-0.52090d0,-0.52270d0,-0.52460d0,-0.52660d0,-0.52860d0,-0.53080d0,&
-0.53320d0,-0.53570d0,-0.53850d0,-0.54150d0,-0.54470d0,-0.54850d0,-0.55270d0,-0.55760d0,-0.56360d0,-0.57120d0,-0.58280d0,&
-0.60870d0,-0.60460d0,-0.60110d0,-0.59620d0,-0.59040d0,-0.58430d0,-0.57800d0,-0.57160d0,-0.56520d0,-0.55890d0,-0.55270d0,&
-0.54670d0,-0.54080d0,-0.53520d0,-0.52960d0,-0.52430d0,-0.51910d0,-0.51410d0,-0.50930d0,-0.50450d0,-0.49990d0,-0.49550d0,&
-0.49110d0,-0.48690d0,-0.48270d0,-0.47870d0,-0.47480d0,-0.47090d0,-0.46720d0,-0.46350d0,-0.45990d0,-0.45640d0,-0.45300d0,&
-0.44960d0,-0.44630d0,-0.44310d0,-0.43990d0,-0.43680d0,-0.43370d0,-0.43070d0,-0.42780d0,-0.42490d0,-0.42210d0,-0.41930d0,&
-0.41650d0,-0.41380d0,-0.41120d0,-0.40860d0/)
contains
  !--------------------------------------------------------------------
  subroutine WIEN2kData
  use param,only : compound,W2kData,RhoPot,Elem,id_at,idW,pi,facmy,rydb
  use maths,only : poisson,interpol_egl3,simpson_s
  use space,only : rx,dx,nx,nxx,nxz,nieq,rmt,ia
  implicit none
  character(len=80) :: txt,fclmsum,fvcoul
  integer :: i,j,k,la,lc,lv,iw,h1(nieq)
  real(8) :: y00,beta,xs,myx,zz(nieq),r1(nieq),a1(nieq)
  real(8),allocatable :: rho4pir2(:,:),a2(:,:)
  real(8),parameter :: thrd=1.d0/3.d0
  Y00=1.d0/sqrt(4.d0*pi)

  la=len_trim(W2kData)
  fclmsum=W2kData(1:la)//'.clmsum_ASA'; lc=len_trim(fclmsum)
  write(61,'(a)') fclmsum(1:lc)
  fvcoul=W2kData(1:la)//'.vcoul_ASA'; lv=len_trim(fvcoul)
  write(61,'(a)') fvcoul(1:lv)
  !CHECK no. atoms.
  open(50,file=fclmsum(1:lc))
  iw=0
  do
    read(50,'(a)',end=1) txt
    if(txt(4:7)=='ATOM') iw=iw+1
  enddo
1 close(50)
  if(iw/=nieq)then
    write(61,'(a,2i4)')'nieq(clmsum_ASA), nieq(inputX) =',iw,nieq
    write(61,'(a)')'when different, check structure; stop'; stop
  endif
  open(50,file=fclmsum(1:lc))
  !starting with WIEN2k labels.
  do i=1,nieq
    do
      read(50,'(a)') txt
      if(txt(4:7)=='ATOM')then
        read(50,'(21x,f4.0,i6,2es20.12)') zz(i),nx(i),r1(i),dx(i)
        exit
      endif
    enddo
  enddo
  close(50)
  nxx=maxval(nx)
  allocate(rx(nxx,nieq),VC(nxx,nieq),VCit(nxx,nieq),Vtot(nxx,nieq))
  allocate(Vxc(nxx,nieq),QvsR(nxx,nieq),myxc(nxx,nieq),rs(nxx,nieq))
  allocate(rho4pir2(nxx,nieq),a2(nxx,nieq),Vtot1(nxx,nieq,ne))
  !
  !RADIAL GRID.
  do i=1,nieq
    do j=1,nx(i)
      rx(j,i)=r1(i)*exp((j-1)*dx(i))
    enddo
  enddo
  !ELECTRON DENSITY indata, spherical part of ASA.
  902 format(2es14.6)
  930 format(3x,4e19.12)
  open(50,file=fclmsum(1:lc))
  do i=1,nieq
    do
      read(50,'(a)') txt
      if(txt(4:7)=='ATOM') exit
    enddo
    do k=1,5; read(50,*); enddo
    !read rho4pir2.
    do k=1,nx(i)/4
      read(50,930) rho4pir2(k*4-3:k*4,i)
    enddo
    read(50,930) rho4pir2(nx(i)/4*4+1:nx(i),i)
  enddo
  close(50)

  !COULOMB POTENTIAL indata, spherical part of ASA.
  open(50,file=fvcoul(1:lv))
  do i=1,nieq
    do
      read(50,'(a)') txt
      if(txt(4:7)=='ATOM') exit
    enddo
    do k=1,5; read(50,*); enddo
    !read vcoul.
    do k=1,nx(i)/4
      read(50,930) VC(k*4-3:k*4,i)
    enddo
    read(50,930) VC(nx(i)/4*4+1:nx(i),i)
    !copy sph. harmonic cvcoul00*r**2 to potential VC*1.
    VC(1:nx(i),i)=VC(1:nx(i),i)*y00/rx(1:nx(i),i)**2
  enddo
  close(50)
  !
  !turning from WIEN2k labels to LEED ones.
  a1=zz;  zz(:)=a1(idW(:))
      do i=1,nieq
      if(zz(i)/=z(i))then
        write(61,'(a,i2,2f5.0)')'i,zz,z=',i,zz(i),z(i),', stop'; stop
      endif
      enddo
  h1=nx;        nx(:)        =h1(idW(:))
  a2=rx;        rx(:,:)      =a2(:,idW(:))
  a2=rho4pir2;  rho4pir2(:,:)=a2(:,idW(:))
  a2=VC;        VC(:,:)      =a2(:,idW(:))
  !
  !APPLICATIONS: nxz, QvsR, rmt, cls.
  do i=1,nieq
    !nx input is defined in WIEN2k/SRC_lapw0/asa_pot_rho.f;
    !nxz is used later in the present program.
    nxz(i)=nx(i)
    !QvsR.
    do j=1,nx(i)
      QvsR(j,i)=simpson_s(rx(1,i),dx(i),nx(i),rho4pir2(1,i),1,j)
    enddo
    !rmt.
    rmt(i)=rx(nx(i),i)
    !cls.
    if(cls(i)/=0.d0)then
      VC(1:nx(i),i)=VC(1:nx(i),i)+cls(i)/rydb !cls(eV), VC(Ry)
    endif
    !write-outs.
    if(RhoPot(1:1)=='y')then
      j=len_trim(compound)
      open(10,file=compound(1:j)//'.RHO4pir2.'//id_at(i),status='replace')
      write(10,902) (rx(j,i),rho4pir2(j,i),j=1,nxz(i))
      open(10,file=compound(1:j)//'.QvsR.'//id_at(i),status='replace')
      write(10,902) (rx(j,i),QvsR(j,i),j=1,nx(i))
      close(10)
    endif
    if(RhoPot(2:2)=='y')then
      j=len_trim(compound)
      open(10,file=compound(1:j)//'.VC.'//id_at(i),status='replace')
      do k=1,nx(i)
        write(10,902) rx(k,i),rx(k,i)*VC(k,i)
      enddo
    endif
    close(10)
  enddo !i
  !logfile.
  write(61,*)
  write(61,906)'atom      ',(i,Elem(int(z(i))),i=1,nieq)
  write(61,908)'Znucleus  ',z(1:nieq)
  write(61,908)'Zelectrons',(VC(1,i)*rx(1,i)*0.5d0,i=1,nieq) !VC in Ry.
  write(61,910)'rZ (B)    ',(rx(nx(i),i),i=1,nieq)
  write(61,906)'NN        ',(ia(2,i),Elem(int(z(ia(2,i)))),i=1,nieq)
  906 format(a/(10(i5,'_',a:)))
  908 format(a/(10f8.2:))
  910 format(a/(10f8.4:))
  !
  !MUFFIN-TIN myxc and rs after HEDIN-LUNDQVIST.
  !rs = radius of sphere accommodating one electron from
  !(4pi/3)*rs**3 = 1/charge density.
  !Ref.: L.Hedin and B.I.Lundqvist, J.Phys.C4,2064(1971).
  do i=1,nieq
    do j=1,nx(i)
      rs(j,i)=(3.d0/(rho4pir2(j,i)/rx(j,i)**2))**thrd
      myx=facmy/rs(j,i)
      xs=rs(j,i)/21.d0
      beta=1.d0+0.7734d0*xs*log(1.d0+1.d0/xs)
      myxc(j,i)=-beta*myx
    enddo
  enddo
  return
  end subroutine WIEN2kData
!-----------------------------------------------------------------------
subroutine FreeAtomData
use param,only : Elem
use maths,only : simpson_s
use space,only : rx,dx,nx,nxx,nieq
implicit none
 character :: atm*2
integer :: i,ir,j
real(8) :: zel,rmin,rmax,dum
!CHARGE DENSITY data, first run.
do ir=1,nieq
  atm=Elem(int(z(ir)))
  j=len_trim(atm)
  open(50,file='chgden'//atm(1:j))
  read(50,*)
  !exponential grid from nx,rmin,rmax.
  read(50,*) z(ir),nx(ir),rmin,rmax
  close(50)
enddo
nxx=maxval(nx)
allocate(rx(nxx,nieq),atrho4pir2(nxx,nieq),QvsR(nxx,nieq))
allocate(VC(nxx,nieq),VCit(nxx,nieq),Vxc(nxx,nieq),Vtot(nxx,nieq))
allocate(myxc(nxx,nieq),rs(nxx,nieq))
allocate(Vtot1(nxx,nieq,ne))
!CHARGE DENSITY data, second run.
write(61,*)
do ir=1,nieq
  atm=Elem(int(z(ir)))
  j=len_trim(atm)
  open(50,file='chgden'//atm(1:j))
  read(50,*)
  read(50,*) z(ir),nx(ir),rmin,rmax
  dx(ir)=log(rmax/rmin)/(nx(ir)-1)
  do i=1,nx(ir)
    rx(i,ir)=rmin*exp(dx(ir)*(i-1))
    read(50,*) dum,atrho4pir2(i,ir)
  enddo
  zel=simpson_s(rx(1,ir),dx(ir),nx(ir),atrho4pir2(1,ir),1,nx(ir))
  write(61,'(a,i2,a,f7.4,a,f7.4)')&
    'FreeAtomData: atom ',ir,', zat=',z(ir),', zel=',zel
  if(abs(z(ir)-zel)>1.d-02)then
    write(61,'(a)')'stop'
    stop
  endif
  call rho_center(ir)
  close(50)
enddo
return
end subroutine FreeAtomData
  !-------------------------------------------------------------------
  !Atom centre correction
  !rho(1:16) calculated by E.L.Shirley's atomic program
  !contains minor discontiuities. A spline calculation on rho(17:30) is
  !used for smoothing rho(1:16) by extrapolation.
  !The smoothed rho can be studied on plots of the files a4pirho,
  !a4pirhop, and a4pirhob, see subroutine superpose.
  subroutine rho_center(ir)
  use maths,only : spcoef,svalue
  use space,only : rx
  implicit none
  integer :: i,j,ir,i1,i2,n,last_interval
  real(8) :: x(14),f(14),fp(14),fb(14),ft(14),fi,xi,dum
  i1=17; i2=30; n=i2-i1+1
  j=0
  do i=i1,i2
    j=j+1
    x(j)=i
    f(j)=atrho4pir2(i,ir)/(rx(i,ir)*rx(i,ir))
  enddo
  call spcoef(n,x,f,1.d+33,1.d+33,fp,fb,ft,'rho_center')
  last_interval=1
  do i=1,i1-1
    xi=i
    call svalue(n,x,f,fp,fb,ft,xi,fi,dum,last_interval,'rho_center')
    atrho4pir2(i,ir)=fi*rx(i,ir)*rx(i,ir)
  enddo
  return
  end subroutine rho_center
!----------------------------------------------------------------------
!SUPERPOSITION calculates muffin-tin potentials using Mattheiss's
!superposition method.
!(1) Lattice: A 3D structure is generated from a given unit cell.
!    A surface slab of n layers is modeled by a unit cell whose depth
!    is equal to n-1 interlayer separations plus a "large" distance.
!    Hence, the atoms situated near the surface of the slab get
!    negligible overlap with the atoms next to the bulk cut-off.
!(2) Superposition: The charge density of an atom A in the unit cell
!    and the charge density contributions from all neighbours of A are
!    added. The same procedure for the Hartree potential.
!(3) Spherically symmetric approximation: A's superposition density and
!    superposition potential are considered as expanded in multipoles
!    around A. The zeroth multipole is kept.
!(4) The outer radius of an atom in the material is determined by the
!    requirement that the atom is neutral. This is the so-called Norman
!    atomic radius.

!Ref.: L.F. Mattheiss, Phys. Rev. 133, A1399 (1964).
!Ref.: T.L. Loucks, Augmented plane wave method (Benjamin, New York,
!      1967).
!Ref.: J.G. Norman, Mol. Phys. 31, 1191 (1975).
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package,
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/Activity.
subroutine superpose
use param,only : compound,RhoPot,Elem,pi,id_at,facmy !,rydb,bohr,radfac
use maths,only : poisson
use space,only : nieq,nx,nxx,rx,dx,nxz,ad,ia
implicit none
integer :: ir,i,imax(nieq),j
real(8) :: atrho4pi(nxx,nieq),rho4pi(nxx,nieq),&
  rho4pir2(nxx,nieq),pi2,fpi,ufpi,myx,xs,beta,a,rMTmax
real(8),parameter :: thrd=1.d0/3.d0
 pi2=pi*pi; fpi=4.d0*pi; ufpi=1.d0/fpi
 !
 !Superposition.
 !Crystal atom:  radial grid of length imax(ir).
 !Neighbor atom: grid of full length nx(ir).
 !# neighbors:   nshell in subroutine NeighborShells.
 !
 !rMTmax < NN distance ad(2,ir), see Refs. in NeighborSum.
 do ir=1,nieq
   rMTmax=ad(2,ir)*0.9d0
   imax(ir)=nint((log(rMTmax/rx(1,ir)))/dx(ir)+1.d0)
   imax(ir)=min(imax(ir),nx(ir))
 enddo
 do ir=1,nieq
   do i=1,nx(ir)
     atrho4pi(i,ir)=atrho4pir2(i,ir)/rx(i,ir)**2
   enddo
 enddo
 call NeighborSum(imax,atrho4pi,rho4pi)
 !
 !Charge and potential versus radius.
 do ir=1,nieq
   do i=1,imax(ir)
     rho4pir2(i,ir)=rho4pi(i,ir)*rx(i,ir)*rx(i,ir)
   enddo
   call poisson(rx(1,ir),dx(ir),imax(ir),z(ir),rho4pir2(1,ir),&  !superposition
                QvsR(1,ir),VC(1,ir))
   if(QvsR(imax(ir),ir)<z(ir))then
     write(61,900)'superposition: ir,imax,QvsR(imax),z=',&
       ir,imax(ir),QvsR(imax(ir),ir),z(ir)
     nxz(ir)=imax(ir)
   else
     do i=1,imax(ir)
       if(QvsR(i,ir)>=z(ir))then
         nxz(ir)=i
         exit
       endif
     enddo
   endif
   900 format(a,i2,i6,2f7.2)
   !The integration constant of Poisson's equation is unknown in
   !Mattheiss's potential model. What can we do?
   !
   !Suggestion:
   !Setting the integration constant of VC(:,:) to maxVC(:,:)=0.
   !
   !This is a HEURISTIC renormalization based on the assumption that
   !all potentials 1:nieq terminate at the same level near the
   !interstitial potential.
   a=maxval(VC(1:nxz(ir),ir))
   VC(1:nxz(ir),ir)=VC(1:nxz(ir),ir)-a
   !
   !Core level shifts:
   !the researcher's adjustments to the atomic potentials (superposition).
   if(cls(ir)/=0.d0) VC(1:nxz(ir),ir)=VC(1:nxz(ir),ir)+cls(ir)
 enddo
 !
 !POTENTIAL AFTER VANHOVE (not updated).
 !call NeighborShells
 !call Poisson(atrho4pir2,atpot,nxz)  !superposition
 !call NeighborSum(imax,atpot,VC)
 !
 if(RhoPot(1:1)=='y')then
   j=len_trim(compound)
   do ir=1,nieq
     open(10,file=compound(1:j)//'.RHO4pir2.'//id_at(ir),status='replace')
     write(10,902) (rx(i,ir),rho4pir2(i,ir),i=1,nxz(ir))
     open(10,file=compound(1:j)//'.QvsR.'//id_at(ir),status='replace')
     write(10,902) (rx(i,ir),QvsR(i,ir),i=1,nxz(ir))
   enddo
   close(10)
 endif
 if(RhoPot(2:2)=='y')then
   j=len_trim(compound)
   do ir=1,nieq
     open(10,file=compound(1:j)//'.VC.'//id_at(ir),status='replace')
     do i=1,nxz(ir)
       if(VC(i,ir)>-10.d0) write(10,902) rx(i,ir),VC(i,ir)
     enddo
   enddo
   close(10)
 endif
 902 format(2es14.6)
  write(61,*)
  write(61,'(a)')'superposition'
  write(61,906)'atom      ',(i,Elem(int(z(i))),i=1,nieq)
  write(61,908)'Znucleus  ',z(1:nieq)
  write(61,908)'Zelectrons',(VC(1,i)*rx(1,i)*0.5d0,i=1,nieq)
  write(61,910)'rZ (B)    ',(rx(nxz(i),i),i=1,nieq)
  write(61,906)'NN        ',(ia(2,i),Elem(int(z(ia(2,i)))),i=1,nieq)
  906 format(a/(10(i5,'_',a:)))
  908 format(a/(10f8.2:))
  910 format(a/(10f8.4:))
 !
 !MUFFIN-TIN myxc and rs after HEDIN-LUNDQVIST.
 !rs = radius of sphere accommodating one electron, rs satisfying
 !(4pi/3)*rs**3 = inverse charge density.
 !Ref.: L.Hedin and B.I.Lundqvist, J.Phys.C4,2064(1971).
 do ir=1,nieq
   do i=1,nxz(ir)
     rs(i,ir)=(3.d0/rho4pi(i,ir))**thrd
     myx=facmy/rs(i,ir)
     xs=rs(i,ir)/21.d0
     beta=1.d0+0.7734d0*xs*log(1.d0+1.d0/xs)
     myxc(i,ir)=-beta*myx
   enddo
 enddo
 return
end subroutine superpose
!----------------------------------------------------------------------
!NEIGHBOR_SUM superposes charge-density (or potential) contributions
!from neighbor atoms. The present code updates the routine 'sumax' of
!the following program packages,
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package,
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/activity.
!
!     ncon  = number of shells included,
!     ia(j) = atomic type in j'th shell,
!     na(j) = number of atoms in j'th shell,
!     ad(j) = distance to j'th shell,
!     imax  = summation up to rx(imax).
subroutine NeighborSum(imax,atrho4pi,rho4pi)
!Author: J.Rundgren.
!Ref.: T.L. Loucks, Augmented plane wave method (Benjamin, New York,
!      1967), eqs. 3.22,3.26,3.28.
!Ref.: A. Messiah, Quantum mechanics (Dover,1999), Vol. 1, Appendix B,
!      formula (99), Green function of Laplace operator, Legendre
!      polynomial expansion, first term.
use maths,only : integral_egt
use space,only : nieq,nx,nxx,rx,dx,ad,na,ia,ncon
implicit none
integer :: ja,ic,i,ir,imax(nieq),jmax
real(8) :: rho4pi(nxx,nieq),atrho4pi(nxx,nieq),x1,x2,a
!Loewdin's spherical harmonics expansion term for L,M=0,0 (Loucks).
do ir=1,nieq
  jmax=imax(ir)
  rho4pi(1:jmax,ir)=atrho4pi(1:jmax,ir) !Catom.
  do ja=2,ncon(ir)
    ic=ia(ja,ir)
    do i=1,jmax !i in Catom
      x1=min(abs(rx(i,ir)-ad(ja,ir)),rx(nx(ic),ic)) !x1,nx in Natom
      x2=min(    rx(i,ir)+ad(ja,ir), rx(nx(ic),ic)) !x2,nx in Natom
      if(x1<x2)then
        a=integral_egt&
          (x1,x2,2,nx(ic),rx(1,ic),dx(ic),atrho4pi(1,ic)) !Natom
        rho4pi(i,ir)=rho4pi(i,ir)+&                       !Catom.
          a*na(ja,ir)*0.5d0/(ad(ja,ir)*rx(i,ir))          !Catom.
      endif
    enddo
  enddo
enddo
return
end subroutine NeighborSum
  !--------------------------------------------------------------------
  real(8) function V0interstitial(Etmp,xtmp)
  !Ref.: L.Hedin and B.I.Lundqvist, J.Phys.C 4, 2064 (1971), eq. (3.6).
  use param,only : fpi3,facmy,prnt !,rydb,TAinp
  use maths,only : interpol_egl3
  use space,only : nieq,neq,volXC,nxz
  implicit none
  real(8),intent(in) :: Etmp,xtmp(nieq)
  integer :: i
  real(8),parameter :: thrd=1.d0/3.d0
  real(8) :: myxc,vMT,myx,xs,beta,pmem,zUC,eMT
  !
  !rsis interstitial.
  zUC=0.d0; eMT=0.d0; vMT=0.d0
  do i=1,nieq
    zUC=zUC+neq(i)*z(i)
    eMT=eMT+neq(i)*interpol_egl3&
              (xtmp(i),nxz(i),rx(1,i),dx(i),QvsR(1,i),'V0interstitial')
    vMT=vMT+neq(i)*fpi3*xtmp(i)**3
  enddo
  rsis=(zUC-eMT)/(volXC-vMT)
  rsis=(fpi3*rsis)**(-thrd)
  if(prnt)then
    write(61,898)'rmt',xtmp(:)
    write(61,900)'zUC,eMT    =',zUC,eMT
    write(61,900)'volXC,volMT=',volXC,vMT
    write(61,900)'      rsis =',rsis
  endif
  898 format(a,4es18.10)
  900 format(a,2es18.10)
  !
  !V0 interstitial from rs and myxc value.
  myx=facmy/rsis
  xs=rsis/21.d0
  beta=1.d0+0.7734d0*xs*log(1.d0+1.d0/xs)
  myxc=-beta*myx
  pmem=0.d0
  V0interstitial=xc_myxc(Etmp,rsis,myxc,pmem)*myxc
  if(prnt) write(61,900)'V0inter    =',V0interstitial
  return
  end function V0interstitial
  !---------------------------------------------------------------------
  subroutine MTxc(Exc)
  !LocalDensityApproximation: Vxc is generated from E-VC.
  !Ref.: Phys.Rev.B 68, 125405 (2003).
  use space,only : nxx,nxz,nieq
  implicit none
  integer :: ir,i
  real(8) :: f(nxx),pmem,Exc
  f=0.d0
  do ir=1,nieq
    pmem=0.d0
    do i=1,nxz(ir)
      f(i)=xc_myxc(Exc-VCit(i,ir),rs(i,ir),myxc(i,ir),pmem)
    enddo
    call xc_smooth(f,nxz(ir))
    do i=1,nxz(ir)
      Vxc(i,ir)=f(i)*myxc(i,ir)
    enddo
  enddo
  return
  end subroutine MTxc
  !-------------------------------------------------------------------
  real(8) function xc_myxc(Etmp,rs,myxc,pmem)
  !XCPOT is excited-state exchange-correlation potential.
  !Ref.: L.Hedin and B.I.Lundqvist, J.Phys.C 4,2064 (1971).
  !Ref.: R.E.Watson, J.F.Herbst, L.Hodges, B.I.Lundqvist, and
  !      J.W.Wilkins, Phys.Rev.B13, 1463 (1976), Appendix B.
  !Ref.: J.Neve, J.Rundgren, and P.Westrin, J.Phys.C 15, 4391 (1982).
  !pF=fermi momentum,
  !p =electron momentum in units of pF, used in self_energy table,
  !se(p,rs)=self energy in units of pF**2,
  !myxc    =ground-state exchange-correlation potential,
  !myxc*(se(p,rs)/se(1,rs))=excited-state exch.-corr. potential.
  !Hedin-Lundqvist equation:
  !  (p*pF)**2+myxc*(se(p,rs)/se(1,rs)=Etmp+pF**2+myxc
  !below rewritten in the form
  !  f=(p-1)*(p+1)+[myxc*(se(p,rs)/se(1,rs)-myxc-Etmp]/pF**2
  implicit none
  integer :: j
  real(8)  :: se1,Etmp,rs,myxc,pF,p1,p2,pFm2,pmem
  real(8),parameter :: pi9_4thrd=1.919158292677504d0
  pF=pi9_4thrd/rs
  pFm2=1.d0/(pF*pF)
  se1=self_energy(1.d0,rs)
  !rough approximation.
  if(pmem==0.d0)then
    p1=1.d0+Etmp*pFm2
  else
    p1=1.d0+(Etmp+(1.d0-self_energy(pmem,rs)/se1)*myxc)*pFm2
  endif
  if(p1<0.d0)then
    p1=0.d0
    goto 999
  endif
  p1=sqrt(p1)
  !successive approximation.
  do j=1,10
    p2=1.d0+(Etmp+(1.d0-self_energy(p1,rs)/se1)*myxc)*pFm2
    if(p2<0.d0)then
      p1=0.d0
      goto 999
    endif
    p2=sqrt(p2)
    if(abs(p1-p2)<1.d-08) exit
    p1=p2
  enddo
  999 continue
  pmem=p1
  xc_myxc=self_energy(p1,rs)/se1
  return
  end function xc_myxc
  !---------------------------------------------------------------------
  real(8) function self_energy(ps,rs)
  !VXCDATA extracts data from the Vxc tables supplied by Hedin-Lundqvist
  !and Sernelius.
  implicit none
  real(8) :: ps,rs
  if(rs>sr(nsr))then
    write(61,*)'VXCDATA warning: rs>6. is put rs=6.'
    rs=sr(nsr)
  endif
  if(ps<=sp(nsp))then
    self_energy=interpol_tl(ps,rs)
  else
    self_energy=extrapol_tl(ps,rs)
  endif
  return
  end function self_energy
  !---------------------------------------------------------------------
  real(8) function extrapol_tl(ps,rs)
  !EXTRAPOL_TL=inverse-wavenumber extrapolation from given Vxc tables.
  implicit none
  real(8),parameter :: ps1=2.80d0,ps2=3.00d0
  real(8) :: ps,rs,v1,v2,g,b,a
  v1=interpol_tl(ps1,rs)
  v2=interpol_tl(ps2,rs)      !value at table margin.
  g=abs(v2*(ps2-ps1)/(v2-v1)) !gradient at table margin.
  b=ps2*(g-ps2)
  a=v2*sqrt(ps2*g)
  extrapol_tl=a/sqrt(ps*ps+b)
  return
  end function extrapol_tl
  !---------------------------------------------------------------------
  real(8) function interpol_tl(ps,rs)
  !INTERPOL_TL=triangulated linear interpolation in given Vxc tables.
  !Ref.: L.F.Shampine, R.C.Allen, and S.Pruess, Fundamentals of
  !      numerical computing (Wiley,1997).
  implicit none
  integer :: k,k1,k2,l1,l2,l
  real(8)  :: ps,rs,a1,a2,a3,p21,r21,s,ta,a,b
  !rectangle containing ps,rs.
  l1=1; l2=nsr
  do
    l=(l1+l2)/2; if(sr(l)>rs)then; l2=l; else; l1=l; endif
    if(l2-l1==1)exit
  enddo
  k1=1; k2=nsp
  do
    k=(k1+k2)/2; if(sp(k)>ps)then; k2=k; else; k1=k; endif
    if(k2-k1==1)exit
  enddo
  !triangles: rectangle with diagonal.
  p21=sp(k2)-sp(k1); r21=sr(l2)-sr(l1); s=r21/p21; ta=r21*p21
  if(rs-sr(l1)>=s*(ps-sp(k1)))then
   !upper triangle.
    a1=phi(ps,rs,sp(k2),sr(l2),sp(k1),sr(l2),ta)
    a2=phi(ps,rs,sp(k1),sr(l2),sp(k1),sr(l1),ta)
    a3=phi(ps,rs,sp(k1),sr(l1),sp(k2),sr(l2),ta)
    a=a1*sdat(k1,l1)*sr(l1)+a2*sdat(k2,l2)*sr(l2)+a3*sdat(k1,l2)*sr(l2)
  else
   !lower triangle.
    a1=phi(ps,rs,sp(k2),sr(l1),sp(k2),sr(l2),ta)
    a2=phi(ps,rs,sp(k2),sr(l2),sp(k1),sr(l1),ta)
    a3=phi(ps,rs,sp(k1),sr(l1),sp(k2),sr(l1),ta)
    a=a1*sdat(k1,l1)*sr(l1)+a2*sdat(k2,l1)*sr(l1)+a3*sdat(k2,l2)*sr(l2)
  endif
  !triangles: rectangle with antidiagonal.
  s=-s
  if(rs-sr(l2)>=s*(ps-sp(k1)))then
   !upper triangle.
    a1=phi(ps,rs,sp(k2),sr(l1),sp(k2),sr(l2),ta)
    a2=phi(ps,rs,sp(k2),sr(l2),sp(k1),sr(l2),ta)
    a3=phi(ps,rs,sp(k1),sr(l2),sp(k2),sr(l1),ta)
    b=a1*sdat(k1,l2)*sr(l2)+a2*sdat(k2,l1)*sr(l1)+a3*sdat(k2,l2)*sr(l2)
  else
   !lower triangle.
    a1=phi(ps,rs,sp(k1),sr(l1),sp(k2),sr(l1),ta)
    a2=phi(ps,rs,sp(k2),sr(l1),sp(k1),sr(l2),ta)
    a3=phi(ps,rs,sp(k1),sr(l2),sp(k1),sr(l1),ta)
    b=a1*sdat(k1,l2)*sr(l2)+a2*sdat(k1,l1)*sr(l1)+a3*sdat(k2,l1)*sr(l1)
  endif
  interpol_tl=0.5d0*(a+b)
  return
  end function interpol_tl
  !---------------------------------------------------------------------
  real(8) function phi(x,y,xa,ya,xb,yb,ta)
  real(8) :: x,y,xa,ya,xb,yb,ta
  phi=(xa*yb-xb*ya+(ya-yb)*x+(xb-xa)*y)/ta
  end function phi
  !---------------------------------------------------------------------
  subroutine xc_smooth(f,n)
 !binomial smoothing of second order, with coefficients 1,2,1;
 !when executed many times, it converges to the normal distribution.
  implicit none
  integer :: i,j,n
  real(8) :: f(n),g(n)
  real(8),parameter :: thrd=1.d0/3.d0
  do j=1,16
   !f goes to g.
    g(1)=(f(2)+f(1)+f(1))*thrd
    do i=2,n-1
      g(i)=(f(i-1)+f(i+1)+f(i)+f(i))*0.25d0
    enddo
    g(n)=(f(n-1)+f(n)+f(n))*thrd
   !g goes back to f.
    f(1)=(g(2)+g(1)+g(1))*thrd
    do i=2,n-1
      f(i)=(g(i-1)+g(i+1)+g(i)+g(i))*0.25d0
    enddo
    f(n)=(g(n-1)+g(n)+g(n))*thrd
  enddo
  return
  end subroutine xc_smooth
  !--------------------------------------------------------------------
  subroutine errV0(X,fitness)
  !ERR_V0 is mean square error of a four-parameter fit to the energy-
  !dependent inner potential.
  implicit none
  integer :: i
  real(8) :: x(4),fitness,errv
  errv=0.d0
  do i=1,ne
    errv=errv+( v0ev(i)-max( x(4), x(1)+x(2)/sqrt(eev(i)+x(3))) )**2
  enddo
  fitness=sqrt(errv/ne)
  return
  end subroutine errV0
  !--------------------------------------------------------------------
  subroutine errV01(X,fitness)
  !ERR_V0 is mean square error of a four-parameter fit to the energy-
  !dependent inner potential.
  implicit none
  integer :: i,m
  real(8) :: X(3),fitness,errv,x4
  x4=minval(v0ev)
  m=minloc(v0ev,1)
  errv=0.d0
  do i=m+1,ne
    errv=errv+( v0ev(i)-max(X(1)+X(2)/sqrt(eev(i)+X(3)),x4) )**2
  enddo
  fitness=sqrt(errv/(ne-m))
  return
  end subroutine errV01
  !--------------------------------------------------------------------
  subroutine errV02(X,fitness)
  !ERR_V0 is mean square error of a four-parameter fit to the energy-
  !dependent inner potential.
  implicit none
  integer :: i,m
  real(8) :: X(3),fitness,errv,x5,tmp
  m=minloc(v0ev,1)
  x5=v0ev(1)
  errv=0.d0
  do i=1,m
    tmp=(eev(i)/eev(m))**2
    tmp=v0ev(i)-(x5+tmp*(X(1)+tmp*(X(2)+tmp*X(3))))
    errv=errv+tmp**2
  enddo
  fitness=sqrt(errv/m)
  return
  end subroutine errV02
 !--------------------------------------------------------------------
  !radial Dirac and Schroedinger waves.
  !Ref.: M.E.Rose, Relativistic Electron Theory (Wiley, 1961).
  subroutine waveq(x,y,dydx)
  use param,only : cm2
  use space,only : rx,dx,nxz
  implicit none
  integer :: i
  real(8)  :: y(2),dydx(2),x,p,vmv0x,qbx,emv
  real(8),parameter :: a6=1.d0/6.d0
  !lagrangian 4-point interpolation.
  p=log(x/rx(1,iatom))/dx(iatom)+1.d0
  i=max(int(p),2); i=min(i,nxz(iatom)-2)
  p=p-i
  vmv0x=   -a6*p*(p-1.d0)*(p-2.d0)*Vtot(i-1,iatom)&
        +0.5d0*(p*p-1.d0)*(p-2.d0)*Vtot(i  ,iatom)&
        -0.5d0*p*(p+1.d0)*(p-2.d0)*Vtot(i+1,iatom)&
                  +a6*p*(p*p-1.d0)*Vtot(i+2,iatom)
  emv=emv0-vmv0x/x
  qbx=kappa/x
  dydx(1)=-y(1)*qbx+y(2)*(1.d0+cm2*emv)
  dydx(2)= y(2)*qbx-y(1)*emv
  return
  end subroutine waveq
end module enrgy

!==============================================================================
module calMT
implicit none
contains
  !--------------------------------------------------------------------
  subroutine MTradii(ie,Exc)
  use param,only : SelMTrad
  use maths,only : interpol_egl3
  use space,only : rx,dx,nxz,rmt,ad,ia,nieq
  use enrgy,only : MTxc,V0interstitial,Vtot,VC,V0,Vxc
  implicit none
  integer :: ie,n,i,j,k,ir(nieq)
  integer,parameter :: n2=129
  real(8) :: Exc,Vrmt,x1,x2,xn,d,dudx,dvdx,a,r(nieq,nieq),&
    x(n2),u(n2),v(n2),xk
  !
  if(SelMTrad=='M')then
    !OPTION: NN midpoint radii.
    !NB, a midpoint radius is generally not unique: in the NN distance
    !table a particular atom can be connected with two or more NN atoms.
    !A mean radius is attributed to each particular atom.
    if(ie==1)then
      ir=0
      do i=1,nieq
        a=ad(2,i)*0.5
        !minimum fix when 1:nxz does not overlap midpoint.
        ir(i)=ir(i)+1; r(i,ir(i))=min(a,rx(nxz(i),i))
        j=ia(2,i)
        ir(j)=ir(j)+1; r(j,ir(j))=min(a,rx(nxz(j),j))
      enddo
      write(61,'(3x,a)')'i,rmt(i) | mean[i-NNi distances]'
      do i=1,nieq
        rmt(i)=sum(r(i,1:ir(i)))/ir(i)
        write(61,921) i,rmt(i),'  |',(2.d0*r(i,j),j=1,ir(i))
      enddo
      write(61,922)'gNN (B) ',(ad(2,i)-rmt(i)-rmt(ia(2,i)),i=1,nieq)
      write(61,'(a)')"gNN<0 when rMT's overlap"
      921 format(i3,f8.4,a,32f8.4)
      922 format(a,8f9.4:/(8x,8f9.4:))
    endif
    !
    !Total potential.
    call MTxc(Exc)
    V0=V0interstitial(Exc,rmt)
    do i=1,nieq
      Vtot(1:nxz(i),i)=VC(1:nxz(i),i)+Vxc(1:nxz(i),i)
      Vrmt=&
        interpol_egl3(rmt(i),nxz(i),rx(1,i),dx(i),Vtot(1,i),'MTr1')
      Vtot(1:nxz(i),i)=Vtot(1:nxz(i),i)-(Vrmt-V0)
    enddo
  elseif(SelMTrad=='I')then
    !OPTION: NN potential intersection radii.
    !NB, a midpoint radius is generally not unique: in the NN distance
    !table a particular atom can be connected with two or more NN atoms.
    !A mean radius is attributed to each particular atom.
    call MTxc(Exc)
    ir=0
    do i=1,nieq
      j=ia(2,i)
      if(ad(2,i)-rx(nxz(i),i)-rx(nxz(j),j)>=0.d0)then
        !nonoverlapping potentials.
        ir(i)=ir(i)+1; r(i,ir(i))=rx(nxz(i),i)
        ir(j)=ir(j)+1; r(j,ir(j))=rx(nxz(j),j)
      else
        !overlapping potentials.
        !linear grid in interval [x1,x2].
        x2=rx(nxz(i),i)
        x1=ad(2,i)-rx(nxz(j),j)
        d=(x2-x1)/(n2-1)
        Vtot(1:nxz(i),i)=VC(1:nxz(i),i)+Vxc(1:nxz(i),i)
        Vtot(1:nxz(j),j)=VC(1:nxz(j),j)+Vxc(1:nxz(j),j)
        do n=1,n2
          x(n)=x1+d*(n-1)
          u(n)=interpol_egl3(x(n),nxz(i),rx(1,i),dx(i),Vtot(1,i),'MTr2')
          xn=ad(2,i)-x(n)
          v(n)=interpol_egl3(xn  ,nxz(j),rx(1,j),dx(j),Vtot(1,j),'MTr3')
        enddo
        k=minloc(abs(u-v),1)
        dudx=(u(k+1)-u(k))/(x(k+1)-x(k))
        dvdx=(v(k+1)-v(k))/(x(k+1)-x(k))
        xk=x(k)+(v(k)-u(k))/(dudx-dvdx)
        !
        ir(i)=ir(i)+1; r(i,ir(i))=xk
        ir(j)=ir(j)+1; r(j,ir(j))=ad(2,i)-xk
      endif
    enddo
    write(61,'(3x,a)')'i,rmt(i) | mean[i-NNi distances]'
    do i=1,nieq
      rmt(i)=sum(r(i,1:ir(i)))/ir(i)
      write(61,921) i,rmt(i),'  |',(2.d0*r(i,j),j=1,ir(i))
    enddo
    write(61,922)'gNN (B) ',(ad(2,i)-rmt(i)-rmt(ia(2,i)),i=1,nieq)
    write(61,'(a)')"gNN<0 when rMT's overlap"
    !
    !Total potential.
    V0=V0interstitial(Exc,rmt)
    do i=1,nieq
      Vtot(1:nxz(i),i)=VC(1:nxz(i),i)+Vxc(1:nxz(i),i)
      Vrmt=&
        interpol_egl3(rmt(i),nxz(i),rx(1,i),dx(i),Vtot(1,i),'MTr6')
      Vtot(1:nxz(i),i)=Vtot(1:nxz(i),i)-(Vrmt-V0)
    enddo
  endif
  return
  end subroutine MTradii
  !----------------------------------------------------------------------
  ! MT_OPTIMISATION calculates muffin-tin radii and interstitial potential.
  !
  ! Atomic MT floor is the potential in the shell between the MT radius
  !     radius and the charge radius.
  ! Interstitial potential is an average of atomic floors, each
  !     inequivalent atom weighted with its occupation in the structure.
  ! The MT radii are optimised in such a way that
  !     (MT potential at MT radius) - (interstitial potential) = constant.
  ! Potential steps are eliminated by a common shift of MT potentials.
  subroutine MTopt(Einc,Exc)
  use param,only : SelVMadl,Elem,rydb,prnt,compound,RhoPot,id_at !,SelVxc,bohr
  use maths,only : DE,Dim_XC,NP,NPfac,nfeval,itval,bestval,interpol_egl3
  use maths,only : dfdx_egl3
  use space,only : nxz,nxx,nieq,rmt,rmt0,rx,dx,ad,ia,rmtn,rmtx,rmt1
  use enrgy,only : ie,ne,MTxc,V0interstitial,VC,VCit,Vtot,QvsR,E_DE,Vxc,z
  use enrgy,only : wrMadl,dVeV,rerr,verr,mVstp,Vrmt,V0,rsis,V01,Vtot1
  implicit none
  integer :: i,j,jc,iter,rl
  real(8) :: VM(nieq),Einc,Exc,a,gNN(nieq),fitness
  real(8) :: XCmin(nieq),XCmax(nieq),bestmem_XC(nieq),qmt(nieq)
  !
  prnt=.false.
  if(ie==1.and.SelVMadl=='y')then
    call ProcRoySoc_ElectroStaticMatrix
  endif
  !
  !monatomic crystal.
  iter=0
  if(nieq==1)then
    call MTxc(Exc)
    Vtot(1:nxz(1),1)=VC(1:nxz(1),1)+Vxc(1:nxz(1),1)
    rmt(1)=ad(2,1)*0.5d0
    write(61,922)'rMT (B) ',rmt(1)
    V0=V0interstitial(Exc,rmt(1))
    write(61,924) V0*rydb,rsis
    a=interpol_egl3(rmt(1),nxz(1),rx(1,1),dx(1),Vtot(1,1),'MTo1')
    Vtot(1:nxz(1),1)=Vtot(1:nxz(1),1)-(a-V0)
    return
  endif
  924 format(/'V0  (eV)  ',f6.2,', rs(B)=',f0.2)

  !polyatomic crystal.
      if(ie==1)then
        VCit=VC+0.5d0
        Vxc=0.d0
        Vtot1=0.d0
      endif

DO ITER=1,6
  write(61,'(i2,78("-"))') ITER; write(*,'(i4)') ITER

      !LocalDensityApproximation: Vxc generated from Exc-VC, see sub MTxc.
      call MTxc(Exc)
      do i=1,nieq
        Vtot(1:nxz(i),i)=VCit(1:nxz(i),i)+Vxc(1:nxz(i),i)
      enddo

  !Max MT radius for DE optimizer
  do i=1,nieq; rmt0(i)=rx(nxz(i),i); enddo
  write(61,922)'rDE (B) ',rmt0(1:nieq)
  !
  if(SelVMadl=='y')then
    write(61,'(a)')'DO DIFFERENTIAL EVOLUTION ON rMT & V0 & Vxc & VM.'
    wrMadl=.false.
  else
    write(61,'(a)')'DO DIFFERENTIAL EVOLUTION ON rMT & V0 & Vxc.'
  endif
  Dim_XC=nieq
  XCmin(1:Dim_XC)=0.5d0
  XCmax(1:Dim_XC)=rmt0(1:nieq)
  do i=1,nieq
    if(rmtn(i)>1.d-04) XCmin(i)=rmtn(i)
    if(rmtx(i)>1.d-04) XCmax(i)=rmtx(i)
  enddo
  if(ie==1.and.iter==1)then
    write(61,922)'rDE min ',XCmin(1:nieq)
    write(61,922)'rDE max ',XCmax(1:nieq)
  endif
  NP=NPfac*Dim_XC
  E_DE=Einc
  call DE(errMT,XCmin,XCmax,bestmem_XC)
  write(61,910)'calls,iteff=',nfeval,itval,'fitness=',bestval
  910 format(a,i9,i6,4x,a,es10.3)
  if(bestval>=1.d+33)then
    write(61,'(a)')'DE failed, stop. Check the NN distances table.'
    stop
  endif
  !
  !errMT is recalculated from rmt determined by DE, because the
  !last fitness value is generally not the lowest value.
  prnt=.false.

      !DO best rmt after DE optimization.
      !errMT calculates  rerr(B),verr(eV),Vrmt,Vstp,mVstp,V0,rsis
      rmt(1:nieq)=bestmem_XC(1:Dim_XC)
      call errMT(bestmem_XC,fitness)
      !ENDDO best rmt.

  prnt=.false.
  !Madelung print-outs.
  if(SelVMadl=='y')then
    write(61,'(2x,a)')'DO MADELUNG POTENTIAL.'
    wrMadl=.true.
    call ProcRoySoc_ElectroStaticPotential(bestmem_XC,VM) !making print-outs.
    write(61,'(2x,a)')'ENDDO MADELUNG POTENTIAL.'
  endif
  write(61,'(a)')'ENDDO DIFFERENTIAL EVOLUTION.'
  !
  write(61,'(a)')'values at r=rMT:'
  write(61,921)'atom   ',(i,Elem(int(z(i))),i=1,nieq)
  write(61,922)'rMT (B)',rmt(1:nieq)
  do i=1,nieq
    gNN(i)=ad(2,i)-rmt(i)-rmt(ia(2,i))
    qmt(i)=interpol_egl3(rmt(i),nxz(i),rx(1,i),dx(i),QvsR(1,i),'MToq')
  enddo
  write(61,922)'gNN (B) ',(gNN(i),i=1,nieq)
  write(61,909)'rms[gNN] (B)=',rerr,', rs (B)=',rsis
  write(61,906)'qMT (#e) ',qmt(1:nieq)
  write(61,906)'dV/drMT (eV/B)',&
    (dfdx_egl3(rmt(i),nxz(i),rx(1,i),dx(i),Vtot(1,i),'MTod')*rydb,i=1,nieq)
  write(61,906)'Vtot@rmt(:)-V0 (eV)',(Vrmt(1:nieq)-V0)*rydb
  write(61,909)'in eV: V0=',V0*rydb,', mean[Vtot@rmt(:)-V0]=',mVstp*rydb,&
    ', rms[Vtot@rmt(:)-V0]=',verr
  write(61,909) '|mean|<=',dVeV
  921 format(a/(10(i5,'_',a:)))
  922 format(a/(10f8.4:))
  906 format(a/(10f8.2:))
  909 format(5(a,f0.3))

      !VCit shifted by mVstp.
      do i=1,nieq
        VCit(1:nxz(i),i)=VCit(1:nxz(i),i)-mVstp
      enddo
      if(abs(mVstp)*rydb<dVeV)then
        V01(ie)=V0
        rmt1(1:nieq,ie)=rmt(1:nieq)
        exit
      endif
  write(61,'(a)')'Next, do another iteration.'
ENDDO !ITER
  write(61,'(/a)')'Finish with no-standing-waves approximation:'
  write(61,'( a)')'shift of MT potential to continuity with V0.'
      do i=1,nieq
        Vtot(1:nxz(i),i)=Vtot(1:nxz(i),i)-Vrmt(i)
        Vtot1(1:nxz(i),i,ie)=Vtot(1:nxz(i),i)
      enddo

  !Vtot and Vxc print-out.
  if(RhoPot(2:2)=='y'.and.ie==1)then
    jc=len_trim(compound)
    do i=1,nieq
      open(10,file=compound(1:jc)//'.Vtot.'//id_at(i),status='replace')
      open(11,file=compound(1:jc)//'.Vxc.'//id_at(i),status='replace')
      do j=1,nxz(i)
        if(rx(j,i)<=rmt(i).and.Vtot(j,i)>-20d0) &
        write(10,'(2es14.6)') rx(j,i),Vtot(j,i)
        write(11,'(2es14.6)') rx(j,i),Vxc(j,i)
      enddo
      close(10); close(11)
    enddo
  endif
  return

  !unformatted store.
  if(ie==ne)then
    write(61,'(80("-")/a)')'Unformatted store: rmt, V0, Vtot'
    j=len_trim(compound)
    inquire(iolength=rl) rmt1 !B.
    write(61,'(a,i0,1x,i0)')'rmt: iolength*4= & nieq*nev*8= ',&
                                  rl*4,nieq*ne*8
    open(62,file=compound(1:j)//'.tm_rmt',form='unformatted',recl=rl*4)
    write(62) rmt1
    close(62)
    !
    inquire(iolength=rl) V01 !Ry.
    write(61,'(a,i0,1x,i0)')'V0: iolength*4= & nev*8= ',&
                                 rl*4,ne*8
    open(62,file=compound(1:j)//'.tm_V0',form='unformatted',recl=rl*4)
    write(62) V01
    close(62)
    !
    inquire(iolength=rl) Vtot1 !Ry.
    write(61,'(a,i0,1x,i0)')'Vtot: iolength*4= & n0*nieq*nev*8= ',&
                                   rl*4,nxx*nieq*ne*8
    open(62,file=compound(1:j)//'.tm_Vtot',form='unformatted',recl=rl*4)
    write(62) Vtot1
    close(62)
  endif
  return
  end subroutine MTopt
  !-------------------------------------------------------------------
  subroutine errMT(X,fitness)
  !errMT measures potential steps at MT radius between interstitial
  !potential and MT potentials.
  use param,only : SelVxc,SelVMadl,rydb
  use maths,only : Dim_XC,interpol_egl3
  use space,only : neq,nieq,nlatp,rx,dx,nxz,ad,ia
  use enrgy,only : E_DE,Vtot,VCit,Vxc,V0,Vstp,mVstp,Vrmt,rerr,verr
  use enrgy,only : rMTweight,V0interstitial
  implicit none
  integer :: ir
  real(8) :: harvest,fitness
  real(8) :: X(1:Dim_XC),VM(nieq),dmt2,gNN(nieq),Exc
  real(8),parameter :: Excn=60.d0/13.60569172d0
  intrinsic random_number
  rerr=0.d0; verr=0.d0
  !Overlapping MT spheres discarded.
  do ir=1,nieq
    gNN(ir)=ad(2,ir)-X(ir)-X(ia(2,ir))
    if(gNN(ir)<0.d0)then
      call random_number(harvest)
      fitness=(1.d0+harvest)*1.d+33
      return
    endif
  enddo
  !
  !interstitial potential.
  if(SelVxc=='E')then; Exc=E_DE; elseif(SelVxc=='C')then; Exc=Excn; endif
  V0=V0interstitial(Exc,X)
  !
  !Madelung potentials.
  if(SelVMadl=='y')then
    call ProcRoySoc_ElectroStaticPotential(X,VM) !***VM in Ry.***
    do ir=1,nieq
    Vtot(1:nxz(ir),ir)=VCit(1:nxz(ir),ir)+Vxc(1:nxz(ir),ir)+VM(ir)
    enddo
  endif
  !
  !Vrmt(:) and Vstp(:) contain core level shifts with respect to the
  !VC(:) given in FreeAtomData and WIEN2kData.
  dmt2=0.d0
  do ir=1,nieq
    dmt2=dmt2+neq(ir)*gNN(ir)**2
    Vrmt(ir)=&
      interpol_egl3(X(ir),nxz(ir),rx(1,ir),dx(ir),Vtot(1,ir),'errMT')
    Vstp(ir)=Vrmt(ir)-V0
  enddo
  !rerr=rms[interMT separation] in Bohr.
  rerr=sqrt(dmt2/nlatp)
  !
  !standard deviation: verr=sd[potential steps at rMT] (in Rydbergs).
  !sd -> 0 implies Vstp(i) -> a=constant for i=1:nieq.
  mVstp=sum(neq(1:nieq)*Vstp(1:nieq))/nlatp
  verr=sqrt(sum(neq(1:nieq)*(Vstp(1:nieq)-mVstp)**2)/nlatp)
  !verr in eV.
  verr=verr*rydb
  !
  fitness=sqrt(rerr**2*rMTweight+verr**2*(1.d0-rMTweight))
  return
  end subroutine errMT
  !---------------------------------------------------------------------
  !ELECTROSTATIC POTENTIAL (so-called MADELUNG potential)
  !calculated from the total set of MT charges in the crystal.
  !
  !For LEED we make an ad-hoc approximation to semi-infinite space.
  !
  subroutine ProcRoySoc_ElectroStaticPotential(X,VM)
  use maths,only : interpol_egl3
  use space,only : nieq,neq,nlatp,nxz,rx,dx
  use enrgy,only : QvsR,z,ewii,dpii,wrMadl
  implicit none
  integer :: i
  real(8) :: X(nieq),VM(nieq),zmt(nieq),qewii(nieq),qdpii(nieq),&
    a,zmtm(nieq)

  !Ref.: Proc.Roy.Soc. A373,27(1980), eqs.(3.23) and (3.24).
  !MT charges zmt.
  do i=1,nieq
    a=interpol_egl3(X(i),nxz(i),rx(1,i),dx(i),QvsR(1,i),'Madelung')
    zmt(i)=a-z(i)
  enddo
  !MT charges relative to mean zmt.
  a=sum(neq(1:nieq)*zmt(1:nieq))/nlatp  !nlatp=sum(nieq)
  zmtm=zmt-a
  do i=1,nieq
    qewii(i)=sum(zmtm(:)*ewii(i,:)) !(H)
    qdpii(i)=-sum(zmtm(:)*dpii(i,:)) !(H)
    VM(i)=qewii(i)+qdpii(i)
  enddo
  !VM shifted to zero mean value.
  !a factor =2.0 for VM in Rydberg atomic units;
  !a factor =0.5 for VM in semi-infinite space (ad-hoc approx.).
  a=sum(neq*VM)/nlatp  !nlatp=sum(neq)
  VM=VM-a
  if(wrMadl)then
    write(61,'(a)')'zMT'
    write(61,902) zmt
    write(61,'(a)')'zMT-mean'
    write(61,902) zmtm
    write(61,'(a)')'Ewald potential sum_j[zmt*ewii](H):'
    write(61,902) qewii
    write(61,'(a)')'Dipole-moment potential -sum_j[zmt*dpii](H):'
    write(61,902) qdpii
    write(61,'(a)')'VM-mean(H)'
    write(61,902) VM
  endif
  902 format(10f8.4:)
  return
  end subroutine ProcRoySoc_ElectroStaticPotential
  !---------------------------------------------------------------------
  !ELECTROSTATIC MATRIX (so-called MADELUNG matrix)
  !for the total set of MT's in the crystal.
  !Output:
  !  Ewald matrix ewii(nieq:nieq),
  !  dipole-moment matrix dpii(nieq:nieq).
  !Ref.: S.W.DELEEUW et al., Proc.Roy.Soc.London A373,27(1980).
  !Ref.: E.R.Smith, Proc.Roy.Soc.London A375,475(1981).
  !Ref.: N.Karasawa & W.A.Goddard, J.Phys.Chem. 93,7320(1989).
  !Ref.: J.Rundgren, this work (2011).
  !
  !***routine in Hartree atomic units***.
  subroutine ProcRoySoc_ElectroStaticMatrix
  use param,only : pi,fpi3,compound,SelVCoul
  use maths,only : interpol_egl3
  use space,only : nieq,neq,ieq,nlatp,rcX,rk,volXC
  use enrgy,only : ewii,dpii,z
  implicit none
!  real(8),parameter :: vacgap=30.d0 !in Bohr
  real(8),parameter :: thrd=1.d0/3.d0,ewerr=1.d-013
  integer,parameter :: np2=20,nbeta=4
  !
  character :: fl*80
  integer :: i,i1,i2,i3,j,n,np,npc,nij,n1,n2,k,nc(np2),l,m
  real(8) :: rc(3,3),gc(3,3),p(3),dpm(3),tbsqpi,p2,b2,erfc,&
    facg,tpi,a,radAT,r,sc,rr,rg,amin,co,dfac,JMS
  real(8) :: zmt(nieq),ql(nlatp),d(3,nlatp,nlatp),&
    ac(np2),beta(nbeta),wr(nlatp,nlatp),wg(nlatp,nlatp),&
    ewll(nlatp,nlatp),ewli(nlatp,nieq),ewiik(nieq,nieq,nbeta),&
    dpll(nlatp,nlatp),dpli(nlatp,nieq),ve(nieq),vd(nieq)
  logical :: writ,ex
  allocate(ewii(nieq,nieq),dpii(nieq,nieq))
  rc=rcX
  !
  write(61,'(a)')'ProcRoySoc_ElectroStaticMatrix:'
  !Reciprocal lattice vectors in units of 1/volXC.
  a=1.d0/volXC
  gc(1,1)=(rc(2,2)*rc(3,3)-rc(3,2)*rc(2,3))*a
  gc(2,1)=(rc(3,2)*rc(1,3)-rc(1,2)*rc(3,3))*a
  gc(3,1)=(rc(1,2)*rc(2,3)-rc(2,2)*rc(1,3))*a
  !
  gc(1,2)=(rc(2,3)*rc(3,1)-rc(3,3)*rc(2,1))*a
  gc(2,2)=(rc(3,3)*rc(1,1)-rc(1,3)*rc(3,1))*a
  gc(3,2)=(rc(1,3)*rc(2,1)-rc(2,3)*rc(1,1))*a
  !
  gc(1,3)=(rc(2,1)*rc(3,2)-rc(3,1)*rc(2,2))*a
  gc(2,3)=(rc(3,1)*rc(1,2)-rc(1,1)*rc(3,2))*a
  gc(3,3)=(rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2))*a
  write(61,'(4x,a)')'Real vectors for slab with vacuum (B)'
  do j=1,3
    write(61,'(4x,3f12.5)') (rc(i,j),i=1,3)
  enddo
  write(61,'(4x,a)')'Reciprocal vectors (1/B)'
  do j=1,3
    write(61,'(4x,3f12.5)') (gc(i,j),i=1,3)
  enddo
  write(61,'(4x,a)')'Real and recip. vectors i,j=1:3 orthonormal?'
  do i=1,3
    do j=1,3
      write(61,'(i6,i2,f9.6)')&
        i,j,rc(1,i)*gc(1,j)+rc(2,i)*gc(2,j)+rc(3,i)*gc(3,j)
    enddo
    if(i==1.or.i==2) write(61,*)
  enddo
  !MT charges zmt fetched from previous S, A, or W run.
  j=len_trim(compound)
  fl='~/phs/'//compound(1:j)//'/'//SelVCoul//'/qMT'
  j=len_trim(fl)
  inquire(file=fl(1:j),exist=ex)
  if(ex)then
    open(50,file=fl(1:j))
    do i=1,nieq
      read(50,*) a
      zmt(i)=a-z(i)
    enddo
    close(50)
  else
    write(61,'(a)')'stop, file '//fl(1:j)//' not available.'; stop
  endif
  write(61,'(a)')'MT charge derived from electron scattering at VC+Vxc,'
  write(61,901)'zMT(#e) ',zmt(1:nieq)
  !Charge neutrality condition of Ewald summation:
  a=sum(neq(1:nieq)*zmt(1:nieq))/nlatp  !nlatp=sum(nieq)
  zmt(1:nieq)=zmt(1:nieq)-a
  write(61,'(a)')'MT charge with respect to global mean MT charge,'
  write(61,901)'zMT-mean',zmt(1:nieq)
  !Subscripts in unit-charge Madelung matrix.
  n1=0
  do i=1,nieq
    ieq(i)=n1+1
    n2=n1+neq(i)
    n1=n2
  enddo
  900 format(a)
  910 format(10i8:)
  write(61,900)'Blocks ii of symmetrically equivalent atoms i:'
  write(61,900)'      ii       i'
  do i=1,nieq
    write(61,910) (j,j=ieq(i),ieq(i)+neq(i)-1)
  enddo
  do i=1,nieq
    ql(ieq(i):ieq(i)+neq(i)-1)=zmt(i)
  enddo
  write(61,900)'MT charge with respect to global mean MT charge, check,'
  i=1
    write(61,901)'qMT(#e) ',(ql(j),j=ieq(i),ieq(i)+neq(i)-1)
  do i=2,nieq
    write(61,911)           (ql(j),j=ieq(i),ieq(i)+neq(i)-1)
  enddo
  901 format(a/(10f8.4:))
  911 format(10f8.4:)
  !
  !LATTICE MEASURES.
  radAT=(volXC/nlatp/fpi3)**thrd
  rr=volXC**thrd
  rg=1.d0/rr
  dfac=pi/(3.d0*volXC)
  do l=1,nlatp
    do m=1,nlatp
      d(:,l,m)=rk(:,l)-rk(:,m)
      dpll(l,m)=dfac*(d(1,l,m)**2+d(2,l,m)**2+d(3,l,m)**2)
    enddo
  enddo
  !
  !do REAL LATTICE [Proc.Roy.Soc.A 375,475(1981)].
  !***Hartree atomic units.***
  write(61,'(a)')"Checking Ewald sum's invariance with respect to beta."
  nij=nlatp*(nlatp+1)/2
  do k=1,nbeta
  beta(k)=1.d0/(radAT*2**(k-1))
  write(61,'(a,i0,a,f0.4)')'  k=',k,', beta=',beta(k)
  tbsqpi=2.d0*beta(k)/sqrt(pi)
  nc=0; ac=0.d0
  do j=1,nlatp
    do i=j,nlatp  !triangular matrix, no factor 1/2.
      do np=2,np2
        npc=np
        wr(j,i)=0.d0
        amin=1.d+33
        do i1=-np,np
        do i2=-np,np
        do i3=-np,np
          if(i==j.and.i1==0.and.i2==0.and.i3==0)cycle
          p(:)=i1*rc(:,1)+i2*rc(:,2)+i3*rc(:,3)+d(:,j,i)
          r=sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
          if(r>np*rr)cycle !spherical summation.
          a=erfc(beta(k)*r)/r; if(amin>a) amin=a
          wr(j,i)=wr(j,i)+a
          !for error estimate.
          nc(np)=nc(np)+1
          ac(np)=ac(np)+abs(ql(j)*ql(i))
        enddo !i3
        enddo !i2
        enddo !i1
        if(amin<ewerr)exit
      enddo !np
      if(i/=j) wr(i,j)=wr(j,i) !square matrix.
    enddo !i
  enddo !j
  do j=1,nlatp
    wr(j,j)=wr(j,j)-tbsqpi
  enddo
  !enddo REAL LATTICE.
  a=nij
  write(61,903)'np,n,np*rr',npc,nint(nc(npc)/a),npc*rr
  903 format(4x,a,i3,i10,2es10.2)
  !
  !do RECIPROCAL LATTICE [Proc.Roy.Soc.A 375,475(1981)].
  !***Hartree atomic units.***
  tpi=2.d0*pi
  facg=1.d0/(pi*volXC)
  b2=(pi/beta(k))**2
  nc=0; ac=0.d0
  do j=1,nlatp
    do i=j,nlatp  !triangular matrix, no factor 1/2.
      do np=2,np2
        npc=np
        amin=1.d+33
        wg(j,i)=0.d0
        do i1=-np,np
        do i2=-np,np
        do i3=-np,np
          if(i1==0.and.i2==0.and.i3==0)cycle
          p(:)=i1*gc(:,1)+i2*gc(:,2)+i3*gc(:,3)
          p2=p(1)*p(1)+p(2)*p(2)+p(3)*p(3)
          r=sqrt(p2)
          if(r>np*rg)cycle !spherical summation.
          a=exp(-b2*p2)/p2; if(amin>a) amin=a
          sc=p(1)*d(1,j,i)+p(2)*d(2,j,i)+p(3)*d(3,j,i)
          co=cos(tpi*sc)
          wg(j,i)=wg(j,i)+a*co
          !for error estimate.
          nc(npc)=nc(npc)+1
          ac(npc)=ac(npc)+abs(ql(j)*ql(i)*co)
        enddo !i3
        enddo !i2
        enddo !i1
        wg(j,i)=facg*wg(j,i)
        if(amin<ewerr)exit
      enddo !np
      if(i/=j) wg(i,j)=wg(j,i)  !square matrix.
    enddo !i
  enddo !j
  !enddo RECIPROCAL LATTICE.
  a=nij
  write(61,903)'np,n,np*rg',npc,nint(nc(npc)/a),npc*rg
  ewll=wr+wg
  !
  !EWALD MATRIX nlatp:nlatp.
  writ=.true.
  if(k==nbeta-1.and.writ)then
    open(40,file='VM.ewll')
    write(40,900)'ewll(1:nlatp,1:nlatp)'
    do i=1,nlatp
      write(40,940) (ewll(i,j),j=1,nlatp)
    enddo
    close(40)
  endif
  940 format(32f7.3:)
  !
  !EWALD MATRIX nlatp:nieq.
    do l=1,nlatp
      do i=1,nieq
        ewli(l,i)=sum(ewll(l,ieq(i):ieq(i)+neq(i)-1))
      enddo
    enddo
  if(k==nbeta-1.and.writ)then
    open(40,file='VM.ewli')
    write(40,900)'ewli(1:nlatp,1:nieq)'
    do l=1,nlatp
      write(40,940) (ewli(l,i),i=1,nieq)
    enddo
    close(40)
  endif
  !
  !EWALD MATRIX nieq:nieq.
  !Equal probability for electron's hitting any of {I|q_J=q_i} sites,
  !therefore, division by neq(ii).
  do j=1,nieq
    do i=1,nieq
      ewiik(i,j,k)=sum(ewli(ieq(i):ieq(i)+neq(i)-1,j))/neq(i)
    enddo
  enddo
  if(k==nbeta-1.and.writ)then
    open(40,file='VM.ewii')
    write(40,900)'ewii(1:nieq,1:nieq)'
    do i=1,nieq
      write(40,940) (ewiik(i,j,k),j=1,nieq)
    enddo
    close(40)
  endif
  !Ewald electron potential 1:nieq.
  do i=1,nieq
    ve(i)=sum(zmt(1:nieq)*ewiik(i,1:nieq,k))
  enddo
  write(61,'(a)')'Ewald potential ve=sum[zmt*ewii] (H):'
  write(61,911) ve(1:nieq)
  enddo !k
  ewii(:,:)=ewiik(:,:,nbeta-1)
  !
  !DIPOLE MATRIX nlatp:nlatp.
  !Ref.: Proc.Roy.Soc.A373,27(1980), eq. (3.23).
  !Ref.: Proc.Roy.Soc.A375,475(1981), eqs. (1.2) and (1.9).
  !***Hartree atomic units.***
  if(writ)then
    open(40,file='VM.dpll')
    write(40,'(a)')'dpll(1:nlatp,1:nlatp)'
    do l=1,nlatp
      write(40,940) (dpll(l,m),m=1,nlatp)
    enddo
    close(40)
  endif
  !
  !DIPOLE MATRIX nlatp:nieq.
  do l=1,nlatp
    do i=1,nieq
      dpli(l,i)=sum(dpll(l,ieq(i):ieq(i)+neq(i)-1))
    enddo
  enddo
  if(writ)then
    open(40,file='VM.dpli')
    write(40,'(a)')'dpli(1:nlatp,1:nieq)'
    do l=1,nlatp
      write(40,940) (dpli(l,i),i=1,nieq)
    enddo
    close(40)
  endif
  !DIPOLE MATRIX nieq:nieq.
  !Equal probability for electron's hitting any of {I|q_J=q_i} sites,
  !hence, division by neq(ii).
  do j=1,nieq
    do i=1,nieq
      dpii(i,j)=sum(dpli(ieq(i):ieq(i)+neq(i)-1,j))/neq(i)
    enddo
  enddo
  if(writ)then
    open(40,file='VM.dpii')
    write(40,'(a)')'dpii(1:nieq,1:nieq)'
    do i=1,nieq
      write(40,940) (dpii(i,j),j=1,nieq)
    enddo
    close(40)
  endif
  !Dipole electron potential 1:nieq.
  !v components are constant when dpm=0; eq. zero after shift.
  !Ref.: Proc.Roy.Soc.A373,27(1980), eq.(3.23).
  do i=1,nieq
    vd(i)=-sum(zmt(1:nieq)*dpii(i,1:nieq))
  enddo
  write(61,900)'Dipole-moment potential vd=-sum[e*zmt(:)*dpii(i,:)] (H):'
  write(61,911) vd(1:nieq)
  write(61,900)'Madelung potential VM=(ve+vd)-mean (H):'
  ve=ve+vd
  a=sum(neq*ve)/nlatp
  ve=ve-a
  write(61,911) ve(1:nieq)
  !DIPOLE MOMENT PER LATTICE POINT.
  a=0.d0
  do l=1,nlatp
    do m=1,nlatp
      a=a-ql(l)*ql(m)*dpll(l,m)
    enddo
  enddo
  a=a/nlatp**2
  do n=1,3
    dpm(n)=sum(ql(1:nlatp)*rk(n,1:nlatp))/nlatp
  enddo
  !checking J(M,S) [Proc.Roy.Soc.A 375,475(1981), eq. (3.19)].
  JMS=sum(dpm**2)*dfac*2.d0
  if(abs(a-JMS)>1.d-08)then
    write(61,'(a)')'stop, check dipole moment code'
    write(*,*)'a,JMS',a,JMS
    stop
  endif
  write(61,904)'dpm= ',dpm,'   J(M,S)=',JMS
  904 format(a,3f9.4,a,es9.2)
  return
  end subroutine ProcRoySoc_ElectroStaticMatrix
end module calMT

!==============================================================================
module calPS
implicit none
contains
  !----------------------------------------------------------------------
  !PartialWaveMethod
  !calculates electron-atom scattering phase shifts on an energy scale
  !referred to an energy-dependent inner potential.
  !
  ! Vfe = fast-electron level  = zero point of the energy scale,
  ! E = Eprimary-Vfe
  ! V0(E) = Vinterstitial-Vfe  = energy-dependent inner potential,
  ! q=sqrt(E-V0(E)) = electron wave number in the crystal,
  ! psr,ps = phaseshift versus q and orbital quantum number l,
  ! psrm,psrp = phaseshift versus q and relativistic quantum number kappa.
  subroutine PartialWaveMethod
  use param,only : SelVxc,selMTrad,SelOutp,SpinPS,PSvsE,BmSel,Elem,bohr,rydb
  use param,only : cm2,wave_start,eta_cutoff,relerr1,relerr2,fpi3,compound
  use maths,only : sph_bessel_j,sph_bessel_y,dfdx_egl3
  use space,only : nieq,neq,volXC,rx,rmt,nxz,ia,ad
  use enrgy,only : ie,ne,eev,Vtot,V0,z,emv0,q,lmax,sdat,sdat1,nsp,nsr
  use enrgy,only : rerr,verr
  use calMT,only : MTopt,MTradii
  implicit none
  character :: fil*80
  integer :: i,j,ir,l,lx,lxq(ne,nieq)
  real(8) :: V0tab(ne),avrmt(nieq),v0c(9),bj(0:lmax+1),by(0:lmax+1)
  real(8) :: psr(ne,0:lmax,nieq),psrm(ne,0:lmax,nieq),psrp(ne,0:lmax,nieq)
  real(8) :: qr,volmt,bohr3,Einc,Exc,avverr,avrerr
!  real(8),parameter :: thrd=1.d0/3.d0
  real(8),parameter :: Excn=60.d0/13.60569172d0
  bohr3=bohr**3

!ODE parameters.
write(61,'(/a,i3)')'SCATTERING: ne=',ne
write(61,'(a,2es8.1)')'wave_start,ps_cutoff=',wave_start,eta_cutoff
write(61,'(a,2es8.1)')'relerr1,2=',relerr1,relerr2
if(SelOutp=='d')then
  open(20,file='LvsKRstart')
  open(21,file='LMAXvsKRmt',status='replace')
endif
!
!SELF-ENERGY.
allocate(sdat(nsp,nsr))
sdat=reshape(sdat1,shape=(/nsp,nsr/))
!
!loop over energy grid:
lx=0; avrmt=0.d0; avrerr=0.d0; avverr=0.d0
psr=0.d0; psrm=0.d0; psrp=0.d0
do ie=1,ne
  write(*,'(i3,1x,f6.1)') ie,eev(ie)
  write(61,'(20("====")/a,i3,1x,f6.1)') "ie,E=",ie,eev(ie)
  Einc=eev(ie)/rydb
  if(SelVxc=='E')then; Exc=Einc; elseif(SelVxc=='C')then; Exc=Excn; endif
  !
  if(selMTrad=='D')then
    !DE optimized radii.
    call MTopt(Einc,Exc)
  else
    !'DM'=midpoint or 'DI'=intersection radii.
    call MTradii(ie,Exc)
  endif
  avrmt=avrmt+rmt !sum over 1:nieq.
  avrerr=avrerr+rerr
  avverr=avverr+verr
  V0tab(ie)=V0
  !
  !calculation of phase shifts, Vtot->total potential times radius.
  do ir=1,nieq
    Vtot(1:nxz(ir),ir)=Vtot(1:nxz(ir),ir)*rx(1:nxz(ir),ir)
  enddo
  emv0=Einc-V0
  emv0=max(emv0,1.d-02)
  q=sqrt(emv0*(1.d0+emv0*cm2))
  do ir=1,nieq
    qr=q*rmt(ir)
    call sph_bessel_j(lmax+1,qr,bj)
    call sph_bessel_y(lmax+1,qr,by)
    call PScalc('-',ir,lx,bj,by,psrm)
    call PScalc('+',ir,lx,bj,by,psrp)
    lxq(ie,ir)=lx
  enddo
enddo !ie
 close(10); close(20); close(21)

avrmt=avrmt/ne
avrerr=avrerr/ne
avverr=avverr/ne
volmt=fpi3*sum(neq(1:nieq)*avrmt(1:nieq)**3) !corresponding to volXC.
write(61,'(20("====")/"mean rMT over E:")')
write(61,921) (ia(1,i),Elem(int(z(i))),i=1,nieq)
write(61,923) (rx(nxz(i),i),i=1,nieq)
write(61,922) (avrmt(i),i=1,nieq)
write(61,924) (( ad(2,i)-avrmt(i)-avrmt(ia(2,i)) ),i=1,nieq)
write(61,925) (ia(2,i),Elem(int(z(ia(2,i)))),i=1,nieq)
write(61,926) avrerr,avverr
write(61,928) volmt/volXC
write(61,'(20("----"))')
921 format('Atom    '/(10(i5,'_',a:)))
923 format('rZ  (B) '/(10f8.4:))
922 format('rMT (B) '/(10f8.4:))
924 format('gNN (B) '/(10f8.4:))
925 format('NN      '/(10(i5,'_',a:)))
926 format('rms[gNN] (B)=',f0.3,', rms[Vtot@rmt] (eV)=',f0.3)
928 format('MT SpaceFill=',f0.3)
  !deallocate(Vtot)
  !
  !Jumps of pi removed from the phaseshift versus energy curves.
  call PSnojump(psrm,psrp)
  !spin averaged phase shift psr.
  do ir=1,nieq
    psr(1:ne,0,ir)=psrm(1:ne,0,ir)
    do l=1,lmax
      do i=1,ne
        psr(i,l,ir)=((l+1)*psrm(i,l,ir)+l*psrp(i,l,ir))/(l+l+1)
      enddo
    enddo
  enddo
  if(PSvsE=='y')then
    i=len_trim(compound)
    do ir=1,nieq
      if(BmSel(ir)/='y')cycle
      fil=compound(1:i)//'.PScrv.su.'; j=len_trim(fil)
      write(fil,'(2(a,i0))') fil(1:j),nint(z(ir)),'.',ir; j=len_trim(fil)
      open(62,file=fil(1:j),status='replace')
      fil=compound(1:i)//'.PScrv.sd.'; j=len_trim(fil)
      write(fil,'(2(a,i0))') fil(1:j),nint(z(ir)),'.',ir; j=len_trim(fil)
      open(63,file=fil(1:j),status='replace')
      fil=compound(1:i)//'.PScrv.sl.'; j=len_trim(fil)
      write(fil,'(2(a,i0))') fil(1:j),nint(z(ir)),'.',ir; j=len_trim(fil)
      open(64,file=fil(1:j),status='replace')
      do j=1,ne
        emv0=eev(j)-V0tab(j)*rydb
        write(62,'(8es14.6)') emv0,(psrm(j,l,ir),l=0,lmax)
        write(63,'(8es14.6)') emv0,(psrp(j,l,ir),l=0,lmax)
        write(64,'(8es14.6)') emv0,( psr(j,l,ir),l=0,lmax)
      enddo
      close(62); close(63); close(64)
    enddo
  endif
  !
  call InnerPotential(V0tab,v0c)
  if(SelOutp=='p')then
    call PStab(psr,'r',V0tab)
    if(SpinPS=='y')then
      call PStab(psrm,'-',V0tab)
      call PStab(psrp,'+',V0tab)
    endif
  endif
  if(SelOutp=='s')then
    call sigma(lxq,psrm,psrp,V0tab)
  endif
  return
  end subroutine PartialWaveMethod
  !----------------------------------------------------------------------------
  !PScalc calculates phase shifts and radial wave functions.
  subroutine PScalc(tag,ir,lx,bj,by,delta)
  use param,only : SelOutp,UvsR,bohr,cm2,node,eta_cutoff,id_at,relerr0,relerr1
  use maths,only : sph_bessel_j
  use space,only : nieq,rx,rmt
  use enrgy,only : ie,ne,kappa,emv0,q,lmax
  implicit none
  logical :: w_xec
  character :: idl*3,tag*1
  integer,parameter :: lsta=5
  integer :: ir,l,lmin,i,ista,j,lx,nps
  real(8)  :: delta(ne,0:lmax,nieq),bj(0:lmax+1),by(0:lmax+1),y(2),&
    ysta(2),or(node),ou(2,node),t_end,qr,w_max(2),ps,a,bjl,bjlp,&
    byl,bylp,c,s,g,gnrm,gnrmp,gp,bsl(0:lmax+1)

qr=q*rmt(ir)
select case('SelOutp')
 case('d')
  write(61,'("tag.atom=",a,", KRMT=",f6.2)') tag//id_at(ir),qr
  write(61,'("  l    ps        no")')
endselect
!
if(tag=='-')then; lmin=0; else; lmin=1; endif
ista=1
do l=lmin,lmax
  if(tag=='-')then; kappa=-l-1; else; kappa=l; endif
  call StartODE(ir,l,lsta,ista,ysta)
  !
  !LvsKRstart: insensitive to energy, hence recorded for ie=ne.
  if(SelOutp=='d'.and.tag/='+'.and.ie==ne.and.ir==1)then
    write(20,904) q*rx(ista,ir),l
  endif
  !
  !PS integrated.
  t_end=rmt(ir)
  w_xec=.true.; w_max=(/1.d0,1.d0/)
  call WaveODE&
    (ir,l,kappa,ista,ysta,t_end,y,qr,relerr0,w_xec,w_max,bj,by,ps,nps)
  w_xec=.false.
  call WaveODE&
    (ir,l,kappa,ista,ysta,t_end,y,qr,relerr1,w_xec,w_max,bj,by,ps,nps)
  delta(ie,l,ir)=ps
  !
  !print wave.
  if(UvsR=='y')then
    if(ie==ne.and.tag/='+')then
      !normalisation of u1(l,r) to r*j(l,q*r) at muffin-tin radius.
      g=y(1)/rmt(ir)
      gp=(-(kappa+1)*y(1)/rmt(ir)+&
         (1.d0+cm2*emv0)*y(2))/rmt(ir)
      bjl=bj(l); bjlp=(l/qr)*bjl-bj(l+1)
      byl=by(l); bylp=(l/qr)*byl-by(l+1)
      c=cos(ps); s=sin(ps)
      gnrm=c*bjl-s*byl; gnrmp=q*(c*bjlp-s*bylp)
      if(l<lsta.and.abs(gnrmp)>=abs(q*gnrm))then
        a=gnrmp/gp
      else
        a=gnrm/g
      endif
      write(idl,'(i0)') idl
      j=len_trim(idl)
      open(10,file='UvsR'//idl(1:j)//id_at(ir))
      do i=1,nps
        call sph_bessel_j(l,q*or(i),bsl)
        write(10,940) or(i)*bohr,ou(1,i)*a/or(i),bsl(l)
      enddo
      close(10)
    endif
  endif
  940 format(4es14.6)

  !PS termination.
  if(SelOutp=='d') write(61,902) tag,ir,l,ps,nps
  if(l>=6) then
    lx=lmax
    if(abs(ps)<eta_cutoff)then
      lx=l-1
      if(SelOutp=='d'.and.tag/='+'.and.ir==1)then
        write(21,904) qr,lx !LMAXvsKRmt
      endif
      exit
    endif
  endif
enddo
  return
  902   format(a,2i5,es12.4,i5)
  904   format(e15.8,i4)
  end subroutine PScalc
  !--------------------------------------------------------------------
  subroutine WaveODE&
    (ir,l,kappa,ista,ysta,t_end,y,qr,relerr,w_xec,w_max,bj,by,ps,nps)
  use param,only : cm2,node
  use maths,only : ODE
  use space,only : rx,rmt
  use enrgy,only : lmax,iatom,waveq,emv0
  implicit none
  logical :: w_xec
  integer,parameter :: neqn=2,nwork=100+21*neqn
  integer :: ir,l,kappa,ista,i,flag,iwork(5),nps
  real(8)  :: t,y(2),ysta(2),tout,t_end,qr,or(node),ou(2,node),&
    work(nwork),relerr,abserr,ps,a,b,h_start,w_max(2),&
    bj(0:lmax+1),by(0:lmax+1)
  t=rx(ista,ir); tout=t_end
  y=ysta
  h_start=0.d0
  abserr=10.d0*epsilon(1.d0)*max(w_max(1),w_max(2))
  nps=0; flag=1
  iatom=ir
  call ODE(waveq,neqn,y,t,tout,relerr,abserr,flag,nwork,work,iwork,&
           h_start,or,ou,nps,node)
  if(flag>=3)then
    write(61,*)'WaveODE: l,flag=',l,flag; stop
  endif
  if(w_xec)then
    w_max=0.d0
    do i=1,nps
      w_max(1)=max(w_max(1),abs(ou(1,i)))
      w_max(2)=max(w_max(2),abs(ou(2,i)))
    enddo
  endif
  a=-(kappa+l+1)*y(1)+rmt(ir)*(1.d0+cm2*emv0)*y(2)
  b=qr*y(1)
  ps=atan((a*bj(l)+b*bj(l+1))/(a*by(l)+b*by(l+1)))
  return
  end subroutine WaveODE
  !----------------------------------------------------------------------
  !StartODE initiates the wavefunctions at r=t1, output is y=ysta(1:2);
  !the wavefunctions are y(1)=u1 and y(2)=c*u2 in the notation u1,u2
  !of Rose;
  !the large Dirac wavefunction component is G=u1/r=y(1)/r;
  !in the limit of c=infinity, G is the Schroedinger wavefunction.
  subroutine StartODE(ir,l,lsta,ista,ysta)
  use param,only : alfa,c_light,cm2,wave_start
  use maths,only : sph_bessel_j
  use space,only : rx,nxz
  use enrgy,only : Vtot,kappa,emv0,q,lmax
  implicit none
  integer :: ir,l,lsta,i,ista,j
  real(8)  :: t1,ysta(2),bsl(0:lmax+1),alfaz,gamma,z0,a0,b0,v0,&
    gg,kg,ee

!potential V(r(1:2)) approximated by -2*z0/r(1:2)+v0 .
z0=-0.5d0*(rx(2,ir)*Vtot(1,ir)-rx(1,ir)*Vtot(2,ir))/(rx(2,ir)-rx(1,ir))
v0=(Vtot(2,ir)-Vtot(1,ir))/(rx(2,ir)-rx(1,ir))
!magnitude of j(l,kr) for small values of kr.
!PScalc sets ista=1 corresponding to l=lmin (lmin=0 or 1).
do i=ista,nxz(ir)
  call sph_bessel_j(lmax+1,q*rx(i,ir),bsl)
  if(abs(rx(i,ir)*bsl(l))>wave_start)then
    j=i
    exit
  endif
enddo
ista=max(j,1)
if(l<lsta)then
    !initial condition of Dirac equation.
    !Ref.: M.E.Rose, Relativistic Electron Theory (Wiley, 1961).
    alfaz=alfa*z0
    gamma=sqrt(kappa*kappa-alfaz*alfaz)
    if(kappa>0)then
      a0=alfaz
      b0=kappa+gamma
    else
      a0=kappa-gamma
      b0=alfaz
    endif
    gg=gamma+gamma+1.d0
    kg=kappa/gg
    gg=gamma/gg
    ee=1.d0+(emv0-v0)*cm2
    ysta(1)=rx(1,ir)*(a0+((1.d0-kg)*ee-gg)*c_light*b0*rx(1,ir))
    ysta(2)=rx(1,ir)*(b0+((1.d0+kg)*ee-gg)*c_light*a0*rx(1,ir))*c_light
else
  !initial condition when l>=lsta:
  !u1(l,r) is similar to r*j(l,kr) and approximately zero over a radial
  !interval [0,t1]. The integration of u1 is stable for r>r1, because
  !the linearly independent companion wavefunction to u1 is similar to
  !r*y(l,kr) and quickly vanishing with increasing r. Ref.: this work.
  t1=rx(ista,ir)
  ysta(1)=t1*bsl(l)
  ysta(2)=((kappa+l+1)*bsl(l)-q*t1*bsl(l+1))/&
             (1.d0+cm2*(emv0-Vtot(ista,ir)/t1))
  !where Vtot is tot.pot.*radius.
endif
  return
  end subroutine StartODE
  !----------------------------------------------------------------------
  subroutine PSnojump(psrm,psrp)
  use param,only : pi
  use space,only : nieq
  use enrgy,only : ne,lmax
  implicit none
  integer :: ir,l,npi,i
  real(8)  :: psrm(ne,0:lmax,nieq),psrp(ne,0:lmax,nieq),del(ne),dif,pih
pih=pi*0.5d0
do ir=1,nieq
  do l=0,lmax
    !smoothing psrm. Ref.: W. Moritz.
    del(1:ne)=psrm(1:ne,l,ir)
    npi=0
    do i=2,ne
      dif=del(i)-del(i-1)
      if(dif>=pih)then
        npi=npi-1
      elseif(dif<=-pih)then
        npi=npi+1
      endif
      psrm(i,l,ir)=del(i)+pi*npi
    enddo
    !smoothing psrp.
    del(1:ne)=psrp(1:ne,l,ir)
    npi=0
    do i=2,ne
      dif=del(i)-del(i-1)
      if(dif>=pih)then
        npi=npi-1
      elseif(dif<=-pih)then
        npi=npi+1
      endif
      psrp(i,l,ir)=del(i)+pi*npi
    enddo
    if(l>0)then
      !both PS curves running close together.
      if(abs(psrm(1,l,ir)-psrp(1,l,ir))>pih)then
        if(psrm(1,l,ir)<=0.d0) psrp(1:ne,l,ir)=psrp(1:ne,l,ir)-pi
        if(psrm(1,l,ir)> 0.d0) psrp(1:ne,l,ir)=psrp(1:ne,l,ir)+pi
      endif
    endif
  enddo !l
enddo !ir
  return
  end subroutine PSnojump
  !----------------------------------------------------------------------
  !INNER_POTENTIAL approximates the inner potential by the expression
  ! V0=a+b/sqrt(e+c).
  !Optimum parameters a,b,c are found by Differential Evolution Method.
  subroutine InnerPotential(V0tab,v0c)
  use param,only : compound,SelVxc,rydb
  use maths,only : DE,Dim_XC,NP,method,itermax,strategy,refresh
  use maths,only : F_XC,F_CR,VTR,CR_XC,itval
  use enrgy,only : ne,eev,errV01,errV02,v0ev
  implicit none
  integer :: i,m,jc,jf
  real(8) :: V0tab(ne),v0c(9),a,dev1,dev2,v0tmp(ne),tmp
  real(8),allocatable :: XCmin(:),XCmax(:),bestmem_XC(:)
  character :: fil*80
 allocate(v0ev(ne))
 v0ev=V0tab*rydb
 jc=len_trim(compound)
 fil=compound(1:jc)//'.V0vsE'
 jf=len_trim(fil)
 open(62,file=fil(1:jf),status='replace')
 !
 ! constant inner potential.
 if(SelVxc=='C')then
   write(61,'(a,f6.2)')'constant inner potential, V0vsE=',v0ev(1)
   v0c=0.d0;  v0c(1)=v0ev(1)
   write(62,890) eev(1),v0ev(1)
   write(62,890) eev(ne),v0ev(1)
   close(62)
   return
 endif
 !
 !ENERGY RANGE BELOW E=V0min.
 Dim_XC=3;  NP=10*Dim_XC
 allocate(XCmin(3),XCmax(3),bestmem_XC(3))
 XCmin=(/-3.0d0,-3.0d0,-3.0d0/)
 XCmax=(/ 3.0d0, 3.0d0, 3.0d0/)
 VTR=1.d-03
 F_XC=0.5d0;  CR_XC=0.8d0
 method=(/0,1,0/)
 strategy=2; F_CR=0.8d0
 itermax=1000; refresh=2000
 call DE(errV02,XCmin,XCmax,bestmem_XC)
 v0c(6)=bestmem_XC(1)
 v0c(7)=bestmem_XC(2)
 v0c(8)=bestmem_XC(3)
 m=minloc(v0ev,1)
 v0c(9)=eev(m)
 a=minval(v0ev); v0tmp=a; v0c(4)=a; v0c(5)=v0ev(1)
 dev2=0.d0
 do i=1,m
   tmp=(eev(i)/eev(m))**2
   v0tmp(i)=v0c(5)+tmp*(v0c(6)+tmp*(v0c(7)+tmp*v0c(8)))
   write(62,890) eev(i),v0ev(i),v0tmp(i)
   dev2=dev2+(v0ev(i)-v0tmp(i))**2
 enddo
 dev2=sqrt(dev2/m)
 !
 !ENERGY RANGE ABOVE E=V0min.
 Dim_XC=3;  NP=10*Dim_XC
 XCmin= (/-0.2d0,-100.d0, 0.d0/)
 XCmax= (/ 0.2d0, -10.d0,40.d0/)
 VTR=1.d-03
 F_XC=0.5d0;  CR_XC=0.8d0
 method=(/0,1,0/)
 strategy=2; F_CR=0.8d0
 itermax=1000; refresh=2000
 call DE(errV01,XCmin,XCmax,bestmem_XC)
 v0c(1:3)=bestmem_XC(1:3)
 a=minval(v0ev); v0tmp=a; v0c(4)=a
 dev1=0.d0
 do i=m+1,ne
   v0tmp(i)=max(v0c(1)+v0c(2)/sqrt(eev(i)+v0c(3)),v0c(4))
   write(62,890) eev(i),v0ev(i),v0tmp(i)
   dev1=dev1+(v0ev(i)-v0tmp(i))**2
 enddo
 close(62)
 dev1=sqrt(dev1/(ne-m))
 890 format(3es14.6)
 !
 !do V0vsE_analytic.
 open (63,file=fil(1:jf)//'_analytic',status='replace')
 write(63,'(a,f5.2)')'Em=',v0c(9)
 write(63,'(a)')'E<Em: V0vsE=c5+c6*(E/Em)**2+c7*(E/Em)**4+c8*(E/Em)**6, where'
 write(63,'(a,f6.2,3f8.4)')'ci=',v0c(5:8)
 write(63,'(a)')'E>Em: V0vsE=max[c1+c2/sqrt(E+c3),c4], where'
 write(63,'(a,f6.2,3f8.2)')'ci=',v0c(1:4)
 write(63,'(a)')'V0vsE and V0vsE_analytic coplotted by "xmgrace -nxy V0vsE".'
 close(63)
 !enddo V0vsE_analytic.
 !
 !interstitial potential variation.
 write(61,902) ne,v0ev(1)-v0ev(ne)
 write(61,904) itval,dev1
 write(*,902) ne,v0ev(1)-v0ev(ne)
 write(*,904) itval,dev1
 902 format('V0vsE(1)-V0vsE(',i3,')=',f6.2,' eV')
 904 format('iter,misfit= ',i4,f6.2,' eV')
  return
  end subroutine InnerPotential
  !----------------------------------------------------------------------
  !PStab transforms the set of phase shifts delta(l,E-V0) onto
  !a uniform E-V0 grid using spline interpolation.
  subroutine PStab(delta,tag,V0tab)
  use param,only : BmSel,rydb,id_at,SelVCoul,SelVxc,SelVMadl,SelMTrad,compound
  use maths,only : spcoef,svalue
  use space,only : nieq
  use enrgy,only : ne,eev,lmax,v0ev
  implicit none
  character :: tag*1,txt*100,txu*80,pid*80
  integer :: i,ir,l,jc,jt,jp,ju
  real(8) :: delta(ne,0:lmax,nieq),V0tab(ne),emv0(ne)
  !energy grid inside crystal.
  emv0=eev-V0tab*rydb
  !
  !phase shift table PS.
  jc=len_trim(compound)
  txt='-VCoul '//SelVCoul//' -Vxc '//SelVxc//' -VMadl '//SelVMadl//' -MTrad '//SelMTrad
  jt=len_trim(txt)
  if(SelVxc=='E')then
    txu=' Ecryst=Einc-V0vsE(Einc)'; ju=len_trim(txu)
  elseif(SelVxc=='C')then
    write(txu,'(a,f5.2,a)')' Ecryst=Einc',v0ev(1),' eV';  ju=len_trim(txu)
  endif
  do ir=1,nieq
    if(BmSel(ir)/='y')cycle
    if(tag=='+')then
      pid=compound(1:jc)//'.PS.su.'//id_at(ir); jp=len_trim(pid)
    elseif(tag=='-')then
      pid=compound(1:jc)//'.PS.sd.'//id_at(ir); jp=len_trim(pid)
    else
      pid=compound(1:jc)//'.PS.sl.'//id_at(ir); jp=len_trim(pid)
    endif
    open(40,file=pid(1:jp),status='replace')
    write(40,900) lmax+1,txt(1:jt)//': '//txu(1:ju)
    do i=1,ne
      write(40,910) emv0(i),(delta(i,l,ir),l=0,lmax)
    enddo
    close(40)
  enddo
  return
  900 format(i2,1x,a)
  910 format(f8.4,20f10.6)
  end subroutine PStab
  !----------------------------------------------------------------------
  subroutine sigma(lxq,psrm,psrp,V0tab)
  !SIGMA calculates differential cross section and total cross section in
  !atomic units (Bohr radius squared).
  use param,only : Omega,Theta,Sherman,CrossSection,BmSel,pi,rydb,cm2,id_at
  use space,only : nieq
  use enrgy,only : ne,eev,lmax
  implicit none
  character :: id_e*7
  integer,parameter :: ntheta=360
  integer :: i,j,l,ir,itheta,lxq(ne,nieq)
  real(8) :: psrm(ne,0:lmax,nieq),&
    psrp(ne,0:lmax,nieq),facr(ne),s(ne),dsdo(0:ntheta,ne),&
    shf(0:ntheta,ne),V0tab(ne),dtheta,emv0,tpi,fpi,frac,fg,rf,rg
  real(8) :: p(0:ntheta,0:lmax,0:1)
  complex(8) :: cf(0:lmax),cg(0:lmax)
  complex(8) :: f,g,cm,cp
  complex(8),parameter :: ci=(0.d0,1.d0)

!prefactor |1/ik|**2 accompanying exp(i*ps)*sin(ps).
tpi=2.d0*pi; fpi=4.d0*pi; dtheta=pi/ntheta
do i=1,ne
  emv0=eev(i)/rydb-V0tab(i)
  facr(i)=fpi/(emv0*(1.d0+emv0*cm2))
enddo
 call legendre_polynomial
do ir=1,nieq
  if(BmSel(ir)/='y')cycle
  do i=1,ne
    s(i)=0.d0
    do l=lxq(i,ir),0,-1 !sum from small to large phase shifts.
      cm=exp(ci*psrm(i,l,ir))*sin(psrm(i,l,ir))
      cp=exp(ci*psrp(i,l,ir))*sin(psrp(i,l,ir))
      cf(l)=(l+1)*cm+l*cp
      cg(l)=cp-cm
      rf=real(cf(l)*conjg(cf(l)),8)
      rg=real(cg(l)*conjg(cg(l)),8)
      s(i)=s(i)+(rf+(l*l+l)*rg)/(l+l+1)
    enddo
    s(i)=s(i)*facr(i)
    !
    do itheta=0,ntheta
      f=cmplx(0.d0,0.d0,8)
      g=cmplx(0.d0,0.d0,8)
      do l=lxq(i,ir),0,-1 !sum from small to large phase shifts.
        f=f+cf(l)*p(itheta,l,0)
        g=g+cg(l)*p(itheta,l,1)
      enddo
      fg=real(f*conjg(f)+g*conjg(g),8)
      dsdo(itheta,i)=fg*facr(i)
      shf(itheta,i)=real(ci*(f*conjg(g)-conjg(f)*g)/fg,8)
    enddo !itheta
  enddo !i
  !
  if(CrossSection=='y')then
    open(10,file='sigma'//id_at(ir),status='replace')
    do i=1,ne
      write(10,910) eev(i)-V0tab(i)*rydb,s(i)
    enddo
    close(10)
  endif
  !
  frac=180.d0/ntheta
  do i=1,ne
    write(id_e,'(i0)') nint(eev(i))
    j=len_trim(id_e)
    if(Omega=='y')then
      open(10,file='dsdo'//id_e(1:j)//id_at(ir),status='replace')
      do itheta=0,ntheta
        write(10,900) itheta*frac,log10(dsdo(itheta,i))
      enddo
      close(10)
    endif
    if(Theta=='y')then
      open(10,file='dsdt'//id_e(1:j)//id_at(ir),status='replace')
      do itheta=0,ntheta
        write(10,900)itheta*frac,dsdo(itheta,i)*sin(itheta*dtheta)*tpi
      enddo
      close(10)
    endif
    if(Sherman=='y')then
      open(10,file='shf'//id_e(1:j)//id_at(ir),status='replace')
      do itheta=0,ntheta
        write(10,900) itheta*frac,shf(itheta,i)
      enddo
      close(10)
    endif
  enddo !i
enddo !ir
  900 format(f8.3,20es14.6)
  910 format(f8.3,es14.6)
  return
  end subroutine sigma
  !-------------------------------------------------------------------
  !LEGENDRE_POLYNOMIAL
  !dtheta      =pi/ntheta,
  !p(theta,l,m)=legendre function of degree l, and order m.
  !qp(l,0)     =legendre polynomial of degree l for a given angle,
  !qp(l,1)     =associated legendre polynomial of degree l and order 1.
  !
  !Test in quatric precision:
  !The recurrence relations of the legendre polynomials are used.
  !The stability the recurrence was tested by a calculation from
  !l=0 to 200 and back to 0. The maximum disagreement for l=0:200 and
  !theta=1:179 degrees was found to be 4.E-30. One has reason to
  !believe that all digits are correct in the double precision output.
  !Ref.: Abramowitz and Stegun, page XIII and Sec. 8.5.3.
  subroutine legendre_polynomial
  use enrgy,only : lmax
  implicit none
  integer :: i,l,ll
  integer,parameter :: ntheta=360
  real(8) :: x,c,s,dtheta
  real(8) :: p(0:ntheta,0:lmax,0:1)
  real(8),allocatable :: qp(:,:)
  dtheta=acos(-1.d0)/ntheta
  !theta equal to 0 and pi.
  p(0,0:lmax,0)=1.d0
  do l=0,lmax,2
    p(ntheta,l,0)=1.d0
  enddo
  do l=1,lmax,2
    p(ntheta,l,0)=-1.d0
  enddo
  p(0,0:lmax,1)=0.d0; p(ntheta,0:lmax,1)=0.d0
  !theta between 0 and pi.
  allocate(qp(0:lmax,0:1))
  do i=1,ntheta-1
    x=dtheta*i
    c=cos(x);    s=sin(x)
    qp(0,0)=1.d0; qp(0,1)=0.d0
    qp(1,0)=c;    qp(1,1)=s
    ll=1
    do l=2,lmax
      ll=ll+2
      qp(l,0)=(ll*c*qp(l-1,0)-(l-1)*qp(l-2,0))/l
      qp(l,1)=(ll*c*qp(l-1,1)- l   *qp(l-2,1))/(l-1)
    enddo
    p(i,0:lmax,0)=qp(0:lmax,0)
    p(i,0:lmax,1)=qp(0:lmax,1)
  enddo
  return
  end subroutine legendre_polynomial
end module calPS

!==============================================================================
 subroutine EEASiSSS(input_file)

 !'Elastic Electron-Atom Scattering in Solids and Solid Surfaces',
 !author: John Rundgren, KTH Royal Institute of Technology, Stockholm,
 !Sweden, jru@kth.se.
 use param,only : SelVCoul
 use enrgy,only : WIEN2kData,FreeAtomData,superpose
 use calPS,only : PartialWaveMethod
 implicit none
 integer :: time(8)
 character(len=255) :: input_file
!f2py character(len=255), optional, intent(in) :: input_file = "inputX"
 intrinsic date_and_time
 !
 open(61,file='logfile',status='replace')
 call date_and_time(values=time)
 write(61,990) time(1:3), time(5:7)
 990 format('EEASiSSS time: ',i4,'/',i0,'/',i0,'  ',i0,':',i0,':',i0)
 !990 format('EEASiSSS time: ',i4,'/',i2,'/',i2,'  ',i2,':',i2,':',i2)
 write(61,'(a)')'READ inputX:'
 !SETUP.
 call ModesOfCalc(input_file)
 call Structure(input_file)
 call NeighborShells
 call EnergyGrid
 !EXECUTION.
 if(SelVCoul=='W')then
   call WIEN2kData
 elseif(SelVCoul=='S')then
   call FreeAtomData
   call superpose
 endif
 call PartialWaveMethod
 !
 call date_and_time(values=time)
 write(61,990) time(1:3),time(5:7)

 return

 end subroutine
!---------------------------------------------------------------------
  subroutine ModesOfCalc(input_file)
  use param,only : SelVxc,SelMTrad,SelVMadl,SelOutp,SpinPS,Omega,Theta
  use param,only : Sherman,CrossSection,RhoPot,PSvsE,UvsR,rydb,SelVCoul !,ATAinp
  use maths,only : method,itermax,strategy,refresh,F_XC,F_CR
  use maths,only : VTR,CR_XC,NPfac
  use enrgy,only : lmax,ee1,ee2,dee,rMTweight,dVeV
  implicit none
  character :: txt*7
  character(len=255), intent(in) :: input_file
  integer :: nes
  real(8) :: es1,des
  open(5,file=trim(input_file))
  do
    read(5,*) txt
    if(txt=='OPTIONS') exit
  enddo
  read(5,*) SelVCoul
    write(61,'(2a)')'SelVCoul==',SelVCoul
    if(.not.(SelVCoul=='S'.or.SelVCoul=='W'))goto 3
  read(5,*) SelVxc
    write(61,'(2a)')'SelVxc==',SelVxc
    if(.not.(SelVxc=='E'.or.SelVxc=='C'))goto 3
  read(5,*) SelVMadl
    write(61,'(2a)')'SelVMadl==',SelVMadl
    if(.not.(SelVMadl=='y'.or.SelVMadl=='n'))goto 3
  read(5,*) SelMTrad
    write(61,'(2a)')'SelMTrad==',SelMTrad
    if(.not.(SelMTrad=='D'.or.SelMTrad=='M'.or.SelMTrad=='I'))goto 3
  read(5,*) SelOutp
    write(61,'(2a)')'SelOutp==',SelOutp
    if(.not.(SelOutp=='p'.or.SelOutp=='s'.or.SelOutp=='d'))goto 3
  read(5,*) SpinPS
    write(61,'(2a)')'SpinPS==',SpinPS
    if(.not.(SpinPS=='y'.or.SpinPS=='n'))goto 3
  read(5,*) Omega
    write(61,'(2a)')'Omega==',Omega
    if(.not.(Omega=='y'.or.Omega=='n'))goto 3
  read(5,*) Theta
    write(61,'(2a)')'Theta==',Theta
    if(.not.(Theta=='y'.or.Theta=='n'))goto 3
  read(5,*) Sherman
    write(61,'(2a)')'Sherman==',Sherman
    if(.not.(Sherman=='y'.or.Sherman=='n'))goto 3
  read(5,*) CrossSection
    write(61,'(2a)')'CrossSection==',CrossSection
    if(.not.(CrossSection=='y'.or.CrossSection=='n'))goto 3
  read(5,*) RhoPot
    write(61,'(2a)')'RhoPot==',RhoPot
    if(.not.(RhoPot(1:1)=='y'.or.Rhopot(1:1)=='n'.or.&
             RhoPot(2:2)=='y'.or.RhoPot(2:2)=='n'))goto 3
  read(5,*) PSvsE
    write(61,'(2a)')'PSvsE==',PSvsE
    if(.not.(PSvsE=='y'.or.PSvsE=='n'))goto 3
  read(5,*) UvsR
    write(61,'(2a)')'UvsR==',UvsR
    if(.not.(UvsR=='y'.or.UvsR=='n'))goto 3
  goto 4
  3 write(61,'(a)')'stop, incorrect input character'; stop
  4 continue

  !energy range.
  read(5,*) ee1,ee2,dee
    write(61,'(a,3f6.0)')'ee1,ee2,dee=',ee1,ee2,dee
  read(5,*) es1,des,nes
    write(61,'(a,2f6.0,i6)')'es1,des,nes=',es1,des,nes

  !no. phase shifts.
  read(5,*) lmax
    write(61,'(a,i4)')'lmax=',lmax
  read(5,*) rMTweight
    write(61,'(a,f5.2)')'rMTweight=',rMTweight
  read(5,*) dVeV
    write(61,'(a,f5.2)')'dVeV=',dVeV
  read(5,*)
    write(61,'(a)')'DIFFERENTIAL EVOLUTION'
  read(5,*) VTR ! input in eV.
    write(61,'(a,es9.2)')'VTR(eV)=',VTR
    VTR=VTR/rydb ! execution in Rydbergs.
  read(5,*) NPfac,F_XC,CR_XC
    write(61,'(a,i2,2f6.2)')'NPfac,F_XC, CR_XC=',NPfac,F_XC,CR_XC
  read(5,*) method
    write(61,'(a,3i2)')'method=',method
  read(5,*) strategy,F_CR
    write(61,'(a,i1,f6.2)')'strategy, F_CR=',strategy,F_CR
  read(5,*) itermax,refresh
    write(61,'(a,2i10)')'itermax, refresh=',itermax,refresh
  close(5)
  return
  end subroutine ModesOfCalc
  !---------------------------------------------------------------------
  subroutine Structure(input_file)
  use param,only : compound,W2kData,id_at,BmSel,SelBm,Nsel,nshell,&
    pi,fpi3,facmy,bohr,rydb,idW
  use space,only : rcX,rk,ad,ia,na,ncon,nx,dx,nxz,nieq,neq,ieq,&
    nlatp,rmt,rmt0,rmtn,rmtx,volUC
  use enrgy,only : z,cls,Vrmt,Vstp
  implicit none
  character :: txt*80
  character(len=255), intent(in) :: input_file
  integer :: i,ir,j,k
  integer,allocatable :: idZ(:),idL(:)
  real(8) :: UOL,rcc1,rcc2,rcc3,rc(3,3)
  real(8),allocatable :: rk_tmp(:,:)
  real(8),parameter :: thrd=1.d0/3.d0
  pi=acos(-1.d0); fpi3=pi/0.75d0; facmy=(18.d0/(pi*pi))**thrd
 !
 open(5,file=trim(input_file))
 read(5,*) txt
    write(61,'(a)') txt(1:10)
 read(5,'(a)') compound; j=len_trim(compound)
    compound=adjustl(compound(1:j)); j=len_trim(compound)
    write(61,'(a)') compound(1:j)
 read(5,'(a)') W2kData; j=len_trim(W2kData)
    write(61,'(a)') W2kData(1:j)
    do i=j,1,-1
      if(W2kData(i:i)=='/')then
        k=i
        exit
      endif
    enddo
 !
 !STRUCTURE data.
 !UOL = UnitOfLength RelativeToResearcher'sUnits ForBohrInside.
 read(5,*) UOL
 write(61,200) UOL
 200 format('UOL=',f10.6)
 !
 !rcX(i,j)=i'th coordinate of the j'th unit cell vector (Researcher'sUnits).
 write(61,'(a)')'UC vectors (B)'
 read(5,*) rcX(:,1)
 read(5,*) rcX(:,2)
 read(5,*) rcX(:,3)
 rcX=rcX*UOL
 write(61,'(3f9.4)') rcX(:,1)
 write(61,'(3f9.4)') rcX(:,2)
 write(61,'(3f9.4)') rcX(:,3)
 !
 !nieq    =number of ineqivalent atoms in unit cell,
 !neq(ir) =number of equivalent atoms of type ir,
 !z(ir)   =atomic number of type ir,
 !nlatp   =# lattice points in unit cell.
 read(5,*) nieq
 write(61,'(a,i0)')'nieq ',nieq
 allocate(z(nieq),cls(nieq),neq(nieq),nx(nieq),dx(nieq),id_at(nieq))
 allocate(BmSel(nieq),SelBm(nieq),nxz(nieq),rmt(nieq),rmt0(nieq))
 allocate(rmtn(nieq),rmtx(nieq),ieq(nieq),Vrmt(nieq),Vstp(nieq))
 allocate(ad(nshell,nieq),na(nshell,nieq),ia(nshell,nieq),ncon(nieq))
 allocate(rk_tmp(3,1000),idZ(nieq),idL(nieq),idW(nieq))
 nlatp=0; Nsel=0
 !rmtn, rmtx in Bohr.
 do ir=1,nieq
   read(5,*) neq(ir),idZ(ir),idL(ir),idW(ir),rmtn(ir),rmtx(ir),cls(ir),BmSel(ir)
   z(ir)=real(idZ(ir),8)
   write(id_at(ir),'(i0,a,i0)') idZ(ir),'.',idL(ir)
   !j=len_trim(id_at(ir))
   rmtn(ir)=rmtn(ir)*UOL; rmtx(ir)=rmtx(ir)*UOL; cls(ir)=cls(ir)/rydb
   write(61,'(a)')'type neq  id_at idW rMTmin(B) rMTmax(B) cls(Ry)  Bmsel '
   write(61,203)&
   ir,neq(ir),id_at(ir)(1:5),idW(ir),rmtn(ir),rmtx(ir),cls(ir),Bmsel(ir)
   203 format(i3,i5,2x,a,i4,f8.4,2f10.4,3x,a)
   !
   !PS selection.
   if(.not.(BmSel(ir)=='y'.or.BmSel(ir)=='n'))then
     write(61,'(a)')'BmSel not "y" or "n", stop';  stop
   endif
   if(BmSel(ir)=='y')then
     Nsel=Nsel+1
     SelBm(Nsel)=ir
   endif
   !
   !inequivalent lattice points.
   do j=1,neq(ir)
     nlatp=nlatp+1
     if(nlatp>1000)then
       write(61,'(a)')'nlatp>dimension 1000, stop'; stop
     endif
     read(5,*) rk_tmp(1:3,nlatp)
     rk_tmp(1:3,nlatp)=rk_tmp(1:3,nlatp)*UOL
     write(61,205) rk_tmp(1:3,nlatp)*bohr,rk_tmp(1:3,nlatp)
   enddo
 enddo
 205 format(3f9.4,' (A)'/3f9.4,' (B)')
 !
 write(61,'(a,i3)')'#atoms in UC=',nlatp
 allocate(rk(1:3,nlatp))
 do j=1,nlatp
   rk(1:3,j)=rk_tmp(1:3,j)
 enddo
 deallocate(rk_tmp)
 !
!volUC = volume of unit cell, where "unit cell" signifies:
!for bulk, crystallographic unit cell;
!for slab, crystallographic unit cell inclusive vacuum slab;
!rcc1 and rcc2 || to the surface.
rc=rcX
rcc1=rc(2,1)*rc(3,2)-rc(3,1)*rc(2,2)
rcc2=rc(3,1)*rc(1,2)-rc(1,1)*rc(3,2)
rcc3=rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2)
volUC=abs(rc(1,3)*rcc1+rc(2,3)*rcc2+rc(3,3)*rcc3)
write(61,'(a,f8.2)')'volUC (B**3)=',volUC
return
end subroutine Structure
!----------------------------------------------------------------------
!NEIGHBOR_SHELLS is nearest neighbor data for atoms in a crystal
!structure.
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package,
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/Activity.

!New in this code is ARMAX(1:nieq). These radii sort out the necessary
!and sufficient lattice points for getting neighbor contributions of
!a given accuracy to charge density and potential.
!The main part of the subroutine is a fortran90 editon of the Ref.

!input:
!rc(i,j) =the i'th coordinate of the j'th axis of the unit cell,
!rk(i,j) =the i'th coordinate of the j'th atom in the unit cell,
!neq(ir) =the number of equivalent atoms of type ir.
!output:
!ia(j,ir)=the type of atoms in the j'th neighbor shell,
!na(j,ir)=the number of atoms in the j'th shell,
!ad(j,ir)=the radius of the j'th shell,
!ncon(ir)=no. of ir type shells, limited by ARMAX(ir) and NSHELL.
subroutine NeighborShells
use param,only : bohr,Elem,nshell
use space,only : nieq,neq,rcX,rk,ad,ia,na,ncon,volUC,volXC
use enrgy,only : z
implicit none
integer :: j,jr,jjr,k,kr,ic,iic,jjc,jc,ir,nc,kc,i,jx,jy,jz
real(8) :: rc(3,3),rj(3),dr,r,a,b
real(8),parameter :: rmax=15.d0
!
rc=rcX
!search over adjacent unit cells to include nshell nearest neighbors.
ia=0
na=0
ad=1.d+33
do jx=-nshell,nshell
  do jy=-nshell,nshell
    do jz=-nshell,nshell
      do j=1,3
        rj(j)=jx*rc(j,1)+jy*rc(j,2)+jz*rc(j,3)
      enddo
      !rj is current unit cell origin. for each atom in this unit cell,
      !find displacement r from kr-type atom in basic unit cell.
      j=0
      do jr=1,nieq
        do jjr=1,neq(jr)
          j=j+1
          k=1
          do kr=1,nieq
            r=sqrt((rj(1)+rk(1,j)-rk(1,k))**2+&
                   (rj(2)+rk(2,j)-rk(2,k))**2+&
                   (rj(3)+rk(3,j)-rk(3,k))**2 )
            !fixed limit for r.
            if(r>rmax)then
              k=k+neq(kr)
              cycle
            endif
            !compare r with nearest neighbor distances already found.
            ic=0
            do
              ic=ic+1
              if(ic>nshell)then
                k=k+neq(kr)
                exit
              endif
              dr=r-ad(ic,kr); if(abs(dr)<1.d-03) dr=0.0d0
              if(dr>0.d0)then
                cycle
              elseif(dr==0.d0)then
                if(ia(ic,kr)/=jr)then
                  cycle
                else
                  na(ic,kr)=na(ic,kr)+1
                  k=k+neq(kr)
                  exit
                endif
              else
                if(ic<nshell)then
                  iic=ic+1
                  do jjc=iic,nshell
                    jc=nshell+iic-jjc
                    ia(jc,kr)=ia(jc-1,kr)
                    na(jc,kr)=na(jc-1,kr)
                    ad(jc,kr)=ad(jc-1,kr)
                  enddo
                endif
                ia(ic,kr)=jr
                na(ic,kr)=1
                ad(ic,kr)=r
                k=k+neq(kr)
                exit
              endif
            enddo !ic
          enddo !kr
        enddo !jjr
      enddo !jr
    enddo !jz
  enddo !jy
enddo !jx
do ir=1,nieq
  ncon(ir)=0
  do ic=1,nshell
    if(na(ic,ir)>0) ncon(ir)=ncon(ir)+1
  enddo
enddo
!display.
write(61,'(/a)')'NN shells'
do ir=1,nieq
  j=0
  do i=2,ncon(ir)
    j=j+na(i,ir)  !j counts the neighbors about atom ir.
  enddo
  write(61,200) ir,Elem(int(z(ir))),ncon(ir)-1,j
  nc=ncon(ir)
  ic=(nc-1)/8+1
  kc=0
  do i=1,ic
    jc=kc+1
    kc=min(nc,kc+8)
    write(61,202) (ad(j,ir)*bohr,j=jc,kc)
    write(61,204) (ad(j,ir),j=jc,kc)
    write(61,206) (na(j,ir),j=jc,kc)
    write(61,208) (ia(j,ir),Elem(int(z(ia(j,ir)))),j=jc,kc)
  enddo
enddo
write(61,*)
write(61,'(a)')'NN distances (B).'
write(61,212) (i,i=1,nieq)
write(61,214) (ad(2,i),i=1,nieq)
!
!calculating volXC.
!dz of volUC.
a=rcX(3,3)-rcX(3,1)
!dz of volXC.
b=maxval(rk(3,:))-minval(rk(3,:))+0.5d0*(ad(2,1)+ad(2,nieq))
volXC=volUC*(b/a)
write(61,'(/a,2f8.2/)')'volUC,volXC (B**3)=',volUC,volXC
return
200 format(/'   atom   ',i2,'_',a,i6,' neighbor shells,',i3,' neighbors')
202 format('dist.(A)',10f8.4:)
204 format('dist.(B)',10f8.4:)
206 format('  number ',10(3x,i2,3x:))
208 format('      NN  ',10(i2,'_',a,3x:))
212 format(10(2x,i2,4x:))
214 format(10(f8.4:))
end subroutine NeighborShells
!---------------------------------------------------------------------
  subroutine EnergyGrid
  !ENERGY GRID can be uniform or non-uniform (by keeping or remowing !'s).
  use space,only : nieq,rmt1
  use enrgy,only : ne,eev,ee1,ee2,dee,V01
  implicit none
  integer :: je
  real(8) :: E,egrid(1000)
  je=0; E=ee1-dee
  1001  continue
        if(E>=ee2) goto 1002
        je=je+1
        if    (E>=600.d0)then; E=E+80.d0;
        !elseif(E>=400.d0)then; E=E+40.d0;
        !elseif(E>=200.d0)then; E=E+20.d0;
        !elseif(E>= 80.d0)then; E=E+10.d0;
        !elseif(E>= 30.d0)then; E=E+ 2.d0;
        else; E=E+dee;
        endif
        E=min(E,ee2); egrid(je)=E
        goto 1001
  1002  continue
  ne=je
  allocate(eev(ne),V01(ne),rmt1(nieq,ne))
  eev(1:ne)=egrid(1:ne)
  return
  end subroutine EnergyGrid
