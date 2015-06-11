!      program hartfock
! Author: Eric L. Shirley, Optical Technology Division, NIST,
!         Gaithersburg, MD 20899-8441 U.S.A.
!         eric.shirley@nist.gov,
! August, 2000.
!
! ref.: E.L. Shirley, PhD thesis, University of Illinois (1991),
!       unpublished.
! ref.: S. Kotochigova, Z.H. Levine, E.L. Shirley, M.D. Stiles,
!       and C.W. Clark, Phys. Rev. A 55, 191 (1997); 
!       ibid. A 56, 5191 (1997).
!
! The author appreciates acknowledgment in publications, by citing the
! above two references when the code is first mentioned, e.g. 
! " ... using an atomic program [*,*], we have ... "
! 
! This program is a code originally written by Eric Shirley during his 
! thesis work at the University of Illinois at Urbana-Champaign.  
! Since then, the code has undergone many minor revisions, refinements,
! improvements in numerical methods used, and so forth. 
!
! As a research tool, the program is provided with the understanding 
! that it may be imperfect. However, many tests argue for its accuracy
! as follows:
!
!  [1] The above work by Kotochigova et al.
!  [2] Results cited in Shirley's thesis, as well the back-to-back 
!      papers, E.L. Shirley and R.M. Martin, Phys. Rev. B 47, 15404 
!      and 15413 (1993).
!  [3] Comparison against results published by Charlotte Froese Fischer,
!      The Hartree-Fock method for atoms: a numerical approach (Wiley,
!      New York, 1977).
!  [4] Band structure results obtained with pseudopotentials calculated
!      in the code.
!  [5] A comparison of Hartree-Fock results to OPM results by 
!      M.D. Stiles (unpublished).
!
! And there are other tests, too.  However, tests are not guarantees,
! as experience shows.
!
! Eric Shirley is committed to helping people use the code by providing
! a reasonable level of advice when special problems arise, or when the
! validity of results is in question (in the ! interest of all parties
! involved).
!----------------------------------------------------------------------
      module dprec
      implicit none
      !double precision definition,
      !http://fortranwiki.org/fortran/show/Real+precision
      integer,parameter :: dp=selected_real_kind(15,307)
      end module dprec
!----------------------------------------------------------------------
      module paths
      contains

      function normalize_windows_path(filepath)
        implicit none
        character(len=255), intent(in) :: filepath
        character(len=255) :: normalize_windows_path
        integer :: i

        normalize_windows_path = filepath
        do i=1,len_trim(filepath)
          if (filepath(i:i) == "/") then
            normalize_windows_path(i:i) = "\\"
          end if
        end do

        return
      end function

      function get_unix_path(filepath)
        implicit none
        character(len=255), intent(in) :: filepath
        character(len=255) :: get_unix_path
        integer :: i

        get_unix_path = filepath
        do i=1,len_trim(filepath)
          if (filepath(i:i) == "\\") then
            get_unix_path(i:i) = "/"
          end if
        end do

        return
      end function

      end module paths
!----------------------------------------------------------------------
      module PARAMS
      use dprec
      implicit none
      integer, parameter  :: iorbs = 36
      integer, parameter  :: lmax = 4
      integer, parameter  :: ihmax = 7
      integer, parameter  :: nrmax = 17010
      integer, parameter  :: ntmax = 10
      integer, parameter  :: npmax = 60
      real(dp) :: rmifac, rmafac, x137   
      !LD: parameter statements separate to real due to limitations of f2py
      parameter (rmifac = 1.d-06)        !0.00000001d0 !JR
      parameter (rmafac = 200.d0)        !800.d0       !JR
      parameter (x137 = 137.0359895d0)
      logical, parameter  :: verb = .true.
      character(len=255)  :: outdir,atdir,outfile
      integer io
      integer logunit                    !LD: file unit for log file (may be stdout)
      integer, parameter :: inpunit = 21 !LD: file unit for atorb input
      end module PARAMS
!----------------------------------------------------------------------
      !program hartfock
      subroutine hartfock(input_file,log_file,output_dir,atom_dir)
      use dprec
      use paths
      use PARAMS
      use iso_fortran_env ! requires compiler adhering to Fortran 2003 standard
      implicit none
      character(len=255):: input_file, log_file, output_dir, atom_dir
      integer, allocatable :: no(:),nl(:),nm(:),is(:)
      integer, allocatable :: njrc(:),ilp(:)
      real(dp), allocatable :: r(:),dr(:),r2(:)
      real(dp), allocatable :: ev(:),occ(:),xnj(:)
      real(dp), allocatable :: ek(:),phe(:,:),orb(:,:)
      real(dp), allocatable :: vctab(:,:)
      real(dp), allocatable :: vi(:,:)
      real(dp), allocatable :: cq(:)
      real(dp), allocatable :: vold(:),vnew(:), wgts(:)
      integer :: i,j,k,nel,nr,nst,iuflag,iu,ir,ixflag,vtry,isuse,nwgt
      integer ::ncore, imul( 0 : 2 )
      integer, parameter :: nwgmx = 10
      real(dp), allocatable :: mlt( : )
      real(dp) :: rel,alfa,etot,dl,zorig,xntot,rmin,rmax
      real(dp) :: rlast, potn, rad
      logical  :: done,ex
      character(len=1) :: ichar,txt
      character(len=3) :: mode
      character(len=10) :: atom
      character(len=255) :: env

      write(outdir,'(255(" "))')
      outdir(1:len_trim(output_dir)) = output_dir
      outdir = adjustl(outdir)

      if (len_trim(log_file) > 0) then
        logunit=61
        open(logunit,file=trim(log_file),status='replace')
        write(logunit,'(2a)') trim(outdir)//adjustl(trim(log_file)),' is logfile.'
      else
        !logunit=6 ! stdout (on ifort & gfortran) => no need to open unit
        logunit=output_unit ! portable stdout but needs Fortran 2003 compliant compiler
        !write(logunit,*) 'Logging to stdout:'
      end if

      !Check atom_dir is not empty
      write(atdir,'(255(" "))')
      if (len_trim(atom_dir) > 0) then
        atdir(1:len_trim(atom_dir)) = atom_dir
      else
        !Try $ATDIR environment variable
        call getenv('ATDIR',env)
        if (len_trim(env) > 0) then
            atdir(1:len_trim(env)) = env
        else
            !Fail-safe is current working directory
            atdir = '.'
        end if
      end if
      atdir = adjustl(atdir)

      ! make atlib directory if needed
      call getenv('SystemDrive',env)
      if (len_trim(env) > 0) then
        call system('mkdir '//trim(normalize_windows_path(atdir)))
      else
        call system('mkdir -p '//trim(atdir))
      end if

      open(inpunit,file=trim(adjustl(input_file)),status='old')


      allocate( no(iorbs),nl(iorbs),nm(iorbs),is(iorbs))
      allocate( njrc(4),ilp(iorbs),r(nrmax),dr(nrmax),r2(nrmax))
      allocate( ev(iorbs),occ(iorbs),xnj(iorbs))
      allocate( ek(iorbs),phe(nrmax,iorbs),orb(nrmax,iorbs))
      allocate( vi(nrmax,7),cq(nrmax),vctab(nrmax,0:3))
      allocate( vold(iorbs),vnew(iorbs), wgts( nwgmx ) )

      rel = 0.d0
      nst = 2
      read(inpunit,'(a)') ichar

!====================================
      do while ( ichar .ne. 'q' )
!====================================
      select case( ichar )
      case ( 'b' )
        call bachelet(vi,r,njrc,nr)
      case ( 'X' )
        read (inpunit,'(2x,1i4,2x,1a3)') nwgt,mode
        if (nwgt.gt.nwgmx) then
          write (6,*) 'bad nwgt'
          stop
        end if
        read (inpunit,*) (wgts(i),i=1,nwgt)
        close(unit=77)
        call grandop(nwgt,wgts,mode)
        open(unit=77,file='tmp',form='formatted',status='unknown')
        rewind 77
      case ( 'd' )
        read (inpunit,*) rel
      case ( 'x' )
        read (inpunit,*) alfa
      case ( 'a' )
        call abinitio(etot,rel,alfa,nr,r,&
                     dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,&
                     no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,&
                     vold,vnew,vtry,isuse)
      case ( 'i' )
!        if the chgden//atom file already exists, the Hartree-Fock
!        calculation is not repeated.
         read(inpunit,'(a)') atom                              !JR
         io=len_trim(atom)                                     !JR
         outfile=trim(outdir)//'/'//'chgden'//atom(1:io)       !LD: included output directory
         io=io+len_trim(outdir)+7                              !LD
         inquire(file=outfile(1:io),exist=ex)                  !JR
         if(ex)then                                            !JR
           write(logunit,'(a/a)') outfile(1:io),' exists, no recalculation'
           do k=1,1000                                         !JR
             read(inpunit,'(a)') txt                           !JR
             if(txt=='w') goto 900                             !JR
           enddo                                               !JR
         else                                                  !JR
           write(logunit,'(2a)') outfile(1:io),' calculation starts'!JR 
         endif                                                 !JR
         call initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,xntot,nel)
         do j=0,3
           do i=1,nr
             vctab(i,j)=0.d0
           end do
         end do
900      continue !JR
      case ( 's' )
        read(inpunit, * ) iu,k,ncore
        allocate( mlt( k ) )
        do i = 0, 2
          imul( i ) = 1
        end do
        i = 1
        do while ( i .le. ncore )
          imul( nl( i ) ) = - imul( nl( i ) )
          i = i + 1
        end do
        read(inpunit, * ) (ilp(j),j=1,k)
        do j = 1, k
          mlt( j ) = imul( nl( ilp( j ) ) )
        end do
        rewind iu
        open(unit=iu,file='orbvu',form='formatted',status='unknown')
        rewind iu
        do i=1,nr
          write (iu,'(1x,10f25.15)') &
           r(i),(mlt(j)*phe(i,ilp(j))/r(i),j=1,k)
        end do
        close(unit=iu)
        deallocate( mlt )
      case ( 'W' )
        call melwriter( nr, r, dr, nel, nl, nrmax, phe, ev )
      case ( 'w' )
        ixflag=1
        iu=-1
        ir=0
        call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,&
                   zorig,xntot,ixflag,nel,&
                   no,nl,xnj,is,ev,ek,occ,njrc,vi,cq,phe,orb)
      case ( 'r' )
        iu=-1
        ir=1
        call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,&
                   zorig,xntot,ixflag,nel,&
                   no,nl,xnj,is,ev,ek,occ,njrc,vi,cq,phe,orb)
        call setgrid(nr,rmin,rmax,r,dr,r2,dl)
      case ( 'u' )
        write(6,*) 'please enter iuflag. (0=u, 1=su, 2=r).'
        read(inpunit,*) iuflag
      case ( 'p' )
        call pseudo(etot,rel,alfa,nr,rmin,&
                   rmax,r,dr,r2,dl,&
                   phe,orb,njrc,vi,cq,zorig,xntot,nel,&
                   no,nl,nm,xnj,ev,occ,is,iuflag,vctab,&
                   vold,vnew,vtry,isuse)
      case ( 'C' )
        call vpscpp( nr, nrmax, r, r2, vi )
      case ( 'v' )
        call realspace(nr,nrmax,r,vi)
      case ( 'V' )
        call fourier(nr,r,dr,vi)
      case ( 'g' )
        call ppopt(nr,r,nrmax,vi,cq,nel,phe,occ)
        read(inpunit, * ) rad
        open( unit=99, file='precpw',&
             form='formatted', status='unknown' )
        rewind 99
        write ( 99, * ) nel, rad
        do i = 1, nel
          write ( 99, * ) nl( i )
          rlast = r( 1 )
          j = 2         
          done = .false.
          do while ( .not. done )
            if ( r( j ) .gt. rlast + 0.01d0 ) then
              potn = vi( j, 2 * nl( i ) + 1 ) + orb( j, i )
              write ( 99, '(3(2x,1e15.8))' ) r( j ), phe( j, i ), potn
              rlast = r( j )
              done = ( rlast .gt. rad )
            end if
            j = j + 1
          end do
        end do
      case ( 'k' )
        call mkkbfile(nel,nl,xnj,ev,dl,nr,r,r2,cq,&
                     nrmax,vi,orb,zorig,njrc)
      case ( 'F' )
        call sigfit(zorig,nrmax,nr,r,dr,nel,nl,ev,phe)
      case ( 'f' )
        call sigfit2(nel,nl,ev,nr,nrmax,dr,r,phe)
      case ( 'G' )
        call leadbehv(nel,nl,nr,nrmax,phe,r)
      case ( 'c' )
        call mkvctab( nr, nrmax, vctab, r, r2 )
      case ( 'z' )
        call corepot( nr, nrmax, vi )
      end select

!====================================
        read(inpunit,'(a)',end=999) ichar  !JR
      end do
!====================================
999   stop                           !JR
      end subroutine hartfock
      !end program hartfock
!-------------------------------------------------------------------
      subroutine sigfit( z, ldim, nr, r, dr, nlev, llev, eps, phi )
      use dprec
      use PARAMS
      implicit none
      integer :: ldim, nr, nlev
      real(dp) :: z
      integer :: llev( nlev )
      real(dp) :: r( nr ), dr( nr ), eps( nlev )
      real(dp) :: phi( ldim, nlev )

      integer i, j, k, n
      real(dp) corpol, expt, rl, rh, sd, sl, sh, rc, sc

      real(dp), external :: sig

      k = int(z + 0.1d0)  !JR
      read  (inpunit,*) corpol,n
      do i=1,n
        read  (inpunit,*) j,expt,rl,rh
        sd=abs(abs(expt)-abs(eps(j)))
        sl=sig(corpol,rl,nr,r,dr,phi(1,j))
        sh=sig(corpol,rh,nr,r,dr,phi(1,j))
        if (((sd-sh)*(sd-sl).ge.0.d0).and.&
           (abs(rl-rh).gt.0.0001d0)      ) then
          rc=-1.d0
        else
  324     rc=rl+(rh-rl)/2.d0
          sc=sig(corpol,rc,nr,r,dr,phi(1,j))
          write (6,'(1x,8f9.4)') rl,sl,rh,sh,rc,sc,sd,eps(j)-sc
          if (sc.gt.sd) then
            rl=rc
            sl=sc
          else
            rh=rc
            sh=sc
          end if
          if ((abs(sc-sd).gt.0.000001d0).and.&
             (abs(rl-rh).gt.0.0001d0)       ) go to 324
        end if
        write (6,'(1x,2i4,1f10.4)') k,llev(j),rc
        write (16,'(1x,2i4,1f10.4)') k,llev(j),rc
      end do

      return
      end subroutine sigfit
!----------------------------------------------------------------------
      subroutine sigfit2( nlev, llev, eps, nr, ldim, dr, r, phi )
      use dprec
      use PARAMS
      implicit none
      integer nlev, nr, ldim
      integer llev( nlev )
      real(dp) eps( nlev ), dr( nr ), r( nr )
      real(dp) phi( ldim, nlev )
      integer i, k, l, n
      real(dp) lim
      real(dp) tmp, eav, e1, e2, corpol
      real(dp) sc, sd, sl, sh, rc, rl, rh
      character(len=10) :: limunit, eunit
!
      read  (inpunit,*) corpol
      read  (inpunit,*) limunit, lim
      read  (inpunit,*) i, n
      l = llev( i )
      if (n.eq.1) then
        read (inpunit,*) eunit, eav
      else
        read (inpunit,*) eunit, e1, e2
        eav=(e1*real(l,dp)+e2*real(l+1,dp))/real(2*l+1,dp)
      end if
      call convert( limunit, lim )
      call convert( eunit, eav )
      eav=abs(lim)-abs(eav)
      sd=abs(abs(eav)-abs(eps(i)))
      rl= 0.d0
      rh=10.d0
      sl= 0.d0
      sh= 0.d0
 300  if (sl*sh.le.0.00000001d0) then
        rc=rl+(rh-rl)/2.d0
      else
        rc=rl+(rh-rl)*(sd-sl)/(sh-sl)
      end if
      sc=0.d0
      do k=1,nr
        tmp=(1.d0-exp(-(r(k)/rc)**2))**2
        tmp=corpol/(2.d0*r(k)**4)*tmp*tmp
        sc=sc+dr(k)*phi(k,i)*phi(k,i)*tmp
      end do
      if (sc.gt.sd) then
        rl=rc
        sl=sc
      else
        rh=rc
        sh=sc
      end if
      if (abs(sc-sd).gt.0.000001d0) go to 300
      write (6,'(1x,3f14.4)') rc,sc,sc*27.2114d0
      return
      end subroutine sigfit2
!----------------------------------------------------------------------
      subroutine convert( unit, value )
      use dprec
      implicit none
      character(len=10), intent(inout) :: unit
      real(dp), intent(inout) :: value
      real(dp) :: mult

      mult = -1.d0
      if ( unit .eq. 'hartree' ) mult = 1.d0
      if ( unit .eq. 'Ryd' ) mult = 0.5d0
      if ( unit .eq. 'eV' ) mult = 1.d0 / 27.2114d0
      if ( unit .eq. 'invcm' ) mult = 0.000123985d0 / 27.2114d0
      if ( mult .lt. 0.d0 ) stop 'unsupported unit'
      value = value * mult
      unit = 'hartree' 
      return
      end subroutine convert
!----------------------------------------------------------------------
      subroutine leadbehv( nlev, llev, nr, ldim, phi, r )
      use dprec
      implicit none
      integer  :: nlev, nr, ldim
      integer  :: llev( nlev )
      real(dp) :: phi( ldim, nlev ), r( nr )
      integer  :: i, j, pow
      real(dp) :: pref, rtar, term, pi
      parameter (pi = 3.1415926535897932384d0)
      
      do i=1,nlev
        write (6,'(2x,2i5)') i,llev(i)
        pow = llev( i ) + 1
        pref = sqrt( real( 2 * llev( i ) + 1,dp ) / ( 4.d0 * pi ) )
        rtar = 0.d0
        j = 1
        do while ( r( j ) .lt. 0.1d0 )
          if ( r( j ) .ge. rtar ) then
            term = pref * phi( j, i ) / r( j ) ** pow
            write ( 6, '(2x,2f20.10)' ) r( j ), term
            rtar = rtar + 0.01d0
          end if
          j = j + 1
        end do
      end do
      return
      end subroutine leadbehv
!----------------------------------------------------------------------
      subroutine mkvctab( nr, ldim, vctab, r, rsqd )
      use dprec
      use params
      implicit none
      integer :: nr, ldim
      integer :: k
      real(dp) :: vctab( ldim, 0 : 2 ), r( nr ), rsqd( nr )
      real(dp) :: corpol, rs, rp, rd
      real(dp) :: pref, rssqd, rpsqd, rdsqd, fssqd, fpsqd, fdsqd

      read (inpunit,*) corpol,rs,rp,rd
      rssqd = rs ** 2
      rpsqd = rp ** 2
      rdsqd = rd ** 2
      do k=1,nr
        pref = - 0.5d0 * corpol / r( k ) ** 4
        fssqd=(1.d0-exp(-rsqd(k)/rssqd))**4
        fpsqd=(1.d0-exp(-rsqd(k)/rpsqd))**4
        fdsqd=(1.d0-exp(-rsqd(k)/rdsqd))**4
        vctab(k,0)= pref * fssqd
        vctab(k,1)= pref * fpsqd
        vctab(k,2)= pref * fdsqd
      end do
      return
      end subroutine mkvctab
!----------------------------------------------------------------------
      subroutine bachelet(vi,r,njrc,nr)
      use dprec
      use PARAMS
      implicit none
      integer njrc(4),i,ii,nr,nop,k,l
      real(dp) vi(nrmax,7),vcore(nrmax),r(nrmax),errfunc
      real(dp) c(6),aa(6),q(6,6),s(6,6),a(3)
      real(dp) zv,a1,a2,c1,c2,ra1,ra2,pi,rpi,rr
!
      write (6,*) 'PLEASE ENTER ZV,A1,A2,C1,C2 FOR CORE POTENTIAL.'
      read (inpunit,*) zv,a1,a2,c1,c2
      ra1=sqrt(a1)
      ra2=sqrt(a2)
      do 25 i=1,nr
      vcore (i)=-zv/r(i)*(c1*errfunc(ra1*r(i)) &
                        +c2*errfunc(ra2*r(i)))
 25   continue
      write (6,*) 'YOU MAY NOW ENTER WHAT YOU WANT TO ENTER.'
      write (6,*) '1=S,2=PSO,3=PAV,4=DSO,5=DAV,6=FSO,7=FAV,0=STOP.'
 35   read (inpunit,*) nop
      if (nop.eq.0) return
      write (6,*) 'PLEASE ENTER A1,A2,A3,C1,C2,C3,C4,C5,C6.'
      read (inpunit,*) a(1),a(2),a(3),c(1),c(2),c(3),c(4),c(5),c(6)
      pi=4.d0*atan(1.d0)
      rpi=sqrt(pi)
      do 42 i=1,3
        do 40 k=1,3
          s(i  ,k  )=0.2500d0*rpi/(a(i)+a(k))**1.5d0
          s(i+3,k  )=0.3750d0*rpi/(a(i)+a(k))**2.5d0
          s(i  ,k+3)=0.3750d0*rpi/(a(i)+a(k))**2.5d0
          s(i+3,k+3)=0.9375d0*rpi/(a(i)+a(k))**3.5d0
 40     continue
 42   continue
      do 100  l=1  ,6
        do 45   i=l+1,6
          q(i,l)=0.d0
 45     continue
        do 75 i=1,l-1
          q(i,l)=s(i,l)
          do 55 k=1,i-1
            q(i,l)=q(i,l)-q(k,i)*q(k,l)
 55       continue
          q(i,l)=q(i,l)/q(i,i)
 75     continue
        q(l,l)=s(l,l)
        do 85 k=1,l-1
          q(l,l)=q(l,l)-q(k,i)*q(k,i)
 85     continue
        q(l,l)=sqrt(q(l,l))
 100  continue
      call patt(q)
      do 120 i=1,6
        aa(i)=0.d0
        do 110 l=1,6
          aa(i)=aa(i)+q(i,l)*c(l)
 110    continue
 120  continue
      if (2*(nop/2).eq.nop) then
        do 235 i=1,nr
          rr=r(i)*r(i)
          vi(i,nop)=-(aa(1)+rr*aa(4))*exp(-a(1)*rr)&
                   -(aa(2)+rr*aa(5))*exp(-a(2)*rr)&
                   -(aa(3)+rr*aa(6))*exp(-a(3)*rr)
 235    continue
      else
        ii=(nop+1)/2
        njrc(ii)=1
        do 240 i=1,nr
          rr=r(i)*r(i)
          vi(i,nop)=-(aa(1)+rr*aa(4))*exp(-a(1)*rr)&
                   -(aa(2)+rr*aa(5))*exp(-a(2)*rr)&
                   -(aa(3)+rr*aa(6))*exp(-a(3)*rr)&
                   +vcore(i)
 240    continue
      endif
      goto 35
      end subroutine bachelet
!----------------------------------------------------------------------
      SUBROUTINE PATT(Q)
      use dprec
      implicit none
      integer :: ichk,i,j,isum,k,l,m
      real(dp) :: Q(6,6),QI(6,6),DELTA(6,6)

      ICHK=0
      DO 100 I=1,6
      DO 100 J=1,6
      QI(I,J)=0.D0
  100 CONTINUE

      QI(1,1)=1.D0/Q(1,1)
      QI(1,2)=-Q(1,2)/(Q(1,1)*Q(2,2))
      QI(1,3)=(Q(1,2)*Q(2,3)-Q(1,3)*Q(2,2))/(Q(1,1)*Q(2,2)*Q(3,3))
      QI(1,4)=(-Q(1,2)*Q(2,3)*Q(3,4)+Q(1,2)*Q(2,4)*Q(3,3)+&
      Q(1,3)*Q(2,2)*Q(3,4)-Q(1,4)*Q(2,2)*Q(3,3))/&
      (Q(1,1)*Q(2,2)*Q(3,3)*Q(4,4))
      QI(1,5)=(Q(1,2)*Q(2,3)*Q(3,4)*Q(4,5)-Q(1,2)*Q(2,3)*Q(3,5)*Q(4,4)-&
      Q(1,2)*Q(2,4)*Q(3,3)*Q(4,5)+Q(1,2)*Q(2,5)*Q(3,3)*Q(4,4)-&
      Q(1,3)*Q(2,2)*Q(3,4)*Q(4,5)+Q(1,3)*Q(2,2)*Q(3,5)*Q(4,4)+&
      Q(1,4)*Q(2,2)*Q(3,3)*Q(4,5)-Q(1,5)*Q(2,2)*Q(3,3)*Q(4,4))/&
      (Q(1,1)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5))
      QI(1,6)=(-Q(1,2)*Q(2,3)*Q(3,4)*Q(4,5)*Q(5,6)+&
      Q(1,2)*Q(2,3)*Q(3,4)*Q(4,6)*Q(5,5)+&
      Q(1,2)*Q(2,3)*Q(3,5)*Q(4,4)*Q(5,6)-&
      Q(1,2)*Q(2,3)*Q(3,6)*Q(4,4)*Q(5,5)+&
      Q(1,2)*Q(2,4)*Q(3,3)*Q(4,5)*Q(5,6)-&
      Q(1,2)*Q(2,4)*Q(3,3)*Q(4,6)*Q(5,5)-&
      Q(1,2)*Q(2,5)*Q(3,3)*Q(4,4)*Q(5,6)+&
      Q(1,2)*Q(2,6)*Q(3,3)*Q(4,4)*Q(5,5)+&
      Q(1,3)*Q(2,2)*Q(3,4)*Q(4,5)*Q(5,6)-&
      Q(1,3)*Q(2,2)*Q(3,4)*Q(4,6)*Q(5,5)-&
      Q(1,3)*Q(2,2)*Q(3,5)*Q(4,4)*Q(5,6)+&
      Q(1,3)*Q(2,2)*Q(3,6)*Q(4,4)*Q(5,5)-&
      Q(1,4)*Q(2,2)*Q(3,3)*Q(4,5)*Q(5,6)+&
      Q(1,4)*Q(2,2)*Q(3,3)*Q(4,6)*Q(5,5)+&
      Q(1,5)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,6)-&
      Q(1,6)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5))/&
      (Q(1,1)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5)*Q(6,6))
      QI(2,2)=1.D0/Q(2,2)
      QI(2,3)=-Q(2,3)/(Q(2,2)*Q(3,3))
      QI(2,4)=(Q(2,3)*Q(3,4)-Q(2,4)*Q(3,3))/(Q(2,2)*Q(3,3)*Q(4,4))
      QI(2,5)=(-Q(2,3)*Q(3,4)*Q(4,5)+Q(2,3)*Q(3,5)*Q(4,4)+&
      Q(2,4)*Q(3,3)*Q(4,5)-Q(2,5)*Q(3,3)*Q(4,4))/&
      (Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5))
      QI(2,6)=(Q(2,3)*Q(3,4)*Q(4,5)*Q(5,6)-Q(2,3)*Q(3,4)*Q(4,6)*Q(5,5)-&
      Q(2,3)*Q(3,5)*Q(4,4)*Q(5,6)+Q(2,3)*Q(3,6)*Q(4,4)*Q(5,5)-&
      Q(2,4)*Q(3,3)*Q(4,5)*Q(5,6)+Q(2,4)*Q(3,3)*Q(4,6)*Q(5,5)+&
      Q(2,5)*Q(3,3)*Q(4,4)*Q(5,6)-Q(2,6)*Q(3,3)*Q(4,4)*Q(5,5))/&
      (Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5)*Q(6,6))
      QI(3,3)=1.D0/Q(3,3)
      QI(3,4)=-Q(3,4)/(Q(3,3)*Q(4,4))
      QI(3,5)=(Q(3,4)*Q(4,5)-Q(3,5)*Q(4,4))/(Q(3,3)*Q(4,4)*Q(5,5))
      QI(3,6)=(-Q(3,4)*Q(4,5)*Q(5,6)+Q(3,4)*Q(4,6)*Q(5,5)+&
      Q(3,5)*Q(4,4)*Q(5,6)-Q(3,6)*Q(4,4)*Q(5,5))/&
      (Q(3,3)*Q(4,4)*Q(5,5)*Q(6,6))
      QI(4,4)=1.D0/Q(4,4)
      QI(4,5)=-Q(4,5)/(Q(4,4)*Q(5,5))
      QI(4,6)=(Q(4,5)*Q(5,6)-Q(4,6)*Q(5,5))/(Q(4,4)*Q(5,5)*Q(6,6))
      QI(5,5)=1.D0/Q(5,5)
      QI(5,6)=-Q(5,6)/(Q(5,5)*Q(6,6))
      QI(6,6)=1.D0/Q(6,6)
! CHECK INVERSE
      WRITE(9,1000)
 1000 FORMAT(/4X,'QUALITY OF ANALYTIC INVERSION BY PATTNAIK ET AL'/)
      DO 150 I=1,6
      DO 200 J=1,6
      DELTA(I,J)=0.D0
      IF (I.EQ.J) DELTA(I,J)=-1.D0

      DO 300 ISUM=1,6
      DELTA(I,J)=DELTA(I,J)+QI(I,ISUM)*Q(ISUM,J)
  300 CONTINUE

  200 CONTINUE

      WRITE(9,2000) (DELTA(I,K),K=1,6)
 2000 FORMAT(4X,6D18.7)

  150 CONTINUE

! TRANSFER INVERSE MATRIX ONTO ORIGINAL MATRIX
      DO 500 L=1,6
      DO 500 M=1,6
      Q(L,M)=QI(L,M)
  500 CONTINUE

      RETURN
      END subroutine patt
!----------------------------------------------------------------------
      subroutine chops(nr,r,dl,rpot,rtot,etot,orb,p,q)
      use dprec
      use PARAMS
      implicit none
      integer nr,istart,ii
      real(dp) r(nr),dl,rpot,rtot,etot,dlr,sum
      real(dp) orb(nr),p(nr),q(nr),addin(nrmax)
!
      dlr=dl/45.d0
      do istart=1,4
        ii=istart
        do while (ii+4.le.nr)
          ii=ii+4
        end do
        ii=ii-4
        sum=0.d0
        do while (ii.ge.istart)
          sum=sum+dlr*(14.d0*(p(ii+0)*p(ii+0)+p(ii+4)*p(ii+4))+&
                      64.d0*(p(ii+1)*p(ii+1)+p(ii+3)*p(ii+3))+&
                      24.d0*(p(ii+2)*p(ii+2)                ) )
          addin(ii)=sum
          ii=ii-4
        end do
      end do
      do istart=5,8
        sum=0.d0
        do ii=istart,nr,4
          sum=sum+dlr*(14.d0*(p(ii-0)*p(ii-0)*r(ii-0)+&
                             p(ii-4)*p(ii-4)*r(ii-4) )+&
                      64.d0*(p(ii-1)*p(ii-1)*r(ii-1)+&
                             p(ii-3)*p(ii-3)*r(ii-3) )+&
                      24.d0*(p(ii-2)*p(ii-2)*r(ii-2) ) )
          addin(ii)=addin(ii)+sum/r(ii)
        end do
      end do
      do istart=1,4
        do ii=istart,nr-4,4
          orb(ii)=orb(ii)+rpot*addin(ii)
          if (istart.eq.1) then
          etot=etot+rtot*dlr/2.d0*(&
           +addin(ii+0)*q(ii+0)*q(ii+0)*14.d0*r(ii+0)& 
           +addin(ii+1)*q(ii+1)*q(ii+1)*64.d0*r(ii+1)& 
           +addin(ii+2)*q(ii+2)*q(ii+2)*24.d0*r(ii+2)& 
           +addin(ii+3)*q(ii+3)*q(ii+3)*64.d0*r(ii+3)& 
           +addin(ii+4)*q(ii+4)*q(ii+4)*14.d0*r(ii+4))
          end if
        end do
      end do
      do ii=nr-9,nr
        orb(ii)=orb(nr-10)*r(nr-10)/r(ii)
      end do
      return
      end subroutine chops
!----------------------------------------------------------------------
      subroutine corepot( n, lead, vi )
      use dprec
      implicit none
      integer :: n, lead
      real(dp) :: vi( lead, 7 )
      integer :: i, j
      real(dp) :: dum, vadd
!
      open( unit=99, file='Real.Input',&
           form='formatted', status='unknown' )
      rewind 99
      do i = 1, n
        read ( 99, * ) dum, vadd
        do j = 1, 7
          vi( i, j ) = vi( i, j ) + vadd
        end do
      end do
      close( unit=99 )
      return
      end subroutine corepot
!     program erftest
!     implicit none
!     real(dp)   errfunc
!     external errfunc
!10   read (5,*) x
!     if (x.lt.0.d0) x=sqrt(-x*4.d0*atan(1.d0))
!     write (6,22) errfunc(x),1.d0-errfunc(x)
!22   format(1x,2f30.20)
!     goto 10
!     end
!----------------------------------------------------------------------
      function errfunc(ss)
      use dprec
      implicit none
      integer :: isi, n, k
      real(dp) :: errfunc, ss, s2, sn, pih, temp
!
      if (ss.gt. 8.d0) then
        errfunc=1.d0
        return
      endif
      pih=sqrt(4.d0*atan(1.d0))
      s2=ss*ss
      temp=0.d0
      isi=1
      if (ss.lt.3.65d0) then
        n=1
        k=0
        sn=ss
 100    temp=temp+sn/real(isi*n,dp)
        n=n+2
        k=k+1
        sn=sn*s2/real(k,dp)
        isi=-isi
        if ((n-10).lt.int(100.d0*s2)) goto 100
        temp=temp*2.d0/pih
      else
        n=0
        sn=2.d0*ss
 200    temp=temp+real(isi,dp)/sn
        n=n+1
        sn=-sn*2.d0*s2/real(2*n-1,dp)
        if ((2*n-1).lt.2.d0*s2) goto 200
        temp=1.d0-temp*2.d0/pih*exp(-s2)
      endif
      errfunc=temp
      return
      end function errfunc
!----------------------------------------------------------------------
      subroutine elener(l,xkappa,xj,zeff,e,phi,v,&
                       xm1,xm2,nr,r,r2,dl,rel,isuse)
      use dprec
      use PARAMS
      implicit none
      real(dp) xkappa,xj,zeff,e,dl,rel,plead
      integer l,nr
      real(dp) phi(nrmax),v(nrmax),xm1(nrmax),xm2(nrmax)
      real(dp) r(nrmax),r2(nrmax)
      integer isuse
      real(dp) xnorm,aa,xo
      integer j, istop, ief, nn

      call getplead(l,xj,rel,xkappa,plead,zeff)
      isuse = nr - 10
      do j = 1, nr
        if ( r( j ) .le. 5.d0 ) isuse = j
      end do

      istop=isuse
      call intego(e,l,xkappa,1000,nn,istop,ief,xo,phi,zeff,v,xm1,&
                 xm2,nr,r,r2,dl,rel,plead)
      do j = istop + 1, nr
        phi(j)=0.d0
      end do

      if (abs(abs(xj)-real(l,dp)).gt.0.25d0)&
       call augment(e,l,xj,phi,v,nr,r,dl,rel)
      do j = nr - 5, nr
        phi( j ) = 0.d0
      end do

      aa=0.d0
      do j=1,nr-4,4
        aa=aa+phi(j+0)*phi(j+0)*dl*r(j+0)*14.d0/45.d0&
            +phi(j+1)*phi(j+1)*dl*r(j+1)*64.d0/45.d0&
            +phi(j+2)*phi(j+2)*dl*r(j+2)*24.d0/45.d0&
            +phi(j+3)*phi(j+3)*dl*r(j+3)*64.d0/45.d0&
            +phi(j+4)*phi(j+4)*dl*r(j+4)*14.d0/45.d0
      end do
      xnorm=1.d0/sqrt(aa)
      do j=1,nr
        phi(j)=phi(j)*xnorm
      end do
      return
      end subroutine elener
!----------------------------------------------------------------------
      subroutine getpot(etot,rel,alfa,dl,nr,dr,r,r2,phe,&
                       ratio,orb,occ,is,nel,nl,nm,no,xnj,xnum,&
                       iuflag,cq,vold,vnew)
      use dprec
      use PARAMS
      implicit none
      integer  :: nr,nel,iuflag
      real(dp) :: etot,rel,alfa,dl,ratio,xnum
      real(dp) :: fa

      integer  :: is(iorbs),nl(iorbs),nm(iorbs),no(iorbs)
      real(dp) :: dr(nrmax),r(nrmax),r2(nrmax)
      real(dp) :: phe(nrmax,iorbs),orb(nrmax,iorbs)
      real(dp) :: occ(iorbs),xnj(iorbs)
      real(dp) :: cq(nrmax)
      real(dp) :: valtab( nrmax ), ualtab( nrmax )
      real(dp) :: cg(0:6,0:6,0:12,-6:6,-6:6),pin(0:8,0:8,0:16)
      real(dp) :: vold(iorbs),vnew(iorbs)

      logical  :: df, spherical
      integer  :: bcount
      real(dp) :: bwgt

      real(dp) :: sum, chk, tmp, exsum
      real(dp) :: etkin,etnuc,etcou,etlda,zvalue
      common /parts/etkin,etnuc,etcou,etlda,zvalue
      save /parts/

      integer  :: i,j,k
      integer  :: li,mi,lj,mj,la,lmn,lmx,jstart
      real(dp) :: ratcom,ri,rj,rc,coeff,ccg,col
      real(dp) :: etemp,xnum2,etni

      call clebschgordan(nel,iorbs,nl,cg)
      call getillls(pin)

      df = (  abs(alfa)  .gt.  0.00000000009d0  )
      ratcom=1.d0-ratio
      do i=1,nel
        vold(i)=0.d0
        do k=1,nr
          vold(i)=vold(i)+orb(k,i)*phe(k,i)*phe(k,i)*dr(k)
          orb(k,i)=ratcom*orb(k,i)
        end do
      end do

!  zero out parts of total energy which are to be computed
      etcou=0.d0
      etnuc=0.d0
      etlda=0.d0

!  here we do the hartree term
!==============================
      do i=1,nel
      etni=0.d0
!==============================

      li=nl (i)
      mi=nm (i)

!  part of electron-nucleus term ... 

      fa=zvalue*occ(i)*dl/45.d0
      etni=0.d0
      do j=1,nr-4,4
        etnuc=etnuc-fa*&
         ((phe(j+0,i)*phe(j+0,i)+phe(j+4,i)*phe(j+4,i))*14.d0&
         +(phe(j+1,i)*phe(j+1,i)+phe(j+3,i)*phe(j+3,i))*64.d0&
         +(phe(j+2,i)*phe(j+2,i)                      )*24.d0)
        etni=etni+dl/45.d0*&
         ((phe(j+0,i)*phe(j+0,i)+phe(j+4,i)*phe(j+4,i))*14.d0&
         +(phe(j+1,i)*phe(j+1,i)+phe(j+3,i)*phe(j+3,i))*64.d0&
         +(phe(j+2,i)*phe(j+2,i)                      )*24.d0)
      end do

!==============================
      jstart=i+1
      if ((xnj(i).lt.0.d0).or.(occ(i).gt.1.d0).or.df) jstart=i
      do j=jstart,nel
!==============================

      lj=nl (j)
      mj=nm (j)

!========================
      if (abs(occ(i))>1.d-12.or.abs(occ(j))>1.d-12) then  !JR
!========================
      spherical=((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.&
                (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0).or.df)
      spherical=((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.&
                (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0))
                           lmx=2*nl(i)
      if (nl(i).gt.nl(j))  lmx=2*nl(j)
      if (spherical) lmx=0
!     lmx=0
!============
      do la=lmx,0,-2
!============
      coeff=real((li+li+1)*(lj+lj+1),dp)/real((la+la+1)**2,dp)*&
         cg(li,li,la,mi,-mi)*cg(lj,lj,la,mj,-mj)*&
         cg(li,li,la,0 , 0 )*cg(lj,lj,la,0 , 0 )
      if (mi+mj.ne.2*((mi+mj)/2)) coeff=-coeff
      call getrs(coeff,i,j,occ(i),occ(j),ratio,ri,rj,rc)
      call mkvaltab( nr, r, dl, phe( 1, i ), phe( 1, i ),& 
                    valtab, la )
      call mkvaltab( nr, r, dl, phe( 1, j ), phe( 1, j ),& 
                    ualtab, la )
      do k=1,nr
        orb(k,j)=orb(k,j)+valtab( k )*ri
        orb(k,i)=orb(k,i)+ualtab( k )*rj
      end do
      etemp=etot
      bcount = 0
      bwgt = 14.d0 / 45.d0
      do k = 1, nr
        bwgt = bwgt * rc * dl * r( k ) * 0.5d0
        etot = etot + bwgt * ualtab( k ) *& 
                     phe( k, i ) * phe( k, i ) +&
                     bwgt * valtab( k ) *&
                     phe( k, j ) * phe( k, j )
        bcount = bcount + 1
        if ( bcount .eq. 4 ) then 
          bcount = 0
          bwgt = 28.d0 / 45.d0
        else
          bwgt = 64.d0 / 45.d0
          if ( bcount .eq. 2 ) bwgt = 24.d0 / 45.d0
        end if
      end do
      etcou=etcou+(etot-etemp)
!============
      end do
!============
!========================
      end if
!========================
!==============================
      end do
!==============================
!==============================
      end do
!==============================

!     if ( ratio .gt. 0.99d0 ) then !JR
      if ( ratio .gt. 1.d+33 ) then !JR
      open(unit=66,file='radpot',form='formatted',status='unknown')
      rewind 66
      chk = 0.d0
      sum = 0.d0
      do k = 1, nr
        tmp = 0.d0
        do i = 1, nel
          tmp = tmp + phe( k, i ) * phe( k, i ) * occ( i )
        end do
        sum = sum + tmp * dr( k ) / r( k )
        chk = chk + tmp * dr( k )
      end do
      chk = 0.d0
      do k = 1, nr
        tmp = 0.d0
        do i = 1, nel
          tmp = tmp + phe( k, i ) * phe( k, i ) * occ( i )
        end do
        chk = chk + tmp * dr( k ) * 0.5
        sum = sum - tmp * dr( k ) / r( k ) * 0.5
        write ( 66, '(3f15.10)' ) r( k ), sum + chk / r( k )
        chk = chk + tmp * dr( k ) * 0.5
        sum = sum - tmp * dr( k ) / r( k ) * 0.5
      end do
      close(unit=66)
      end if

!     if ( ratio .gt. 0.99d0 ) then !JR
      if ( ratio .gt. 1.d+33 ) then !JR
        open(unit=66,file='hapot',form='formatted',status='unknown')
        rewind 66
        do i = 1, nel
          do k = 1, nr
            write ( 66, '(2x,4i5,2(2x,1e20.12))' ) &
           i, nel, k, nr, r( k ), orb( k, i )
          end do
        end do
        close( unit=66 )
      end if

!  here we do exchange and correlation

!===========================================================
      if (df) then
!===========================================================
        etemp=etot
        call dft(dl,rel,alfa,nr,nrmax,nel,nl,xnj,is,occ,r2,&
                cq,phe,orb,ratio,etot,x137)
        etlda=etlda+(etot-etemp)
!===========================================================
      else
!===========================================================
        xnum2=xnum*xnum
!==============================
        do i=1,nel
!==============================
          li=nl (i)
          mi=nm (i)
                                                    jstart=i+1
          if ((xnj(i).lt.0.d0).or.(occ(i).gt.1.d0)) jstart=i
!==============================
          do j=jstart,nel
!==============================
            lj=nl (j)
            mj=nm (j)
            if (abs(occ(i))>1.d-12.or.abs(occ(j))>1.d-12) then  !JR
              spherical= .not. .true.
              if (occ(i).gt.1.d0) spherical=.true.
              if (occ(j).gt.1.d0) spherical=.true.
              if (xnj(i).lt.0.d0) spherical=.true.
              if (xnj(j).lt.0.d0) spherical=.true.
              if ((is(i).eq.is(j)).or.spherical) then
                lmx=li+lj
                lmn=iabs(mi-mj)
                if (spherical) lmn=0
                do la=lmx,lmn,-2
                  if (spherical) then
                    coeff=pin(li,lj,la)/4.d0
                  else
                  col=real((li+li+1)*(lj+lj+1),dp)/real((la+la+1)**2,dp)
                    ccg=cg(li,lj,la,-mi,mj)*cg(li,lj,la,0,0)
                    coeff=col*ccg*ccg
                  end if
                  call getrs(coeff,i,j,occ(i),occ(j),ratio,ri,rj,rc)
                  call mkvaltab( nr, r, dl, phe( 1, i ), phe( 1, j ),& 
                                valtab, la )
                  exsum = 0.d0
                  bcount = 0
                  bwgt = 14.d0 / 45.d0
                  do k = 1, nr
                    etot = etot -& 
                      rc * valtab( k ) * dl * r( k ) * bwgt *&
                      phe( k, i ) * phe( k, j )
                    exsum = exsum -& 
                      rc * valtab( k ) * dl * r( k ) * bwgt *&
                      phe( k, i ) * phe( k, j )
                    bcount = bcount + 1
                    if ( bcount .eq. 4 ) then 
                      bcount = 0
                      bwgt = 28.d0 / 45.d0
                    else
                      bwgt = 64.d0 / 45.d0
                      if ( bcount .eq. 2 ) bwgt = 24.d0 / 45.d0
                    end if
                  end do
!CD                write ( 6, * ) i, j, exsum / ( occ( i ) * occ( j ) )
                  call fockin( nr, orb( 1, i ), valtab,&
                              phe( 1, i ), phe( 1, j ), xnum, rj )
                  call fockin( nr, orb( 1, j ), valtab,&
                              phe( 1, j ), phe( 1, i ), xnum, ri )
                end do
              end if
            end if
!==============================
          end do
!==============================
!==============================
        end do
!==============================
!===========================================================
      end if
!===========================================================

!     if ( ratio .gt. 0.99d0 ) then !JR
      if ( ratio .gt. 1.d+33 ) then !JR
        open(unit=66,file='hfpot',form='formatted',status='unknown')
        rewind 66
        do i = 1, nel
          do k = 1, nr
            write ( 66, '(2x,4i5,2(2x,1e20.12))' )& 
           i, nel, k, nr, r( k ), orb( k, i )
          end do
        end do
        close( unit=66 )
      end if

!  where we might do restriction

      if (iuflag.ne.0) call rest(nrmax,nr,nel,no,nl,is,iuflag,occ,orb)

!  figure total kinetic energy

      etkin=etot-(etlda+etcou+etnuc)

!  figure out the first-order ptbn. th'y effects on new
!  eigenvalue for an orbital

      do i=1,nel
        vnew(i)=0.d0
        do k=1,nr
          vnew(i)=vnew(i)+orb(k,i)*phe(k,i)*phe(k,i)*dr(k)
        end do
      end do

      return
      end subroutine getpot
!----------------------------------------------------------------------
      subroutine getrs(coeff,i,j,oi,oj,ratio,ri,rj,rc)
      use dprec
      implicit none
      real(dp) coeff,oi,oj,ratio,ri,rj,rc
      integer i,j
!
      if (i.eq.j) coeff=coeff*0.5
      ri=coeff*oi*ratio
      rj=coeff*oj*ratio
      rc=coeff*oi*oj
      return
      end subroutine getrs
!----------------------------------------------------------------------
      subroutine dft(dl,rel,al,nr,mr,ne,l,j,s,o,r2,cq,ph,or,&
                    ra,et,x137)
      use dprec
      implicit none
      integer nr,mr,ne,l(ne),s(ne),i,k,ii
      real(dp) j(ne),o(ne),r2(mr),cq(mr)
      real(dp) ph(mr,ne),or(mr,ne)
      real(dp) rel,al,ra,et,dl,den,occ,pr,xn,fx,fc,bfac
      real(dp) xn1,ux1,uc1,uxc1,xn2,ux2,uc2,uxc2,nex,ec,x137

      pr=0.0001d0
      if (al.gt.0.d0) then
        fx=1.0d0
        fc=1.0d0
      else
        fx=1.5d0*abs(al)
        fc=0.0d0
      end if

      ii=0
      do i=1,nr
        xn =0.d0
        xn1=0.d0
        xn2=0.d0
!       xn1=0.00000000001d0
!       xn2=0.00000000001d0
        do k=1,ne
          occ=o(k)
          den=occ*ph(i,k)*ph(i,k)
          if ((j(k).lt.-pr).or.(occ.gt.real(2*l(k)+1,dp)+pr)) then
            xn=xn+den
          else
            if (s(k).eq.1) then
              xn1=xn1+den
            else
              xn2=xn2+den
            end if
          end if
        end do
        xn=xn+cq(i)
        xn1=xn1+0.5d0*xn
        xn2=xn2+0.5d0*xn
        if ((xn1+xn2)/r2(i).lt.1.d-30) then
          nex=0.d0
          ec=0.d0
          ux1=0.d0 
          ux2=0.d0 
          uc1=0.d0 
          uc2=0.d0 
        else
          call exchcorr(rel,r2(i),xn1,xn2,nex,ec,ux1,ux2,uc1,uc2,x137)
        end if

        if (ii.eq.0) bfac=14.d0/45.d0
        if ((ii.eq.1).or.(ii.eq.3)) bfac=64.d0/45.d0
        if (ii.eq.2) bfac=24.d0/45.d0
        if (ii.eq.4) bfac=28.d0/45.d0
        et=et+dl*sqrt(r2(i))*bfac*(fc*ec*(xn1+xn2)+fx*nex)
        uxc1 = ra * ( fx*ux1 + fc*uc1 )
        uxc2 = ra * ( fx*ux2 + fc*uc2 )

        do k=1,ne
          occ=o(k)
          if ((j(k).lt.-pr).or.(occ.gt.real(2*l(k)+1,dp)+pr)) then
            or(i,k)=or(i,k)+0.5d0*(uxc1+uxc2)
          else
            if (s(k).eq.1) then
              or(i,k)=or(i,k)+uxc1
            end if
            if (s(k).eq.2) then
              or(i,k)=or(i,k)+uxc2
            end if
          end if
        end do

        ii=ii+1
        if (ii.eq.5) ii=1
      end do

      return
      end subroutine dft
!----------------------------------------------------------------------
      subroutine rest(nrmax,nr,nel,no,nl,is,iu,occ,orb)
      use dprec
      implicit none
      integer nr,nel,no(nel),nl(nel),is(nel),nrmax,iu,i,ii,jj,icond,k
      real(dp) occ(nel),orb(nrmax,nel),orba,div
      logical nsame,lsame,ssame
!
      jj=1
 8960 ii=jj
 8965 if (ii.lt.nel) then
        nsame=(no(jj).eq.no(ii+1))
        lsame=(nl(jj).eq.nl(ii+1))
        ssame=(is(jj).eq.is(ii+1))
                                                     icond=0
        if ((iu.eq.2).and.nsame.and.lsame          ) icond=1
        if ((iu.eq.1).and.nsame.and.lsame.and.ssame) icond=1
        if (icond.eq.1) then
          ii=ii+1
          go to 8965
        end if
      end if
      div=0.d0
      do k=jj,ii
        div=div+occ(k)
      end do
      if (div.gt.0.000001d0) then
        div=1.d0/div
        do i=1,nr
          !JR: ifort compiler fails with this piece of code.
          !orba=0.d0
          !do k=jj,ii
          !  orba=orba+orb(i,k)*occ(k)
          !end do
          !orba=orba*div
          !do k=jj,ii
          !  orb(i,k)=orba
          !end do
          orba=sum(orb(i,jj:ii)*occ(jj:ii))*div
          orb(i,jj:ii)=orba
        end do
      end if
      if (ii.ne.nel) then
        jj=ii+1
        go to 8960
      end if
      return
      end subroutine rest
!-----------------------------------------------------------------------
!  subr gets Clebsch-Gordan coefficients, in the form of 
!  cg(l1,l2,L,m1,m2) = <l1,m1;l2,m2|L,m1+m2>, according to Rose's 
!  'Elementary Theory of Angular Momentum', p. 39, Wigner's formula.
!  those coefficients listed are only those for which l1.ge.l2.
!  coefficients known to be zero because of either the L or M 
!  selection rules are not computed, and should not be sought.

      subroutine clebschgordan(nel,nelmx,nl,cg)
      use dprec
      implicit none
      integer nel,nelmx,nl(nelmx)
      integer i,lmx,l1,l2,l3,m1,m2,m3,lmin,numin,numax,nu

      real(dp) cg(0:6,0:6,0:12,-6:6,-6:6),si(0:32),fa(0:32)
      real(dp) pref,sum

      lmx=0
      do i=1,nel
        if (nl(i).gt.lmx) lmx=nl(i)
      end do
      si(0)=1.d0
      fa(0)=1.d0
      do i=1,32
        si(i)=-si(i-1)
        fa(i)=real(i,dp)*fa(i-1)
      end do

!===============================
      do l1=0,lmx
        do l2=0,l1
          do m1=-l1,l1
            do m2=-l2,l2
!===============================

      m3=m1+m2
                            lmin=iabs(l1-l2)
      if (lmin.lt.iabs(m3)) lmin=iabs(m3)

      do l3=lmin,l1+l2
        pref=real(2*l3+1,dp)
        pref=pref*fa(l3+l1-l2)/fa(l1+l2+l3+1)
        pref=pref*fa(l3-l1+l2)/fa(l1-m1)
        pref=pref*fa(l1+l2-l3)/fa(l1+m1)
        pref=pref*fa(l3+m3)/fa(l2-m2)
        pref=pref*fa(l3-m3)/fa(l2+m2)
        pref=sqrt(pref)
        sum=0.d0
                              numax=l3-l1+l2
        if ((l3+m3).lt.numax) numax=l3+m3
                               numin=0
        if (l1-l2-m3.lt.numin) numin=-(l1-l2-m3)
        do nu=numin,numax
          sum=sum+&
          (si(nu+l2+m2)/fa(nu))*fa(l2+l3+m1-nu)*fa(l1-m1+nu)&
          /fa(l3-l1+l2-nu)/fa(l3+m3-nu)/fa(nu+l1-l2-m3)
        end do
        cg(l1,l2,l3,m1,m2)=pref*sum
        cg(l2,l1,l3,m2,m1)=si(l1+l2+l3)*pref*sum
      end do

!===============================
            end do
          end do
        end do
      end do
!===============================
      return
      end subroutine clebschgordan
!----------------------------------------------------------------------
      subroutine fockin( nr, pot, val, ph1, ph2, num, fact )
      use dprec
      implicit none
      integer nr, k, kk
      logical past
      real(dp) pot( nr ), val( nr ), ph1( nr ), ph2( nr )
      real(dp) num, fact, xx
!
      past = .false.
      k = 2
      do while ( .not. past )
        past = ( abs( ph1( k ) / ph1( k - 1 ) ) .lt. 1.d0 )
        k = k + 1
        if ( k .gt. nr ) stop ' bad fockin '
      end do
      do kk = 1, k
        if (abs(ph1( kk ) * ph2( kk ))>1.d-12) then  !JR
          xx= ph2( kk ) / ph1( kk )
          pot( kk ) = pot( kk ) - val( kk ) * fact * xx
        end if
      end do
      do kk = k + 1, nr
        if (abs(ph1( kk ) * ph2( kk ))>1.d-12) then  !JR
          xx= ph2( kk ) / ph1( kk )
          if ( abs( xx ) .gt. num ) xx = num * num / xx
          pot( kk ) = pot( kk ) - val( kk ) * fact * xx
        end if
      end do
      return
      end subroutine fockin
!----------------------------------------------------------------------
      subroutine grandop(nwgt,wgts,mode)
      use dprec
      implicit none
      integer nwgt,iz,i,j,nel,ibig,ilit
      real(dp) zorig,xz
      integer no(1000),nl(1000),is(1000)
      real(dp) xnj(1000),occ(1000),ev(1000),wgts(nwgt)
      real(dp) other(5),tmp(1000)
      character*1 dumcha,first,second
      character*3 mode,fn3,front
      character*4 lets,fn4
      character*10 digs
      logical four
!
      open(unit=77,file='tmp',form='formatted',status='unknown')
      rewind 77
      do j=1,5
        other(j)=0.d0
      end do
      do i=1,nwgt
        read (77,*) zorig,nel
        iz=int(zorig)  !JR
        read (77,*) (tmp(j),j=1,5)
        do j=1,5
          other(j)=other(j)+wgts(i)*tmp(j)
        end do
        read (77,*) (no(j),nl(j),xnj(j),is(j),occ(j),tmp(j),j=1,nel)
        do j=1,nel
          ev(j)=ev(j)+wgts(i)*tmp(j)
        end do
      end do
      close(unit=77)
      xz=0.d0
      do j=1,nel
        xz=xz+occ(j)
      end do
      ibig=iz/10
      ilit=iz-10*ibig
      ibig=ibig+1
      ilit=ilit+1
      digs='0123456789'
      lets='spdf'
      fn4(1:1)=digs(ibig:ibig)
      fn4(2:2)=digs(ilit:ilit)
      fn3(1:1)=digs(ibig:ibig)
      fn3(2:2)=digs(ilit:ilit)
      open(unit=98,file='occinfo',form='formatted',status='unknown')
      rewind 98
      do j=1,iz-1
        read (98,'(1a1)') dumcha
      end do
      read (98,'(8x,2a1)') first,second
      close(unit=98)
      if (second.eq.' ') then
        four=.false.
      else
        four=.true.
      end if
      fn3(3:3)=first
      fn4(3:3)=first
      fn4(4:4)=second
      if (four) then
        open(unit=99,file=fn4,form='formatted',status='unknown')
        write (6,*) fn4
      else
        open(unit=99,file=fn3,form='formatted',status='unknown')
        write (6,*) fn3
      end if
      rewind 99
      write (99,'(1a7,1f16.6)') 'Etot  =',other(1)
      write (99,'(1a7,1f16.6)') 'Ekin  =',other(2)
      write (99,'(1a7,1f16.6)') 'Ecoul =',other(3)
      write (99,'(1a7,1f16.6)') 'Eenuc =',other(4)
      write (99,'(1a7,1f16.6)') 'Exc   =',other(5)
      do j=1,nel
        front(1:1)=digs(no(j)+1:no(j)+1)
        front(2:2)=lets(nl(j)+1:nl(j)+1)
        if ((mode.eq.'lda').or.(mode.eq.'scr')) then
          front(3:3)=' '
        else
          if (mode.eq.'lsd') then
            if (is(j).eq.1) then
              front(3:3)='D'
            else
              front(3:3)='u'
            end if
          else
            if (abs(xnj(j)).gt.real(nl(j),dp)) then
              front(3:3)='P'
            else
              front(3:3)='M'
            end if
          end if
        end if
        write (99,'(1a3,1f17.6)') front,ev(j)
      end do
      close(unit=99)
      return
      end subroutine grandop
!----------------------------------------------------------------------
      subroutine gropen(icr)
      implicit none
      integer :: icr
      character(len=6) :: s6
      character(len=7) :: s7, u7
      character(len=9) :: f9
      parameter(u7='unknown',f9='formatted')
      call encrypt(icr,s6,s7)
      if (icr.lt.10) open(unit=icr,file=s6,form=f9,status=u7)
      if (icr.ge.10) open(unit=icr,file=s7,form=f9,status=u7)
      rewind icr
      return
      end subroutine gropen
!----------------------------------------------------------------------
      subroutine encrypt(iu,s6,s7)
      implicit none
      integer :: iu
      character(len=6) :: s6
      character(len=7) :: s7
      if ((iu.lt.1) .or.(iu.gt.99)) stop 'bad iu!!!'
      open(unit=99,file='toupee',form='formatted',status='unknown')
      rewind 99
      if (iu.lt.10) write (99,61) 'fort.',iu
      if (iu.ge.10) write (99,71) 'fort.',iu
      rewind 99
      if (iu.lt.10) read  (99,62) s6
      if (iu.ge.10) read  (99,72) s7
      close(unit=99)
      return
   61 format(1x,1a5,1i1)
   62 format(1x,1a6)
   71 format(1x,1a5,1i2)
   72 format(1x,1a7)
      end subroutine encrypt
!------------------------------------------------------------------------
!  exchange correlation routine, by Ceperley-Alder, as parametrized by
!  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
!  between the unpolarized and polarized gas for the correlation part.
      subroutine cacorr(nst,rel,rr,rh1,rh2,ex,ec,ux1,ux2,uc1,uc2)
      use dprec
      use PARAMS
      implicit none
      integer :: nst
      real(dp) :: trd,ft,rh,rh1,rh2,pi,fp,xl,xn,xn1,xn2,rr,ex,ec,ux1,ux2
      real(dp) :: uc1,uc2,rs,zeta,exchfactor,befactor,fe1,fe2,fu1,fu2
      real(dp) :: ex1,ex2,beta,beta1,beta2,b2,eta,rootr,denom,ecu,ucu,ecp
      real(dp) :: ucp,xlr,au,bu,cu,du,ap,bp,cp,tp,rel,f,dfdz,gamma,rlr

      trd=1.d0/3.d0
      ft=4.d0/3.d0
      rh=rh1+rh2

!  if one spin type, average polarization

      if (nst.eq.1) then
        rh1=rh/2.d0
        rh2=rh/2.d0
      endif

!  get the n's, and the rs.

      pi=3.14159265358979d0
      fp=4.d0*pi
      xn1=rh1/(rr*fp)
      xn2=rh2/(rr*fp)
      xn=xn1+xn2

!  effect cutoff, to avoid overflow

      if ((nst.eq.3).or.(xn.lt.0.00000001d0)) then

        ex=0.d0
        ec=0.d0
        ux1=0.d0
        ux2=0.d0
        uc1=0.d0      
        uc2=0.d0

      else
        rs=(3.d0/(fp*xn))**trd
        zeta=(xn1-xn2)/xn
!       exchfactor=-0.930525736d0
        exchfactor=-1.5d0*(0.75d0/pi)**trd

        befactor=(9.d0*pi/4.d0)**trd/x137
        if (abs(xn1)<1.d-12) then  !JR
          fe1=1.d0
          fu1=1.d0
          ex1=0.d0
          ux1=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=sqrt(1.d0+b2)
          xl=log(beta+eta)
          fe1=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu1=-0.5d0+1.5d0*xl/beta/eta
          ex1=exchfactor*xn1**trd
          ux1=4.d0*ex1/3.d0
        endif
        if (abs(xn2)<1.d-12) then
          fe2=1.d0
          fu2=1.d0
          ex2=0.d0
          ux2=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=sqrt(1.d0+b2)
          xl=log(beta+eta)
          fe2=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu2=-0.5d0+1.5d0*xl/beta/eta
          ex2=exchfactor*xn2**trd
          ux2=4.d0*ex2/3.d0
        endif
!  these next lines do the Ceperley-Alder correlation
        if (rs.ge.1.d0) then

          rootr=sqrt(rs)

          gamma=-0.1423d0
          beta1=1.0529d0
          beta2=0.3334d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecu=gamma/denom
          ucu=ecu*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

          gamma=-0.0843d0
          beta1=1.3981d0
          beta2=0.2611d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecp=gamma/denom
          ucp=ecp*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

        else

          xlr=log(rs)
          rlr=rs*xlr

          au= 0.0311d0
          bu=-0.048d0
          cu= 0.002d0
          du=-0.0116d0
          ecu=au*xlr+bu+cu*rlr+du*rs
          ucu=au*xlr+(bu-au/3.d0)+2.d0/3.d0*cu*rlr+(2.d0*du-cu)*rs/3.d0

          ap= 0.01555d0
          bp=-0.0269d0
          cp= 0.0007d0
          tp=-0.0048d0
          ecp=ap*xlr+bp+cp*rlr+tp*rs
          ucp=ap*xlr+(bp-ap/3.d0)+2.d0/3.d0*cp*rlr+(2.d0*tp-cp)*rs/3.d0

        endif

!  if we are nonrel, turn off the MacDonald-Vosko correction.

        if (abs(rel)<1.d-12) then  !JR
          fe1=1.d0
          fu1=1.d0
          fe2=1.d0
          fu2=1.d0
        endif

!  interpolate the correlation energies.

        denom=2.d0**ft-2.d0
        f=((1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0)/denom
        dfdz=ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
        ec=ecu+f*(ecp-ecu)
        uc1=ucu+f*(ucp-ucu)+(ecp-ecu)*(1.d0-zeta)*dfdz
        uc2=ucu+f*(ucp-ucu)-(ecp-ecu)*(1.d0+zeta)*dfdz        

!  get the final functional and potential.

        ex=(xn1*fe1*ex1+xn2*fe2*ex2)/xn
        ux1=fu1*ux1
        ux2=fu2*ux2
        uc1=uc1
        uc2=uc2
      endif
      return
      end subroutine cacorr
!---------------------------------------------------------------------
      subroutine hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,&
                       zorig,xntot,ixflag,nel,&
                       no,nl,xnj,is,ev,ek,occ,njrc,vi,cq,phe,orb)
      use dprec
      use PARAMS
      implicit none
      integer :: iu,ir,nst,nr,ixflag,nel
      real(dp) :: etot,rel,rmin,rmax,zorig,xntot,r1,rj

      integer :: no(iorbs),nl(iorbs),is(iorbs),njrc(4)
      real(dp) :: xnj(iorbs),ev(iorbs),ek(iorbs),occ(iorbs)
      real(dp) :: vi(nrmax,7),cq(nrmax)
      real(dp) :: phe(nrmax,iorbs),orb(nrmax,iorbs)
      real(dp) :: rho(nrmax)
      integer i,j,k,iu1

      if (iu.lt.0) then
        !'w' has iu=-1.
!       write (6,*) 'enter filename ...'
!JR        read (inpunit,52) filename
!JR 52     format (a)
        iu1=1
!JR        open (unit=iu1,status='unknown',file=filename)
!JR        open(iu1,file=outfile(1:io),status='unknown')   !JR
      else
        iu1=iu
        call gropen(iu1)
      end if
      rewind iu1

      if (ir.eq.0) then
        !'w' has ir=0.
!JR     write ( 6, * ) ' nst = ', nst
!JR     write (iu1,102)  etot,nst,rel,nr,rmin,rmax,zorig,xntot,nel 
      else
        read  (iu1,102)  etot,nst,rel,nr,rmin,rmax,zorig,xntot,nel 
      end if
 102  format (f15.6,i2,f4.1,i5,d15.8,d15.8,f6.1,f12.8,i3)

      if (ir.eq.0) then
        !'w' has ir=0, ixflag=1.
!JR     write (iu1,112)  ixflag
      else
        read  (iu1,112)  ixflag
      end if
 112  format (i4)

      if (ir.eq.0) then
        !'w' has ir=0.
!JR     write (iu1,152) (no(j),nl(j),xnj(j),is(j),ev(j),ek(j),occ(j),&
!JR                     j=1,nel)
      else
        read  (iu1,152) (no(j),nl(j),xnj(j),is(j),ev(j),ek(j),occ(j),&
                        j=1,nel)
      end if
 152  format (i3,i2,f4.1,i2,f12.6,f12.6,f12.8)

      if (ir.eq.0) then
        !'w' has ir=0.
!JR     write (iu1,202) (njrc(i),i=1,4)
      else
        read  (iu1,202) (njrc(i),i=1,4)
      end if
 202  format (4i5)

      if (njrc(1)+njrc(2)+njrc(3)+njrc(4).ne.0) then
        if (ir.eq.0) then
          write (iu1,252) ((vi(k,i),k=1,nr),i=1,7)
          write (iu1,252) (cq(k),k=1,nr)
        else
          read  (iu1,252) ((vi(k,i),k=1,nr),i=1,7)
          read  (iu1,252) (cq(k),k=1,nr)
        end if
      end if
 252  format (4d18.11)

      if (ixflag.eq.0) then
        if (ir.eq.0) then
          write (iu1,252) ((phe(j,i),j=1,nr),i=1,nel)
        else
          read  (iu1,252) ((phe(j,i),j=1,nr),i=1,nel)
        end if
      else
        if (ir.eq.0) then
          !'w' has ir=0, ixflag=1.
          do j=1,nr
            rho(j)=0.d0
            do i=1,nel
              rho(j)=rho(j)+phe(j,i)*phe(j,i)*occ(i)
            enddo
          enddo
          open(iu1,file=outfile(1:io),status='unknown')              !JR
          write(iu1,'(a)') &                                         !JR
     &   'line2: Z,NR,R(1),R(NR) /line3,4,...: R(N),RHO*4*PI*R(N)**2'!JR
          r1=rmin*(rmax/rmin)**(dfloat(1)/dfloat(nr))                !JR
          write(iu1,54) zorig,nr,r1,rmax                             !JR
          do j=1,nr                                                  !JR
            rj=rmin*(rmax/rmin)**(dfloat(j)/dfloat(nr))              !JR
            write(iu1,56) rj,rho(j)                                  !JR
          enddo                                                      !JR
!         write (iu1,252) ((phe(j,i),orb(j,i),j=1,nr),i=1,nel)
        else
          read  (iu1,252) ((phe(j,i),orb(j,i),j=1,nr),i=1,nel)
        end if
      end if
!54    format(f5.2,i6,1p,2e22.14)                                    !JR
!56    format(1p,e9.3,e22.14)                                        !JR
 54   format(f5.2,i6,2es22.14)                                       !JR
 56   format(es10.3,es22.14)                                         !JR

      close (unit=iu1)
      return
      end subroutine hfdisk
!----------------------------------------------------------------------
      subroutine getillls(pin)
      use dprec
      implicit none
      real(dp) :: pin(0:8,0:8,0:16)
      real(dp) :: fa(0:40),si(0:40)
      integer :: i,ia,ib,ic,l,m,n,ll,mm,nn
      real(dp) :: xi,xf,af,bf,cf

      fa(0)=1.d0
      si(0)=1.d0
      do i=1,32
        fa(i)=real(i,dp)*fa(i-1)
        si(i)=-si(i-1)
      end do
      do l=0,8
        do m=0,8
          do n=m+l,0,-2
            xi=0.d0
            xf=2.d0/2.d0**real(n+l+m,dp)
            nn=int(0.0001d0+0.5d0*(real(n,dp)+1.d0))  !JR
            mm=int(0.0001d0+0.5d0*(real(m,dp)+1.d0))  !JR
            ll=int(0.0001d0+0.5d0*(real(l,dp)+1.d0))  !JR
            do ia=nn,n
              af=si(ia)*fa(ia+ia)/fa(ia)/fa(n-ia)/fa(ia+ia-n)
              do ib=ll,l
                bf=si(ib)*fa(ib+ib)/fa(ib)/fa(l-ib)/fa(ib+ib-l)
                do ic=mm,m
                  cf=si(ic)*fa(ic+ic)/fa(ic)/fa(m-ic)/fa(ic+ic-m)
                  xi=xi+xf*af*bf*cf/real(ia*2+ib*2+ic*2-n-l-m+1,dp)
                end do
              end do
            end do
            pin(l,m,n)=xi
          end do
        end do
      end do
      return
      end subroutine getillls
!----------------------------------------------------------------------
      subroutine abinitio(etot,rel,alfa,nr,r,dr,r2,dl,&
                         phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,&
                         ev,occ,is,ek,orb,iuflag,cq,&
                         vold,vnew,vtry,isuse)
      use dprec
      use PARAMS
      implicit none
      integer :: nr,nel,iuflag,mode
      real(dp) :: etot,rel,alfa,dl,zorig,xntot
      integer :: njrc(4),no(iorbs),nl(iorbs),nm(iorbs),is(iorbs)
      real(dp) :: r(nrmax),dr(nrmax),r2(nrmax)
      real(dp) :: phe(nrmax,iorbs),vi(nrmax,7),xnj(iorbs)
      real(dp) :: ev(iorbs),occ(iorbs),ek(iorbs)
      real(dp) :: orb(nrmax,iorbs),cq(nrmax)
      real(dp) :: rp(nrmax,0:15)
      real(dp) :: etkin,etnuc,etcou,etlda,zvalue
      common/parts/etkin,etnuc,etcou,etlda,zvalue
      save/parts/
      integer :: vtry,isuse
      real(dp) :: vold(iorbs),vnew(iorbs)
      integer :: i, j, k, nfc, nj
      real(dp) :: ratio, etol, xnum, eerror

      integer :: ifil
      real(dp) :: wgt, rval
      real(dp), allocatable :: vpert( : ), tab( : )

!  this will be good for going up to and including l=3...

      zvalue=zorig
      do i=0,7
        do k=1,nr
          rp(k,i)=r( k ) ** i
        end do
      end do

!  read in nfc, nel.  refer to the documentation for their meanings.

      read(inpunit, * ) nfc, nel, ratio, etol, xnum, ifil
      allocate( vpert( nr ) )
      do i = 1, nr
        vpert( i ) = 0.d0
      end do
      if ( ifil .gt. 0 ) then
        allocate( tab( ifil ) )
        read(inpunit, * ) wgt
        open( unit=99, file='vvalence', form='formatted',&
             status='unknown' )
        rewind 99
        do i = 1, nr
          read ( 99, * ) rval, tab
          vpert( i ) = tab( ifil )
        end do
        close( unit=99 )
        deallocate( tab )
        ifil = 1
      end if

!  get quantum numbers for levels, Hartree charge, init ev's...

      xntot=0.d0
      do i=nfc+1,nel
        read (inpunit,*) no(i),nl(i),nm(i),xnj(i),is(i),occ(i)
        ev(i)=0.d0
        if ( no( i ) .le. 0 ) read(inpunit, * ) ev( i )
        xntot=xntot+occ(i)
        do j=1,nr
          orb(j,i)=0.d0
        end do
        do j=1,nr
          phe(j,i)=0.d0
        end do
      end do
      if (njrc(1).eq.0) then
        do j=1,nr
          cq(j)=0.d0
        end do
      end if
!-------------------------------------
      mode = 1
      vtry=0
 110  continue
!-------------------------------------
      call atsolve(etot,rel,alfa,eerror,nfc,nr,r,dr,r2,dl,phe,&
                  njrc,vi,zorig,nel,no,nl,nm,xnj,&
                  ev,occ,is,ek,ratio,orb,&
                  xnum,iuflag,cq,vold,vnew,vtry,isuse,&
                  ifil, wgt, vpert )
      eerror=eerror*(1.d0-ratio)/(ratio*ratio)
!  give status report during the job
      if (verb) then
        write (logunit,'(3f16.8,1i5)') eerror,etot,ratio,mode
      endif
!JR      open(unit=99,file='aprog',form='formatted',status='unknown')
!JR      rewind 99
!JR      write (99,'(1x,3f16.8,1i5)') eerror,etot,ratio,mode
!JR      close(unit=99)
!  loop if we needed to ...
!-------------------------------------
      vtry=1
      if ( mode .eq. 1 ) then
        if (eerror.gt.etol) go to 110
      end if
      ratio = 1.d0
      mode = mode + 1
      if ( mode .lt. 4 ) go to 110
!-------------------------------------
!  write out info about the atom.
      do i=1,nel
        nj=int(xnj(i)+xnj(i))  !JR
        write (logunit,'(1x,2i4,i2,i4,a2,i4,f10.4,f18.6)') &
       no(i),nl(i),nm(i),nj,'/2',is(i),occ(i),ev(i)
      end do
      write (logunit,'(1x,1a16,2es22.14)') &
           'total energy =  ',etot,etot*27.2114d0
           
!JR      open(unit=99,file='chgsum',form='formatted',status='unknown')
!JR      rewind 99
!JR      do i=1,nr
!JR        tmp=0.d0
!JR        do j=1,nel
!JR          tmp=tmp+occ(j)*phe(i,j)*phe(i,j)
!JR        end do
!JR        tot = tmp + cq( i )
!JR        n = tot / ( 4.d0 * 3.1415926535898d0 * r( i ) ** 2 )
!JR        write (99,'(1p,e8.2,3e18.10)') r(i), tmp, tot, n
!JR      end do
!JR      close(unit=99)
!JR      write (77,'(1x,1f10.0,1i8)') zorig,nel
!JR      write (77,'(10x,1f20.10)') etot,etkin,etcou,etnuc,etlda

!LD       write(77,'(a,f9.0,i4,a)')'Z,nel=',zorig,nel, &                 !JR
!LD       ', energy in Hartree, Rydberg, and eV'                         !JR
!LD       write(77,'(a,f15.6)') 'Etot =',etot &                          !JR
!LD                            ,'Etkin=',etkin &                         !JR
!LD                            ,'Etcou=',etcou &                         !JR
!LD                            ,'Etnuc=',etnuc &                         !JR
!LD                            ,'Etlda=',etlda                           !JR
!LD       do i=1,nel
!LD         write(77,'(1x,2i4,f6.1,i4,f6.1,2f14.6,f14.4)') &             !JR
!LD        no(i),nl(i),xnj(i),is(i),occ(i),ev(i),ev(i)*2.d0, &           !JR
!LD        ev(i)*27.21138344d0                                           !JR
!LD       end do

      vtry=0
      deallocate( vpert )
      return
      end subroutine abinitio
!----------------------------------------------------------------------
      subroutine atsolve(etot,rel,alfa,eerror,nfc,nr,r,dr,r2,dl,phe,&
                        njrc,vi,zorig,nel,no,nl,nm,xnj,&
                        ev,occ,is,ek,ratio,orb,&
                        xnum,iuflag,cq,vold,vnew,vtry,isuse,&
                        ifil, wgt, vpert )
      use dprec
      use PARAMS
      implicit none
      integer :: nfc, nr
      real(dp) :: etot, rel, alfa, eerror, dl
      real(dp) :: r( nr ), dr( nr ), r2( nr )
      real(dp) :: phe( nrmax, iorbs )

      integer :: njrc( 4 ), nel
      integer :: no( nel ), nl( nel ), nm( nel )
      real(dp) :: vi( nrmax, 7 ), zorig
      real(dp) :: xnj( nel ) 

      integer :: is( nel )
      real(dp) :: ratio
      real(dp) :: ev( nel ), occ( nel ), ek( nel )
      real(dp) :: orb( nrmax, nel )
   
      integer :: iuflag, vtry, isuse
      real(dp) :: xnum
      real(dp) :: cq( nr ), vold( nel ), vnew( nel )

      integer :: ifil
      real(dp) :: wgt, vpert( nr )

      integer :: i, j, ll, idoflag
      real(dp) :: zeff,xkappa,evi,ekk,dq
      real(dp), allocatable :: xm1( : ), xm2( : ), v( : )
      allocate( xm1( nr ), xm2( nr ), v( nr ) )

      eerror=0.d0
      etot=0.d0
      do i=1,nel
        if (i.gt.nfc) then
          idoflag=1
          call setqmm(orb(1,i),nl(i),xnj(i),idoflag,v,zeff,&
                     zorig,rel,nr,r,dl,xm1,xm2,njrc,vi)
          if ( ifil .eq. 1 ) then
            do j = 1, nr
              v( j ) = v( j ) + wgt * vpert( j )
            end do
          end if
          xkappa=-1
          if (abs(xnj(i)).gt.real(nl(i),dp)+0.25d0) xkappa=-nl(i)-1
          if (abs(xnj(i)).lt.real(nl(i),dp)-0.25d0) xkappa= nl(i)
          if (vtry.eq.1) evi=ev(i)+vnew(i)-vold(i)
          isuse=njrc(nl(i)+1)
          if ( no( i ) .gt. 0 ) then
            call elsolve(i,no(i),nl(i),xkappa,xnj(i),zorig,zeff,&
                        evi,phe(1,i),v,xm1,xm2,nr,r,r2,dl,rel,&
                        vtry,isuse)
            if (occ(i).gt.0.d0) then
              if (abs(ev(i)-evi).gt.eerror) eerror=abs(ev(i)-evi)
            end if
            ev(i)=evi
          else
            evi = ev( i ) 
            call elener(nl(i),xkappa,xnj(i),zeff,&
                       evi,phe(1,i),v,xm1,xm2,nr,r,r2,dl,rel,isuse)
          end if
          ll=2
          ekk=0.d0
          do j=nr,1,-1
            dq=phe(j,i)*phe(j,i)
            ekk=ekk+(evi-orb(j,i))*r(j)*dq*real(ll,dp)
            ll=6-ll
          end do
          ek(i)=dl*ekk/3.d0
        end if 
        etot=etot+ek(i)*occ(i)
      end do
      call getpot(etot,rel,alfa,dl,nr,dr,r,r2,phe,ratio,orb,&
                 occ,is,nel,nl,nm,no,xnj,xnum,iuflag,cq,&
                 vold,vnew)

      deallocate( xm1, xm2, v )
      return
      end subroutine atsolve
!-----------------------------------------------------------------------
      subroutine setqmm(orb,l,xj,idoflag,v,zef,&
                       zorig,rel,nr,r,dl,xm1,xm2,njrc,vi)
      use dprec
      use PARAMS
      implicit none
      integer :: l,idoflag,nr
      real(dp) :: xj,zef,zorig,rel,dl
      integer :: njrc(4)
      real(dp) :: v(nr),r(nr),orb(nr)
      real(dp) :: xm1(nr),xm2(nr),vi(nrmax,7)
      
      integer :: j,lp,lpx,lp2,l1,l2,lu,ij
      real(dp) :: alpha,aa,a2,zaa,za2,d1,d2,w1,w2
      real(dp) :: dvdl,ddvdll,dvdr,ddvdrr

      real(dp), allocatable :: vu( : )
      real(dp), allocatable :: o1( : ), o2( : ), o3( : )
      real(dp), allocatable :: vief( : )

      allocate( vu( nr ), o1( nr ), o2( nr ), o3( nr ), vief( nr ) )

! setting up constants which might be used
      alpha=rel/x137
      aa=alpha*alpha
      a2=0.5d0*aa
      lp=l+1
                          lpx=lp
      if (lp.gt.4)        lpx=4
                          lp2=l+l+1
      if (lp2.gt.7)       lp2=7
                          zef=zorig
      if (njrc(lpx).ne.0) zef=0
      zaa=zef*aa
      za2=zef*a2
! we carry on only if idoflag is not zero
      if (idoflag.eq.0) return
      d1=0.5d0/dl
      d2=1.d0/(dl*dl)
      do j=1,nr
        o1(j)=1.d0/r(j)
        o2(j)=o1(j)*o1(j)
        o3(j)=o1(j)*o2(j)
      end do
! below, first fork=full potential, second fork=pseudopotential
!-----------------------------
      if (njrc(lpx).eq.0) then
!-----------------------------
      if (idoflag.ne.1) go to 50
      do j=1,nr
        v(j)=-zef*o1(j)+orb(j)
      end do
   50 do j=1,nr
        vu(j)=orb(j)
      end do
!-----------------------------
      else
!-----------------------------
      if (idoflag.ne.1) go to 70
! insertion A
      lu=0
      ij=int(2.1*(abs(xj)-real(l,dp)))  !JR
      if (l.eq.0) then
        lu=1
      else
        if (ij.lt.0) lu=2*l
        if (ij.gt.0) lu=2*l+1
      end if
! insertion B
      if (lu.gt.0) then
        do j=1,nr
          vief(j)=vi(j,lu)
        end do
      else
        l1=2*l
        l2=2*l+1
        w1=real(l,dp)/real(2*l+1,dp)
        w2=real(l+1,dp)/real(2*l+1,dp)
        do j=1,nr
          vief(j)=w1*vi(j,l1)+w2*vi(j,l2)
        end do
      end if
      do j=1,nr
        v(j)=vief(j)+orb(j)
      end do
   70 do j=1,nr
        vu(j)=v(j)
      end do
!-----------------------------
      end if
!-----------------------------
!  following indept of full potential versus pseudopential
      do j=3,nr-2
        dvdl=(8.d0*(vu(j+1)-vu(j-1))-(vu(j+2)-vu(j-2)))/(12.d0*dl)
        ddvdll=(16.d0*(vu(j+1)+vu(j-1))-(vu(j+2)+vu(j-2))&
              -30.d0*vu(j)) / (12.d0*dl*dl)
        dvdr=dvdl*o1(j)
        ddvdrr=(ddvdll-dvdl)*o2(j)
        xm1(j)=-a2*dvdr-za2*o2(j)
        xm2(j)=-a2*ddvdrr+zaa*o3(j)
      end do
!  inner end point
      xm1(1)=xm1(3)+za2*(o2(3)-o2(1))
      xm2(1)=xm2(3)-zaa*(o3(3)-o3(1))
      xm1(2)=xm1(3)+za2*(o2(3)-o2(2))
      xm2(2)=xm2(3)-zaa*(o3(3)-o3(2))
!  outer end point
      xm1(nr-1)=xm1(nr-2)
      xm2(nr-1)=xm2(nr-2)
      xm1(nr)=xm1(nr-1)
      xm2(nr)=xm2(nr-1)
      deallocate( vu, o1, o2, o3, vief )
      return
      end subroutine setqmm
!----------------------------------------------------------------------
      subroutine augment(e,l,xj,phi,v,nr,r,dl,rel)
      use dprec
      use PARAMS
      implicit none
      integer :: l,nr
      real(dp) :: e,xj,dl,rel
      real(dp) :: phi(nrmax),v(nrmax),r(nrmax)
      integer :: j
      real(dp) :: c,cc,c2,xkappa,g0,ga,gb,gc,gg,f0
      real(dp) :: phi2(nrmax)
!
      c=x137*rel
      cc=c*c
      c2=cc+cc
      xkappa=-1
      if (abs(xj).gt.real(l,dp)+0.25d0) xkappa=-l-1
      if (abs(xj).lt.real(l,dp)-0.25d0) xkappa= l
      do j=4,nr-3
        if (abs(phi(j))>1.d-12) then  !JR
          g0=phi(j)
          ga=(phi(j+1)-phi(j-1))
          gb=(phi(j+2)-phi(j-2))/2.d0
          gc=(phi(j+3)-phi(j-3))/3.d0
          gg=((1.5d0*ga-0.6d0*gb+0.1d0*gc)/(2.d0*dl)+xkappa*g0)/r(j)
          f0=c*gg/(e-v(j)+c2)
          phi2(j)=sqrt(g0*g0+f0*f0)
          if (g0.lt.0.d0) phi2(j)=-phi2(j)
        else
          phi2(j)=phi(j)
        end if
      end do
      do j=1,3
        phi2(j)=phi(j)*phi(4)/phi2(4)
      end do
      do j=1,nr
        phi(j)=phi2(j)
      end do
      return
      end subroutine augment
!----------------------------------------------------------------------
      subroutine initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,xntot,nel)
      use dprec
      use PARAMS
      implicit none
      integer :: nr,nel
      real(dp) :: zorig,rmin,rmax,dl,xntot
      integer :: njrc(4)
      real(dp) :: r(nrmax),dr(nrmax),r2(nrmax)
      integer :: j

      read (inpunit,*) zorig,nr
      rmin=rmifac/zorig         !from params.f,
      rmax=rmafac/sqrt(zorig)  !from params.f
      write(logunit,'(a,2es10.2)')'rmin,rmax=',rmin,rmax
      rmin=1.d-10   !JR
      rmax=60.d0    !JR
      call setgrid(nr,rmin,rmax,r,dr,r2,dl)
      do j=1,4
        njrc(j)=0
      end do
      xntot=0.d0
      nel=0
      return
      end subroutine initiali
!----------------------------------------------------------------------
      subroutine setgrid(nr,rmin,rmax,r,dr,r2,dl)
      use dprec
      use PARAMS
      implicit none
      integer :: nr
      real(dp) :: rmin,rmax,dl
      real(dp) :: r(nrmax),dr(nrmax),r2(nrmax)
      integer :: i
      real(dp) :: rat, xrat, xr1

      rat=rmax/rmin
      dl=log(rat)/real(nr,dp)
      xrat=exp(dl)
      xr1=sqrt(xrat)-sqrt(1.d0/xrat)
      do i=1,nr
        r(i)=rmin*xrat**real(i,dp)
        dr(i)=r(i)*xr1
        r2(i)=r(i)*r(i)
      end do
      return
      end subroutine setgrid
!----------------------------------------------------------------------
      subroutine melwriter( nr, r, dr, nlev, llev, ldim, phi, eps )
      use dprec
      use params
      implicit none
      integer :: nr, nlev, ldim
      integer :: llev( nlev )
      real(dp) :: r( nr ), dr( nr ), phi( ldim, nlev ), eps( nlev )

      integer :: ic, ivl, ivh, nrtab, iadd
      real(dp) :: tmp, drtab

      integer :: j, jj, jrad, iv, idiff, iunit
      real(dp) :: m, rlast

      read(inpunit, * ) ic, ivl, ivh, tmp, nrtab, drtab, iadd
      jrad = 1
      do while ( r( jrad ) .lt. tmp )
        jrad = jrad + 1
      end do
      iunit = 36 + llev( ic )
      do iv = ivl, ivh
        tmp = phi( jrad, iv ) / abs( phi( jrad, iv ) )
        write ( 6, '(2x,1f10.6,1f5.1)' ) eps( iv ), tmp
        idiff = ( llev( iv ) - llev( ic ) ) ** 2
        m = 0.d0
        if ( idiff .eq. 1 ) then
          open( unit=39, file='rmelconv',&
               form='formatted', status='unknown' )
          rewind 39
          do j = 1, nr
            m = m + dr( j ) * r( j ) * phi( j, ic ) * phi( j, iv )
            write ( 39, * ) ' rmel ', r( j ),  m * tmp
          end do
          m = m * tmp
          close( unit=39 )
        end if
        write ( iunit, '(1a1,2i5,2f20.10)' )&
         '#', llev( iv ), nrtab, m, eps( iv )
        rlast = 0.d0
        j = 1
        jj = 0
        do while ( jj .lt. nrtab )
          if ( r( j ) .gt. rlast + drtab ) then
            write ( iunit, '(1x,1i5,2f20.10)' )&
             iv + iadd, r( j ), phi( j, iv ) * tmp
            rlast = r( j )
            jj = jj + 1
          end if
          j = j + 1
        end do
      end do
      return
      end subroutine melwriter
!----------------------------------------------------------------------
      subroutine ppopt( nr, r, ldim, vps, cq, nlev, phi, occ )
      use dprec
      use params
      implicit none
      integer :: nr, ldim, nlev
      real(dp) :: r( nr ), vps( ldim, 7 ), cq( nr )
      real(dp) :: phi( ldim, nlev ), occ( nlev )
      integer :: iu, k, j
      real(dp) :: tmp
      character(len=11) :: jive

      iu = 91
      open(unit=iu,file='fort.91',&
           form='formatted', status='unknown' )
      rewind iu
      read(inpunit,'(1a11)') jive
      write(iu,'(1a11)') jive
      write(iu,*) 3,nr,abs(r(nr-10)*vps(nr-10,1))
      do k=1,nr
        write(iu,*) r(k)
      end do
      write(iu,*) 0
      do k=1,nr
        write(iu,*) vps(k,1)
      end do
      write(iu,*) 1
      do k=1,nr
        write(iu,*) vps(k,3)
      end do
      write(iu,*) 2
      do k=1,nr
        write(iu,*) vps(k,5)
      end do
      do k=1,nr
        tmp=0.d0
        do j=1,nlev
          tmp=tmp+phi(k,j)*phi(k,j)*occ(j)
        end do
        write( iu, '(2(2x,1d22.15))' ) cq(k),tmp
      end do
      return
      end subroutine ppopt
!----------------------------------------------------------------------
      subroutine realspace( nr, nrmax, r, vps )
      use dprec
      implicit none
      integer :: nr, nrmax
      real(dp) :: r( nr ), vps( nrmax, 7 )
      integer :: k
      character(len=7), parameter :: u7 = 'unknown'
      character(len=9), parameter :: f9 = 'formatted'

      open(unit=10,file='fort.10',form=f9,status=u7)
      open(unit=11,file='fort.11',form=f9,status=u7)
      open(unit=12,file='fort.12',form=f9,status=u7)
      rewind 10
      rewind 11
      rewind 12
      do k=1,nr
        write (10,'(1x,3f25.15)') r(k), vps(k,1)*r(k), vps(k,1)
        write (11,'(1x,3f25.15)') r(k), vps(k,3)*r(k), vps(k,3)
        write (12,'(1x,3f25.15)') r(k), vps(k,5)*r(k), vps(k,5)
      end do
      close(unit=10)
      close(unit=11)
      close(unit=12)
      return
      end subroutine realspace
!----------------------------------------------------------------------
      subroutine vpscpp( nr, ldim, r, rsqd, vps )
      use dprec
      use PARAMS
      implicit none
      integer :: nr, ldim
      real(dp) :: r( nr ), rsqd( nr ), vps( ldim, 7 )
      integer :: k
      real(dp) :: corpol, rs, rp, rd, tmp, rb, term, fs, fp, fd

      read (inpunit,*) corpol,rs,rp,rd,tmp
      do k=1,nr
        if (tmp.ge.0.9999d0) then
          fs=(1.d0-exp(-rsqd(k)/(rs*rs)))**tmp
          fp=(1.d0-exp(-rsqd(k)/(rp*rp)))**tmp
          fd=(1.d0-exp(-rsqd(k)/(rd*rd)))**tmp
        else
          if (tmp.ge.-1.5d0) then
            rb=r(k)/rs
            fs=1.d0-exp(-rb)*(1.d0+rb*(1.d0+rb*(0.5d0+0.125d0*rb)))
            rb=r(k)/rp
            fp=1.d0-exp(-rb)*(1.d0+rb*(1.d0+rb*(0.5d0+0.125d0*rb)))
            rb=r(k)/rd
            fd=1.d0-exp(-rb)*(1.d0+rb*(1.d0+rb*(0.5d0+0.125d0*rb)))
          else
            fs=rsqd(k)/(rsqd(k)+rs*rs)
            fp=rsqd(k)/(rsqd(k)+rp*rp)
            fd=rsqd(k)/(rsqd(k)+rd*rd)
          end if
        end if
        term=-0.5d0*corpol/rsqd(k)**2
        vps(k,1)=vps(k,1)+term*fs*fs
        vps(k,2)=vps(k,2)+term*fp*fp
        vps(k,3)=vps(k,3)+term*fp*fp
        vps(k,4)=vps(k,4)+term*fd*fd
        vps(k,5)=vps(k,5)+term*fd*fd
      end do
      return
      end subroutine vpscpp
!-----------------------------------------------------------------------
      subroutine elsolve(i,n,l,xkappa,xj,zorig,zeff,e,phi,v,&
                        xm1,xm2,nr,r,r2,dl,rel,vtry,isuse)
      use dprec
      use PARAMS
      implicit none
      real(dp) :: xkappa,xj,zorig,zeff,e,dl,rel,plead
      integer  :: i,n,l,nr
      real(dp) :: phi(nrmax),v(nrmax),xm1(nrmax),xm2(nrmax)
      real(dp) :: phis(nrmax)
      real(dp) :: r(nrmax),r2(nrmax)

      integer :: vtry,isuse
      logical :: clear

      real(dp) :: el,eh,etol,xnorm,aa,eadj,xi,xo,x0, mult
      integer  :: j, jj, istop, ief, nn

      clear=.false.
      el=-1.2d0*zorig*zorig/real(n*n,dp)-60.d0
      eh=0.d0
      etol=0.00000000001d0
      call getplead(l,xj,rel,xkappa,plead,zeff)

      if (isuse.ne.0) go to 145

      jj=0
      if (vtry.eq.1) go to 165
 155  e=(el+eh)/2.d0
 165  continue
      jj=jj+1
      if (jj.gt.20) go to 145


!  here we do it two ways

      ief = 1
      do while ( ief .ne. 0 )
        ief = 0
        istop = 0
!       write ( 6, * ) el, eh, e
        call intego(e,l,xkappa,n,nn,istop,ief,xo,phi,zeff,v,xm1,&
                   xm2,nr,r,r2,dl,rel,plead)
        if ( ief .eq. - 1 ) then
          el = e
          e = 0.5d0 * ( el + eh )
        end if
        if ( ief .eq. + 1 ) then
          eh = e
          e = 0.5d0 * ( el + eh )
        end if
      end do
      mult = phi( istop )
      do j=1,istop
        phis(j)=phi(j)/phi(istop)
      end do
!     write ( 6, * ) el, eh, e
      call integi(e,l,xkappa,istop,xi,phi,zeff,v,xm1,&
                 xm2,nr,r,r2,dl,rel)
      do j=nr-2,istop,-1
        phi(j)=phi(j)/phi(istop)
      end do
      do j=1,istop
        phi(j)=phis(j)
      end do
      do j = 1, nr
        phi( j ) = phi( j ) * mult
      end do
      aa=0.d0
      do j=1,nr-4,4
        aa=aa+phi(j+0)*phi(j+0)*dl*r(j+0)*14.d0/45.d0&
            +phi(j+1)*phi(j+1)*dl*r(j+1)*64.d0/45.d0&
            +phi(j+2)*phi(j+2)*dl*r(j+2)*24.d0/45.d0&
            +phi(j+3)*phi(j+3)*dl*r(j+3)*64.d0/45.d0&
            +phi(j+4)*phi(j+4)*dl*r(j+4)*14.d0/45.d0
      end do
      xnorm=1.d0/sqrt(aa)
      do j=1,nr
        phi(j)=phi(j)*xnorm
      end do
      if (nn.gt.n-l-1) then
        eh=e
        go to 155
      end if
      if (nn.lt.n-l-1) then
        el=e
        go to 155
      end if
      eadj=e+1.01d0*(xo-xi)/2.d0*phi(istop)*phi(istop)
      if (eadj.lt.el) eadj=el
      if (eadj.gt.eh) eadj=eh
      if (eadj.lt.e) then
        eh=e
        e=eadj
      else
        el=e
        e=eadj
      end if
      if (abs(eh-el).gt.etol) go to 165
      if (.not. clear) then
        clear=.true.
        go to 165
      end if

      go to 200

!  here we do it one way

  145 e=0.5d0*(el+eh)
      istop=0
      call integ(e,l,xkappa,n,nn,istop,ief,x0,phi,zeff,v,xm1,&
                xm2,nr,r,r2,dl,rel,plead)
      if (nn.lt.n-l-1) ief=-1
      if (ief.ne.1) then
        el=e
        if (el.gt.-0.001d0) then
          write (6,*) ' mixing too strong ... ', i
          stop
        end if
      end if
      if (ief.ne.-1) eh=e
      if (abs(eh-el).gt.etol) go to 145

  200 continue
      if (abs(abs(xj)-real(l,dp)).gt.0.25d0)&
       call augment(e,l,xj,phi,v,nr,r,dl,rel)
      do j = nr - 5, nr
        phi( j ) = 0.d0
      end do

      aa=0.d0
      do j=1,nr-4,4
        aa=aa+phi(j+0)*phi(j+0)*dl*r(j+0)*14.d0/45.d0&
            +phi(j+1)*phi(j+1)*dl*r(j+1)*64.d0/45.d0&
            +phi(j+2)*phi(j+2)*dl*r(j+2)*24.d0/45.d0&
            +phi(j+3)*phi(j+3)*dl*r(j+3)*64.d0/45.d0&
            +phi(j+4)*phi(j+4)*dl*r(j+4)*14.d0/45.d0
      end do
      xnorm=1.d0/sqrt(aa)
      do j=1,nr
        phi(j)=phi(j)*xnorm
      end do
      return
      end subroutine elsolve
!---------------------------------------------------------------------
      subroutine integ(e,l,xkappa,n,nn,istop,ief,x0,phi,&
                      z,v,xm1,xm2,nr,r,r2,dl,rel,plead)
      use dprec
      use PARAMS
      implicit none
      real(dp) :: e,xkappa,x0,z,dl,rel,plead
      integer  :: l,n,nn,istop,ief,nr
      real(dp) :: phi(nrmax),v(nrmax),xm1(nrmax),xm2(nrmax)
      real(dp) :: r(nrmax),r2(nrmax)

      real(dp) :: dl2,dl5,c,alpha,a2,za2,xl,xlp,xl2,xl4
      real(dp) :: ss,rtest,ss2,xm,tm,xmx,t,xm0,tmp,maxval
      real(dp) :: xk0,dk0,xk2,dk2,p0,p1,p2
      integer  :: i,j,nnideal,is0,il

      il=1
      maxval=10000000000.d0
      dl2=dl*dl/12.d0
      dl5=10.d0*dl2
      c=x137
      alpha=rel/c
      za2=z*z*alpha*alpha
      a2=alpha*alpha/2.d0
      xl=l
      xlp=l+1
      xl2=0.5d0+xl
      xl4=xl2*xl2

!  set leading power, with 1/2 for Desclaux Numerov.

      if (abs(rel)<1.d-12) then  !JR
        ss=xlp
      else
        rtest=xkappa*xkappa-za2
        if (rtest.lt.0.d0) then
          write (6,*) 'Z>137 IS TOO BIG.'
          stop
        end if  
        ss=sqrt(rtest)
      end if
      ss=plead
      ss2=ss-0.5d0

!  we shall set ief to -1 if energy is too low, +1 if too high.

      ief=0

!  see Desclaux for origin of equations. here, get two points...

      t=e-v(1)
      xm0=1.d0+a2*t
      tm=xm0+xm0
      xmx=xm1(1)/xm0
      xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
      dk0=1.d0+dl2*xk0
      p0=dk0
      phi(1)=p0*sqrt(xm0*r(1))/dk0

      t=e-v(2)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(2)/xm
      xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
      dk2=1.d0+dl2*xk2
      p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*sqrt(xm0/xm)
      phi(2)=p1*sqrt(xm*r(2))/dk2

!  if istop is set, stop there, if zero, it will be ctp.

      is0=istop
      if (istop.eq.0) then
        do j=nr-1,2,-1
          if (e.gt.v(j)) go to 15
        end do
        ief=-1
        return
 15     istop=j
      end if

!  initialize number of nodes, and determine the ideal number.

      nn=0
      nnideal=n-l-1

!  integ out.  count nodes, stop along way if there are too many.

      do i=3,istop+2
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
          if (nn.gt.nnideal) then
            ief=1
            return
          end if
        end if
        p0=p1
        p1=p2
      end do
      if (istop.gt.0) call getxval(phi(istop-2),dl,r(istop),x0)
      if (is0.ne.0) then
        tmp=abs(phi(istop)) ! HERE CD JIVE was made dabs!!!!
        do i=1,istop+2
          phi(i)=phi(i)/tmp
        end do
        return
      end if
      do i=istop+3,nr-1
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        if (p2/p1.gt.1.d0) then
          ief=-1
          return
        end if
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
          if (nn.gt.nnideal) then
            ief=1
            return
          end if
        end if
        p0=p1
        p1=p2
      end do
      return
      end subroutine integ
!---------------------------------------------------------------------
      subroutine intego(e,l,xkappa,n,nn,istop,ief,x0,phi,z,v,xm1,&
                      xm2,nr,r,r2,dl,rel,plead)
      use dprec
      use PARAMS
      implicit none
      real(dp) :: e,xkappa,x0,z,dl,rel,plead
      integer  :: l,n,nn,istop,ief,nr
      real(dp) :: phi(nrmax),v(nrmax),xm1(nrmax),xm2(nrmax)
      real(dp) :: r(nrmax),r2(nrmax)

      real(dp) :: dl2,dl5,c,alpha,za2,a2,xl,xlp,xl2,xl4
      real(dp) :: ss,rtest,ss2,t,xm0,xm,tm,xmx
      real(dp) :: xk0,xk2,dk0,dk2,p0,p1,p2,maxval
      integer  :: i, nnideal, il, count

      real(dp) poteff
      logical cont, ready

      maxval=10000000000.d0
      il=1
      dl2=dl*dl/12.d0
      dl5=10.d0*dl2
      c=x137
      alpha=rel/c
      za2=z*z*alpha*alpha
      a2=alpha*alpha/2.d0
      xl=l
      xlp=l+1
      xl2=0.5d0+xl
      xl4=xl2*xl2

!  set leading power, with 1/2 for Desclaux Numerov.

      if (abs(rel)<1.d-12) then  !JR
        ss=xlp
      else
        rtest=xkappa*xkappa-za2
        if (rtest.lt.0.d0) then
          write (6,*) 'Z>137 IS TOO BIG.'
          stop
        end if  
        ss=sqrt(rtest)
      end if
      ss=plead
      ss2=ss-0.5d0

!  we shall set ief to -1 if energy is too low, +1 if too high.

      ief=0

!  see Desclaux for origin of equations. here, get two points...

      t=e-v(1)
      xm0=1.d0+a2*t
      tm=xm0+xm0
      xmx=xm1(1)/xm0
      xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
      dk0=1.d0+dl2*xk0
      p0=dk0
      phi(1)=p0*sqrt(xm0*r(1))/dk0

      t=e-v(2)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(2)/xm
      xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
      dk2=1.d0+dl2*xk2
      p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*sqrt(xm0/xm)
      phi(2)=p1*sqrt(xm*r(2))/dk2

!  initialize number of nodes, and determine the ideal number.

      nn=0
      nnideal=n-l-1

!  integ out.  count nodes, stop along way if there are too many.

      cont = .true.
      ready = .false.
      count = 0
      i = 3 
      do while ( cont .and. ( ief .eq. 0 ) )
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        if ( dk2 .lt. 0.d0 ) then
!CD        write ( 6, * ) 'dk2 trouble ...' 
          ief = - 1
        end if
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
        end if
        if ( nn .gt. nnideal ) then
          if ( ief .eq. 0 ) ief = 1
        end if
        if ( nn .eq. nnideal ) then
          if ( p2 * p1 .gt. 0.d0 ) then
            if ( p2 / p1 .lt. 1.d0 ) ready = .true.
          end if
        end if
        if ( ( istop .eq. 0 ) .and. ( ready ) ) then
          poteff = v( i ) + xl * xlp / ( 2.d0 * r2( i ) )
          if ( e .lt. poteff ) then
            istop = i - 2
            cont = .false.
          end if
        end if
        if ( ( istop .ne. 0 ) .and. ( i .eq. istop + 2 ) ) then
          cont = .false.
        end if
        p0 = p1
        p1 = p2
        i = i + 1
        if ( i .gt. nr ) then
          if ( ief .eq. 0 ) ief = -1
        end if
      end do
      if ( ief .eq. 0 ) then
         call getxval( phi( istop - 2 ), dl, r( istop ), x0 )
      end if
      return
      end subroutine intego
!---------------------------------------------------------------------
      subroutine integi(e,l,xkappa,istop,x0,phi,z,v,xm1,&
                      xm2,nr,r,r2,dl,rel)
      use dprec
      use PARAMS
      implicit none
      real(dp) :: phi(nrmax),v(nrmax),xm1(nrmax),xm2(nrmax)
      real(dp) :: r(nrmax),r2(nrmax)
      real(dp) :: dl2,dl5,c,alpha,za2,a2,xl,xlp,xl2,xl4
      real(dp) :: p0,p1,p2,t,xm,tm,xmx
      real(dp) :: xk2,dk2,x0,e,rel,dl,z,xkappa,verybig
      integer  :: i,l,istop,nr,ih

      verybig=1.d+10 !10000000000.d0
      ih=nr-2
      dl2=dl*dl/12.d0
      dl5=10.d0*dl2
      c=x137
      alpha=rel/c
      za2=z*z*alpha*alpha
      a2=alpha*alpha/2.d0
      xl=l
      xlp=l+1
      xl2=0.5d0+xl
      xl4=xl2*xl2

! set up bogus tail start...

      p0=1.d0
      p1=1.d0
      do i=istop-2,nr
        phi(i)=0.d0
      end do
      i=nr-2
      t=e-v(i)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(i)/xm
      xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
      dk2=1.d0+dl2*xk2

!  integ out.  count nodes, stop along way if there are too many.

      do i=nr-3,istop-2,-1
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.verybig)then
          !ifort fails when temper is in-line optimized, use option -O0.
          call temper(i,ih,phi(i),p0,p1,p2)
        endif
        p0=p1
        p1=p2
      end do
      call getxval(phi(istop-2),dl,r(istop),x0)
      return
      end subroutine integi
!-----------------------------------------
      subroutine getxval(phi,dl,rctr,x0)
      use dprec
      implicit none
      real(dp), intent(in)  :: phi(-2:2), dl, rctr
      real(dp), intent(out) :: x0
      real(dp) :: psip2, psip1, psip
!
      psip2=phi(2)-phi(-2)
      psip1=phi(1)-phi(-1)
      psip=(8.d0*psip1-psip2)/(12.d0*dl*rctr)
      x0=psip/phi(0)
      return
      end subroutine getxval
!---------------------------
      subroutine getplead(l,xnj,rel,xkappa,plead,zeff)
      use dprec
      use PARAMS
      implicit none
      integer  :: l
      real(dp) :: xnj,rel,xkappa,plead,zeff,zzaa
!
      if (rel.lt.0.5d0) then
        plead=l+1
      else
        zzaa=(zeff/x137)*(zeff/x137)
        if (abs(abs(xnj)-real(l,dp)).gt.0.25d0) then
          plead=sqrt(xkappa*xkappa-zzaa)
        else
          plead=0.5d0*(1.d0-zeff+ &
           sqrt((zeff-1.d0)*(zeff-1.d0)+ &
                 4.d0*(real(l*(l+1),dp)+zeff-zzaa)))
        end if
      end if
      return
      end subroutine getplead
!---------------------------------------------------------------------
      subroutine temper( i1, i2, phi, p0, p1, p2 )
      use dprec
      implicit none
      integer  :: i1, i2, i
      real(dp), intent(inout) :: phi( i1 : i2 ), p0, p1, p2
      real(dp) :: pdiv2

      pdiv2 = abs( p2 )
      do i = i1, i2
        phi( i ) = phi( i ) / pdiv2
      end do
      p0 = p0 / pdiv2
      p1 = p1 / pdiv2
      p2 = p2 / pdiv2
      return
      end subroutine temper
!----------------------------------------------------------------------
      subroutine pseudo(etot,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,&
                       phe,orb,njrc,vi,cq,zorig,xntot,nel,&
                       no,nl,nm,xnj,ev,occ,is,iuflag,vctab,&
                       vold,vnew,vtry,isuse)
      use dprec
      use PARAMS
      implicit none
      real(dp), dimension(nrmax) :: r, dr, r2, xm1, xm2
      real(dp) :: vi(nrmax,7),phe(nrmax,iorbs)
      real(dp) :: rpower(nrmax,0:15),orb(nrmax,iorbs)
      real(dp) :: vctab(nrmax,0:3),cq(nrmax)

      integer  :: njrc(4),mjrc(4)
      integer  :: no(iorbs),nl(iorbs),nm(iorbs),is(iorbs)
      real(dp) :: xnj(iorbs),ev(iorbs),occ(iorbs)
      real(dp) :: vq(nrmax)

      integer  :: vtry,isuse
      real(dp) :: vold(iorbs),vnew(iorbs)
      real(dp) :: pi,tmp,aa,bb
      integer  :: i,j,idoflag,np,lfc,nr,im,k,nel,lu,in,iuflag,ij,l
      real(dp) :: corpol,rnorm,zuse,zorig,ratio,dl,cor0,corp,corm,f,fp
      real(dp) :: rhs,xl,xh,br,diff,ruse,xkappa,xntot,zeff,rel
      real(dp) :: rmin,rmax,xi,xnum,etot,alfa,asym

      pi=4.d0*atan(1.d0)

      do i=1,4
        if (njrc(i).gt.0) stop 'cannot repseudize as of now'
        mjrc(i)=0
      end do

      read (inpunit,*) np,corpol,rnorm
      zuse=zorig
      do i=1,np-1
        zuse=zuse-occ(i)
      end do
      read (inpunit,*) lfc,ratio
      do i=1,nr
        cq(i)=0.d0
        vq(i)=0.d0
      end do
      if (lfc.ne.0) then
        im=0
        do i=nr,1,-1
          do k= 1,np-1
            cq(i)=cq(i)+phe(i,k)*phe(i,k)*occ(k)
          end do
          do k=np,nel
            vq(i)=vq(i)+phe(i,k)*phe(i,k)*occ(k)
          end do
          if ((im.eq.0).and.((cq(i)*ratio).gt.vq(i))) im=i
        end do
        if (ratio.lt.0) im=int(log(-ratio/r(1))/dl)
        write (6,*)   im 
        write (6,*) r(im)
        cor0=cq(im  )/r(im  )
        corp=cq(im+1)/r(im+1)
        corm=cq(im-1)/r(im-1)
        f =cor0
        fp=(corp-corm)/(2.d0*dl*r(im))
        rhs=r(im)*fp/f
        if (rhs.gt.0.d0) then
          xl=0.0000000d0
        else
          xl=0.5d0*pi
        end if
        xh=xl+0.5d0*pi
   20   br=(xl+xh)/2.d0
        diff=tan(br)-(br/rhs)
        if (verb) write (6,22) br,diff
   22   format(1x,2f20.10)
        if (diff.ge.0.d0) xh=br
        if (diff.le.0.d0) xl=br
        if (abs(xh-xl).ge.0.0000001d0) go to 20
        bb=br/r(im)
        aa=f/dsin(bb*r(im))
        write (6,32) 'lfc a=',aa,'  b=',bb
   32   format(1x,1a6,1f10.6,1a4,1f10.6)
        open(unit=99,file='corchg',form='formatted',status='unknown')
        rewind 99
        do j=1,nr
          tmp=cq(j)
          if (j.lt.im) cq(j)=r(j)*aa*dsin(bb*r(j))
          write (99,'(1x,3f20.10)') r(j),tmp,cq(j)
        end do
        close(unit=99)
      end if
      ruse=0.d0
      xkappa=-1.d0
   42 format(1x,1a2,1i1,1a7,1f4.1)
   44 format(1x,1a13,1i5,3x,1f4.1,3x,1f14.6)
      xntot=0.d0
      do i=np,nel
        write (6,42) 'l=',nl(i),' ... j=',abs(xnj(i))
        lu=2*nl(i)+1
        if (abs(xnj(i))+0.25d0.lt.real(nl(i),dp)) lu=2*nl(i)
        do j=1,nr
          orb(j,i)=orb(j,i)+vctab(j,nl(i))
        end do
        idoflag=1
        call setqmm(orb(1,i),nl(i),xnj(i),idoflag,vi(1,lu),zeff,&
                   zorig,rel,nr,r,dl,xm1,xm2,mjrc,vi)
        call zout(nr,orb(1,i))
        call pseudize(orb(1,i),ev(i),nl(i),xnj(i),no(i),njrc,zeff,&
                     vi(1,lu),xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
        no(i)=nl(i)+1
        vtry=1
        isuse=0
        call elsolve(i,no(i),nl(i),xkappa,xnj(i),&
                    zorig,zeff,ev(i),phe(1,i),vi(1,lu),&
                    xm1,xm2,nr,r,r2,dl,ruse,vtry,isuse)
        write (6,44) 'solution ... ',nl(i),xnj(i),ev(i)
        in=i-(np-1)
        do j=1,nr
          phe(j,in)=phe(j,i)
        end do
        no (in)=no (i)
        nl (in)=nl (i)
        nm (in)=nm (i)
        xnj(in)=xnj(i)
        is (in)=is (i)
        ev (in)=ev (i)
        occ(in)=occ(i)
        xntot=xntot+occ(in)
      end do

      nel=nel-(np-1)
      do i=0,7
        xi=i
        do k=1,nr
          rpower(k,i)=r(k)**xi
        end do
      end do
      xnum=10000.d0
      ratio=1.d0
      call getpot(etot,rel,alfa,dl,nr,dr,r,r2,phe,&
                 ratio,orb,occ,is,nel,nl,nm,no,xnj,xnum,&
                 iuflag,cq,vold,vnew)
      do i=1,nel
        lu=2*nl(i)+1
        if (abs(xnj(i))+0.25d0.lt.real(nl(i),dp)) lu=2*nl(i)
        do j=1,nr
          vi(j,lu)=vi(j,lu)-orb(j,i)
        end do
        ij=int(2.1d0*(abs(xnj(i))-real(nl(i),dp)))  !JR
        if ((nl(i).gt.0).and.(ij.eq.0)) then
          do j=1,nr
            vi(j,lu-1)=vi(j,lu)
          end do
        end if
        vi(1,lu)=vi(2,lu)
      end do
      rel=0.d0
      do k=1,nr
        if (r(k).gt.rnorm) then
          asym=-zuse/r(k)-corpol/(2.d0*r(k)**4.d0)
          do l=1,7
            vi(k,l)=asym
          end do
        end if
      end do
      return
      end subroutine pseudo
!----------------------------------------------------------------------
      subroutine parabreg(f,fp,fpp,rf,vf)
      use dprec
      implicit none
      real(dp), intent(in)  :: rf(3), vf(3)
      real(dp), intent(out) :: f, fp, fpp
      real(dp) :: r21, r32, v21, v32

      f=vf(2)
      r21=rf(2)-rf(1)
      r32=rf(3)-rf(2)
      v21=vf(2)-vf(1)
      v32=vf(3)-vf(2)
      fp=(v21+v32)/(r21+r32)
      fpp=2.d0*(v32/r32-v21/r21)/(r21+r32)
      return
      end subroutine parabreg
!----------------------------------------------------------------------
      function hb(x,factor)
      use dprec
      implicit none
      real(dp), intent(in) :: x, factor
      real(dp) :: hb
!
      if (x.gt.3.d0) then
        hb=0.d0
      else
        hb=0.01d0**((dsinh(x/factor)/1.1752d0)**2.d0)
      end if
      return
      end function hb
!----------------------------------------------------------------------
      subroutine fitx0(orb,rcut,njrc,e,l,xj,n,jrt,xideal,phi,&
                      zeff,v,xm1,xm2,nr,r,dr,r2,dl,rel,factor,xkappa)
      use dprec
      use PARAMS
      implicit none
      integer  :: l,n,jrt,nr
      real(dp) :: rcut,e,xj,xideal,zeff,dl,rel,factor,xkappa
      real(dp) :: plead
      integer  :: njrc(4)
      real(dp) :: orb(nrmax),phi(nrmax),v(nrmax)
      real(dp) :: xm1(nrmax),xm2(nrmax)
      real(dp) :: r(nrmax),dr(nrmax),r2(nrmax)

      integer  :: idoflag,nn,ief,ii
      real(dp) :: vl,vh,dum1,dum2(nrmax,7),xactual,xla,xerror,dxdla
      real(dp) :: vmaybe,tmp 
      real(dp) :: hb
      external hb

      vl=-1000000.d0
      vh=+1000000.d0
 115  idoflag=2
      call setqmm(orb,l,xj,idoflag,v,zeff,&
                 dum1,rel,nr,r,dl,xm1,xm2,njrc,dum2)
      call getplead(l,xj,rel,xkappa,plead,zeff)
      call integ(e,l,xkappa,n,nn,jrt,ief,xactual,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,rel,plead)
      if (nn.ne.0) then
        vl=v(1)
        xla=1.d0
      else
        if (xactual.gt.xideal) then
          vh=v(1)
        else
          vl=v(1)
        end if
        xerror=xideal-xactual
        if (abs(xerror).lt.0.000000001d0) return
        dxdla=0.d0
        do ii=1,jrt
          tmp=r(ii)/rcut
          dxdla=dxdla+dr(ii)*phi(ii)*phi(ii)*hb(tmp,factor)
        end do
        dxdla=2.d0*dxdla/(phi(jrt)*phi(jrt))
        xla=xerror/dxdla
      end if
      vmaybe=v(1)+xla
      if ((vmaybe.gt.vh).or.(vmaybe.lt.vl)) xla=(vl+vh)/2.d0-v(1)
      do ii=1,jrt-1
         v(ii)=v(ii)+xla*hb(r(ii)/rcut,factor)
      end do
      go to 115
      end subroutine fitx0
!----------------------------------------------------------------------
      subroutine pseudize(orb,ev,l,xj,n,njrc,zeff,v,&
                         xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
      use dprec
      use PARAMS
      implicit none
      integer  :: l,n,nr
      real(dp) :: ev,xj,zeff,rmin,rmax,dl,rel
      integer  :: njrc(4)
      real(dp) :: orb(nrmax),v(nrmax),xm1(nrmax),xm2(nrmax)
      real(dp) :: r(nrmax),dr(nrmax),r2(nrmax)

      real(dp) :: phi(nrmax),rf(3),vf(3)
      real(dp) :: phi0(nrmax),yl(nrmax),vraw(nrmax)

      integer  :: ii,jj,j,k,jrc,jrt,lp,nn,istop,ief
      real(dp) :: xkappa,xdummy,rcut,factor,rtest
      real(dp) :: switch,ee,de,xp,xm,dqddel
      real(dp) :: ruse,dvdl,ddvdll,dldr,ddldrr,rr
      real(dp) :: v0,v1,v2,b0,b2,b4,xi0,xi1,xi2
      real(dp) :: f,fp,fpp,psi,psip,psipp,ph2,quant,deltal
      real(dp) :: c0,x0,xn0,c00,x00,xn00,plead
      real(dp) :: elsa,elsb,elssnh,elscsh,derr1,derr2

      real(dp) :: hb
      external hb

      lp=l+1
      xkappa=-1
      if (abs(xj).gt.real(l,dp)+0.25d0) xkappa=-l-1
      if (abs(xj).lt.real(l,dp)-0.25d0) xkappa= l
      write (6,*) xkappa
      istop=nr-10
      call getplead(l,xj,rel,xkappa,plead,zeff)
      call integ(ev,l,xkappa,n,nn,istop,ief,xdummy,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,rel,plead)
      read (inpunit,*) rcut,factor

!  if rcut is negative, use it as fraction from outermost node
!  to outermost maximum for computing rcut ...
!  if there are no nodes, the effective node becomes the origin
!
      if (rcut.lt.0.d0) then
        j=1
        do ii=1,n-l-1
          do while (phi(j+1)*phi(j).gt.0.d0)
            j=j+1
          end do 
        end do
        k=j+1
        do while (phi(k+1)/phi(k).gt.1.d0)
          k=k+1
        end do
        rcut=r(j)+abs(rcut)*(r(k)-r(j))
        write (6,22) k,r(k)
        write (6,22) j,r(j)
   22   format(1x,1i5,1f10.4)
      end if
      jrc=int(1.d0+real(nr-1,dp)*log(rcut /rmin)/log(rmax/rmin)) !JR
      rcut=r(jrc)
      rtest=2.d0*rcut
      jrt=int(1.d0+real(nr-1,dp)*log(rtest/rmin)/log(rmax/rmin)) !JR
      njrc(l+1)=jrt
      rtest=r(jrt)
      switch=phi(jrt)/abs(phi(jrt))
      write (6,92) 'rcutoff = ',rcut ,'  jrc = ',jrc
      write (6,92) 'rtest   = ',rtest,'  jrt = ',jrt
 92   format (1x,1a10,1f8.4,1a8,1i5)
      call integ(ev,l,xkappa,n,nn,jrt,ief,x00,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,rel,plead)
      do ii=1,jrt
        phi(ii)=phi(ii)/phi(jrt)
      end do
      xn00=0.d0
      do ii=1,jrt-1
        xn00=xn00+dr( ii)*phi( ii)*phi( ii)
      end do
      xn00=xn00+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
      de=0.0001d0
      ee=ev+de/2.d0
      call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,rel,plead)
      ee=ev-de/2.d0
      call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,rel,plead)
      c00=(xm-xp)/(2.d0*de)
      write (6,*) 'dx/de, x, and norm ... '
      if (verb) write (6,94) c00,x00,xn00
 94   format (1x,3d15.8)
      ruse=0.d0
      v0=v(jrc)
      dvdl  =(8.d0*(v(jrc+1)-v(jrc-1))-(v(jrc+2)-v(jrc-2)))&
              /(12.d0*dl)
      ddvdll=(16.d0*(v(jrc+1)+v(jrc-1))&
      -30.d0*v(jrc)-v(jrc+2)-v(jrc-2))&
                  /(12.d0*dl*dl)
      dldr=1.d0/r(jrc)
      ddldrr=-1.d0/r2(jrc)
      v1=dvdl*dldr
      v2=dvdl*ddldrr+ddvdll*dldr*dldr
      b4=(v2*rcut-v1)/(8.d0*rcut**3.d0)
      b2=(v1-4.d0*b4*rcut**3.d0)/(2.d0*rcut)
      b0=v0-b4*rcut**4.d0-b2*rcut**2.d0
      do ii=1,jrc
        rr=r(ii)
        v(ii)=b0+b2*rr**2.d0+b4*rr**4.d0
      end do
      xkappa=-1
      if (abs(xj).gt.real(l,dp)+0.25d0) xkappa=-l-1
      if (abs(xj).lt.real(l,dp)-0.25d0) xkappa= l
      call fitx0(orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,&
                xm1,xm2,nr,r,dr,r2,dl,ruse,factor,xkappa)
      do ii=1,jrt
        phi0(ii)=phi(ii)
        vraw(ii)=v(ii)
      end do
      xi0=0.d0
      xi1=0.d0
      xi2=0.d0
      do ii=1,jrt
        f=hb(r(ii)/rcut,factor)
        ph2=dr(ii)*phi0(ii)*phi0(ii)
        xi0=xi0+ph2
        if (ii.le.jrt) then
          xi1=xi1+ph2*f
          xi2=xi2+ph2*f*f
        end if
      end do
      ph2=phi0(jrt)*phi0(jrt)
      xi0=xi0/ph2
      xi1=xi1/ph2
      xi2=xi2/ph2
      quant=xi1*xi1+xi2*(c00-xi0)
      if (quant.gt.0.d0) then
        deltal=(sqrt(xi1*xi1+xi2*(c00-xi0))-xi1)/xi2
      else
        deltal=(c00-xi0)/(2.d0*xi1)
      end if
      write (6,222) 'deltal = ',deltal
 222  format (1x,1a9,1f11.8)
 225  continue
      do ii=1,jrt
        yl (ii)=phi0(ii)*hb(r(ii)/rcut,factor)
        phi(ii)=phi0(ii)+deltal*yl(ii)
        if (phi(ii).lt.0.d0) then
          write (6,*) 'big trouble!!! cross axis!!!'
          stop
        end if
      end do
      do ii=2,jrt-1
        if (abs(phi(ii))<1.d-12.or.abs(yl(ii))<1.d-12) go to 1170  !JR
        jj=ii
        if (ii.eq.1) jj=2
        do j=jj-1,jj+1
          rf(2+j-jj)=r(j)
          vf(2+j-jj)=hb(r(j)/rcut,factor)
        end do
        call parabreg(f,fp,fpp,rf,vf)
        do j=jj-1,jj+1
          vf(2+j-jj)=phi0(j)
        end do
        call parabreg(psi,psip,psipp,rf,vf)
        elsa=log(0.01d0)/(1.1752d0**2)
        elsb=1.d0/(factor*rcut)
        elssnh=dsinh(elsb*r(ii))
        elscsh=dcosh(elsb*r(ii))
        derr1=2.d0*elsa*elsb*elssnh*elscsh
        derr2=2.d0*elsa*elsb*elsb*(elscsh*elscsh+elssnh*elssnh)&
             +4.d0*elsa*elsa*elsb*elsb*elssnh*elssnh*elscsh*elscsh
        if ((r(ii)/rcut).gt.3.d0) then
          derr1=0.d0
          derr2=0.d0
        end if
        v(ii)=vraw(ii)+&
         (1.d0-phi0(ii)/phi(ii))*(2.d0*psip/psi*derr1+derr2)/2.d0
      end do
      v(1)=v(2)
      v(nr)=v(nr-1)*r(nr)/r(nr-1)
 1170 call fitx0(orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,&
                xm1,xm2,nr,r,dr,r2,dl,ruse,factor,xkappa)
      call getplead(l,xj,ruse,xkappa,plead,zeff)
      call integ(ev,l,xkappa,n,nn,jrt,ief,x0,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,ruse,plead)
      do ii=1,jrt
        phi(ii)=phi(ii)/phi(jrt)
      end do
      xn0=0.d0
      do ii=1,jrt-1
        xn0=xn0+dr(ii)*phi(ii)*phi(ii)
      end do
      xn0=xn0+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
      de=0.0001d0
      ee=ev+de/2.d0
      call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,ruse,plead)
      ee=ev-de/2.d0
      call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,&
                xm1,xm2,nr,r,r2,dl,ruse,plead)
      c0=(xm-xp)/(2.d0*de)
      if (verb) write (6,94)  c0,x0,xn0
      if (abs(c0-c00).ge.0.000000001d0) then
        dqddel=2.*(xi1+deltal*xi2)
        deltal=deltal+(c00-c0)/dqddel
        go to 225
      end if

      write (6,*) 'ncpp achieved !!!'
      return
      end subroutine pseudize
!----------------------------------------------------------------------
      subroutine zout(n,x)
      use dprec
      implicit none
      integer  :: n,i
      real(dp), intent(out) :: x(n)

      do i=1,n
        x(i)=0.d0
      end do
      return
      end subroutine zout
!----------------------------------------------------------------------
      subroutine fourier(nr,r,dr,vi)
      use dprec
      use PARAMS
      implicit none
      integer :: nr
      real(dp) :: r(nrmax),dr(nrmax),vi(nrmax,7)
      real(dp) :: a(nrmax),v1(nrmax)

      character(len=7), parameter :: u7='unknown'
      character(len=9), parameter :: f9='formatted'
      integer  :: l,lp2,i,ii
      real(dp) :: dl,dl1,dl2,al,ar,q,vq

      dl=log(r(2)/r(1))
      dl1=12.d0*dl
      dl2=12.d0*dl*dl

      do l=0,2
        lp2=l+l+1
        do i=1,nr
          a(i)=r(i)*vi(i,lp2)
        end do
        do i=3,nr-2
          al =(8.d0*(a(i+1)-a(i-1))-(a(i+2)-a(i-2)))/dl1
          ar =al/r(i)
          v1(i)=ar
        end do
        if(l.eq.0)open(unit=20,file='fort.20',form=f9,status=u7)
        if(l.eq.1)open(unit=21,file='fort.21',form=f9,status=u7)
        if(l.eq.2)open(unit=22,file='fort.22',form=f9,status=u7)
        rewind 20+l
        do ii=1,200
          q=real(ii,dp)/10.d0
          vq=0.d0
          do i=3,nr-2
            vq=vq+dr(i)*dcos(q*r(i))*v1(i)
          end do
          write (20+l,*) q,vq
        end do
        close(unit=20+l)
      end do
      return
      end subroutine fourier
!------------------------------------------------------------------------
!  exchange correlation routine, ceperley-alder data, 
!  as parametrized by vosko, wilk and nusair., with macdonald vosko.
      subroutine exchcorr(rel,rr,rh1,rh2,nex,ec,ux1,ux2,uc1,uc2,x137)
      use dprec
      implicit none
      real(dp) :: rel,rr,rh1,rh2,ec,nex,ux1,ux2,uc1,uc2,x137
      real(dp) :: trd,ft,pi,fp,xn1,xn2,xn

      trd = 1.d0 / 3.d0
      ft = 4.d0 / 3.d0
      pi = 3.141592653589793238462643383279d0
      fp = 4.d0 * pi
      xn1=rh1/(rr*fp)
      xn2=rh2/(rr*fp)
      xn=xn1+xn2
      call getx(xn1,xn2,nex,ux1,ux2,trd,ft,pi,x137,rel)
      nex = nex * rr * fp
      call getc(xn1,xn2,xn,ec,uc1,uc2)
      return
      end subroutine exchcorr
!---------------------------------------------
      subroutine getc(xn1,xn2,xn,ec,uc1,uc2)
      use dprec
      implicit none
      real(dp) :: xn1,xn2,ec,uc1,uc2,nec,xn

      if (xn.lt.1.d-10) then
        uc1=0.d0
        uc2=0.d0
        ec=0.d0
      else
        call getc2(xn1,xn2,ec,nec,uc1,uc2)
      end if
      return
      end subroutine getc
!---------------------------------------------
      subroutine getc2(n1,n2,ec,nec,uc1,uc2)
      use dprec
      implicit none
      real(dp) pi,rs,ecp,ecf,ac,fd,zeta,fn,f,fpp0,ft,trd
      real(dp) beta,ec,nec,n1,n2,n,eecp,eecf,eeca,z4
      real(dp) dfdz,dzdn1,dzdn2,uc1,uc2,der1,der2
      real(dp) dera,derp,derf,drsdn,dmlt
      external eeca,eecp,eecf

      parameter(ft=4.d0/3.d0,trd=1.d0/3.d0)
      fd=2.d0**ft-2.d0
      pi = 3.141592653589793238462643383279d0
      fpp0=(8.d0/9.d0)/fd

      n=n1+n2
      zeta=(n1-n2)/n
      dzdn1=  1.d0/n-zeta/n
      dzdn2= -1.d0/n-zeta/n
      z4=zeta*zeta*zeta*zeta
      fn=(1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0
      f=fn/fd
      dfdz=(ft*(1.d0+zeta)**trd-ft*(1.d0-zeta)**trd)/fd
      rs=(3.d0/(4.d0*pi*n))**trd
      drsdn=-trd*rs/n
      dmlt=drsdn*0.5d0/sqrt(rs)

      ecp = eecp(rs,derp)
      ecf = eecf(rs,derf)
      ac  = eeca(rs,dera)
      derp=derp*dmlt
      derf=derf*dmlt
      dera=dera*dmlt

      beta=fpp0*(ecf-ecp)/ac-1.d0
      ec=ecp+ac*(f/fpp0)*(1.d0+beta*z4)

      der1=derp+dera*(f/fpp0)       *(1.d0+beta*z4)&
         +ac  *dfdz*dzdn1/fpp0*(1.d0+beta*z4)&
         +ac*(f/fpp0)*fpp0*z4*((derf-derp)-dera*(ecf-ecp)/ac)/ac&
         +ac  *(f/fpp0)       *beta*zeta*zeta*zeta*4.d0*dzdn1
      der2=derp+dera*(f/fpp0)       *(1.d0+beta*z4)&
         +ac  *dfdz*dzdn2/fpp0*(1.d0+beta*z4)&
         +ac*(f/fpp0)*fpp0*z4*((derf-derp)-dera*(ecf-ecp)/ac)/ac&
         +ac  *(f/fpp0)       *beta*zeta*zeta*zeta*4.d0*dzdn2

      nec=n*ec
      uc1=ec+n*der1
      uc2=ec+n*der2
      return
      end subroutine getc2
!----------------------------------------------------------------------
      function eecp(rs,der)
      use dprec
      implicit none
      real(dp) :: rs,eecp,a,b,c,x0,vwngen,der
      external vwngen

      parameter(a=0.0621814d0,b=3.72744d0,c=12.9352d0,x0=-0.10498d0)
      eecp=vwngen(rs,a,b,c,x0,der)
      return
      end function eecp
!----------------------------------------------------------------------
      function eecf(rs,der)
      use dprec
      implicit none
      real(dp) :: a, b, c, x0
      real(dp) :: rs, eecf, vwngen, der
      external vwngen
      parameter (a=0.0310907d0,b=7.06042d0,c=18.0578d0,x0=-0.32500d0)
      
      eecf=vwngen(rs,a,b,c,x0,der)
      return
      end function eecf
!----------------------------------------------------------------------
      function eeca(rs,der)
      use dprec
      implicit none
      real(dp) :: b, c, x0, pi
      real(dp) :: rs, eeca, a, vwngen, der
      external vwngen
      parameter (b=1.13107d0,c=13.0045d0,x0=-0.00475840d0)
      parameter (pi=3.141592653589793238462643383279d0)
      
      a=-1.d0/(3.d0*pi*pi)
      eeca=vwngen(rs,a,b,c,x0,der)
      return
      end function eeca
!----------------------------------------------------------------------
      function vwngen(rs,a,b,c,x0,vwnder)
      use dprec
      implicit none
      real(dp) :: rs,a,b,c,x0,t1,t2,t3,x
      real(dp) :: q,t,t3a,t3b,xfcn,vwngen,vwnder
      real(dp) :: tder,tder1,tder2,tder3,tder3a,tder3b
      real(dp) :: xx,xx0,xd,xd0
      external xfcn
!
      x=sqrt(rs)
      q=sqrt(4.d0*c-b*b) 

      xx=xfcn(x,b,c)
      xd=2.d0*x+b
      xx0=xfcn(x0,b,c)
      xd0=2.d0*x0+b

      t=atan(q/(2.d0*x+b))
      tder=-2.d0*q/(q*q+(2.d0*x+b)*(2.d0*x+b))

      t1    = log(x*x/xx)
      tder1 = 2.d0/x-xd/xx

      t2    = t    * 2.d0*b/q
      tder2 = tder * 2.d0*b/q

      t3a=log((x-x0)*(x-x0)/xx)
      tder3a=2.d0/(x-x0)-xd/xx

      t3b    = t    * 2.d0*(2.d0*x0+b)/q
      tder3b = tder * 2.d0*(2.d0*x0+b)/q

      t3   =-b*x0/xx0*(t3a+t3b)
      tder3=-b*x0/xx0*(tder3a+tder3b)

      vwngen=0.5d0*a*(t1+t2+t3)
      vwnder=0.5d0*a*(tder1+tder2+tder3)
      return
      end function vwngen
!----------------------------------------------------------------------
      function xfcn(x,b,c)
      use dprec
      implicit none
      real(dp) :: x, b, c, xfcn
      
      xfcn=c+x*(b+x)
      return
      end function xfcn
!----------------------------------------------------------------------
      subroutine getx(xn1,xn2,nex,ux1,ux2,trd,ft,pi,x137,rel)
      use dprec
      implicit none
      real(dp) :: xn1,xn2,nex,ux1,ux2,trd,ft,pi,x137,rel
      real(dp) :: fe1,fu1,fe2,fu2,ee,bb,ex1,ex2
      real(dp) :: b, e, d

      ee=-1.5d0*(0.75d0/pi)**trd
      bb=(9.d0*pi/4.d0)**(1.d0/3.d0)/x137

      fe1 = 1.d0
      fu1 = 1.d0
      ex1 = ee * xn1 ** trd
      ux1 = ex1 * ft
      if ( rel .gt. 0.5d0 ) then
        d = bb * ( 8.d0 * pi / 3.d0 ) ** trd
        b = d * xn1 ** trd
        e = sqrt( 1.d0 + b ** 2 )
        if ( b .gt. 0.0001d0 ) then
          fe1 =  1.0d0 - 1.5d0 * ( b*e - log(b+e) ) ** 2 / b ** 4
          fu1 = -0.5d0 + 1.5d0 * log( b + e ) / ( b * e )
        else
          fe1 = 1.d0 - 2.d0 * b ** 2 / 3.d0 + 2.d0 * b ** 4 / 5.d0 &
               - 48.d0 * b ** 6 / 175.d0 
          fu1 = -0.5d0 + 1.5d0 * ( 1.d0 - 2 * b ** 2 / 3.d0 + &
                                  8.d0 * b ** 4 / 15.d0 - &
                                  16.d0 * b ** 6 / 35.d0 )
        end if
      end if

      fe2 = 1.d0
      fu2 = 1.d0
      ex2 = ee * xn2 ** trd
      ux2 = ex2 * ft
      if ( rel .gt. 0.5d0 ) then
        d = bb * ( 8.d0 * pi / 3.d0 ) ** trd
        b = d * xn2 ** trd
        e = sqrt( 1.d0 + b ** 2 )
        if ( b .gt. 0.0001d0 ) then
          fe2 =  1.0d0 - 1.5d0 * ( b*e - log(b+e) ) ** 2 / b ** 4
          fu2 = -0.5d0 + 1.5d0 * log( b + e ) / ( b * e )
        else
          fe2 = 1.d0 - 2.d0 * b ** 2 / 3.d0 + 2.d0 * b ** 4 / 5.d0 &
               - 48.d0 * b ** 6 / 175.d0 
          fu2 = -0.5d0 + 1.5d0 * ( 1.d0 - 2 * b ** 2 / 3.d0 + &
                                  8.d0 * b ** 4 / 15.d0 - &
                                  16.d0 * b ** 6 / 35.d0 )
        end if
      end if

!  get final functional, potential.

      nex = xn1 * fe1 * ex1 + xn2 * fe2 * ex2
      ux1 = fu1 * ux1
      ux2 = fu2 * ux2
      return
      end subroutine getx
!----------------------------------------------------------------------
      subroutine mkkbfile( nel, nl, xnj, ev, dl, nr, r, rsqd, cq,& 
                          ldim, vi, orb, zorig, njrc )
      use dprec
      use PARAMS
      implicit none
      integer  :: nel, nr, ldim
      real(dp) :: dl, zorig
      integer  :: nl( nel ), njrc( 4 )
      real(dp) :: xnj( nel ), ev( nel ), r( nr ), cq( nr ), rsqd(nr)
      real(dp) :: vi( ldim, 5 ), orb( ldim, nel )

      integer, allocatable  :: twoj( : )
      real(dp), allocatable :: vps( : ), xm1( : ), xm2( : )
      real(dp), allocatable :: phi( : , : ), psi( : )

      integer  :: i, j, k, n, iu, flag, neff, jrt, ief, nn
      real(dp) :: e, de, el, eh, rstop, kappa, zeff, rel, plead, x

      allocate( twoj( nel ), vps( nr ), xm1( nr ), xm2( nr ) )
      allocate( phi( nr, -1 : 1 ), psi( nr ) )

      read (inpunit,*) de,rstop,el,eh,n
      jrt=int(1.0000001d0+log(rstop/r(1))/dl)  !JR
      do i=1,nel
        twoj(i)=int(0.0001d0+2.d0*abs(xnj(i)))  !JR
      end do
      call gropen(13)
      write (13,'(1x,2i5,1f20.10)') nel
      write (13,'(1x,2i5,1f20.10)') (nl(i),twoj(i),ev(i),i=1,nel)
      write (13,'(1x,2i5,1f20.10)') jrt,nr,de
      write (13,'(1x,4f20.15)') r(1),r(2)
      do i=1,nr
        write (13,'(1x,4(1x,1d19.12))') (vi(i,k),k=1,5),cq(i)
        write (13,'(1x,4(1x,1d19.12))') (orb(i,j),j=1,nel)
      end do
      do k=1,nel
        iu=40+k
        call gropen(iu)
        kappa=-1.d0
        zeff=0.d0
        rel=0.d0
        flag=1
        neff=1000
        call setqmm(orb(1,k),nl(k),xnj(k),flag,vps,zeff,&
                   zorig,rel,nr,r,dl,xm1,xm2,njrc,vi)
        call getplead(nl(k),xnj(k),rel,kappa,plead,zeff)
        do i=-1,1
          do j=1,nr
            phi(j,i)=0.d0
          end do
          e=ev(k)+de*real(i,dp)
          call integ(e,nl(k),kappa,neff,nn,jrt,ief,x,phi(1,i),&
                    zeff,vps,xm1,xm2,nr,r,rsqd,dl,rel,plead)
        end do
        do j=1,nr
          write (13,'(1x,1i5,3(2x,1d19.12))') j,(phi(j,i),i=-1,1)
        end do
        call getplead(nl(k),xnj(k),rel,kappa,plead,zeff)
        do i=0,n
          e=el+(eh-el)*real(i,dp)/real(n,dp)
          call integ(e,nl(k),kappa,neff,nn,jrt,ief,x,psi,&
                    zeff,vps,xm1,xm2,nr,r,rsqd,dl,rel,plead)
          write (iu,'(1x,2f20.10)') x, e
        end do
      end do
      return
      end subroutine mkkbfile
!----------------------------------------------------------------------
      subroutine mkvaltab( nr, r, dl, phi1, phi2, v, l )
      use dprec
      implicit none
      integer  :: nr, l
      real(dp) :: dl, r( nr ), phi1( nr ), phi2( nr ), v( nr )
      integer  :: i, j
      real(dp) :: s, pref
      real(dp), allocatable :: rnum( : ), rden( : )
      real(dp), allocatable :: fint( : ), fext( : )
      real(dp), allocatable :: int( : ), ext( : )

      allocate( rnum( nr ), rden( nr ), fint( nr ), fext( nr ) )
      allocate( int( nr ), ext( nr ) )
      pref = dl / 45.d0
      do i = 1, nr
        int( i ) = 0.d0
        ext( i ) = 0.d0
      end do
      do i = 1, nr
        rnum( i ) = r( i ) ** l
        rden( i ) = r( i ) ** ( l + 1 )
        if (abs(rden( i ))>1.d-12) rden( i ) = 1.d0 / rden( i )  !JR
        fint( i ) = pref * rnum( i ) * r( i ) * phi1( i ) * phi2( i )
        fext( i ) = pref * rden( i ) * r( i ) * phi1( i ) * phi2( i )
      end do
      do i = 5, 8
        s = 0.d0
        do j = i, nr, 4
          s = s + 14.d0 * ( fint( j - 4 ) + fint( j ) ) + &
                 64.d0 * ( fint( j - 3 ) + fint( j - 1 ) ) + &
                 24.d0 * fint( j - 2 )
          int( j ) = s
        end do
      end do
      do i = nr - 7, nr - 4
        s = 0.d0
        do j = i, 1, -4
          s = s + 14.d0 * ( fext( j + 4 ) + fext( j ) ) + &
                 64.d0 * ( fext( j + 3 ) + fext( j + 1 ) ) + &
                 24.d0 * fext( j + 2 )
          ext( j ) = s
        end do
      end do
      do i = 1, 4
        ext( i ) = ext( 5 )
      end do
      do i = nr - 3, nr
        int( i ) = int( nr - 4 )
      end do
      do i = 1, nr
        v( i ) = rnum( i ) * ext( i ) + rden( i ) * int( i )
      end do
      deallocate( rnum, rden, fint, fext, int, ext )
      return
      end subroutine mkvaltab
!----------------------------------------------------------------------
      function sig(corpol,ru,nr,r,dr,phi)
      use dprec
      implicit none
      integer  :: nr,i
      real(dp) :: r(nr),dr(nr),phi(nr),corpol,ru
      real(dp) :: zero,one,half,four,sig,xiu,pref,rsub,f,vcpp

      parameter(zero=0.,one=1.,half=0.5,four=4.)
      sig=zero
      xiu=one/ru
      pref=half*corpol
      do 100 i=1,nr
        rsub=r(i)*xiu
        f=(one-exp(-rsub*rsub))
        vcpp=pref*f*f*f*f/(r(i)**four)
        sig=sig+dr(i)*phi(i)*phi(i)*vcpp 
  100 continue
      return
      end function sig
!----------------------------------------------------------------------
      subroutine sgetc(xn1,xn2,ec,nec)
      use dprec
      implicit none
      real(dp) :: nec,xn,xn1,xn2,pi,rs,z,x,ecp,ecf,eca,vcp,vcf,vca,ecud
      real(dp) :: vcd,vcu,ec
!     ceperly-alder (vw)
!     
!     The Vosko-Wilk-Nusair parameterization is used.
!     See Can. J.  Phys. 58, 1200 (1980) and
!     Phys. Rev. B 22, 3812 (1980)
!
!     The vwn* subroutines are courtesy of Mark Stiles

      xn=xn1+xn2
      pi=4.d0*atan(1.d0)
      rs=(3.d0/(4.d0*pi*xn))**(1.d0/3.d0)
      z=(xn1-xn2)/xn
      x = sqrt(rs)
      call vwncop(x,ecp,vcp)
      call vwncof(x,ecf,vcf)
      call vwncoa(x,eca,vca)
!     
      call vwnmix(z,ecp,ecf,eca,vcp,vcf,vca,ecud,vcd,vcu)
      ec=ecud
      nec=xn*ec
      return
      end subroutine sgetc
!----------------------------------------------------------------------
      subroutine vwncop(x,ec,vc)
!     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
!     correlation energy for the paramagnetic limit
!     Can. J.  Phys. 58, 1200 (1980)
!     Phys. Rev. B 22, 3812 (1980)
!     on input x is the square root of rs
      use dprec
      implicit none
      real(dp) :: x2,x,xox
      real(dp) :: taninq,xxb,ec,vc
      real(dp) :: one,two,three,third,a,b,c,q,x0,c1,c2,c3
      
      parameter (one=1.0d0)
      parameter (two=2.0d0)
      parameter (three=3.0d0)
      parameter (third=one/three)
      parameter (a=0.0310907d0)
      parameter (b=3.72744d0)
      parameter (c=12.9352d0)
      parameter (q=6.15199081975908d0)
      parameter (x0=-0.10498d0)
      parameter (c1=two * b / q )
      parameter (c2=two * ( b + two * x0 ) / q )
      parameter (c3=b * x0 / ( c + b * x0 + x0**2 ))
      
      x2= x*x
      xox= x2 + b * x + c
      taninq= atan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= log( x2 / xox ) &
          + c1 * taninq &
          - c3 * ( log( (x-x0)**2 / xox ) &
                 + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
      return
      end subroutine vwncop
!----------------------------------------------------------------------
      subroutine vwncof(x,ec,vc)
!     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
!     correlation energy for the ferromagnetic limit
!     Can. J.  Phys. 58, 1200 (1980)
!     Phys. Rev. B 22, 3812 (1980)
!     on input x is the square root of rs
      use dprec
      implicit none
      real(dp) :: x2,x,xox
      real(dp) :: taninq,xxb,ec,vc
      real(dp) :: one,two,three,third,a,b,c,q,x0,c1,c2,c3

      parameter (one=1.0d0)
      parameter (two=2.0d0)
      parameter (three=3.0d0)
      parameter (third=one/three)
      parameter (a=.01554535d0)
      parameter (b=7.06042d0)
      parameter (c=18.0578d0)
      parameter (q=4.73092690956011d0)
      parameter (x0=-.325d0)
      parameter (c1=two * b / q )
      parameter (c2=two * ( b + two * x0 ) / q )
      parameter (c3=b * x0 / ( c + b * x0 + x0**2 ))
      
      x2= x*x
      xox= x2 + b * x + c
      taninq= atan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= log( x2 / xox ) &
          + c1 * taninq &
          - c3 * ( log( (x-x0)**2 / xox ) &
                 + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
      return
      end subroutine vwncof
!----------------------------------------------------------------------
      subroutine vwncoa(x,ec,vc)
!     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
!     correlation contribution to the spin stiffness in the 
!     paramagnetic limit
!     Can. J.  Phys. 58, 1200 (1980)
!     Phys. Rev. B 22, 3812 (1980)
!     on input x is the square root of rs
      use dprec
      implicit none
      real(dp) :: x2,x,xox
      real(dp) :: taninq,xxb,ec,vc
      real(dp) :: one,two,three,third,a,b,c,q,x0,c1,c2,c3

      parameter (one=1.0d0)
      parameter (two=2.0d0)
      parameter (three=3.0d0)
      parameter (third=one/three)
      parameter (a=-.01688685d0)
      parameter (b=1.13107d0)
      parameter (c=13.0045d0)
      parameter (q=7.12310891781812d0)
      parameter (x0=-.0047584d0)
      parameter (c1=two * b / q )
      parameter (c2=two * ( b + two * x0 ) / q )
      parameter (c3=b * x0 / ( c + b * x0 + x0**2 ))
      
      x2= x*x
      xox= x2 + b * x + c
      taninq= atan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= log( x2 / xox ) &
          + c1 * taninq &
          - c3 * ( log( (x-x0)**2 / xox ) &
                 + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
      return
      end subroutine vwncoa
!--------------------------------------------------------
      subroutine vwnmix(pol,ecp,ecf,eca,vcp,vcf,vca,ec,vcd,vcu)
!     mixing between paramagnetic limits and ferromagnetic limits
!     Vosko - Wilk - Nusair
!     Can. J.  Phys. 58, 1200 (1980)
!     Phys. Rev. B 22, 3812 (1980)
      use dprec
      implicit none
      real(dp) :: fup,pol
      real(dp) :: ec,vc,fdn,fupth,fdnth,fmix,dfmix,pol3,pol4,fmpol4,dfmp4
      real(dp) :: ecp,ecf,eca,vcpol,vcp,vcf,vca,vcd,vcu
      real(dp) :: one,two,three,four,third,fthird,cmix1,cfppin

      parameter (one=1.0d0)
      parameter (two=2.0d0)
      parameter (three=3.0d0)
      parameter (four=4.0d0)
      parameter (third=one/three)
      parameter (fthird=four*third)
      parameter (cmix1=1.92366105093154d0)
      parameter (cfppin=0.5848223622634647d0)
      
      fup= one + pol
      fdn= one - pol
      fupth= fup**third
      fdnth= fdn**third
      fmix= ( fup*fupth + fdn*fdnth - two ) * cmix1
      dfmix= fthird * ( fupth - fdnth ) * cmix1

      pol3= pol**3
      pol4= pol3 * pol
      fmpol4= fmix * pol4
      dfmp4= four * pol3 * fmix + pol4 * dfmix

      ec=    ecp * (  one - fmpol4 ) &
          + ecf *          fmpol4 &
          + eca * ( fmix - fmpol4 ) * cfppin
      vcpol= &
          - ecp * dfmp4 &
          + ecf * dfmp4 &
          + eca * ( dfmix - dfmp4 ) * cfppin
      vc=  - pol * vcpol &
          + vcp * (  one - fmpol4 ) &
          + vcf *          fmpol4 &
          + vca * ( fmix - fmpol4 ) * cfppin

      vcd = (vc + vcpol)
      vcu = (vc - vcpol)
      return
      end subroutine vwnmix

