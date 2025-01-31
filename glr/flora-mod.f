      module ArraySizes 
        implicit none
      
        ! Declare main parameters
        integer :: ix = 150        ! Number of axial grid points
        integer :: jx = 50         ! Number of radial (psi) grid points
        integer :: npob = 200      ! Number of p(B) points
        integer :: ksimp = 23      ! Number of points for Simpson's Rule integration
        integer :: nmax = 5        ! Upper limit of time steps for problem
        integer :: nps = 100       ! (Explanation needed for nps)
        integer :: neng = 500      ! (Explanation needed for neng)
        
        ! Secondary array sizes
        integer :: kxx = 0         ! Secondary array size, kxx = (ix-2)*(jx-2)
        integer :: kbw = 0         ! Secondary array size, kbw = min(ix-1,jx-1)
        integer :: lda = 0         ! Secondary array size, lda = 3*kbw+1
      
        namelist /flora_ArraySizes/ nmax

      end module ArraySizes
      
      module Splines
        implicit none
      
        ! Define the size of the arrays
        integer :: mspin
        integer :: mspout
      
        ! Declare the arrays
        integer, allocatable :: lsp(:)      ! Integer array of size mspout
        real, allocatable :: dsp(:)         ! Real array of size mspout
        real, allocatable :: csp1(:)        ! Real array of size mspin
        real, allocatable :: csp2(:)        ! Real array of size mspin
        real, allocatable :: csp3(:)        ! Real array of size mspin
        real, allocatable :: csp4(:)        ! Real array of size mspin

        contains

          ! Allocate subroutine
          subroutine allocate_arrays_Splines(mspin_in, mspout_in)
            integer, intent(in) :: mspin_in, mspout_in
        
            ! Assign values from input arguments
            mspin = mspin_in
            mspout = mspout_in
        
            ! Allocate arrays based on the given dimensions
            allocate(lsp(mspout))
            allocate(dsp(mspout))
            allocate(csp1(mspin))
            allocate(csp2(mspin))
            allocate(csp3(mspin))
            allocate(csp4(mspin))
        
          end subroutine allocate_arrays_Splines
      
      end module Splines
      
      
      module Const_1
        implicit none
      
        ! Declare all the real variables
        real :: startime          ! Initial value from 'call second' (Cray).
        real :: endtime           ! Final value from 'call second' (Cray).
        real :: totalcpu          ! Total timeloop CPU from 'call second' (Cray/Sun)
        real :: t                 ! Time-related variable
        
        real :: gam1
        real :: gam2
        real :: fac1
        real :: fac2
        real :: bias = 0.5       ! Time centering parameter, ==0 fully centered, ==1 fully forward bias.
        
        real :: du = 1.0
        real :: dv = 1.0
        real :: dt = 1.0e-5      ! Time step.
        real :: ex0 = 0.0        ! Initial perturbation coefficient, ==1 random perturbation.
        real :: ex1 = 0.0        ! Initial perturbation coefficient, ==1 for cosine cosine distribution.
        
        real :: fi1 = 1.0        ! Minimum z boundary condition, ==1 for 0 slope, ==-1 for 0 value.
        real :: fizx = 1.0       ! Maximum z boundary condition, ==1 for 0 slope, ==-1 for 0 value.
        real :: fj1 = 1.0        ! Minimum psi boundary condition, ==1 for 0 slope, ==-1 for 0 value.
        real :: fjrx = 1.0       ! Maximum psi boundary condition, ==1 for 0 slope, ==-1 for 0 value.
        
        real :: fpsi = 0.5       ! Grid stretching parameter [0,1].
        real :: fu = 0.5         ! Grid stretching parameter [0,1].
        real :: fv = 0.5         ! Grid stretching parameter [0,1].
        real :: fz = 0.5         ! Grid stretching parameter [0,1].
        
        real :: azm = 1.0
        real :: apsim = 1.0
        
        real :: u0
        real :: v0
        real :: amass = 3.34e-24
        real :: omegexb
        real :: sf6 = 1.0        ! Arbitrary scaling factor on gyroscopic FLR term xxx.
        real :: sf8 = 1.0        ! Arbitrary scaling factor on quasi-elastic FLR term yyy.
        real :: zedge
        
        real :: xu
        real :: xv
        real :: vw
        real :: psiw
        real :: dvin
        real :: dvout

        namelist /flora_Const_1/ bias, dt, ex0, ex1, fi1, fizx, fj1,
     *   fjrx, fpsi, fu, fv, fz, sf6, sf8
      
      end module Const_1
      
      
      module Const_1a
        implicit none
      
        ! Declare the integer variables
        integer :: iscreen = 1           ! Tells how often to write timestep+ to screen.
        integer :: mm = 4                ! Azimuthal mode number.
        integer :: lmax = 2              ! Iteration limit.
        integer :: isw                   ! Not initialized, may need to be set in program.
        integer :: ndiag = 100           ! Number of time steps between diagnostic plots.
        integer :: kplot                 ! Not initialized, may need to be set in program.
        integer :: npm                   ! Not initialized, may need to be set in program.
        integer :: kplotm = 0            ! Index of spatial location of time history plots.
        integer :: kzs = 1               ! Flute mode initialization.
        integer :: n                     ! Time loop counter, initialized as needed.
        integer :: ltt                   ! Not initialized, may need to be set in program.
        integer :: nen = 1               ! Number of time steps between energy checks.
        integer :: lee                   ! Not initialized, may need to be set in program.
      
        namelist /flora_Const_1a/ mm, lmax, ndiag, kzs, nen
      end module Const_1a
      
      
      module Const_2
        implicit none
      
        ! Declare real variables for constants
        real :: psi0rel = 0.50
        real :: z1rel = 0.5
        real :: z2rel = 1.0
        real :: z3rel = 0.625
        real :: z0rel = 0.75
        real :: nslosh
        real :: nsloshks
        real :: bmg = 1.0e4
        real :: ncenter = 1.0e13
        real :: pslosh
        real :: psloshks
        real :: pcenter
        real :: ztrans = 0.4
        real :: ltrans = 0.05
        real :: bm
        real :: ltran
        real :: ztran
        real :: rp1 = 2.0
        real :: rp1ks = 100.0
        real :: bcen
        real :: pring
        real :: epsp = 1.0e-6
        real :: phicen = 0.1
        real :: phiplg = 0.1
        real :: xpot = 1.0
        real :: ypot = 2.0
        real :: wpot = 2.0
        real :: pfudge = 0.0
        real :: pfudgeks = 0.0
        real :: rpx
        real :: rpxks = 0.0
        real :: phice
        real :: phipl
        real :: betslsh = 0.25
        real :: betslsks = 0.0
        real :: betcent = 0.10
        real :: z0
        real :: z1
        real :: z2
        real :: z3
        real :: psi0
        real :: betcene = 0.1
        real :: betslse = 0.1
        real :: psloshe
        real :: pcentee
        real :: bmax
        real :: alsi
        real :: alsiks
        real :: bm1
        real :: bm1ks
        real :: bvoutks
        real :: psls1
        real :: cold = 1.0e-5
        real :: alfcold = 1.0
        real :: p2wide = 0.1
        real :: psi3rel = 1.05
        real :: psi3
        real :: p1max
        real :: bv0
        real :: bv3
        real :: bceng = 1.0e3
        real :: psloshin
        real :: poshinks
        real :: psloshen
        real :: nsloshin = 2.0e13
        real :: noshinks = 0.0
        real :: p3a
        real :: p3b
        real :: p3c
        real :: p3d
        real :: psim
        real :: pe10 = 0.0
        real :: ae1
        real :: be1
        real :: ce1
        real :: de1
        real :: psi0erel = 0.5
        real :: psi0e
        real :: psime
        real :: p2ewide = 0.2
        real :: wp2e
        real :: p1floor = 0.0
        real :: p2flag = 1.0e-3
        real :: p2floor = 0.0
        real :: fring = 0.9
        real :: dip = 1.0
        real :: psistr
        real :: psislp = 0.01
        real :: psihrel = 0.0
        real :: psih
        real :: dpsihrel = 0.0
        real :: dpsih

        namelist /flora_Const_2/ psi0rel, ncenter, ztrans, ltrans, rp1, 
     *    rp1ks, phicen, phiplg, xpot, ypot, pfudge, pfudgeks, rpxks, 
     *    betslsh, betslsks, betcent, betcene, betslse, cold, alfcold, 
     *    p2wide, psi3rel, bceng, nsloshin, noshinks, pe10, psi0erel, 
     *    p2ewide, p2flag, p2floor, fring, dip, psihrel, dpsihrel,
     *    z0rel, z1rel, z2rel, z3rel
      
      end module Const_2
      
      module Const_3
        implicit none
      
        real :: at0, bt0, ct0, cp0
        real :: ap1, ap2, ap3
        real :: at1, at2, at3
        real :: bp1, bp2, bp3
        real :: bt1, bt2, bt3
        real :: cp1, cp2, cp3
        real :: ct1, ct2, ct3
        real :: ct2ks
        real :: ct3ks
        real :: bmx1 = 0.0
        real :: bmx2 = 0.0
        real :: bmx3 = 0.0
        real :: bmn1
        real :: bmn2
        real :: ppas1
        real :: ppas2 = 0.0
        real :: ppas3 = 0.0
        real :: p1trap
        real :: z1min
        real :: z2min
        real :: dpas1 = 0.0
        real :: d1trap = 0.0
        real :: betrap = 0.0
        real :: betpas1 = 0.0
        real :: bvx2
        real :: bvx3
        real :: z2ct

        namelist /flora_Const_3/ bmx1, bmx2, bmx3, ppas2, ppas3, dpas1,
     *        d1trap, betrap, betpas1
      
      end module Const_3
      
      module Const_4
        implicit none
      
        real, dimension(:), allocatable :: psi, z, u, v, dpsi, dz, vpsi, uuz
        real, dimension(:), allocatable :: vpsih, uuzh

        contains

        ! Subroutine to allocate arrays (set array size here)
        subroutine allocate_arrays_Const_4(size)
          integer, intent(in) :: size  ! Example: size of the arrays
      
          ! Allocate arrays based on the given size
          allocate(psi(size))
          allocate(z(size))
          allocate(u(size))
          allocate(v(size))
          allocate(dpsi(size))
          allocate(dz(size))
          allocate(vpsi(size))
          allocate(uuz(size))
          allocate(vpsih(size))
          allocate(uuzh(size))
      
          ! Initialize arrays to zero (optional, based on your needs)
          psi = 0.0
          z = 0.0
          u = 0.0
          v = 0.0
          dpsi = 0.0
          dz = 0.0
          vpsi = 0.0
          uuz = 0.0
          vpsih = 0.0
          uuzh = 0.0
        end subroutine allocate_arrays_Const_4
      
      end module Const_4
      
      
      module Const_5
        implicit none
      
        real, dimension(:), allocatable :: xrtime, time, enpot, tenergy, enkin
        real, dimension(:,:), allocatable :: xrspz, xrsppsi, xflute 
        real, dimension(:), allocatable :: timengy, tenrel, enbend, encurve, enflr

      contains
      
        ! Subroutine to allocate arrays in Const_5
        subroutine allocate_arrays_Const_5(nmax, ix, jx, nps, neng)
          implicit none
      
          integer, intent(in) :: nmax, ix, jx, nps, neng  ! Input dimensions
      
          ! Allocate arrays with specified sizes
          allocate(xrtime(nmax))
          allocate(xrspz(ix, nps))
          allocate(xrsppsi(jx, 2 * nps))
          allocate(time(nmax))
          allocate(xflute(ix, nps))
          allocate(enpot(neng))
          allocate(tenergy(neng))
          allocate(enkin(neng))
          allocate(timengy(neng))
          allocate(tenrel(neng))
          allocate(enbend(neng))
          allocate(encurve(neng))
          allocate(enflr(neng))
      
        end subroutine allocate_arrays_Const_5
      
      end module Const_5
      
      module Const_6
        implicit none
      
        ! Declare scalar constants
        real :: lb = 1.0e9
        real :: rw = 10.0         ! Wall radius [cm], user input
        real :: cee = 3.0e10     ! Another constant

        namelist /flora_Const_6/ rw
      
      end module Const_6
      
      module Const_7
        implicit none
      
        real, dimension(:), allocatable :: h12, hzt0, h34
        real, dimension(:), allocatable :: abp, bbp, cbp, abf, bbf, cbf
        real, dimension(:), allocatable :: hp3, htrans, abq, bbq, cbq
        real, dimension(:), allocatable :: ebp, fbp, gbp
        real, dimension(:), allocatable :: hp0, hpm, hpme, hflr
        real, dimension(:), allocatable :: hzp0, hzp1, hzp2, hzp3
        real, dimension(:), allocatable :: hzt1, hzt2, hzt3, hzt2ks, hzt3ks
        real, dimension(:), allocatable :: dterjb
        real :: betring = 0.0    ! Peak hot electron perpendicular beta, wrt bmn1

        contains

          ! Subroutine to allocate arrays in Const_7
          subroutine allocate_arrays_Const_7(ix, jx)
            implicit none
        
            integer, intent(in) :: ix, jx  ! Input dimensions
        
            ! Allocate arrays with specified sizes
            allocate(h12(ix))
            allocate(hzt0(ix))
            allocate(h34(ix))
            allocate(abp(ix))
            allocate(bbp(ix))
            allocate(cbp(ix))
            allocate(abf(ix))
            allocate(bbf(ix))
            allocate(cbf(ix))
            allocate(hp3(jx))
            allocate(htrans(ix))
            allocate(abq(ix))
            allocate(bbq(ix))
            allocate(cbq(ix))
            allocate(ebp(ix))
            allocate(fbp(ix))
            allocate(gbp(ix))
            allocate(hp0(jx))
            allocate(hpm(jx))
            allocate(hpme(jx))
            allocate(hflr(jx))
            allocate(hzp0(ix))
            allocate(hzp1(ix))
            allocate(hzp2(ix))
            allocate(hzp3(ix))
            allocate(hzt1(ix))
            allocate(hzt2(ix))
            allocate(hzt3(ix))
            allocate(hzt2ks(ix))
            allocate(hzt3ks(ix))
            allocate(dterjb(ix))
        
          end subroutine allocate_arrays_Const_7
      
      end module Const_7
      
      
      module Const_8
        implicit none
      
        real, allocatable, dimension(:) :: bvac, dbvdz, d2bvdz2, d3bvdz3, bint
        real, allocatable, dimension(:) :: dp1dpsi, p1
        real, allocatable, dimension(:,:) :: p1k
        real, allocatable, dimension(:) :: hpk0, deli1, deli2, deli3, deli4
        real, allocatable, dimension(:,:) :: rzz, qubv
        real, allocatable, dimension(:) :: rzz1a, rzz1b, rzz2, rzz3
        real, allocatable, dimension(:) :: rzz1anew, rzz1bnew, rzz2new, rzz3new, rzznew, rzzold
        real, allocatable, dimension(:) :: rd2diff, rnew, rnewd2, rold, roldd2
        real, allocatable, dimension(:) :: icapi2, icapi3, icapi5
        real, allocatable, dimension(:) :: icap2new, icap3new, icap5new
        real, allocatable, dimension(:) :: phi1, phi2
        real, allocatable, dimension(:,:) :: pperp, ppar, dbdpsi
        real, allocatable, dimension(:) :: jperpz, jperpb, jparb, jdfb, jd2fb
        real, allocatable, dimension(:) :: ddbparb, ddbprpb, jperpbpg, jparbpg
        real, allocatable, dimension(:) :: jdfbpg, jd2fbpg, ddbparbpg, ddbprpbpg
        real, allocatable, dimension(:) :: jparz, jd2fz, jdfz, jbtotz
        real, allocatable, dimension(:) :: jmbselfz, ddbprpz, ddbparz
        real :: bvmaxp = 20.0     ! Max B for jperpz
        real :: bvmaxpg = 2.10    ! Max B for jperpz for plug
        real :: bvadj = 1.0       ! Self-B adjustment for btnorm in plug
        real :: alfpg = 0.5       ! Underrelaxation for iteration in plug
        real :: alfks = 0.5       ! Underrelaxation for iteration in ks
        real :: alfsubp = 0.5
        integer :: delsubp = 0     ! Loop index for ks pob2poz
        integer :: isubp = 140
        logical :: dosubp = .true.
        logical :: doflat = .false.
        logical :: doxrlin = .false.
        logical :: doxrtanh = .true.
        logical :: dorzzlb = .false.
        logical :: dowriteb = .false.
        logical :: doplugpob = .false.
        real, allocatable, dimension(:) :: posh, nosh
        real, allocatable, dimension(:) :: bvacnorm
        real, dimension(100) :: bti82
        real, dimension(100) :: btisubp, btisubpm1, btisubpm2, btisubpm3
        real, dimension(100) :: jpzsubp, jpzsubpm1, jpzsubpm2, jpzsubpm3
      
        ! Additional arrays for diagnostic information
        real, allocatable, dimension(:) :: bj2old, pperpold, pperpold25, pparold
        real, allocatable, dimension(:) :: betaold, qubold, qubvold, dbdpsiold
        real, allocatable, dimension(:) :: dflute3, p2, dp2dpsi
        real, allocatable, dimension(:) :: dflute1, dflute2
        real, allocatable, dimension(:) :: flute1, flute2, flute3
        real, allocatable, dimension(:) :: dpost1, dpost2
        real, allocatable, dimension(:) :: dpost1r, dpost2r, dpostjb, dpostjbo
        real :: post1, post2
        real :: post1r, post2r, postjb
        real :: post1sum, post2sum, postjbsum, postjbsumo
        real :: post1cc, post2cc, postjbcc
        real :: post1pg, post2pg, postjbpg
        real :: post1ks, post2ks, postjbks
        real :: p1kssubp, p1ksmax, p2ksmax, pjbksmax, pjbksmaxo
        real :: pnorm
        real :: pperpinjo, btksinjo, btksinj, pparinj, pperpinj
        real, allocatable, dimension(:) :: ptot, beta
        real, allocatable, dimension(:) :: betanorm
        real, allocatable, dimension(:,:) :: p2k, ering
        real, allocatable, dimension(:) :: deli5, del1new, del2new, del3new, del5new
        real, allocatable, dimension(:,:) :: pperps, errprp, errprl, dqubdb
        real, allocatable, dimension(:,:) :: epsi, pperpe
        real, allocatable, dimension(:) :: xxxfreq, omeg1wkb, omeg2wkb
        real, allocatable, dimension(:) :: gamwkb, dflute4, rhoave, xxxave, yyyave
        real, allocatable, dimension(:) :: p2t, dp2dpst, p3, dp3dpsi
        real, allocatable, dimension(:) :: hpkm, hpkme
        real, allocatable, dimension(:) :: droave, droterm
        real :: grow, growmax, xfreqmax, rzzrmax

        namelist /flora_Const_8/  doflat, doxrlin, doxrtanh, dorzzlb, dowriteb

      contains

        ! Subroutine to allocate arrays (you can adjust dimensions as needed)
        subroutine allocate_arrays_Const_8(ix, jx, ksimp, lda)
          implicit none
          integer, intent(in) :: ix, jx, ksimp, lda
      
          ! Allocate the arrays based on the dimensions
          allocate(bvac(ix))
          allocate(dbvdz(ix))
          allocate(d2bvdz2(ix))
          allocate(d3bvdz3(ix))
          allocate(bint(ix))
          allocate(dp1dpsi(jx))
          allocate(p1(jx))
          allocate(p1k(ksimp,jx))
          allocate(hpk0(ksimp))
          allocate(deli1(ksimp))
          allocate(deli2(ksimp))
          allocate(deli3(ksimp))
          allocate(deli4(ksimp))
          allocate(rzz(ix,jx))
          allocate(rzz1a(ix))
          allocate(rzz1b(ix))
          allocate(rzz2(ix))
          allocate(rzz3(ix))
          allocate(rzz1anew(ix))
          allocate(rzz1bnew(ix))
          allocate(rzz2new(ix))
          allocate(rzz3new(ix))
          allocate(rzznew(ix))
          allocate(rzzold(ix))
      
          ! Allocate remaining arrays
          allocate(rd2diff(ix))
          allocate(rnew(ix))
          allocate(rnewd2(ix))
          allocate(rold(ix))
          allocate(roldd2(ix))
          allocate(icapi2(ix))
          allocate(icapi3(ix))
          allocate(icapi5(ix))
          allocate(icap2new(ix))
          allocate(icap3new(ix))
          allocate(icap5new(ix))
          allocate(dbdpsi(ix,jx))
          allocate(dbdpsiold(ix))
          allocate(phi1(ix))
          allocate(phi2(ix))
          allocate(pperp(ix,jx))
          allocate(ppar(ix,jx))
          allocate(jperpz(ix))
          allocate(jperpb(NINPUT))
          allocate(jparb(NINPUT))
          allocate(jdfb(NINPUT))
          allocate(jd2fb(NINPUT))
          allocate(ddbparb(NINPUT))
          allocate(ddbprpb(NINPUT))
      
          ! Continue allocating the remaining arrays
          allocate(jperpbpg(NINPUT))
          allocate(jparbpg(NINPUT))
          allocate(jdfbpg(NINPUT))
          allocate(jd2fbpg(NINPUT))
          allocate(ddbparbpg(NINPUT))
          allocate(ddbprpbpg(NINPUT))
      
          allocate(jparz(ix))
          allocate(jd2fz(ix))
          allocate(jdfz(ix))
          allocate(jbtotz(ix))
          allocate(jmbselfz(ix))
          allocate(ddbprpz(ix))
          allocate(ddbparz(ix))
          allocate(posh(ix))
          allocate(nosh(ix))
          allocate(bvacnorm(ix))

      
          ! Constants with fixed values
          allocate(bj2old(ix))
          allocate(pperpold(ix))
          allocate(pperpold25(ix))
          allocate(pparold(ix))
          allocate(betaold(ix))
          allocate(qubold(ix))
          allocate(qubvold(ix))
      
          ! Arrays for the diagnostic outputs
          allocate(dflute3(ix))
          allocate(qubv(ix,jx))
          allocate(p2(jx))
          allocate(dp2dpsi(jx))
          allocate(dflute1(ix))
          allocate(dflute2(ix))
          allocate(flute1(jx))
          allocate(flute2(jx))
          allocate(flute3(jx))
          allocate(dpost1(ix))
          allocate(dpost2(ix))
          allocate(dpost1r(ix))
          allocate(dpost2r(ix))
          allocate(dpostjb(ix))
          allocate(dpostjbo(ix))
          allocate(ptot(ix))
          allocate(beta(ix))
      
          ! Arrays related to new variables
          allocate(betanorm(ix))
          allocate(p2k(ksimp,jx))
          allocate(deli5(ksimp))
          allocate(del1new(ksimp))
          allocate(del2new(ksimp))
          allocate(del3new(ksimp))
          allocate(del5new(ksimp))
      
          ! Arrays for further diagnostics
          allocate(pperps(ix,jx))
          allocate(errprp(ix,jx))
          allocate(errprl(ix-2,jx-1))
          allocate(dqubdb(ix,jx))
          allocate(xxxfreq(jx))
          allocate(pperpe(ix,jx))
          allocate(epsi(ix,jx))
          allocate(omeg1wkb(jx))
          allocate(omeg2wkb(jx))
          allocate(gamwkb(jx))
          allocate(dflute4(ix))
          allocate(rhoave(jx))
          allocate(xxxave(jx))
          allocate(yyyave(jx))
          allocate(p2t(jx))
          allocate(dp2dpst(jx))
          allocate(p3(jx))
          allocate(dp3dpsi(jx))
          allocate(hpkm(ksimp))
          allocate(hpkme(ksimp))
          allocate(droave(jx))
          allocate(droterm(ix))
          allocate(ering(ix,jx))
      
        end subroutine allocate_arrays_Const_8
      
      end module Const_8
      
      module Btable
        implicit none
      
        integer :: nks = 0  ! Number of points in imported B in KS region
        real, allocatable, dimension(:) :: zks  ! Z positions for imported B(z) in KS region [cm]
        real, allocatable, dimension(:) :: bks  ! B at zks for imported B(z) in KS region [dimensionless]
        
      end module Btable
      
      
      module Coils
        implicit none
      
        integer :: ncoil = 2
        real, dimension(3) :: ass, als, zs, bs
        real :: z1c = 0.0
        real :: z2c = 0.0
        real :: z3c = 0.0

        namelist /flora_coils/ ncoil, ass, als, z1c, z2c, z3c
      
      end module Coils
      
      module Fstor
        implicit none
      
        ! Declare arrays for f1 through f7, g1 through g4, etc.
        real, allocatable, dimension(:,:) :: f1, f2, f3, f4, f5, f7
        real, allocatable, dimension(:,:) :: g1, g2, g3, g4
        real, allocatable, dimension(:,:) :: b, rho, qub, r, yyy, xxx
        
        real, allocatable, dimension(:) :: xioo, xio, xiol
        real, allocatable, dimension(:) :: xroo, xro, xrol
        real, allocatable, dimension(:,:) :: xro_t, xio_t
      
        real, allocatable, dimension(:) :: omstri1, omgbi1, omstri34, omstri66, omstr148
        real, allocatable, dimension(:) :: omebj1, omstrj1
      
        ! Declare scalar constants
        real :: swg1 = 1.0   ! Arbitrary scaling factor on curvature-drive term
        real :: swg2 = 1.0   ! Scaling factor on some line-bending terms
        real :: swg3 = 1.0   ! Scaling factor on other line-bending terms
        real :: swg4 = 1.0   ! Another scaling factor on line-bending terms
      
        real :: exbratio = -1.0   ! omegexb = exbratio*omegstr
        real :: alfrigid = 1.0    ! Larger for faster decay of rigid rotor dpdpsi/p

        contains

          ! Subroutine to allocate arrays in Fstor
          subroutine allocate_arrays_Fstor(ix, jx, kxx, nmax)
            implicit none
        
            integer, intent(in) :: ix, jx, kxx, nmax  ! Input dimensions
        
            ! Allocate arrays with specified sizes
            allocate(f1(ix, jx), f2(ix, jx), f3(ix, jx), f4(ix, jx), f5(ix, jx), f7(ix, jx))
            allocate(g1(ix, jx), g2(ix, jx), g3(ix, jx), g4(ix, jx))
            allocate(b(ix, jx), rho(ix, jx), qub(ix, jx), r(ix, jx), yyy(ix, jx), xxx(ix, jx))
            allocate(xioo(kxx), xio(kxx), xiol(kxx), xroo(kxx), xro(kxx), xrol(kxx))
            allocate(omstri1(jx), omgbi1(jx), omstri34(jx), omstri66(jx), omstr148(jx))
            allocate(omebj1(ix), omstrj1(ix))
            allocate(xro_t(nmax,kxx), xio_t(nmax,kxx))
        
          end subroutine allocate_arrays_Fstor
      
      end module Fstor
      
      module GEnergy
        implicit none
      
        ! Declare arrays as allocatable with unknown size
        real, dimension(:,:), allocatable :: drxi, drxio, drxr, drxro
        real, dimension(:), allocatable :: dxiz2, dxrz2, dum1, dum2
        real, dimension(:,:), allocatable :: psuma, psumb  ! 2D array for psuma
        real, dimension(:), allocatable :: quad, quad1  ! Arrays for quad and quad1
        real, dimension(:), allocatable :: rxi, rxr  ! Arrays for rxi and rxr
        real, dimension(:), allocatable :: xig, xigo, xrg, xrgo  ! Arrays for xig, xigo, xrg, xrgo

        contains

          ! Subroutine to allocate arrays in GEnergy
          subroutine allocate_arrays_GEnergy(ix, jx)
            implicit none
        
            integer, intent(in) :: ix, jx  ! Input dimensions
        
            ! Allocate arrays with specified sizes
            allocate(drxi(ix, jx))
            allocate(drxio(ix, jx))
            allocate(drxr(ix, jx))
            allocate(drxro(ix, jx))
            allocate(dxiz2(ix))
            allocate(dxrz2(ix))
            allocate(dum1(ix))
            allocate(dum2(ix))
            allocate(psuma(jx, 12))
            allocate(psumb(jx, 4))
            allocate(quad(12))
            allocate(quad1(4))
            allocate(rxi(jx))
            allocate(rxr(jx))
            allocate(xig(ix))
            allocate(xigo(ix))
            allocate(xrg(ix))
            allocate(xrgo(ix))
        
          end subroutine allocate_arrays_GEnergy
      
      end module GEnergy
      
      module Gfiducials
        implicit none
      
        integer :: nregions = 0
        integer, dimension(3) :: ibvmin, ibvmax
      
      end module Gfiducials
      
      module Glrexec
        implicit none
        real, dimension(:,:), allocatable :: rhs1b, rhs2b 

        contains

          ! Subroutine to allocate arrays in Glrexec
          subroutine allocate_arrays_Glrexec(kxx)
            implicit none
        
            integer, intent(in) :: kxx  ! Input dimension
        
            ! Allocate arrays with specified sizes
            allocate(rhs1b(kxx, 1))
            allocate(rhs2b(kxx, 1))
        
          end subroutine allocate_arrays_Glrexec
      
      end module Glrexec
      
      
      module Matrix
        implicit none
      
        ! Declare integer and real arrays
        integer, allocatable, dimension(:) :: ipvt  ! Pivot indices for LU decomposition
        real, allocatable, dimension(:,:) :: a1, a2, a3  ! Coefficient matrices
        real, allocatable, dimension(:,:) :: b1  ! Another coefficient matrix
        real, allocatable, dimension(:) :: bc  ! Boundary conditions or vector for right-hand side
        real, allocatable, dimension(:) :: rhs1, rhs2  ! Right-hand side vectors for the system
        real, allocatable, dimension(:,:) :: abar  ! 2D array used in the system (perhaps for LU factorization)

        contains
      
          ! Subroutine to allocate the arrays
          subroutine allocate_arrays_Matrix(kxx, jx, lda)
            implicit none
            integer, intent(in) :: kxx, jx, lda
        
            ! Allocate 1D and 2D arrays with the provided sizes
            allocate(ipvt(kxx))            ! Pivot indices for LU decomposition
            allocate(a1(kxx, 9))           ! Coefficient matrix a1
            allocate(a2(kxx, 9))           ! Coefficient matrix a2
            allocate(a3(kxx, 9))           ! Coefficient matrix a3
            allocate(b1(kxx, 3))           ! Coefficient matrix b1
            allocate(bc(jx))               ! Boundary condition vector bc
            allocate(rhs1(kxx))            ! Right-hand side vector rhs1
            allocate(rhs2(kxx))            ! Right-hand side vector rhs2
            allocate(abar(lda, kxx))       ! 2D array abar
        
          end subroutine allocate_arrays_Matrix
          
        end module Matrix
        
        module TMInput
          implicit none
        
          ! Declare input variables with appropriate types
          integer :: long = 1  ! Switch which sets hot electron z-length as elongated (long==1) or regular (long==0)
          real :: rw1 = 0.0   ! Equal, or slightly less (within one grid cell), to rw, if zero, gets set to rw
          real :: zmax = 0.0  ! Maximum z of the domain
          integer :: info      ! Returned by dgbtrf and dgbtrs
          real :: bsub         ! Value subtracted from bvac
          real :: frbvac = 0.0  ! Allowed range is 0. to 1.-delta
          logical :: dobvsub = .false.  ! If true, subtract frbvac*bvac(ix)
      
          namelist /flora_TMInput/ long, rw1, zmax
        end module TMInput

       module ncInts

          implicit none

          integer(4) :: ncid
          integer(4) :: dim_z, dim_psi, dim_t, dim_kxx
          integer(4) :: var_dt
          integer(4) :: var_mm
          integer(4) :: var_fi1, var_fj1, var_fizx, var_fjrx
          integer(4) :: var_bvac, var_b, var_r, var_rho, var_pperp, var_ppar, var_epsi
          integer(4) :: var_xro, var_xio, var_flute3

          character(len=256) :: input_file = 'in_flora.nml'
          character(len=256) :: output_file = 'flora.nc'

       end module ncInts
      
