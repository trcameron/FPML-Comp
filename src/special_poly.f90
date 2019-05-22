!********************************************************************************
!   SPECIAL_POLY: Test FPML Compensated Polishing on Special Polynomials
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 25 April 2018
!********************************************************************************
! The accuracy of FPML with compensated polishing routine is tested on several
! special polynomial. 
!********************************************************************************
program special_poly
    use fpml_comp
    use mproutines
    implicit none
    ! testing variables
    integer                                     :: deg, j, poly_num
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    integer, parameter                          :: nitmax=30
    logical, dimension(:), allocatable          :: conv
    real(kind=dp), dimension(:), allocatable    :: berr, cond   
    complex(kind=dp), dimension(:), allocatable :: p, roots
    ! NAG variables
    logical, parameter                          :: scal = .false.
    integer                                     :: ifail
    real(kind=dp), allocatable                  :: a(:,:), w(:), z(:,:)
    complex(kind=dp), allocatable               :: zeros(:)
    
    call mpinit
    
    ! Testing: special polynomials
    open(unit=1,file="data_files/special_poly.dat")
    write(1,'(A)') 'Poly No., FPML, FPML-CP, C02AFF'
    do poly_num=27,27
        write(1, '(I10,A)', advance='no') poly_num, ', '
        ! make polynomial
        call make_poly(poly_num, deg, exact_roots, p)
        ! FPML
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        call main(p, deg, roots, berr, cond, conv, nitmax)
        write(1, '(ES15.2, A)', advance='no') max_rel_err(roots,exact_roots,deg), ', '
        deallocate(roots, berr, cond, conv)
        ! FPML-CP
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        call MLPolish(p, deg, roots)
        write(1, '(ES15.2, A)', advance='no') max_rel_err(roots,exact_roots,deg), ', '
        deallocate(roots, berr, cond, conv)
        ! C02AFF
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg), zeros(deg))
        do j=0,deg
            a(1,j) = dble(p(deg+1-j))
            a(2,j) = aimag(p(deg+1-j))
        end do
        ifail = 0
        call c02aff(a,deg,scal,z,w,ifail)
        zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
        write(1, '(ES15.2)') max_rel_err(zeros,exact_roots,deg)
        deallocate(a,w,z,zeros)
        ! deallocate polynomial arrays
        deallocate(exact_roots, p)
    end do
    
    ! close file
    close(1)
    
contains
    !************************************************
    !                       make_poly               *
    !************************************************
    ! Make special polynomial corresponding to poly_num.
    !************************************************
    subroutine make_poly(poly_num, deg, exact_roots, p)
        implicit none
        ! argument variables
        integer, intent(in)                         :: poly_num
        integer, intent(out)                        :: deg
        complex(kind=dp), allocatable, intent(out)  :: exact_roots(:), p(:)
        ! local variables
        integer                                     :: j
        real(kind=dp), parameter                    :: pi = 3.141592653589793_dp
        
        select case (poly_num)
            case(1)
            ! Poly 1: Wilkinson deg 10
            deg = 10
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(2)
            ! Poly 2: Wilkinson deg 15
            deg = 15
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(3)
            ! Poly 3: Wilkinson deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(4)
            ! Poly 4: scaled and shifted Wilkinson deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(dble(-21+2*(j-1))/10d0,0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(5)
            ! Poly 5: reverse Wilkinson deg 10
            deg = 10
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(1d0/dble(j),0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(6)
            ! Poly 6: reverse Wilkinson deg 15
            deg = 15
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(1d0/dble(j),0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(7)
            ! Poly 7: reverse Wilkinson deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(1d0/dble(j),0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(8)
            ! Poly 8: prescribed roots of varying scale deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(2d0**(-10+(j-1)),0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(9)
            ! Poly 9: prescribed roots of varying scale -3 deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots = (/ (cmplx(2d0**(-10+(j-1))-3d0,0,kind=dp), j=1,deg)/)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(10)
            ! Poly 10: Chebyshev polynomial deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            p(1) = 1d0
            p(2) = 0d0
            p(3) = -200d0
            p(4) = 0d0
            p(5) = 6600d0
            p(6) = 0d0
            p(7) = -84480d0
            p(8) = 0d0
            p(9) = 549120d0
            p(10) = 0d0
            p(11) = -2050048d0
            p(12) = 0d0
            p(13) = 4659200d0
            p(14) = 0d0
            p(15) = -6553600d0
            p(16) = 0d0
            p(17) = 5570560d0
            p(18) = 0d0
            p(19) = -2621440d0
            p(20) = 0d0
            p(21) = 524288d0
            exact_roots = (/ (cmplx(cos((2.0_dp*j-1.0_dp)*pi/(2.0_dp*deg)),0,kind=dp), j=1,deg)/)
            case(11)
            ! Poly 11: z^i for i=0,20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            p = (/ (cmplx(1,0,kind=dp), j=1,deg+1)/)
            exact_roots = (/ (cmplx(cos(2*j*pi/21d0),sin(2*j*pi/21d0),kind=dp), j=1,deg)/)
            case(12)
            ! Poly 12: C. Traverso 24 MPSolve
            deg = 24
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-3.52d2, 0d0, kind=dp)
            exact_roots(2) = cmplx(-3.52d2, 0d0, kind=dp)
            exact_roots(3) = cmplx(-2.8371450777d2, -2.9920517772d2, kind=dp)
            exact_roots(4) = cmplx(-2.8371450777d2,  2.9920517772d2, kind=dp)
            exact_roots(5) = cmplx(-2.7867414048d2,  6.1005469197d2, kind=dp)
            exact_roots(6) = cmplx(-2.7867414048d2, -6.1005469197d2, kind=dp)
            exact_roots(7) = cmplx(-2.74892372d2, 0d0, kind=dp)
            exact_roots(8) = cmplx(-2.014171531d2, 0d0, kind=dp)
            exact_roots(9) = cmplx(-1.255366582d2, 0d0, kind=dp)
            exact_roots(10) = cmplx(-9.599999999d1, 0d0, kind=dp)
            exact_roots(11) = cmplx(-8.8692435121d1,  5.5009607430d2, kind=dp)
            exact_roots(12) = cmplx(-8.869243512d1, -5.5009607430d2, kind=dp)
            exact_roots(13) = cmplx(-1.6000000000d1, 0d0, kind=dp)
            exact_roots(14) = cmplx(8.23178509855d1, 0d0, kind=dp)
            exact_roots(15) = cmplx(8.8692435121d1, -5.50096074303d2, kind=dp)
            exact_roots(16) = cmplx(8.8692435121d1,  5.5009607430d2, kind=dp)
            exact_roots(17) = cmplx(1.9293739373d2,  1.60865921259d3, kind=dp)
            exact_roots(18) = cmplx(1.929373937d2, -1.6086592125d3, kind=dp)
            exact_roots(19) = cmplx(2.0141715312d2, 0d0, kind=dp)
            exact_roots(20) = cmplx(2.7489237213d2, 0d0, kind=dp)
            exact_roots(21) = cmplx(7.52d2, 0d0, kind=dp)
            exact_roots(22) = cmplx(7.52d2, 0d0, kind=dp)
            exact_roots(23) = cmplx(9.1106065d2,  1.5722d0, kind=dp)
            exact_roots(24) = cmplx(9.1106065d2, -1.5722d0, kind=dp)
            ! polynomial
            p(1) = -54765291428198020791747503747742749163073958404455022926495744d0
            p(2) = -4052135566767965847649766745769409681058667331648450681896960d0
            p(3) = -31969984081155943263834965670035075493639295858977076674560d0
            p(4) = 575060225471570237690073740639182419333523437771848417280d0
            p(5) = 7337981286595499156409929740830030318565357725459415040d0
            p(6) = 6611223380089859336490797585290455483968982077145088d0
            p(7) = -195514288747757987122118583800597358656801082441728d0
            p(8) = -726907419403715013562762609680450059293446635520d0
            p(9) = 197178719520196724204974332265013056299335680d0
            p(10) = 5968852409133617129605588058090797893943296d0
            p(11) = 16576506891508825500182005531742679597056d0
            p(12) = 23375026506968330494765978581548924928d0
            p(13) = 2206941937668751746514177591607296d0
            p(14) = -75617855277818001758431020580864d0
            p(15) = -204797687173976372829472423936d0
            p(16) = -143150263927579584306872320d0
            p(17) =  20214880144364480233472d0
            p(18) =  453786251090072698880d0
            p(19) =  1265052493274939392d0
            p(20) = -968887355572224d0
            p(21) =  1015406084096d0
            p(22) = -3949133824d0
            p(23) =  3284992d0
            p(24) = -1728d0
            p(25) = 1d0
            case(13)
            ! Poly 13: Mandelbrot 31 MPSolve
            deg = 31
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-1.996376137,0d0,kind=dp)
            exact_roots(2) = cmplx(-1.966773216,0d0,kind=dp)
            exact_roots(3) = cmplx(-1.907280091,0d0,kind=dp)
            exact_roots(4) = cmplx(-1.772892903,0d0,kind=dp)
            exact_roots(5) = cmplx(-1.754877666,0d0,kind=dp)
            exact_roots(6) = cmplx(-1.47601464272,0d0,kind=dp)
            exact_roots(7) = cmplx(-1.284084925525, 4.272688960406d-1,kind=dp)
            exact_roots(8) = cmplx(-1.284084925525,-4.272688960406d-1,kind=dp)
            exact_roots(9) = cmplx(-1.138000666650,-2.403324012620d-1,kind=dp)
            exact_roots(10) = cmplx(-1.138000666650, 2.403324012620d-1,kind=dp)
            exact_roots(11) = cmplx(-1d0,0d0,kind=dp)
            exact_roots(12) = cmplx(-5.968916446451269d-1, 6.629807445770295d-1,kind=dp)
            exact_roots(13) = cmplx(-5.968916446451269d-1,-6.629807445770295d-1,kind=dp)
            exact_roots(14) = cmplx(-2.17526747030511d-1,-1.11445426587329,kind=dp)
            exact_roots(15) = cmplx(-2.17526747030511d-1, 1.11445426587329,kind=dp)
            exact_roots(16) = cmplx(-1.6359826155202d-1, 1.09778064288827,kind=dp)
            exact_roots(17) = cmplx(-1.6359826155202d-1,-1.09778064288827,kind=dp)
            exact_roots(18) = cmplx(-1.225611668766536d-1,-7.4486176661974423d-1,kind=dp)
            exact_roots(19) = cmplx(-1.225611668766536d-1, 7.4486176661974423d-1,kind=dp)
            exact_roots(20) = cmplx(-1.13418655949436d-1,-8.605694725015730d-1,kind=dp)
            exact_roots(21) = cmplx(-1.13418655949436d-1,8.605694725015730d-1,kind=dp)
            exact_roots(22) = cmplx(-1.5570386020902d-2, 1.020497366498289d0,kind=dp)
            exact_roots(23) = cmplx(-1.5570386020902d-2,-1.020497366498289d0,kind=dp)
            exact_roots(24) = cmplx(3.59892739012579001d-1, 6.84762020211812856d-1,kind=dp)
            exact_roots(25) = cmplx(3.59892739012579001d-1,-6.84762020211812856d-1,kind=dp)
            exact_roots(26) = cmplx(3.8900684056977123543d-1,-2.1585065087081910777d-1,kind=dp)
            exact_roots(27) = cmplx(3.8900684056977123543d-1, 2.1585065087081910777d-1,kind=dp)
            exact_roots(28) = cmplx(3.96534570032415023d-1, 6.04181810488988837d-1,kind=dp)
            exact_roots(29) = cmplx(3.96534570032415023d-1,-6.04181810488988837d-1,kind=dp)
            exact_roots(30) = cmplx(4.433256333996235387d-1, 3.729624166628465083d-1,kind=dp)
            exact_roots(31) = cmplx(4.433256333996235387d-1,-3.729624166628465083d-1,kind=dp)
            ! polynomial
            p(1) = 1d0
            p(2) = 1d0
            p(3) = 2d0
            p(4) = 5d0
            p(5) = 14d0
            p(6) = 42d0
            p(7) = 100d0
            p(8) = 221d0
            p(9) = 470d0
            p(10) = 958d0
            p(11) = 1860d0
            p(12)= 3434d0
            p(13) = 6036d0
            p(14) = 10068d0
            p(15) = 15864d0
            p(16) = 23461d0
            p(17) = 32398d0
            p(18) = 41658d0
            p(19) = 49700d0
            p(20) = 54746d0
            p(21) = 55308d0
            p(22) = 50788d0
            p(23) = 41944d0
            p(24) = 30782d0
            p(25) = 19788d0
            p(26) = 10948d0
            p(27) = 5096d0
            p(28) = 1932d0
            p(29) = 568d0
            p(30) = 120d0
            p(31) = 16d0        
            p(32) = 1d0
            case(14)
            ! Poly 14: Mandelbrot 63 MPSolve
            deg = 63
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-1.999095682327018473210d0,0d0,kind=dp)
            exact_roots(2) = cmplx(-1.9918141725491222157325609498622881d0,0d0,kind=dp)
            exact_roots(3) = cmplx(-1.977179587006257387346088520662828616836d0,0d0,kind=dp)
            exact_roots(4) = cmplx(-1.953705894284396245427622199013653238901d0,0d0,kind=dp)
            exact_roots(5) = cmplx(-1.927147709363950262460068188946594278007d0,0d0,kind=dp)
            exact_roots(6) = cmplx(-1.8848035715866817923294780929158396496359d0,0d0,kind=dp)
            exact_roots(7) = cmplx(-1.8323152027512291920848975260425181432293d0,0d0,kind=dp)
            exact_roots(8) = cmplx(-1.76926167027683114607548022863625740038777d0, &
                5.6919500395600315304900187298015859319654d-2,kind=dp)
            exact_roots(9) = cmplx(-1.76926167027683114607548022863625740038777d0, &
                -5.6919500395600315304900187298015859319654d-2,kind=dp)
            exact_roots(10) = cmplx(-1.674066091474787971565296029172325596206403d0,0d0,kind=dp)
            exact_roots(11) = cmplx(-1.5748891397523009698199655524959742837719482d0,0d0,kind=dp)
            exact_roots(12) = cmplx(-1.408446485740072654917577008805998851928020904d0, &
                 -1.36171997304659915684707793608163610038822995d-1,kind=dp)
            exact_roots(13) = cmplx(-1.408446485740072654917577008805998851928020904d0, &
                 1.36171997304659915684707793608163610038822995d-1,kind=dp)
            exact_roots(14) = cmplx(-1.29255806103352208716418470636149411998013630326d0, &
                 4.3819881608663183712973712432734844004535476504d-1,kind=dp)
            exact_roots(15) = cmplx(-1.29255806103352208716418470636149411998013630326d0, &
                 -4.3819881608663183712973712432734844004535476504d-1,kind=dp)
            exact_roots(16) = cmplx(-1.26228728143847254301011194120806575232050489502d0, &
                 4.0810432411269038329016065742601506306041169168d-1,kind=dp)
            exact_roots(17) = cmplx(-1.26228728143847254301011194120806575232050489502d0, &
                 -4.0810432411269038329016065742601506306041169168d-1,kind=dp)
            exact_roots(18) = cmplx(-1.25273588401203794629581100256433997387062287256d0, &
                 -3.4247064788975089386187578687092843396383393805d-1,kind=dp)
            exact_roots(19) = cmplx(-1.25273588401203794629581100256433997387062287256d0, &
                 3.4247064788975089386187578687092843396383393805d-1,kind=dp)
            exact_roots(20) = cmplx(-1.02819385245481759930249745596731843328070508279421d0, &
                 -3.61376517118561592479460832997830315786692639704085d-1,kind=dp)
            exact_roots(21) = cmplx(-1.02819385245481759930249745596731843328070508279421d0, &
                 3.61376517118561592479460832997830315786692639704085d-1,kind=dp)
            exact_roots(22) = cmplx(-6.23532485956252757990016587001026776428072703359878868d-1, &
                 6.81064414225239608090835812686561539088332735217609127d-1,kind=dp)
            exact_roots(23) = cmplx(-6.23532485956252757990016587001026776428072703359878868d-1, &
                 -6.81064414225239608090835812686561539088332735217609127d-1,kind=dp)
            exact_roots(24) = cmplx(-6.2243629504129358796016350694723840189750985673649588591d-1, &
                 4.2487843647562918431157443880525338683545992964599689876d-1,kind=dp)
            exact_roots(25) = cmplx(-6.2243629504129358796016350694723840189750985673649588591d-1, &
                 -4.2487843647562918431157443880525338683545992964599689876d-1,kind=dp)
            exact_roots(26) = cmplx(-5.308278048599427289214772971196026578135170949646890946d-1, &
                 6.682887255592057714440924655647011851367651843270734380d-1,kind=dp)
            exact_roots(27) = cmplx(-5.308278048599427289214772971196026578135170949646890946d-1, &
                 -6.682887255592057714440924655647011851367651843270734380d-1,kind=dp)
            exact_roots(28) = cmplx(-2.72102461488938894219383324518026874585894699621947085d-1, &
                 -8.42364690294128145503155708242929569550778268698265965d-1,kind=dp)
            exact_roots(29) = cmplx(-2.72102461488938894219383324518026874585894699621947085d-1, &
                 8.42364690294128145503155708242929569550778268698265965d-1,kind=dp)
            exact_roots(30) = cmplx(-2.24915951286740054685326255204118310792682454680693d-1, &
                 1.11626015745499183500126825424467009109873946082435d0,kind=dp)
            exact_roots(31) = cmplx(-2.24915951286740054685326255204118310792682454680693d-1, &
                 -1.11626015745499183500126825424467009109873946082435d0,kind=dp)
            exact_roots(32) = cmplx(-2.0728383545566641282413385018667121332401155604017d-1, &
                 1.11748077249496291137377567312207879579746389236127d0,kind=dp)
            exact_roots(33) = cmplx(-2.0728383545566641282413385018667121332401155604017d-1, &
                 -1.11748077249496291137377567312207879579746389236127d0,kind=dp)
            exact_roots(34) = cmplx(-1.7457822113571696945156643266162905020167505710204d-1, &
                 1.07142767145403118922964631021955987671322451961088d0,kind=dp)
            exact_roots(35) = cmplx(-1.7457822113571696945156643266162905020167505710204d-1, &
                 -1.07142767145403118922964631021955987671322451961088d0,kind=dp)
            exact_roots(36) = cmplx(-1.57516053475965356164335109644674141293297577896685d-1, &
                 -1.10900651411360717797175198615475582901468585712356d0,kind=dp) 
            exact_roots(37) = cmplx(-1.57516053475965356164335109644674141293297577896685d-1, &
                 1.10900651411360717797175198615475582901468585712356d0,kind=dp)
            exact_roots(38) = cmplx(-1.274999735463630001995395653459879637298616757217284d-1, &
                 9.874609094894567922074076807929788675642068522522938d-1,kind=dp)
            exact_roots(39) = cmplx(-1.274999735463630001995395653459879637298616757217284d-1, &
                 -9.874609094894567922074076807929788675642068522522938d-1,kind=dp)
            exact_roots(40) = cmplx(-1.42334819203540667677618453136202688358025283954839d-2, &
                 -1.0329147752136441093950134026551104360994260360822540d0,kind=dp)
            exact_roots(41) = cmplx(-1.42334819203540667677618453136202688358025283954839d-2, &
                 1.0329147752136441093950134026551104360994260360822540d0,kind=dp)
            exact_roots(42) = cmplx(-6.98356849626139181796649107548406610452886379651341d-3, &
                 -1.0036038622882895485307049669513531297649273745391915d0,kind=dp)
            exact_roots(43) = cmplx(-6.98356849626139181796649107548406610452886379651341d-3, &
                 1.0036038622882895485307049669513531297649273745391915d0,kind=dp)
            exact_roots(44) = cmplx( 1.4895466603687646529815779208794106185666477731693128d-2, &
                 -8.481487619084165277193311117832376290806619901265058603d-1,kind=dp)
            exact_roots(45) = cmplx( 1.4895466603687646529815779208794106185666477731693128d-2, &
                 8.481487619084165277193311117832376290806619901265058603d-1,kind=dp)
            exact_roots(46) = cmplx( 1.211927861059064863147044434105037593859287800520963579338d-1, &
                 6.1061169221075421167538724415035774824319702690063863369691d-1,kind=dp)
            exact_roots(47) = cmplx( 1.211927861059064863147044434105037593859287800520963579338d-1, &
                 -6.1061169221075421167538724415035774824319702690063863369691d-1,kind=dp)
            exact_roots(48) = cmplx( 3.52482539722363278193253964052161589243593334212239870706d-1, &
                 -6.98337239583330331258141954760484537633150485928512286760d-1,kind=dp)
            exact_roots(49) = cmplx( 3.52482539722363278193253964052161589243593334212239870706d-1, &
                 6.98337239583330331258141954760484537633150485928512286760d-1,kind=dp)
            exact_roots(50) = cmplx( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
                 -1.4474937132163286474711018201298830556966056842762643026975894d-1,kind=dp)
            exact_roots(51) = cmplx( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
                 1.4474937132163286474711018201298830556966056842762643026975894d-1,kind=dp)
            exact_roots(52) = cmplx( 3.76893240379311323690004017968512473363482317941533875341d-1, &
                 6.78568693190448141957540792996773280196881194582788907016d-1,kind=dp)
            exact_roots(53) = cmplx( 3.76893240379311323690004017968512473363482317941533875341d-1, &
                 -6.78568693190448141957540792996773280196881194582788907016d-1,kind=dp)
            exact_roots(54) = cmplx( 3.865391765961580265082930869043677799284877313516569138807d-1, &
                 5.693247113031029032137923571351905081619323911951388853856d-1,kind=dp)
            exact_roots(55) = cmplx( 3.865391765961580265082930869043677799284877313516569138807d-1, &
                 -5.693247113031029032137923571351905081619323911951388853856d-1,kind=dp)
            exact_roots(56) = cmplx( 4.12916024722700479197334566382612257174765142865547121703d-1, &
                 6.148067601433856949545497204007997358291659758563137777616d-1,kind=dp)
            exact_roots(57) = cmplx( 4.12916024722700479197334566382612257174765142865547121703d-1, &
                 -6.148067601433856949545497204007997358291659758563137777616d-1,kind=dp)
            exact_roots(58) = cmplx( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
                 2.267599044353486186978765599716989721202321914603899690444951d-1,kind=dp)
            exact_roots(59) = cmplx( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
                 -2.267599044353486186978765599716989721202321914603899690444951d-1,kind=dp)
            exact_roots(60) = cmplx( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
                 -3.96170128033165002412596877271155937712569079351815707744770d-1,kind=dp)
            exact_roots(61) = cmplx( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
                 3.96170128033165002412596877271155937712569079351815707744770d-1,kind=dp)
            exact_roots(62) = cmplx( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
                 3.47758700883481983632188723200264206004781117755664551397643d-1,kind=dp)
            exact_roots(63) = cmplx( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
                 -3.47758700883481983632188723200264206004781117755664551397643d-1,kind=dp) 
            ! polynomial
            p(1) = 1d0
            p(2) = 1d0
            p(3) = 2d0
            p(4) = 5d0
            p(5) = 14d0
            p(6) = 42d0
            p(7) = 132d0
            p(8) = 365d0
            p(9) = 950d0
            p(10) = 2398d0
            p(11) = 5916d0
            p(12) = 14290d0
            p(13) = 33708d0
            p(14) = 77684d0
            p(15) = 175048d0
            p(16) = 385741d0
            p(17) = 831014d0
            p(18) = 1749654d0
            p(19) = 3598964d0
            p(20) = 7228014d0
            p(21) = 14162220d0
            p(22) = 27049196d0
            p(23) = 50323496d0
            p(24) = 91143114d0
            p(25) = 160617860d0
            p(26) = 275276716d0
            p(27) = 458591432d0
            p(28) = 742179284d0
            p(29) = 1166067016d0
            p(30) = 1777171560d0
            p(31) = 2625062128d0
            p(32) = 3754272037d0
            p(33) = 5193067630d0
            p(34) = 6939692682d0
            p(35) = 8948546308d0
            p(36) = 11120136162d0
            p(37) = 13299362332d0
            p(38) = 15286065700d0
            p(39) = 16859410792d0
            p(40) = 17813777994d0
            p(41) = 17999433372d0
            p(42) = 17357937708d0
            p(43) = 15941684776d0
            p(44) = 13910043524d0
            p(45) = 11500901864d0
            p(46) = 8984070856d0
            p(47) = 6609143792d0
            p(48) = 4562339774d0
            p(49) = 2943492972d0
            p(50) = 1766948340d0
            p(51) = 981900168d0
            p(52) = 502196500d0
            p(53) = 234813592d0
            p(54) = 99582920d0
            p(55) = 37945904d0
            p(56) = 12843980d0
            p(57) = 3807704d0
            p(58) = 971272d0
            p(59) = 208336d0
            p(60) = 36440d0
            p(61) = 4976d0
            p(62) = 496d0
            p(63) = 32d0
            p(64) = 1d0
            case(15)
            ! Poly 15: Jenkins Traub p1(z) with a=1e-8
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(1d-8,0d0,kind=dp)
            exact_roots(2) = cmplx(-1d-8,0d0,kind=dp)
            exact_roots(3) = cmplx(1d0,0d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(16)
            ! Poly 16: Jenkins Traub p1(z) with a=1e-15
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(1d-15,0d0,kind=dp)
            exact_roots(2) = cmplx(-1d-15,0d0,kind=dp)
            exact_roots(3) = cmplx(1d0,0d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(17)
            ! Poly 17: Jenkins Traub p1(z) with a=1e+8
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(1d+8,0d0,kind=dp)
            exact_roots(2) = cmplx(-1d+8,0d0,kind=dp)
            exact_roots(3) = cmplx(1d0,0d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(18)
            ! Poly 18: Jenkins Traub p1(z) with a=1e+15
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(1d+15,0d0,kind=dp)
            exact_roots(2) = cmplx(-1d+15,0d0,kind=dp)
            exact_roots(3) = cmplx(1d0,0d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(19)
            ! Poly 19: Jenkins Traub p3(z) deg 10
            deg = 10
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(1d-1,0d0,kind=dp)
            exact_roots(2) = cmplx(1d-2,0d0,kind=dp)
            exact_roots(3) = cmplx(1d-3,0d0,kind=dp)
            exact_roots(4) = cmplx(1d-4,0d0,kind=dp)
            exact_roots(5) = cmplx(1d-5,0d0,kind=dp)
            exact_roots(6) = cmplx(1d-6,0d0,kind=dp)
            exact_roots(7) = cmplx(1d-7,0d0,kind=dp)
            exact_roots(8) = cmplx(1d-8,0d0,kind=dp)
            exact_roots(9) = cmplx(1d-9,0d0,kind=dp)
            exact_roots(10) = cmplx(1d-10,0d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(20)
            ! Poly 20: Jenkins Traub p3(z) deg 20
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(1d-1,0d0,kind=dp)
            exact_roots(2) = cmplx(1d-2,0d0,kind=dp)
            exact_roots(3) = cmplx(1d-3,0d0,kind=dp)
            exact_roots(4) = cmplx(1d-4,0d0,kind=dp)
            exact_roots(5) = cmplx(1d-5,0d0,kind=dp)
            exact_roots(6) = cmplx(1d-6,0d0,kind=dp)
            exact_roots(7) = cmplx(1d-7,0d0,kind=dp)
            exact_roots(8) = cmplx(1d-8,0d0,kind=dp)
            exact_roots(9) = cmplx(1d-9,0d0,kind=dp)
            exact_roots(10) = cmplx(1d-10,0d0,kind=dp)
            exact_roots(11) = cmplx(1d-11,0d0,kind=dp)
            exact_roots(12) = cmplx(1d-12,0d0,kind=dp)
            exact_roots(13) = cmplx(1d-13,0d0,kind=dp)
            exact_roots(14) = cmplx(1d-14,0d0,kind=dp)
            exact_roots(15) = cmplx(1d-15,0d0,kind=dp)
            exact_roots(16) = cmplx(1d-16,0d0,kind=dp)
            exact_roots(17) = cmplx(1d-17,0d0,kind=dp)
            exact_roots(18) = cmplx(1d-18,0d0,kind=dp)
            exact_roots(19) = cmplx(1d-19,0d0,kind=dp)
            exact_roots(20) = cmplx(1d-20,0d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(21)
            ! Poly 21: Jenkins Traub p4(z)
            deg = 6
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 0.1_dp
            exact_roots(2) = 0.1_dp
            exact_roots(3) = 0.1_dp
            exact_roots(4) = 0.5_dp
            exact_roots(5) = 0.6_dp
            exact_roots(6) = 0.7_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(22)
            ! Poly 22: Jenkins Traub p5(z)
            deg = 10
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 0.1_dp
            exact_roots(2) = 0.1_dp
            exact_roots(3) = 0.1_dp
            exact_roots(4) = 0.1_dp
            exact_roots(5) = 0.2_dp
            exact_roots(6) = 0.2_dp
            exact_roots(7) = 0.2_dp
            exact_roots(8) = 0.3_dp
            exact_roots(9) = 0.3_dp
            exact_roots(10) = 0.4_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(23)
            ! Poly 23: Jenkins Traub p6(z)
            deg = 5
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 0.1_dp
            exact_roots(2) = 1.001_dp
            exact_roots(3) = 0.998_dp
            exact_roots(4) = 1.00002_dp
            exact_roots(5) = 0.99999_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(24)
            ! Poly 24: Jenkins Traub p7(z) with a=0
            deg = 7
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 0.001_dp
            exact_roots(2) = 0.01_dp
            exact_roots(3) = 0.1_dp
            exact_roots(4) = 0.1_dp
            exact_roots(5) = 0.1_dp
            exact_roots(6) = 1.0_dp
            exact_roots(7) = 10.0_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(25)
            ! Poly 25: Jenkins Traub p7(z) with a=10^(-10)
            deg = 7
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 0.001_dp
            exact_roots(2) = 0.01_dp
            exact_roots(3) = 0.1_dp
            exact_roots(4) = cmplx(0.1_dp,1E-10_dp,kind=dp)
            exact_roots(5) = cmplx(0.1_dp,-1E-10_dp,kind=dp)
            exact_roots(6) = 1.0_dp
            exact_roots(7) = 10.0_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(26)
            ! Poly 26: Jenkins Traub p7(z) with a=10^(-6)
            deg = 7
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 0.001_dp
            exact_roots(2) = 0.01_dp
            exact_roots(3) = 0.1_dp
            exact_roots(4) = cmplx(0.1_dp,1E-6_dp,kind=dp)
            exact_roots(5) = cmplx(0.1_dp,-1E-6_dp,kind=dp)
            exact_roots(6) = 1.0_dp
            exact_roots(7) = 10.0_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(27)
            ! Poly 27: Jenkins Traub p8(z)
            deg = 5
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = -1.0_dp
            exact_roots(2) = -1.0_dp
            exact_roots(3) = -1.0_dp
            exact_roots(4) = -1.0_dp
            exact_roots(5) = -1.0_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(28)
            ! Poly 28: Jenkins Traub p9(z)
            deg = 20
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(+1d-2,0d0,kind=dp)
            exact_roots(2) = cmplx(-1d-2,0d0,kind=dp)
            exact_roots(3) = cmplx(0d0,+1d+2,kind=dp)
            exact_roots(4) = cmplx(0d0,-1d+2,kind=dp)
            exact_roots(5) = cmplx(-95.105651629515357212d0,+30.901699437494742410d0,kind=dp)
            exact_roots(6) = cmplx(-95.105651629515357212d0,-30.901699437494742410d0,kind=dp)
            exact_roots(7) = cmplx(-0.0080901699437494742410d0,+0.0058778525229247312917d0,kind=dp)
            exact_roots(8) = cmplx(-0.0080901699437494742410d0,-0.0058778525229247312917d0,kind=dp)
            exact_roots(9) = cmplx(-58.778525229247312917d0,+80.901699437494742410d0,kind=dp)
            exact_roots(10) = cmplx(-58.778525229247312917d0,-80.901699437494742410d0,kind=dp)
            exact_roots(11) = cmplx(-0.0030901699437494742410d0,+0.0095105651629515357212d0,kind=dp)
            exact_roots(12) = cmplx(-0.0030901699437494742410d0,-0.0095105651629515357212d0,kind=dp)    
            exact_roots(13) = cmplx(+95.105651629515357212d0,+30.901699437494742410d0,kind=dp)
            exact_roots(14) = cmplx(+95.105651629515357212d0,-30.901699437494742410d0,kind=dp)
            exact_roots(15) = cmplx(+0.0080901699437494742410d0,+0.0058778525229247312917d0,kind=dp)
            exact_roots(16) = cmplx(+0.0080901699437494742410d0,-0.0058778525229247312917d0,kind=dp)
            exact_roots(17) = cmplx(+58.778525229247312917d0,+80.901699437494742410d0,kind=dp)
            exact_roots(18) = cmplx(+58.778525229247312917d0,-80.901699437494742410d0,kind=dp)
            exact_roots(19) = cmplx(+0.0030901699437494742410d0,+0.0095105651629515357212d0,kind=dp)
            exact_roots(20) = cmplx(+0.0030901699437494742410d0,-0.0095105651629515357212d0,kind=dp)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(29)
            ! Poly 29: Jenkins Traub p10(z) with a=1e+3
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 1E+3_dp
            exact_roots(2) = 1.0_dp
            exact_roots(3) = 1E-3_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(30)
            ! Poly 30: Jenkins Traub p10(z) with a=1e+6
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 1E+6_dp
            exact_roots(2) = 1.0_dp
            exact_roots(3) = 1E-6_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(31)
            ! Poly 31: Jenkins Traub p10(z) with a=1e+9
            deg = 3
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = 1E+9_dp
            exact_roots(2) = 1.0_dp
            exact_roots(3) = 1E-9_dp
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(32)
            ! Poly 32: Jenkins Traub p11(z) with m=15
            deg = 60
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            do j=-14,14
                exact_roots(15+j) = cmplx(cos(j*pi/30.0_dp),sin(j*pi/30.0_dp),kind=dp)
            end do
            do j=15,45
                exact_roots(15+j) = 0.9_dp*cmplx(cos(j*pi/30.0_dp),sin(j*pi/30.0_dp),kind=dp)
            end do
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(33)
            ! Poly 33: Jenkins Traub p11(z) with m=20
            deg = 80
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            do j=-19,19
                exact_roots(20+j) = cmplx(cos(j*pi/40.0_dp),sin(j*pi/40.0_dp),kind=dp)
            end do
            do j=20,60
                exact_roots(20+j) = 0.9_dp*cmplx(cos(j*pi/40.0_dp),sin(j*pi/40.0_dp),kind=dp)
            end do
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(34)
            ! Poly 34: Jenkins Traub p11(z) with m=25
            deg = 100
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            do j=-24,24
                exact_roots(25+j) = cmplx(cos(j*pi/50.0_dp),sin(j*pi/50.0_dp),kind=dp)
            end do
            do j=25,75
                exact_roots(25+j) = 0.9_dp*cmplx(cos(j*pi/50.0_dp),sin(j*pi/50.0_dp),kind=dp)
            end do
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
            case(35)
            ! Poly 35: Cleve's Corner Polynomial
            deg = 6
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(4,0,kind=dp)
            exact_roots(2) = cmplx(2,1,kind=dp)
            exact_roots(3) = cmplx(2,-1,kind=dp)
            exact_roots(4) = cmplx(0.25_dp,0,kind=dp)
            exact_roots(5) = 1/exact_roots(2)
            exact_roots(6) = 1/exact_roots(3)
            ! polynomial
            call RootCoeffCplx(deg,exact_roots,p)
        end select
    end subroutine make_poly
    !************************************************
    !                       max_rel_err             *
    !************************************************
    ! Compute maximum relative error between roots 
    ! and exact roots.
    !************************************************
    function max_rel_err(roots, exact_roots, deg) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:), roots(:)
        ! local variables
        logical                         :: m(deg)
        integer                         :: i, j, k
        real(kind=dp)                   :: delta, err, res, relerr, rmax
        ! intrinsic procedures
        intrinsic                       :: abs, max, huge
        
        ! initialize
        m = .True.
        rmax = huge(1.0_dp)
        res = mu
        ! main loop
        do i=1,deg
            ! find corresponding exact root (yet to be matched)
            err = rmax
            do j=1,deg
                if(m(j)) then
                    delta = abs(exact_roots(j) - roots(i))
                    if(delta < err) then
                        err = delta
                        k = j
                    end if
                end if
            end do
            ! mark corresponding root as matched
            m(k) = .False.
            ! calculate relative error on this root and update max
            ! zero roots give unhelpful relative errors, so revert to absolute error
            if(abs(roots(i))<mu .or. abs(exact_roots(k))<mu) then
                relerr = err
            else
                relerr = err/abs(exact_roots(k))
            end if
            res = max(relerr, res)
        end do
        res = max(res,mu)
    end function max_rel_err
end program special_poly