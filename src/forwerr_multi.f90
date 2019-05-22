!********************************************************************************
!   FORWERR_MULTI: Compare the forward error of FPML, FPML-CP, and C02AFF.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 28 April 2019
!********************************************************************************
!   
!********************************************************************************
program forwerr_multi
    use fpml_comp
    use mproutines
    implicit none
    ! testing variables
    integer                                     :: deg, j, poly_num
    complex(kind=dp), dimension(:), allocatable :: exact_roots
    ! FPML variables
    integer, parameter                          :: nitmax = 30
    logical, allocatable                        :: conv(:)
    real(kind=dp), allocatable                  :: berr(:), cond(:)   
    complex(kind=dp), allocatable               :: p(:), roots(:)
    ! NAG variables
    logical, parameter                          :: scal = .false.
    integer                                     :: ifail
    real(kind=dp), allocatable                  :: a(:,:), w(:), z(:,:)
    complex(kind=dp), allocatable               :: zeros(:)
    
    call mpinit
    
    ! Testing: polynomials with multiple and near multiple roots
    open(unit=1,file="data_files/forwerr_multi.dat")
    write(1,'(A)') 'Poly No., FPML, FPML-P, FPML-CP, C02AFF'
    do poly_num=1,7
        write(1, '(I10,A)', advance='no') poly_num, ', '
        ! make polynomial
        call make_poly(poly_num, deg, exact_roots, p)
        ! FPML
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        call main(p, deg, roots, berr, cond, conv, nitmax)
        write(1, '(ES15.2, A)', advance='no') max_rel_err(roots,exact_roots,deg), ', '
        deallocate(roots, berr, cond, conv)
        ! FPML-P
        allocate(roots(deg))
        call NPolish(p, deg, roots)
        write(1, '(ES15.2, A)', advance='no') max_rel_err(roots,exact_roots,deg), ', '
        deallocate(roots)
        ! FPML-CP
        allocate(roots(deg))
        call MLPolish(p, deg, roots)
        write(1, '(ES15.2, A)', advance='no') max_rel_err(roots,exact_roots,deg), ', '
        deallocate(roots)
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
    ! Make special polynomial corresponding to 
    ! poly_num and deg > 1.
    !************************************************
    subroutine make_poly(poly_num, deg, exact_roots, p)
        implicit none
        ! argument variables
        integer, intent(in)                         :: poly_num
        integer, intent(out)                        :: deg
        complex(kind=dp), allocatable, intent(out)  :: exact_roots(:), p(:)
        ! local variables
        integer                                     :: j
        
        select case (poly_num)
            ! case 1
            case(1)
            ! Mandelbrot Polynomial, theoretical limiting accuracy 2.53E-5
            deg = 63
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-0.199909568232701847321062999922E1_dp,0E0_dp,kind=dp)
            exact_roots(2) = cmplx(-0.199181417254912221573256094986E1_dp,0E0_dp,kind=dp)
            exact_roots(3) = cmplx(-0.197717958700625738734608852066E1_dp,0E0_dp,kind=dp)
            exact_roots(4) = cmplx(-0.195370589428439624542762219901E1_dp,0E0_dp,kind=dp)
            exact_roots(5) = cmplx(-0.192714770936395026246006818895E1_dp,0E0_dp,kind=dp)
            exact_roots(6) = cmplx(-0.188480357158668179232947809292E1_dp,0E0_dp,kind=dp)
            exact_roots(7) = cmplx(-0.183231520275122919208489752604E1_dp,0E0_dp,kind=dp)
            exact_roots(8) = cmplx(-0.176926167027683114607548022864E1_dp,-0.5691950039560031530490018730E-1_dp,kind=dp)
            exact_roots(9) = cmplx(-0.176926167027683114607548022864E1_dp,0.5691950039560031530490018730E-1_dp,kind=dp)
            exact_roots(10) = cmplx(-0.167406609147478797156529602917E1_dp,0E0_dp,kind=dp)
            exact_roots(11) = cmplx(-0.157488913975230096981996555250E1_dp,0E0_dp,kind=dp)
            exact_roots(12) = cmplx(-0.140844648574007265491757700881E1_dp,-0.13617199730465991568470779361E0_dp,kind=dp)
            exact_roots(13) = cmplx(-0.140844648574007265491757700881E1_dp,0.13617199730465991568470779361E0_dp,kind=dp)
            exact_roots(14) = cmplx(-0.1292558061033522087164184706361E1_dp,-0.438198816086631837129737124327E0_dp,kind=dp)
            exact_roots(15) = cmplx(-0.1292558061033522087164184706361E1_dp,0.438198816086631837129737124327E0_dp,kind=dp)
            exact_roots(16) = cmplx(-0.1262287281438472543010111941208E1_dp,-0.408104324112690383290160657426E0_dp,kind=dp)
            exact_roots(17) = cmplx(-0.1262287281438472543010111941208E1_dp,0.408104324112690383290160657426E0_dp,kind=dp)
            exact_roots(18) = cmplx(-0.1252735884012037946295811002564E1_dp,-0.342470647889750893861875786871E0_dp,kind=dp)
            exact_roots(19) = cmplx(-0.1252735884012037946295811002564E1_dp,0.342470647889750893861875786871E0_dp,kind=dp)
            exact_roots(20) = cmplx(-0.1028193852454817599302497455967E1_dp,-0.361376517118561592479460832998E0_dp,kind=dp)
            exact_roots(21) = cmplx(-0.1028193852454817599302497455967E1_dp,0.361376517118561592479460832998E0_dp,kind=dp)
            exact_roots(22) = cmplx(-0.623532485956252757990016587001E0_dp,-0.681064414225239608090835812687E0_dp,kind=dp)
            exact_roots(23) = cmplx(-0.623532485956252757990016587001E0_dp,0.681064414225239608090835812687E0_dp,kind=dp)
            exact_roots(24) = cmplx(-0.622436295041293587960163506947E0_dp,-0.424878436475629184311574438805E0_dp,kind=dp)
            exact_roots(25) = cmplx(-0.622436295041293587960163506947E0_dp,0.424878436475629184311574438805E0_dp,kind=dp)
            exact_roots(26) = cmplx(-0.530827804859942728921477297120E0_dp,-0.668288725559205771444092465565E0_dp,kind=dp)
            exact_roots(27) = cmplx(-0.530827804859942728921477297120E0_dp,0.668288725559205771444092465565E0_dp,kind=dp)
            exact_roots(28) = cmplx(-0.272102461488938894219383324518E0_dp,-0.842364690294128145503155708243E0_dp,kind=dp)
            exact_roots(29) = cmplx(-0.272102461488938894219383324518E0_dp,0.842364690294128145503155708243E0_dp,kind=dp)
            exact_roots(30) = cmplx(-0.224915951286740054685326255204E0_dp,-0.1116260157454991835001268254245E1_dp,kind=dp)
            exact_roots(31) = cmplx(-0.224915951286740054685326255204E0_dp,0.1116260157454991835001268254245E1_dp,kind=dp)
            exact_roots(32) = cmplx(-0.207283835455666412824133850187E0_dp,-0.1117480772494962911373775673122E1_dp,kind=dp)
            exact_roots(33) = cmplx(-0.207283835455666412824133850187E0_dp,0.1117480772494962911373775673122E1_dp,kind=dp)
            exact_roots(34) = cmplx(-0.174578221135716969451566432662E0_dp,-0.1071427671454031189229646310220E1_dp,kind=dp)
            exact_roots(35) = cmplx(-0.174578221135716969451566432662E0_dp,0.1071427671454031189229646310220E1_dp,kind=dp)
            exact_roots(36) = cmplx(-0.157516053475965356164335109645E0_dp,-0.1109006514113607177971751986155E1_dp,kind=dp)
            exact_roots(37) = cmplx(-0.157516053475965356164335109645E0_dp,0.1109006514113607177971751986155E1_dp,kind=dp)
            exact_roots(38) = cmplx(-0.127499973546363000199539565346E0_dp,-0.987460909489456792207407680793E0_dp,kind=dp)
            exact_roots(39) = cmplx(-0.127499973546363000199539565346E0_dp,0.987460909489456792207407680793E0_dp,kind=dp)
            exact_roots(40) = cmplx(-0.14233481920354066767761845314E-1_dp,-0.1032914775213644109395013402655E1_dp,kind=dp)
            exact_roots(41) = cmplx(-0.14233481920354066767761845314E-1_dp,0.1032914775213644109395013402655E1_dp,kind=dp)
            exact_roots(42) = cmplx(-0.6983568496261391817966491075E-2_dp,-0.1003603862288289548530704966951E1_dp,kind=dp)
            exact_roots(43) = cmplx(-0.6983568496261391817966491075E-2_dp,0.1003603862288289548530704966951E1_dp,kind=dp)
            exact_roots(44) = cmplx(0.14895466603687646529815779209E-1_dp,-0.848148761908416527719331111783E0_dp,kind=dp)
            exact_roots(45) = cmplx(0.14895466603687646529815779209E-1_dp,0.848148761908416527719331111783E0_dp,kind=dp)
            exact_roots(46) = cmplx(0.121192786105906486314704443411E0_dp,-0.610611692210754211675387244150E0_dp,kind=dp)
            exact_roots(47) = cmplx(0.121192786105906486314704443411E0_dp,0.610611692210754211675387244150E0_dp,kind=dp)
            exact_roots(48) = cmplx(0.352482539722363278193253964052E0_dp,-0.698337239583330331258141954760E0_dp,kind=dp)
            exact_roots(49) = cmplx(0.352482539722363278193253964052E0_dp,0.698337239583330331258141954760E0_dp,kind=dp)
            exact_roots(50) = cmplx(0.376008681846767559704804317729E0_dp,-0.144749371321632864747110182013E0_dp,kind=dp)
            exact_roots(51) = cmplx(0.376008681846767559704804317729E0_dp,0.144749371321632864747110182013E0_dp,kind=dp)
            exact_roots(52) = cmplx(0.376893240379311323690004017969E0_dp,-0.678568693190448141957540792997E0_dp,kind=dp)
            exact_roots(53) = cmplx(0.376893240379311323690004017969E0_dp,0.678568693190448141957540792997E0_dp,kind=dp)
            exact_roots(54) = cmplx(0.386539176596158026508293086904E0_dp,-0.569324711303102903213792357135E0_dp,kind=dp)
            exact_roots(55) = cmplx(0.386539176596158026508293086904E0_dp,0.569324711303102903213792357135E0_dp,kind=dp)
            exact_roots(56) = cmplx(0.412916024722700479197334566383E0_dp,-0.614806760143385694954549720401E0_dp,kind=dp)
            exact_roots(57) = cmplx(0.412916024722700479197334566383E0_dp,0.614806760143385694954549720401E0_dp,kind=dp)
            exact_roots(58) = cmplx(0.432376192641994507824669648087E0_dp,-0.226759904435348618697876559972E0_dp,kind=dp)
            exact_roots(59) = cmplx(0.432376192641994507824669648087E0_dp,0.226759904435348618697876559972E0_dp,kind=dp)
            exact_roots(60) = cmplx(0.452774498724915493508803077733E0_dp,-0.396170128033165002412596877271E0_dp,kind=dp)
            exact_roots(61) = cmplx(0.452774498724915493508803077733E0_dp,0.396170128033165002412596877271E0_dp,kind=dp)
            exact_roots(62) = cmplx(0.456823285823316651283953236253E0_dp,-0.347758700883481983632188723200E0_dp,kind=dp)
            exact_roots(63) = cmplx(0.456823285823316651283953236253E0_dp,0.347758700883481983632188723200E0_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(1E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(1E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(2E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(5E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(14E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(42E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(132E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(365E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(950E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(2398E0_dp,0E0_dp,kind=dp)
            p(11) = cmplx(5916E0_dp,0E0_dp,kind=dp)
            p(12) = cmplx(14290E0_dp,0E0_dp,kind=dp)
            p(13) = cmplx(33708E0_dp,0E0_dp,kind=dp)
            p(14) = cmplx(77684E0_dp,0E0_dp,kind=dp)
            p(15) = cmplx(175048E0_dp,0E0_dp,kind=dp)
            p(16) = cmplx(385741E0_dp,0E0_dp,kind=dp)
            p(17) = cmplx(831014E0_dp,0E0_dp,kind=dp)
            p(18) = cmplx(1749654E0_dp,0E0_dp,kind=dp)
            p(19) = cmplx(3598964E0_dp,0E0_dp,kind=dp)
            p(20) = cmplx(7228014E0_dp,0E0_dp,kind=dp)
            p(21) = cmplx(14162220E0_dp,0E0_dp,kind=dp)
            p(22) = cmplx(27049196E0_dp,0E0_dp,kind=dp)
            p(23) = cmplx(50323496E0_dp,0E0_dp,kind=dp)
            p(24) = cmplx(91143114E0_dp,0E0_dp,kind=dp)
            p(25) = cmplx(160617860E0_dp,0E0_dp,kind=dp)
            p(26) = cmplx(275276716E0_dp,0E0_dp,kind=dp)
            p(27) = cmplx(458591432E0_dp,0E0_dp,kind=dp)
            p(28) = cmplx(742179284E0_dp,0E0_dp,kind=dp)
            p(29) = cmplx(1166067016E0_dp,0E0_dp,kind=dp)
            p(30) = cmplx(1777171560E0_dp,0E0_dp,kind=dp)
            p(31) = cmplx(2625062128E0_dp,0E0_dp,kind=dp)
            p(32) = cmplx(3754272037E0_dp,0E0_dp,kind=dp)
            p(33) = cmplx(5193067630E0_dp,0E0_dp,kind=dp)
            p(34) = cmplx(6939692682E0_dp,0E0_dp,kind=dp)
            p(35) = cmplx(8948546308E0_dp,0E0_dp,kind=dp)
            p(36) = cmplx(11120136162E0_dp,0E0_dp,kind=dp)
            p(37) = cmplx(13299362332E0_dp,0E0_dp,kind=dp)
            p(38) = cmplx(15286065700E0_dp,0E0_dp,kind=dp)
            p(39) = cmplx(16859410792E0_dp,0E0_dp,kind=dp)
            p(40) = cmplx(17813777994E0_dp,0E0_dp,kind=dp)
            p(41) = cmplx(17999433372E0_dp,0E0_dp,kind=dp)
            p(42) = cmplx(17357937708E0_dp,0E0_dp,kind=dp)
            p(43) = cmplx(15941684776E0_dp,0E0_dp,kind=dp)
            p(44) = cmplx(13910043524E0_dp,0E0_dp,kind=dp)
            p(45) = cmplx(11500901864E0_dp,0E0_dp,kind=dp)
            p(46) = cmplx(8984070856E0_dp,0E0_dp,kind=dp)
            p(47) = cmplx(6609143792E0_dp,0E0_dp,kind=dp)
            p(48) = cmplx(4562339774E0_dp,0E0_dp,kind=dp)
            p(49) = cmplx(2943492972E0_dp,0E0_dp,kind=dp)
            p(50) = cmplx(1766948340E0_dp,0E0_dp,kind=dp)
            p(51) = cmplx(981900168E0_dp,0E0_dp,kind=dp)
            p(52) = cmplx(502196500E0_dp,0E0_dp,kind=dp)
            p(53) = cmplx(234813592E0_dp,0E0_dp,kind=dp)
            p(54) = cmplx(99582920E0_dp,0E0_dp,kind=dp)
            p(55) = cmplx(37945904E0_dp,0E0_dp,kind=dp)
            p(56) = cmplx(12843980E0_dp,0E0_dp,kind=dp)
            p(57) = cmplx(3807704E0_dp,0E0_dp,kind=dp)
            p(58) = cmplx(971272E0_dp,0E0_dp,kind=dp)
            p(59) = cmplx(208336E0_dp,0E0_dp,kind=dp)
            p(60) = cmplx(36440E0_dp,0E0_dp,kind=dp)
            p(61) = cmplx(4976E0_dp,0E0_dp,kind=dp)
            p(62) = cmplx(496E0_dp,0E0_dp,kind=dp)
            p(63) = cmplx(32E0_dp,0E0_dp,kind=dp)
            p(64) = cmplx(1E0_dp,0E0_dp,kind=dp)
            ! case 2
            case(2)
            ! Kameny Polynomial c=10, theoretical limiting accuracy 1.11E-16
            deg = 9
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-0.250708394126194522431084083966E1_dp,0E0_dp,kind=dp)
            exact_roots(2) = cmplx(-0.77475933255040504312683246509E0_dp,-0.239347474286129783412456829950E1_dp,kind=dp)
            exact_roots(3) = cmplx(-0.77475933255040504312683246509E0_dp,0.239347474286129783412456829950E1_dp,kind=dp)
            exact_roots(4) = cmplx(-0.173313478195944654943157024747E0_dp,0E0_dp,kind=dp)
            exact_roots(5) = cmplx(-0.173097223325097844852013748504E0_dp,0E0_dp,kind=dp)
            exact_roots(6) = cmplx(0.173204810760521109164135307321E0_dp,-0.108125591879847889809932772E-3_dp,kind=dp)
            exact_roots(7) = cmplx(0.173204810760521109164135307321E0_dp,0.108125591879847889809932772E-3_dp,kind=dp)
            exact_roots(8) = cmplx(0.202830184318137779601570296422E1_dp,-0.147928157900441456648083759833E1_dp,kind=dp)
            exact_roots(9) = cmplx(0.202830184318137779601570296422E1_dp,0.147928157900441456648083759833E1_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(9E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(-600E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(10000E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(100E0_dp,0E0_dp,kind=dp)
            ! case 3
            case(3)
            ! Kameny Polynomial c=10^3, theoretical limiting accuracy 1.11E-16
            deg = 9
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-0.158489318488962525196907753240E2_dp,0E0_dp,kind=dp)
            exact_roots(2) = cmplx(-0.48975892839992994551324676957E1_dp,-0.150732300552288401309352587876E2_dp,kind=dp)
            exact_roots(3) = cmplx(-0.48975892839992994551324676957E1_dp,0.150732300552288401309352587876E2_dp,kind=dp)
            exact_roots(4) = cmplx(-0.173205080767700380719051028416E-2_dp,0E0_dp,kind=dp)
            exact_roots(5) = cmplx(-0.173205080746075077991838239885E-2_dp,0E0_dp,kind=dp)
            exact_roots(6) = cmplx(0.173205080756887729350044634151E-2_dp,-0.10812651363606394264E-12_dp,kind=dp)
            exact_roots(7) = cmplx(0.173205080756887729350044634151E-2_dp,0.10812651363606394264E-12_dp,kind=dp)
            exact_roots(8) = cmplx(0.128220552084474257149779093577E2_dp,-0.93157684943778791903260859371E1_dp,kind=dp)
            exact_roots(9) = cmplx(0.128220552084474257149779093577E2_dp,0.93157684943778791903260859371E1_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(9E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(-6000000E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(1000000000000E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(1000000E0_dp,0E0_dp,kind=dp)
            ! case 4
            case(4)
            ! Kameny Polynomial c=10^6, theoretical limiting accuracy 8.09E-9
            deg = 9
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-0.251188643150958006331217160138E3_dp,0E0_dp,kind=dp)
            exact_roots(2) = cmplx(-0.77621559527630266049907587459E2_dp,-0.238894595888056617951698213041E3_dp,kind=dp)
            exact_roots(3) = cmplx(-0.77621559527630266049907587459E2_dp,0.238894595888056617951698213041E3_dp,kind=dp)
            exact_roots(4) = cmplx(-0.173205080756887729353086560209E-5_dp,0E0_dp,kind=dp)
            exact_roots(5) = cmplx(-0.173205080756887729352402708092E-5_dp,0E0_dp,kind=dp)
            exact_roots(6) = cmplx(0.173205080756887729352744634151E-5_dp,-0.341926059E-26_dp,kind=dp)
            exact_roots(7) = cmplx(0.173205080756887729352744634151E-5_dp,0.341926059E-26_dp,kind=dp)
            exact_roots(8) = cmplx(0.203215881103109269215516167528E3_dp,-0.147644979987489859882850769633E3_dp,kind=dp)
            exact_roots(9) = cmplx(0.203215881103109269215516167528E3_dp,0.147644979987489859882850769633E3_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(9E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(-6000000000000E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(1000000000000000000000000E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(1000000000000E0_dp,0E0_dp,kind=dp)
            ! case 5
            case(5)
            ! Multiplicity I
            deg = 55
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(2) = cmplx(-0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(3) = cmplx(-0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(4) = cmplx(-0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(5) = cmplx(-0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(6) = cmplx(-0.951677511788434760780566076462E0_dp,-0.156577036126765257092467586019E0_dp,kind=dp)
            exact_roots(7) = cmplx(-0.951677511788434760780566076462E0_dp,0.156577036126765257092467586019E0_dp,kind=dp)
            exact_roots(8) = cmplx(-0.947111794804360643443628545664E0_dp,-0.45987543603302557306898551016E-1_dp,kind=dp)
            exact_roots(9) = cmplx(-0.947111794804360643443628545664E0_dp,0.45987543603302557306898551016E-1_dp,kind=dp)
            exact_roots(10) = cmplx(-0.935147084823051581325738268702E0_dp,-0.276366395363089095701385358741E0_dp,kind=dp)
            exact_roots(11) = cmplx(-0.935147084823051581325738268702E0_dp,0.276366395363089095701385358741E0_dp,kind=dp)
            exact_roots(12) = cmplx(-0.899451950794801943937381293773E0_dp,-0.394568367012592113985865567971E0_dp,kind=dp)
            exact_roots(13) = cmplx(-0.899451950794801943937381293773E0_dp,0.394568367012592113985865567971E0_dp,kind=dp)
            exact_roots(14) = cmplx(-0.846928410223891427281972003960E0_dp,-0.507591758019284469647521636947E0_dp,kind=dp)
            exact_roots(15) = cmplx(-0.846928410223891427281972003960E0_dp,0.507591758019284469647521636947E0_dp,kind=dp)
            exact_roots(16) = cmplx(-0.779195724062158658036386916327E0_dp,-0.613081457400403535591975580892E0_dp,kind=dp)
            exact_roots(17) = cmplx(-0.779195724062158658036386916327E0_dp,0.613081457400403535591975580892E0_dp,kind=dp)
            exact_roots(18) = cmplx(-0.697743392981121831801589116937E0_dp,-0.709081007070975090262670227767E0_dp,kind=dp)
            exact_roots(19) = cmplx(-0.697743392981121831801589116937E0_dp,0.709081007070975090262670227767E0_dp,kind=dp)
            exact_roots(20) = cmplx(-0.604121886825094827201080735988E0_dp,-0.793887314714735149770017023538E0_dp,kind=dp)
            exact_roots(21) = cmplx(-0.604121886825094827201080735988E0_dp,0.793887314714735149770017023538E0_dp,kind=dp)
            exact_roots(22) = cmplx(-0.500000000000000000000000000000E0_dp,-0.866025403784438646763723170753E0_dp,kind=dp)
            exact_roots(23) = cmplx(-0.500000000000000000000000000000E0_dp,0.866025403784438646763723170753E0_dp,kind=dp)
            exact_roots(24) = cmplx(-0.387174132037879582862879778205E0_dp,-0.924251993602053902039550620488E0_dp,kind=dp)
            exact_roots(25) = cmplx(-0.387174132037879582862879778205E0_dp,0.924251993602053902039550620488E0_dp,kind=dp)
            exact_roots(26) = cmplx(-0.267556959428167474972897159502E0_dp,-0.967566373494130224905327435791E0_dp,kind=dp)
            exact_roots(27) = cmplx(-0.267556959428167474972897159502E0_dp,0.967566373494130224905327435791E0_dp,kind=dp)
            exact_roots(28) = cmplx(-0.143155839532593341810009310219E0_dp,-0.995221832724034632237407792227E0_dp,kind=dp)
            exact_roots(29) = cmplx(-0.143155839532593341810009310219E0_dp,0.995221832724034632237407792227E0_dp,kind=dp)
            exact_roots(30) = cmplx(-0.16045517947534302681825863292E-1_dp,-0.1006735041626386010701028433184E1_dp,kind=dp)
            exact_roots(31) = cmplx(-0.16045517947534302681825863292E-1_dp,0.1006735041626386010701028433184E1_dp,kind=dp)
            exact_roots(32) = cmplx(0.111662481142841735790134563805E0_dp,-0.1001892206008032435299692737060E1_dp,kind=dp)
            exact_roots(33) = cmplx(0.111662481142841735790134563805E0_dp,0.1001892206008032435299692737060E1_dp,kind=dp)
            exact_roots(34) = cmplx(0.237852321567861477202899799722E0_dp,-0.980751410679651646035316169672E0_dp,kind=dp)
            exact_roots(35) = cmplx(0.237852321567861477202899799722E0_dp,0.980751410679651646035316169672E0_dp,kind=dp)
            exact_roots(36) = cmplx(0.360437139274099643543556447604E0_dp,-0.943640865401314173256925688139E0_dp,kind=dp)
            exact_roots(37) = cmplx(0.360437139274099643543556447604E0_dp,0.943640865401314173256925688139E0_dp,kind=dp)
            exact_roots(38) = cmplx(0.477392314844102847110474634171E0_dp,-0.891152945345826129924640822366E0_dp,kind=dp)
            exact_roots(39) = cmplx(0.477392314844102847110474634171E0_dp,0.891152945345826129924640822366E0_dp,kind=dp)
            exact_roots(40) = cmplx(0.586788023132224235580014756671E0_dp,-0.824134044230047196220784957068E0_dp,kind=dp)
            exact_roots(41) = cmplx(0.586788023132224235580014756671E0_dp,0.824134044230047196220784957068E0_dp,kind=dp)
            exact_roots(42) = cmplx(0.686820425487486753027867016412E0_dp,-0.743670357958749198637378689329E0_dp,kind=dp)
            exact_roots(43) = cmplx(0.686820425487486753027867016412E0_dp,0.743670357958749198637378689329E0_dp,kind=dp)
            exact_roots(44) = cmplx(0.775840929682842739762490529077E0_dp,-0.651069801113127816898897251309E0_dp,kind=dp)
            exact_roots(45) = cmplx(0.775840929682842739762490529077E0_dp,0.651069801113127816898897251309E0_dp,kind=dp)
            exact_roots(46) = cmplx(0.852382996720766159190958138062E0_dp,-0.547840332639938504859225568105E0_dp,kind=dp)
            exact_roots(47) = cmplx(0.852382996720766159190958138062E0_dp,0.547840332639938504859225568105E0_dp,kind=dp)
            exact_roots(48) = cmplx(0.915186027756796108167454007547E0_dp,-0.435665032553729539118325217757E0_dp,kind=dp)
            exact_roots(49) = cmplx(0.915186027756796108167454007547E0_dp,0.435665032553729539118325217757E0_dp,kind=dp)
            exact_roots(50) = cmplx(0.963215921119370289466626861950E0_dp,-0.316374328936911136376493374596E0_dp,kind=dp)
            exact_roots(51) = cmplx(0.963215921119370289466626861950E0_dp,0.316374328936911136376493374596E0_dp,kind=dp)
            exact_roots(52) = cmplx(0.995681949864471509895340077408E0_dp,-0.191915823987060561907739223463E0_dp,kind=dp)
            exact_roots(53) = cmplx(0.995681949864471509895340077408E0_dp,0.191915823987060561907739223463E0_dp,kind=dp)
            exact_roots(54) = cmplx(0.1012049674656226877398138236602E1_dp,-0.64322209034046437066803422584E-1_dp,kind=dp)
            exact_roots(55) = cmplx(0.1012049674656226877398138236602E1_dp,0.64322209034046437066803422584E-1_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(1E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(6E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(15E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(20E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(15E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(6E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(1E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(11) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(12) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(13) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(14) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(15) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(16) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(17) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(18) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(19) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(20) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(21) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(22) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(23) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(24) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(25) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(26) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(27) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(28) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(29) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(30) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(31) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(32) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(33) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(34) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(35) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(36) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(37) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(38) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(39) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(40) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(41) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(42) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(43) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(44) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(45) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(46) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(47) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(48) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(49) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(50) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(51) = cmplx(1E0_dp,0E0_dp,kind=dp)
            p(52) = cmplx(5E0_dp,0E0_dp,kind=dp)
            p(53) = cmplx(10E0_dp,0E0_dp,kind=dp)
            p(54) = cmplx(10E0_dp,0E0_dp,kind=dp)
            p(55) = cmplx(5E0_dp,0E0_dp,kind=dp)
            p(56) = cmplx(1E0_dp,0E0_dp,kind=dp)
            ! case 6
            case(6)
            ! Multiplicity II
            deg = 62
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(-0.998026728428271561952336806863E0_dp,-0.62790519529313376076178224566E-1_dp,kind=dp)
            exact_roots(2) = cmplx(-0.998026728428271561952336806863E0_dp,0.62790519529313376076178224566E-1_dp,kind=dp)
            exact_roots(3) = cmplx(-0.982287250728688681085641742865E0_dp,-0.187381314585724630542550734447E0_dp,kind=dp)
            exact_roots(4) = cmplx(-0.982287250728688681085641742865E0_dp,0.187381314585724630542550734447E0_dp,kind=dp)
            exact_roots(5) = cmplx(-0.951056516295153572116439333379E0_dp,-0.309016994374947424102293417183E0_dp,kind=dp)
            exact_roots(6) = cmplx(-0.951056516295153572116439333379E0_dp,0.309016994374947424102293417183E0_dp,kind=dp)
            exact_roots(7) = cmplx(-0.904827052466019527713668647933E0_dp,-0.425779291565072648862502445744E0_dp,kind=dp)
            exact_roots(8) = cmplx(-0.904827052466019527713668647933E0_dp,0.425779291565072648862502445744E0_dp,kind=dp)
            exact_roots(9) = cmplx(-0.844327925502015078548558063967E0_dp,-0.535826794978996618271308767868E0_dp,kind=dp)
            exact_roots(10) = cmplx(-0.844327925502015078548558063967E0_dp,0.535826794978996618271308767868E0_dp,kind=dp)
            exact_roots(11) = cmplx(-0.770513242775789230803009636396E0_dp,-0.637423989748689710176712811676E0_dp,kind=dp)
            exact_roots(12) = cmplx(-0.770513242775789230803009636396E0_dp,0.637423989748689710176712811676E0_dp,kind=dp)
            exact_roots(13) = cmplx(-0.684547105928688673732283357621E0_dp,-0.728968627421411523146730319055E0_dp,kind=dp)
            exact_roots(14) = cmplx(-0.684547105928688673732283357621E0_dp,0.728968627421411523146730319055E0_dp,kind=dp)
            exact_roots(15) = cmplx(-0.587785252292473129168705954639E0_dp,-0.809016994374947424102293417183E0_dp,kind=dp)
            exact_roots(16) = cmplx(-0.587785252292473129168705954639E0_dp,0.809016994374947424102293417183E0_dp,kind=dp)
            exact_roots(17) = cmplx(-0.50000000000000000000000000000E0_dp,-0.217944947177033677611849099193E1_dp,kind=dp)
            exact_roots(18) = cmplx(-0.50000000000000000000000000000E0_dp,-0.217944947177033677611849099193E1_dp,kind=dp)
            exact_roots(19) = cmplx(-0.50000000000000000000000000000E0_dp,-0.217944947177033677611849099193E1_dp,kind=dp)
            exact_roots(20) = cmplx(-0.50000000000000000000000000000E0_dp,0.217944947177033677611849099193E1_dp,kind=dp)
            exact_roots(21) = cmplx(-0.50000000000000000000000000000E0_dp,0.217944947177033677611849099193E1_dp,kind=dp)
            exact_roots(22) = cmplx(-0.50000000000000000000000000000E0_dp,0.217944947177033677611849099193E1_dp,kind=dp)
            exact_roots(23) = cmplx(-0.481753674101715274987191502872E0_dp,-0.876306680043863587308115903922E0_dp,kind=dp)
            exact_roots(24) = cmplx(-0.481753674101715274987191502872E0_dp,0.876306680043863587308115903922E0_dp,kind=dp)
            exact_roots(25) = cmplx(-0.368124552684677959156947147493E0_dp,-0.929776485888251403660942556222E0_dp,kind=dp)
            exact_roots(26) = cmplx(-0.368124552684677959156947147493E0_dp,0.929776485888251403660942556222E0_dp,kind=dp)
            exact_roots(27) = cmplx(-0.248689887164854788242283746006E0_dp,-0.968583161128631119490168375465E0_dp,kind=dp)
            exact_roots(28) = cmplx(-0.248689887164854788242283746006E0_dp,0.968583161128631119490168375465E0_dp,kind=dp)
            exact_roots(29) = cmplx(-0.125333233564304245373118759817E0_dp,-0.992114701314477831049793042786E0_dp,kind=dp)
            exact_roots(30) = cmplx(-0.125333233564304245373118759817E0_dp,0.992114701314477831049793042786E0_dp,kind=dp)
            exact_roots(31) = cmplx(0E0_dp,-0.100000000000000000000000000000E1_dp,kind=dp)
            exact_roots(32) = cmplx(0E0_dp,0.100000000000000000000000000000E1_dp,kind=dp)
            exact_roots(33) = cmplx(0.125333233564304245373118759817E0_dp,-0.992114701314477831049793042786E0_dp,kind=dp)
            exact_roots(34) = cmplx(0.125333233564304245373118759817E0_dp,0.992114701314477831049793042786E0_dp,kind=dp)
            exact_roots(35) = cmplx(0.248689887164854788242283746006E0_dp,-0.968583161128631119490168375465E0_dp,kind=dp)
            exact_roots(36) = cmplx(0.248689887164854788242283746006E0_dp,0.968583161128631119490168375465E0_dp,kind=dp)
            exact_roots(37) = cmplx(0.333333333333333333333333333333E0_dp,0E0_dp,kind=dp)
            exact_roots(38) = cmplx(0.333333333333333333333333333333E0_dp,0E0_dp,kind=dp)
            exact_roots(39) = cmplx(0.368124552684677959156947147493E0_dp,-0.929776485888251403660942556222E0_dp,kind=dp)
            exact_roots(40) = cmplx(0.368124552684677959156947147493E0_dp,0.929776485888251403660942556222E0_dp,kind=dp)
            exact_roots(41) = cmplx(0.481753674101715274987191502872E0_dp,-0.876306680043863587308115903922E0_dp,kind=dp)
            exact_roots(42) = cmplx(0.481753674101715274987191502872E0_dp,0.876306680043863587308115903922E0_dp,kind=dp)
            exact_roots(43) = cmplx(0.587785252292473129168705954639E0_dp,-0.809016994374947424102293417183E0_dp,kind=dp)
            exact_roots(44) = cmplx(0.587785252292473129168705954639E0_dp,0.809016994374947424102293417183E0_dp,kind=dp)
            exact_roots(45) = cmplx(0.684547105928688673732283357621E0_dp,-0.728968627421411523146730319055E0_dp,kind=dp)
            exact_roots(46) = cmplx(0.684547105928688673732283357621E0_dp,0.728968627421411523146730319055E0_dp,kind=dp)
            exact_roots(47) = cmplx(0.770513242775789230803009636396E0_dp,-0.637423989748689710176712811676E0_dp,kind=dp)
            exact_roots(48) = cmplx(0.770513242775789230803009636396E0_dp,0.637423989748689710176712811676E0_dp,kind=dp)
            exact_roots(49) = cmplx(0.844327925502015078548558063967E0_dp,-0.535826794978996618271308767868E0_dp,kind=dp)
            exact_roots(50) = cmplx(0.844327925502015078548558063967E0_dp,0.535826794978996618271308767868E0_dp,kind=dp)
            exact_roots(51) = cmplx(0.904827052466019527713668647933E0_dp,-0.425779291565072648862502445744E0_dp,kind=dp)
            exact_roots(52) = cmplx(0.904827052466019527713668647933E0_dp,0.425779291565072648862502445744E0_dp,kind=dp)
            exact_roots(53) = cmplx(0.951056516295153572116439333379E0_dp,-0.309016994374947424102293417183E0_dp,kind=dp)
            exact_roots(54) = cmplx(0.951056516295153572116439333379E0_dp,0.309016994374947424102293417183E0_dp,kind=dp)
            exact_roots(55) = cmplx(0.982287250728688681085641742865E0_dp,-0.187381314585724630542550734447E0_dp,kind=dp)
            exact_roots(56) = cmplx(0.982287250728688681085641742865E0_dp,0.187381314585724630542550734447E0_dp,kind=dp)
            exact_roots(57) = cmplx(0.998026728428271561952336806863E0_dp,-0.62790519529313376076178224566E-1_dp,kind=dp)
            exact_roots(58) = cmplx(0.998026728428271561952336806863E0_dp,0.62790519529313376076178224566E-1_dp,kind=dp)
            exact_roots(59) = cmplx(0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(60) = cmplx(0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(61) = cmplx(0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(62) = cmplx(0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(125E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(-1175E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(4215E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(-7444E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(7393E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(-5133E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(3402E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(-1917E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(741E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(-316E0_dp,0E0_dp,kind=dp)
            p(11) = cmplx(115E0_dp,0E0_dp,kind=dp)
            p(12) = cmplx(-15E0_dp,0E0_dp,kind=dp)
            p(13) = cmplx(9E0_dp,0E0_dp,kind=dp)
            p(14) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(15) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(16) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(17) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(18) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(19) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(20) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(21) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(22) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(23) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(24) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(25) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(26) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(27) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(28) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(29) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(30) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(31) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(32) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(33) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(34) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(35) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(36) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(37) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(38) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(39) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(40) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(41) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(42) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(43) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(44) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(45) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(46) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(47) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(48) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(49) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(50) = cmplx(0E0_dp,0E0_dp,kind=dp)
            p(51) = cmplx(125E0_dp,0E0_dp,kind=dp)
            p(52) = cmplx(-1175E0_dp,0E0_dp,kind=dp)
            p(53) = cmplx(4215E0_dp,0E0_dp,kind=dp)
            p(54) = cmplx(-7444E0_dp,0E0_dp,kind=dp)
            p(55) = cmplx(7393E0_dp,0E0_dp,kind=dp)
            p(56) = cmplx(-5133E0_dp,0E0_dp,kind=dp)
            p(57) = cmplx(3402E0_dp,0E0_dp,kind=dp)
            p(58) = cmplx(-1917E0_dp,0E0_dp,kind=dp)
            p(59) = cmplx(741E0_dp,0E0_dp,kind=dp)
            p(60) = cmplx(-316E0_dp,0E0_dp,kind=dp)
            p(61) = cmplx(115E0_dp,0E0_dp,kind=dp)
            p(62) = cmplx(-15E0_dp,0E0_dp,kind=dp)
            p(63) = cmplx(9E0_dp,0E0_dp,kind=dp)
            ! case 7
            case(7)
            ! Multiplicity III
            deg = 17
            allocate(exact_roots(deg), p(deg+1))
            ! roots
            exact_roots(1) = cmplx(0.100000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(2) = cmplx(0.200000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(3) = cmplx(0.300000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(4) = cmplx(0.400000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(5) = cmplx(0.500000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(6) = cmplx(0.600000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(7) = cmplx(0.700000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(8) = cmplx(0.800000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(9) = cmplx(0.900000000000000000000000000000E1_dp,0E0_dp,kind=dp)
            exact_roots(10) = cmplx(0.100000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(11) = cmplx(0.110000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(12) = cmplx(0.120000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(13) = cmplx(0.130000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(14) = cmplx(0.140000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(15) = cmplx(0.150000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(16) = cmplx(0.150000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            exact_roots(17) = cmplx(0.150000000000000000000000000000E2_dp,0E0_dp,kind=dp)
            ! polynomial
            p(1) = cmplx(-294226732800000E0_dp,0E0_dp,kind=dp)
            p(2) = cmplx(1015541906400000E0_dp,0E0_dp,kind=dp)
            p(3) = cmplx(-1518791527728000E0_dp,0E0_dp,kind=dp)
            p(4) = cmplx(1327137724803600E0_dp,0E0_dp,kind=dp)
            p(5) = cmplx(-766908691489440E0_dp,0E0_dp,kind=dp)
            p(6) = cmplx(313437620164824E0_dp,0E0_dp,kind=dp)
            p(7) = cmplx(-94377698961000E0_dp,0E0_dp,kind=dp)
            p(8) = cmplx(21485772576905E0_dp,0E0_dp,kind=dp)
            p(9) = cmplx(-3758453397270E0_dp,0E0_dp,kind=dp)
            p(10) = cmplx(509681511053E0_dp,0E0_dp,kind=dp)
            p(11) = cmplx(-53726158200E0_dp,0E0_dp,kind=dp)
            p(12) = cmplx(4387265090E0_dp,0E0_dp,kind=dp)
            p(13) = cmplx(-274687140E0_dp,0E0_dp,kind=dp)
            p(14) = cmplx(12932122E0_dp,0E0_dp,kind=dp)
            p(15) = cmplx(-442800E0_dp,0E0_dp,kind=dp)
            p(16) = cmplx(10405E0_dp,0E0_dp,kind=dp)
            p(17) = cmplx(-150E0_dp,0E0_dp,kind=dp)
            p(18) = cmplx(1E0_dp,0E0_dp,kind=dp)
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
end program forwerr_multi