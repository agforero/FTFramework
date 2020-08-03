DOUBLE PRECISION FUNCTION hh(npt, func, ier)

!  Half-Hermite integration, i.e. integration from 0 to infinity of

!                 f(x).exp(-x**2)

!  Arguments:-
!  NPT    = INPUT, no.of points at which the function, f(x), is
!           to be evaluated.   NPT must be between 2 and 20 incl.
!  FUNC   = INPUT, name of the user's function to supply values of f(x).
!           The function may have any name, but only the one argument, x.
!           FUNC - or whatever name you use, must be declared double
!           precision and external in the calling program.
!  IER    = OUTPUT, error indicator
!         = 0 if no error detected
!         = 1 if an illegal value was supplied for NPT.

!  N.B.  The function f(x) must have continuous derivatives of all orders
!  throughout the range of integration.   The abscissa and weights are such
!  that the approximate integral returned is accurate to about 14 decimal
!  digits if n points are used and f(x) is a polynomial of degree (2n-1) or
!  less.   The abscissae and weights supplied give slightly more accurate
!  results than those given in the two references below, which were used as
!  starting values.

!  References:-
!  Steen, M.M., Byrne, G.D. & Gelbard, E.M. (1969).   Gaussian quadrature
!    for the integrals ... .   Mathematics of Computation, v.23, 661-671.
!  Kahaner, D., Tietjen, G. & Beckman, R. (1982). Gaussian-quadrature
!    formulas for ... .   J. Statist. Comput. Simul., v.15, 155-160.

!  Programmer:   Alan Miller
!                CSIRO Division of Mathematics & Statistics
!                Private Bag 10
!                Clayton, Victoria 3169
!  e-mail: amiller @ bigpond.net.au

!  `Translated' from Fortran 77 code - 24 May 1996
!  Latest revision - 21 January 2000

!***********************************************************************

IMPLICIT NONE
INTEGER, INTENT(IN)  :: npt
INTEGER, INTENT(OUT) :: ier
DOUBLE PRECISION     :: func
EXTERNAL func

! Local variables
DOUBLE PRECISION :: xpt(209), wt(209)
INTEGER          :: i, ipos

DATA xpt(1:119)/ &
  0.300193931060839D+00, 0.125242104533372D+01, &                                                !2
  0.190554149798192D+00, 0.848251867544577D+00, 0.179977657841573D+01, &                         !3
  0.133776446996068D+00, 0.624324690187190D+00, 0.134253782564499D+01, 0.226266447701036D+01, &  !4
  0.100242151968216D+00, 0.482813966046201D+00, 0.106094982152572D+01, 0.177972941852026D+01, &  !5
  0.266976035608766D+01, &
  0.786006594130979D-01, 0.386739410270631D+00, 0.866429471682044D+00, 0.146569804966352D+01, &  !6
  0.217270779693900D+01, 0.303682016932287D+01, &
  0.637164846067008D-01, 0.318192018888619D+00, 0.724198989258373D+00, 0.123803559921509D+01, &  !7
  0.183852822027095D+01, 0.253148815132768D+01, 0.337345643012458D+01, &
  0.529786439318511D-01, 0.267398372167765D+00, 0.616302884182400D+00, 0.106424631211622D+01, &  !8
  0.158885586227006D+01, 0.218392115309586D+01, 0.286313388370807D+01, 0.368600716272440D+01, &
  0.449390308011905D-01, 0.228605305560523D+00, 0.532195844331623D+00, 0.927280745338049D+00, &  !9
  0.139292385519585D+01, 0.191884309919739D+01, 0.250624783400570D+01, 0.317269213348120D+01, &
  0.397889886978974D+01, &
  0.387385243256994D-01, 0.198233304012949D+00, 0.465201111814507D+00, 0.816861885591907D+00, &  !10
  0.123454132402774D+01, 0.170679814968865D+01, 0.222994008892444D+01, 0.280910374689825D+01, &
  0.346387241949537D+01, 0.425536180636561D+01, &
  0.338393212317745D-01, 0.173955727710236D+00, 0.410873840972387D+00, 0.726271784259897D+00, &  !11
  0.110386324646491D+01, 0.153229503457537D+01, 0.200578290246814D+01, 0.252435214151921D+01, &
  0.309535170986922D+01, 0.373947860994358D+01, 0.451783596718736D+01, &
  0.298897007696644D-01, 0.154204878265825D+00, 0.366143962974312D+00, 0.650881015845205D+00, &  !12
  0.994366869880792D+00, 0.138589120364956D+01, 0.181884860842823D+01, 0.229084273867285D+01, &
  0.280409679339362D+01, 0.336727070416293D+01, 0.400168347567348D+01, 0.476821628798986D+01, &
  0.266511266223847D-01, 0.137891855469939D+00, 0.328828675158344D+00, 0.587378531473706D+00, &  !13
  0.901480884535392D+00, 0.126129650258238D+01, 0.166003713190544D+01, 0.209410900418749D+01, &
  0.256320702525468D+01, 0.307091234102964D+01, 0.362669201180255D+01, 0.425220740047932D+01, &
  0.500800834323777D+01, &
  0.239567892441074D-01, 0.124240344109914D+00, 0.297338568916482D+00, 0.533329214443288D+00, &  !14
  0.821873189116935D+00, 0.115406707387441D+01, 0.152327479144572D+01, 0.192533821047710D+01, &
  0.235860076664477D+01, 0.282409375484497D+01, 0.332626935870514D+01, 0.387510499098758D+01, &
  0.449243807160692D+01, 0.523843136267529D+01, &
  0.216869426988512D-01, 0.112684196952642D+00, 0.270492620388998D+00, 0.486902289234274D+00, &  !15
  0.753043574230519D+00, 0.106093087200583D+01, 0.140425480937058D+01, 0.177864621852044D+01, &
  0.218170796284551D+01, 0.261306067251355D+01, 0.307461793949750D+01, 0.357140797746735D+01, &
  0.411373591845599D+01, 0.472351289497065D+01, 0.546048877386486D+01/
DATA xpt(120:209)/ &
  0.197536584600773D-01, 0.102802245237917D+00, 0.247397669452455D+00, 0.446696225961683D+00, &  !16
  0.693073720302000D+00, 0.979404170330730D+00, 0.129978932127704D+01, 0.164985424039743D+01, &
  0.202680815216887D+01, 0.242945049160214D+01, 0.285826652854327D+01, 0.331576927503870D+01, &
  0.380737711675590D+01, 0.434360634547017D+01, 0.494637720404839D+01, 0.567501793404192D+01, &
  0.180910833291262D-01, 0.942756133296731D-01, 0.227367679140888D+00, 0.411621543218912D+00, &  !17
  0.640465304258667D+00, 0.907557245160356D+00, 0.120744250530777D+01, 0.153586090376516D+01, &
  0.188983689868465D+01, 0.226768435443910D+01, 0.266903218379336D+01, 0.309496584156224D+01, &
  0.354841560757166D+01, 0.403506890899052D+01, 0.456557652510097D+01, 0.516182632161921D+01, &
  0.588272607689950D+01, &
  0.166490322202372D-01, 0.868590621084087D-01, 0.209868348130364D+00, 0.380820011068889D+00, &  !18
  0.594030720753567D+00, 0.843863021635420D+00, 0.112530737639273D+01, 0.143428771048619D+01, &
  0.176777393564494D+01, 0.212379869512189D+01, 0.250145945877649D+01, 0.290097223563891D+01, &
  0.332384856894201D+01, 0.377331513435006D+01, 0.425524711266144D+01, 0.478037806093043D+01, &
  0.537053657981972D+01, 0.608421686390343D+01, &
  0.153886660136205D-01, 0.803613182824543D-01, 0.194478326810452D+00, 0.353607595856074D+00, &  !19
  0.552816862912666D+00, 0.787094922636796D+00, 0.105186215045210D+01, 0.134326400045299D+01, &
  0.165829830732126D+01, 0.199484680375345D+01, 0.235167370676724D+01, 0.272844153701442D+01, &
  0.312579015345213D+01, 0.354553959918365D+01, 0.399113055980209D+01, 0.446857002866801D+01, &
  0.498863874530620D+01, 0.557308861757198D+01, 0.628001037094314D+01, &
  0.142795096999168D-01, 0.746313003921914D-01, 0.180861563058038D+00, 0.329433356064288D+00, &  !20
  0.516050543061530D+00, 0.736255457580891D+00, 0.985873575037527D+00, 0.126128901610277D+01, &
  0.155957964520966D+01, 0.187856191930298D+01, 0.221679416538765D+01, 0.257357782082628D+01, &
  0.294898974678723D+01, 0.334398137986164D+01, 0.376059993291077D+01, 0.420244260042114D+01, &
  0.467560884779448D+01, 0.519090168697500D+01, 0.576998516556776D+01, 0.647055838706458D+01/
DATA wt(1:119)/ &
  0.640529179684379D+00, 0.245697745768379D+00, &                                                !2
  0.446029770466658D+00, 0.396468266998335D+00, 0.437288879877644D-01, &                         !3
  0.325302999756919D+00, 0.421107101852062D+00, 0.133442500357520D+00, 0.637432348625728D-02, &  !4
  0.248406152028443D+00, 0.392331066652399D+00, 0.211418193076057D+00, 0.332466603513439D-01, &  !5
  0.824853344515628D-03, &
  0.196849675488598D+00, 0.349154201525395D+00, 0.257259520584421D+00, 0.760131375840057D-01, &  !6
  0.685191862513597D-02, 0.984716452019267D-04, &
  0.160609965149261D+00, 0.306319808158099D+00, 0.275527141784906D+00, 0.120630193130784D+00, &  !7
  0.218922863438067D-01, 0.123644672831056D-02, 0.110841575911059D-04, &
  0.134109188453360D+00, 0.268330754472639D+00, 0.275953397988422D+00, 0.157448282618790D+00, &  !8
  0.448141099174629D-01, 0.536793575602533D-02, 0.202063649132411D-03, 0.119259692659534D-05, &
  0.114088970242111D+00, 0.235940791223676D+00, 0.266425473630252D+00, 0.183251679101671D+00, &  !9
  0.713440493066984D-01, 0.139814184155624D-01, 0.116385272078542D-02, 0.305670214897907D-04, &
  0.123790511337534D-06, &
  0.985520975190362D-01, 0.208678066608076D+00, 0.252051688403725D+00, 0.198684340038460D+00, &  !10
  0.971984227601550D-01, 0.270244164355872D-01, 0.380464962250372D-02, 0.228886243045297D-03, &
  0.434534479845945D-05, 0.124773714818325D-07, &
  0.862207055348204D-01, 0.185767318954432D+00, 0.235826124129156D+00, 0.205850326842101D+00, &  !11
  0.119581170616438D+00, 0.431443275887789D-01, 0.886764989495983D-02, 0.927141875111555D-03, &
  0.415719321683689D-04, 0.586857646864747D-06, 0.122714514000882D-08, &
  0.762461467930431D-01, 0.166446068879474D+00, 0.219394898128707D+00, 0.207016508679094D+00, &  !12
  0.137264362796474D+00, 0.605056743489164D-01, 0.165538019564075D-01, 0.258608378835667D-02, &
  0.206237541067489D-03, 0.706650986752706D-05, 0.759131547256598D-07, 0.118195417166772D-09, &
  0.680463904418563D-01, 0.150057211706640D+00, 0.203606639691744D+00, 0.204104355198729D+00, &  !13
  0.150119228251119D+00, 0.774536315632889D-01, 0.264891667292508D-01, 0.562343031211025D-02, &
  0.683241179366847D-03, 0.424853319211805D-04, 0.113557101360564D-05, 0.946453645906370D-08, &
  0.111810461699377D-10, &
  0.612109812196030D-01, 0.136062058638519D+00, 0.188856801742046D+00, 0.198577828793095D+00, &  !14
  0.158617339225760D+00, 0.928167848032935D-01, 0.379316402929447D-01, 0.102563915568667D-01, &
  0.172277191301072D-02, 0.165956353055286D-03, 0.819589395640479D-05, 0.173876626376522D-06, &
  0.114293991639056D-08, 0.104120023691655D-11, &
  0.554433542997126D-01, 0.124027715695603D+00, 0.175290920953716D+00, 0.191488332592396D+00, &  !15
  0.163473809496946D+00, 0.105937660378078D+00, 0.500270400435842D-01, 0.164429780046165D-01, &
  0.357320679306892D-02, 0.482896942477554D-03, 0.374909051856946D-04, 0.149368597328174D-05, &
  0.255270858132649D-07, 0.134217892411487D-09, 0.956229145184876D-13/
DATA wt(120:209)/ &
  0.505246320213779D-01, 0.113608556894151D+00, 0.162921292314545D+00, 0.183562801116246D+00, &  !16
  0.165438637755610D+00, 0.116572490553503D+00, 0.619996960991566D-01, 0.239197096186835D-01, &
  0.640991442405013D-02, 0.113569531068878D-02, 0.125286221329562D-03, 0.795049571962246D-05, &
  0.259000761941506D-06, 0.361154913974278D-08, 0.153767791618984D-10, 0.867420445249463D-14, &
  0.462902182569444D-01, 0.104528587434386D+00, 0.151692844038654D+00, 0.175289523252223D+00, &  !17
  0.165194113560473D+00, 0.124758912731550D+00, 0.732500435683077D-01, 0.322890873617879D-01, &
  0.102921626932112D-01, 0.227638568842155D-02, 0.333037962863638D-03, 0.303707486255998D-04, &
  0.159468239105431D-05, 0.429767400483309D-07, 0.494454281458509D-09, 0.172326283208696D-11, &
  0.778200825411988D-15, &
  0.426142619814402D-01, 0.965666206000141D-01, 0.141519805800333D+00, 0.166987766912260D+00, &  !18
  0.163315996982125D+00, 0.130699466898224D+00, 0.833782286268739D-01, 0.411112305586427D-01, &
  0.151576539874531D-01, 0.403326905747473D-02, 0.744434183403738D-03, 0.909462923551198D-04, &
  0.693223179081282D-05, 0.304417379930170D-06, 0.685706226732795D-08, 0.657362956131057D-10, &
  0.189340654173429D-12, 0.691219900671972D-16, &
  0.393990986682379D-01, 0.895444367628123D-01, 0.132305744000463D+00, 0.158859868228014D+00, &  !19
  0.160268450441027D+00, 0.134675758617491D+00, 0.921641776964749D-01, 0.499738916193989D-01, &
  0.208459843479637D-01, 0.648435446662631D-02, 0.145409393575714D-02, 0.226202869969883D-03, &
  0.233075944706917D-04, 0.149954388953745D-05, 0.555953942746990D-07, 0.105623488322490D-08, &
  0.851137824908803D-11, 0.204347153060500D-13, 0.608418799770205D-17, &
  0.365679216320084D-01, 0.833175344016776D-01, 0.123954178541192D+00, 0.151028580007023D+00, &  !20
  0.156414467004551D+00, 0.136992234955135D+00, 0.995292472472200D-01, 0.585337706979062D-01, &
  0.271357803755605D-01, 0.964473365921540D-02, 0.255171509156881D-02, 0.486448997660902D-03, &
  0.643514517998962D-04, 0.564239148952609D-05, 0.309083652232771D-06, 0.975641105342091D-08, &
  0.157608671739272D-09, 0.107594574673723D-11, 0.216986354627585D-14, 0.531122306167733D-18/

!     Check for permissible value of NPT.

hh = 0.d0
ier = 1
IF(npt < 2 .OR. npt > 20) RETURN

!     IPOS = position of first abscissa & weight to use.

ier = 0
ipos = npt*(npt-1)/2

!     Evaluate approximation = Sum w(i).f(xi)

DO i = 1, npt
  hh = hh + wt(ipos)*func(xpt(ipos))
  ipos = ipos + 1
END DO
RETURN

END FUNCTION hh
