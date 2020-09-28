function [S, T, R]=TrigNodeVect(n)

% TrigNodeVect: Generate unequally spaced nodes on a unit triangle
% 
% Calling Sequences:
% 
%       [S, T]=TrigNodeVect(n)
% 
%       [S, T, R]=TrigNodeVect(n)
% 
% INPUTS:
% 
%       n - The number of nodes on edges. All three edges have the same number.            
%
% OUTPUT: 
%
%      S, T  -  The nodes on the triangle.
%
%      R = 1-S-T;    The third coordinates.
%

switch n
    case 2
        S=[0;1;0];
        T=[0;0;1];
    case 3
        S=[0;0.500000000000000;1;0;0.500000000000000;0];
        T=[0;0;0;0.500000000000000;0.500000000000000;1];
    case 4
        S=[0;0.276393202250021;0.723606797749979;1;0;0.333333333333333;0.723606797749979;0;0.276393202250021;0];
        T=[0;0;0;0;0.276393202250021;0.333333333333333;0.276393202250021;0.723606797749979;0.723606797749979;1];
    case 5
        S=[0;0.172673164646011;0.500000000000000;0.827326835353989;1;0;0.204263221667533;0.591473556664935;0.827326835353989;0;0.204263221667533;0.500000000000000;0;0.172673164646011;0];
        T=[0;0;0;0;0;0.172673164646011;0.204263221667533;0.204263221667533;0.172673164646011;0.500000000000000;0.591473556664935;0.500000000000000;0.827326835353989;0.827326835353989;1];
    case 6
        S=[0;0.117472338035268;0.357384241759677;0.642615758240323;0.882527661964732;1;0;0.133862391058592;0.429424071138550;0.732275217882817;0.882527661964732;0;0.141151857722900;0.429424071138550;0.642615758240323;0;0.133862391058592;0.357384241759677;0;0.117472338035268;0];
        T=[0;0;0;0;0;0;0.117472338035268;0.133862391058592;0.141151857722900;0.133862391058592;0.117472338035268;0.357384241759677;0.429424071138550;0.429424071138550;0.357384241759677;0.642615758240323;0.732275217882817;0.642615758240323;0.882527661964732;0.882527661964732;1];
    case 7
        S=[0;0.0848880518607165;0.265575603264643;0.500000000000000;0.734424396735357;0.915111948139284;1;0;0.0938818899324123;0.312271549365030;0.587914600449680;0.812236220135175;0.915111948139284;0;0.0998138501852897;0.333333333333333;0.587914600449680;0.734424396735357;0;0.0998138501852897;0.312271549365030;0.500000000000000;0;0.0938818899324123;0.265575603264643;0;0.0848880518607165;0];
        T=[0;0;0;0;0;0;0;0.0848880518607165;0.0938818899324123;0.0998138501852897;0.0998138501852897;0.0938818899324123;0.0848880518607165;0.265575603264643;0.312271549365030;0.333333333333333;0.312271549365030;0.265575603264643;0.500000000000000;0.587914600449680;0.587914600449680;0.500000000000000;0.734424396735357;0.812236220135175;0.734424396735357;0.915111948139284;0.915111948139284;1];
    case 8
        S=[0;0.0641299257451967;0.204149909283429;0.395350391048761;0.604649608951239;0.795850090716571;0.935870074254803;1;0;0.0693964244038337;0.233867594559152;0.462489692311687;0.692667217403640;0.861207151192333;0.935870074254803;0;0.0734651880372086;0.254028315852829;0.491943368294341;0.692667217403640;0.795850090716571;0;0.0750206153766251;0.254028315852829;0.462489692311687;0.604649608951239;0;0.0734651880372086;0.233867594559152;0.395350391048761;0;0.0693964244038337;0.204149909283429;0;0.0641299257451967;0];
        T=[0;0;0;0;0;0;0;0;0.0641299257451967;0.0693964244038337;0.0734651880372086;0.0750206153766251;0.0734651880372086;0.0693964244038337;0.0641299257451967;0.204149909283429;0.233867594559152;0.254028315852829;0.254028315852829;0.233867594559152;0.204149909283429;0.395350391048761;0.462489692311687;0.491943368294341;0.462489692311687;0.395350391048761;0.604649608951239;0.692667217403640;0.692667217403640;0.604649608951239;0.795850090716571;0.861207151192333;0.795850090716571;0.935870074254803;0.935870074254803;1];
    case 9
        S=[0;0.0501210022942700;0.161406860244631;0.318441268086911;0.500000000000000;0.681558731913089;0.838593139755369;0.949878997705730;1;0;0.0533863720337145;0.180729238628503;0.366630325707285;0.575663964519859;0.763149661126985;0.893227255932571;0.949878997705730;0;0.0561211002445120;0.196164522084847;0.398904544536864;0.607670955830306;0.763149661126985;0.838593139755369;0;0.0577057097728568;0.202190910926272;0.398904544536864;0.575663964519859;0.681558731913089;0;0.0577057097728568;0.196164522084847;0.366630325707285;0.500000000000000;0;0.0561211002445120;0.180729238628503;0.318441268086911;0;0.0533863720337145;0.161406860244631;0;0.0501210022942700;0];
        T=[0;0;0;0;0;0;0;0;0;0.0501210022942700;0.0533863720337145;0.0561211002445120;0.0577057097728568;0.0577057097728568;0.0561211002445120;0.0533863720337145;0.0501210022942700;0.161406860244631;0.180729238628503;0.196164522084847;0.202190910926272;0.196164522084847;0.180729238628503;0.161406860244631;0.318441268086911;0.366630325707285;0.398904544536864;0.398904544536864;0.366630325707285;0.318441268086911;0.500000000000000;0.575663964519859;0.607670955830306;0.575663964519859;0.500000000000000;0.681558731913089;0.763149661126985;0.763149661126985;0.681558731913089;0.838593139755369;0.893227255932571;0.838593139755369;0.949878997705730;0.949878997705730;1];
    case 10
        S=[0;0.0402330459167706;0.130613067447247;0.261037525094778;0.417360521166807;0.582639478833194;0.738962474905222;0.869386932552753;0.959766954083229;1;0;0.0423571277701257;0.143561040345549;0.295321367983800;0.477008481669465;0.659161505105111;0.812217519633416;0.915285744459749;0.959766954083229;0;0.0442214400210348;0.154779464078739;0.322662471097680;0.515889725239321;0.690441071842523;0.812217519633416;0.869386932552753;0;0.0455171269110895;0.161447803662999;0.333333333333333;0.515889725239321;0.659161505105111;0.738962474905222;0;0.0459830366610701;0.161447803662999;0.322662471097680;0.477008481669465;0.582639478833194;0;0.0455171269110895;0.154779464078739;0.295321367983800;0.417360521166807;0;0.0442214400210348;0.143561040345549;0.261037525094778;0;0.0423571277701257;0.130613067447247;0;0.0402330459167706;0];
        T=[0;0;0;0;0;0;0;0;0;0;0.0402330459167706;0.0423571277701257;0.0442214400210348;0.0455171269110895;0.0459830366610701;0.0455171269110895;0.0442214400210348;0.0423571277701257;0.0402330459167706;0.130613067447247;0.143561040345549;0.154779464078739;0.161447803662999;0.161447803662999;0.154779464078739;0.143561040345549;0.130613067447247;0.261037525094778;0.295321367983800;0.322662471097680;0.333333333333333;0.322662471097680;0.295321367983800;0.261037525094778;0.417360521166807;0.477008481669465;0.515889725239321;0.515889725239321;0.477008481669465;0.417360521166807;0.582639478833194;0.659161505105111;0.690441071842523;0.659161505105111;0.582639478833194;0.738962474905222;0.812217519633416;0.812217519633416;0.738962474905222;0.869386932552753;0.915285744459749;0.869386932552753;0.959766954083229;0.959766954083229;1];
    case 11
        S=[0;0.0329992847959704;0.107758263168428;0.217382336501897;0.352120932206530;0.500000000000000;0.647879067793470;0.782617663498103;0.892241736831572;0.967000715204030;1;0;0.0344373791988935;0.116700407191735;0.242003617815889;0.397822719945324;0.564895016965348;0.721259513704035;0.847561916091014;0.931125241602213;0.967000715204030;0;0.0357376767172513;0.124807520331361;0.263448843250168;0.433646400100095;0.605957336482739;0.750384959337279;0.847561916091014;0.892241736831572;0;0.0367368684800755;0.130593820267093;0.276256593159807;0.447486813680386;0.605957336482739;0.721259513704035;0.782617663498103;0;0.0372822630893281;0.132707199799810;0.276256593159807;0.433646400100095;0.564895016965348;0.647879067793470;0;0.0372822630893281;0.130593820267093;0.263448843250168;0.397822719945324;0.500000000000000;0;0.0367368684800755;0.124807520331361;0.242003617815889;0.352120932206530;0;0.0357376767172513;0.116700407191735;0.217382336501897;0;0.0344373791988935;0.107758263168428;0;0.0329992847959704;0];
        T=[0;0;0;0;0;0;0;0;0;0;0;0.0329992847959704;0.0344373791988935;0.0357376767172513;0.0367368684800755;0.0372822630893281;0.0372822630893281;0.0367368684800755;0.0357376767172513;0.0344373791988935;0.0329992847959704;0.107758263168428;0.116700407191735;0.124807520331361;0.130593820267093;0.132707199799810;0.130593820267093;0.124807520331361;0.116700407191735;0.107758263168428;0.217382336501897;0.242003617815889;0.263448843250168;0.276256593159807;0.276256593159807;0.263448843250168;0.242003617815889;0.217382336501897;0.352120932206530;0.397822719945324;0.433646400100095;0.447486813680386;0.433646400100095;0.397822719945324;0.352120932206530;0.500000000000000;0.564895016965348;0.605957336482739;0.605957336482739;0.564895016965348;0.500000000000000;0.647879067793470;0.721259513704035;0.750384959337279;0.721259513704035;0.647879067793470;0.782617663498103;0.847561916091014;0.847561916091014;0.782617663498103;0.892241736831572;0.931125241602213;0.892241736831572;0.967000715204030;0.967000715204030;1];
    case 12
        S=[0;0.0275503638885589;0.0903603391779967;0.183561923484070;0.300234529517326;0.431723533572536;0.568276466427464;0.699765470482675;0.816438076515930;0.909639660822003;0.972449636111441;1;0;0.0285572825748119;0.0967094313735544;0.201522017527657;0.335060234684335;0.484539594094790;0.634193696883831;0.768232031628810;0.873804401907343;0.942885434850376;0.972449636111441;0;0.0294861667191023;0.102625508984075;0.217955593073041;0.365107394762556;0.525007749287737;0.674753412465733;0.794748982031849;0.873804401907343;0.909639660822003;0;0.0302459508435329;0.107290994461226;0.229783470479469;0.382937072911091;0.540433059041061;0.674753412465733;0.768232031628810;0.816438076515930;0;0.0307460684318346;0.109884855949707;0.234125854177818;0.382937072911091;0.525007749287737;0.634193696883831;0.699765470482675;0;0.0309208118104204;0.109884855949707;0.229783470479469;0.365107394762556;0.484539594094790;0.568276466427464;0;0.0307460684318346;0.107290994461226;0.217955593073041;0.335060234684335;0.431723533572536;0;0.0302459508435329;0.102625508984075;0.201522017527657;0.300234529517326;0;0.0294861667191023;0.0967094313735544;0.183561923484070;0;0.0285572825748119;0.0903603391779967;0;0.0275503638885589;0];
        T=[0;0;0;0;0;0;0;0;0;0;0;0;0.0275503638885589;0.0285572825748119;0.0294861667191023;0.0302459508435329;0.0307460684318346;0.0309208118104204;0.0307460684318346;0.0302459508435329;0.0294861667191023;0.0285572825748119;0.0275503638885589;0.0903603391779967;0.0967094313735544;0.102625508984075;0.107290994461226;0.109884855949707;0.109884855949707;0.107290994461226;0.102625508984075;0.0967094313735544;0.0903603391779967;0.183561923484070;0.201522017527657;0.217955593073041;0.229783470479469;0.234125854177818;0.229783470479469;0.217955593073041;0.201522017527657;0.183561923484070;0.300234529517326;0.335060234684335;0.365107394762556;0.382937072911091;0.382937072911091;0.365107394762556;0.335060234684335;0.300234529517326;0.431723533572536;0.484539594094790;0.525007749287737;0.540433059041061;0.525007749287737;0.484539594094790;0.431723533572536;0.568276466427464;0.634193696883831;0.674753412465733;0.674753412465733;0.634193696883831;0.568276466427464;0.699765470482675;0.768232031628810;0.794748982031849;0.768232031628810;0.699765470482675;0.816438076515930;0.873804401907343;0.873804401907343;0.816438076515930;0.909639660822003;0.942885434850376;0.909639660822003;0.972449636111441;0.972449636111441;1];
    case 13
        S=[0;0.0233450766789181;0.0768262176740638;0.156905765459121;0.258545089454332;0.375356534946880;0.500000000000000;0.624643465053120;0.741454910545668;0.843094234540879;0.923173782325936;0.976654923321082;1;0;0.0240704657839033;0.0814470743610822;0.170234117561123;0.285201872892434;0.417665363109609;0.556358187781008;0.689046102167943;0.804437759433303;0.893803715619238;0.951859068432193;0.976654923321082;0;0.0247492100196793;0.0858290550585865;0.182793864608609;0.309497210699139;0.453581528497662;0.598536238596414;0.727704253858791;0.828341889882827;0.893803715619238;0.923173782325936;0;0.0253281230055739;0.0895018815326008;0.192803566302484;0.326938130282260;0.474649756387801;0.614392867395032;0.727704253858791;0.804437759433303;0.843094234540879;0;0.0257520249396226;0.0919665507044471;0.198412113329939;0.333333333333333;0.474649756387801;0.598536238596414;0.689046102167943;0.741454910545668;0;0.0259764491093830;0.0928369430046755;0.198412113329939;0.326938130282260;0.453581528497662;0.556358187781008;0.624643465053120;0;0.0259764491093830;0.0919665507044471;0.192803566302484;0.309497210699139;0.417665363109609;0.500000000000000;0;0.0257520249396226;0.0895018815326008;0.182793864608609;0.285201872892434;0.375356534946880;0;0.0253281230055739;0.0858290550585865;0.170234117561123;0.258545089454332;0;0.0247492100196793;0.0814470743610822;0.156905765459121;0;0.0240704657839033;0.0768262176740638;0;0.0233450766789181;0];
        T=[0;0;0;0;0;0;0;0;0;0;0;0;0;0.0233450766789181;0.0240704657839033;0.0247492100196793;0.0253281230055739;0.0257520249396226;0.0259764491093830;0.0259764491093830;0.0257520249396226;0.0253281230055739;0.0247492100196793;0.0240704657839033;0.0233450766789181;0.0768262176740638;0.0814470743610822;0.0858290550585865;0.0895018815326008;0.0919665507044471;0.0928369430046755;0.0919665507044471;0.0895018815326008;0.0858290550585865;0.0814470743610822;0.0768262176740638;0.156905765459121;0.170234117561123;0.182793864608609;0.192803566302484;0.198412113329939;0.198412113329939;0.192803566302484;0.182793864608609;0.170234117561123;0.156905765459121;0.258545089454332;0.285201872892434;0.309497210699139;0.326938130282260;0.333333333333333;0.326938130282260;0.309497210699139;0.285201872892434;0.258545089454332;0.375356534946880;0.417665363109609;0.453581528497662;0.474649756387801;0.474649756387801;0.453581528497662;0.417665363109609;0.375356534946880;0.500000000000000;0.556358187781008;0.598536238596414;0.614392867395032;0.598536238596414;0.556358187781008;0.500000000000000;0.624643465053120;0.689046102167943;0.727704253858791;0.727704253858791;0.689046102167943;0.624643465053120;0.741454910545668;0.804437759433303;0.828341889882827;0.804437759433303;0.741454910545668;0.843094234540879;0.893803715619238;0.893803715619238;0.843094234540879;0.923173782325936;0.951859068432193;0.923173782325936;0.976654923321082;0.976654923321082;1];
    case 14
        S=[0;0.0200324773663695;0.0660994730848264;0.135565700454337;0.224680298535677;0.328637993328644;0.441834065558148;0.558165934441852;0.671362006671356;0.775319701464324;0.864434299545663;0.933900526915174;0.979967522633630;1;0;0.0205679541460498;0.0695369450420155;0.145625842484039;0.245264145252982;0.362400529430712;0.488916419233133;0.615508962013507;0.732868123261064;0.832855097873714;0.909388798043818;0.958864091707900;0.979967522633630;0;0.0210742569141664;0.0728353882943967;0.155282336595347;0.264658030462049;0.392839082814729;0.528148578659904;0.657481309412330;0.769004702132488;0.854329223411207;0.909388798043818;0.933900526915174;0;0.0215190596422480;0.0757129612721649;0.163470560769421;0.280122033281361;0.414506450242159;0.550860300719540;0.673058878461158;0.769004702132488;0.832855097873714;0.864434299545663;0;0.0218677314859547;0.0778606601256210;0.169017665999099;0.288792680742725;0.422414638514549;0.550860300719540;0.657481309412330;0.732868123261064;0.775319701464324;0;0.0220905085557809;0.0790123385253670;0.170987099515683;0.288792680742725;0.414506450242159;0.528148578659904;0.615508962013507;0.671362006671356;0;0.0221671615337346;0.0790123385253670;0.169017665999099;0.280122033281361;0.392839082814729;0.488916419233133;0.558165934441852;0;0.0220905085557809;0.0778606601256210;0.163470560769421;0.264658030462049;0.362400529430712;0.441834065558148;0;0.0218677314859547;0.0757129612721649;0.155282336595347;0.245264145252982;0.328637993328644;0;0.0215190596422480;0.0728353882943967;0.145625842484039;0.224680298535677;0;0.0210742569141664;0.0695369450420155;0.135565700454337;0;0.0205679541460498;0.0660994730848264;0;0.0200324773663695;0];
        T=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0200324773663695;0.0205679541460498;0.0210742569141664;0.0215190596422480;0.0218677314859547;0.0220905085557809;0.0221671615337346;0.0220905085557809;0.0218677314859547;0.0215190596422480;0.0210742569141664;0.0205679541460498;0.0200324773663695;0.0660994730848264;0.0695369450420155;0.0728353882943967;0.0757129612721649;0.0778606601256210;0.0790123385253670;0.0790123385253670;0.0778606601256210;0.0757129612721649;0.0728353882943967;0.0695369450420155;0.0660994730848264;0.135565700454337;0.145625842484039;0.155282336595347;0.163470560769421;0.169017665999099;0.170987099515683;0.169017665999099;0.163470560769421;0.155282336595347;0.145625842484039;0.135565700454337;0.224680298535677;0.245264145252982;0.264658030462049;0.280122033281361;0.288792680742725;0.288792680742725;0.280122033281361;0.264658030462049;0.245264145252982;0.224680298535677;0.328637993328644;0.362400529430712;0.392839082814729;0.414506450242159;0.422414638514549;0.414506450242159;0.392839082814729;0.362400529430712;0.328637993328644;0.441834065558148;0.488916419233133;0.528148578659904;0.550860300719540;0.550860300719540;0.528148578659904;0.488916419233133;0.441834065558148;0.558165934441852;0.615508962013507;0.657481309412330;0.673058878461158;0.657481309412330;0.615508962013507;0.558165934441852;0.671362006671356;0.732868123261064;0.769004702132488;0.769004702132488;0.732868123261064;0.671362006671356;0.775319701464324;0.832855097873714;0.854329223411207;0.832855097873714;0.775319701464324;0.864434299545663;0.909388798043818;0.909388798043818;0.864434299545663;0.933900526915174;0.958864091707900;0.933900526915174;0.979967522633630;0.979967522633630;1];
    case 15
        S=[0;0.0173770367480807;0.0574589778885118;0.118240155024092;0.196873397265077;0.289680972643164;0.392323022318103;0.500000000000000;0.607676977681897;0.710319027356836;0.803126602734923;0.881759844975908;0.942541022111488;0.982622963251919;1;0;0.0177807469066201;0.0600660951133402;0.125955724784499;0.212935215244247;0.316682948325886;0.431266348076119;0.549631711042489;0.664320252608101;0.768270051126619;0.855533328086199;0.921768410468651;0.964438506186760;0.982622963251919;0;0.0181654944180083;0.0625884412208126;0.133451164238141;0.228388934084629;0.341951731188402;0.465883722129905;0.590221249377029;0.704954042184125;0.801698045490203;0.874823117558375;0.921768410468652;0.942541022111488;0;0.0185109471293015;0.0648507902716558;0.140068867191844;0.241528798916638;0.361990740761383;0.490254158465614;0.613411467145648;0.719862265616311;0.801698045490204;0.855533328086199;0.881759844975908;0;0.0187947336291338;0.0666570237312470;0.145059733937715;0.250452813535622;0.373187043623925;0.499094372928756;0.613411467145648;0.704954042184125;0.768270051126619;0.803126602734923;0;0.0189967990660133;0.0678270194345691;0.147755100773002;0.253625912752149;0.373187043623925;0.490254158465614;0.590221249377029;0.664320252608101;0.710319027356836;0;0.0191019408813916;0.0682325557401901;0.147755100773002;0.250452813535622;0.361990740761383;0.465883722129905;0.549631711042489;0.607676977681897;0;0.0191019408813916;0.0678270194345691;0.145059733937715;0.241528798916638;0.341951731188402;0.431266348076119;0.500000000000000;0;0.0189967990660133;0.0666570237312470;0.140068867191844;0.228388934084629;0.316682948325886;0.392323022318103;0;0.0187947336291338;0.0648507902716558;0.133451164238141;0.212935215244247;0.289680972643164;0;0.0185109471293015;0.0625884412208126;0.125955724784499;0.196873397265077;0;0.0181654944180083;0.0600660951133402;0.118240155024092;0;0.0177807469066201;0.0574589778885118;0;0.0173770367480807;0];
        T=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0173770367480807;0.0177807469066201;0.0181654944180083;0.0185109471293015;0.0187947336291338;0.0189967990660133;0.0191019408813916;0.0191019408813916;0.0189967990660133;0.0187947336291338;0.0185109471293015;0.0181654944180083;0.0177807469066201;0.0173770367480807;0.0574589778885118;0.0600660951133402;0.0625884412208126;0.0648507902716558;0.0666570237312470;0.0678270194345691;0.0682325557401901;0.0678270194345691;0.0666570237312470;0.0648507902716558;0.0625884412208126;0.0600660951133402;0.0574589778885118;0.118240155024092;0.125955724784499;0.133451164238141;0.140068867191844;0.145059733937715;0.147755100773002;0.147755100773002;0.145059733937715;0.140068867191844;0.133451164238141;0.125955724784499;0.118240155024092;0.196873397265077;0.212935215244247;0.228388934084629;0.241528798916638;0.250452813535622;0.253625912752149;0.250452813535622;0.241528798916638;0.228388934084629;0.212935215244247;0.196873397265077;0.289680972643164;0.316682948325886;0.341951731188402;0.361990740761383;0.373187043623925;0.373187043623925;0.361990740761383;0.341951731188402;0.316682948325886;0.289680972643164;0.392323022318103;0.431266348076119;0.465883722129905;0.490254158465614;0.499094372928756;0.490254158465614;0.465883722129905;0.431266348076119;0.392323022318103;0.500000000000000;0.549631711042489;0.590221249377029;0.613411467145648;0.613411467145648;0.590221249377029;0.549631711042489;0.500000000000000;0.607676977681897;0.664320252608101;0.704954042184125;0.719862265616311;0.704954042184125;0.664320252608101;0.607676977681897;0.710319027356836;0.768270051126619;0.801698045490203;0.801698045490204;0.768270051126619;0.710319027356836;0.803126602734923;0.855533328086199;0.874823117558375;0.855533328086199;0.803126602734923;0.881759844975908;0.921768410468651;0.921768410468652;0.881759844975908;0.942541022111488;0.964438506186760;0.942541022111488;0.982622963251919;0.982622963251919;1];
    case 16
        S=[0;0.0152159768648910;0.0503997334532639;0.103995854069092;0.173805648558753;0.256970289056431;0.350084765549618;0.449336863239025;0.550663136760975;0.649915234450382;0.743029710943569;0.826194351441247;0.896004145930908;0.949600266546736;0.984784023135109;1;0;0.0155259965438665;0.0524112861498674;0.110001241802560;0.186476475894754;0.278678962080471;0.382203680198341;0.491675156854688;0.601184336282628;0.704819625809520;0.797198267897582;0.873904113219911;0.931765437358897;0.968948006912267;0.984784023135109;0;0.0158232764912355;0.0543690070588537;0.115882463138564;0.198834857896559;0.299487597046781;0.411951006533964;0.528742725475426;0.641773724956084;0.743507500235625;0.827957170600959;0.891261985882293;0.931765437358897;0.949600266546736;0;0.0160946449775297;0.0561603662604766;0.121220436514134;0.209792476703622;0.317128638346459;0.435339262304032;0.554529428756521;0.664679107085624;0.757559126971731;0.827957170600959;0.873904113219911;0.896004145930908;0;0.0163252562076647;0.0576576418678164;0.125528416210754;0.218089026619997;0.329085953506098;0.448331903668268;0.563821946760006;0.664679107085624;0.743507500235625;0.797198267897582;0.826194351441247;0;0.0165014121100089;0.0587386779971345;0.128341932897020;0.222582142825634;0.333333333333333;0.448331903668268;0.554529428756521;0.641773724956084;0.704819625809520;0.743029710943569;0;0.0166119835190312;0.0593062679906101;0.129321475391936;0.222582142825634;0.329085953506098;0.435339262304032;0.528742725475426;0.601184336282628;0.649915234450382;0;0.0166496862906236;0.0593062679906101;0.128341932897020;0.218089026619997;0.317128638346459;0.411951006533964;0.491675156854688;0.550663136760975;0;0.0166119835190312;0.0587386779971345;0.125528416210754;0.209792476703622;0.299487597046781;0.382203680198341;0.449336863239025;0;0.0165014121100089;0.0576576418678164;0.121220436514134;0.198834857896559;0.278678962080471;0.350084765549618;0;0.0163252562076647;0.0561603662604766;0.115882463138564;0.186476475894754;0.256970289056431;0;0.0160946449775297;0.0543690070588537;0.110001241802560;0.173805648558753;0;0.0158232764912355;0.0524112861498674;0.103995854069092;0;0.0155259965438665;0.0503997334532639;0;0.0152159768648910;0];
        T=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0152159768648910;0.0155259965438665;0.0158232764912355;0.0160946449775297;0.0163252562076647;0.0165014121100089;0.0166119835190312;0.0166496862906236;0.0166119835190312;0.0165014121100089;0.0163252562076647;0.0160946449775297;0.0158232764912355;0.0155259965438665;0.0152159768648910;0.0503997334532639;0.0524112861498674;0.0543690070588537;0.0561603662604766;0.0576576418678164;0.0587386779971345;0.0593062679906101;0.0593062679906101;0.0587386779971345;0.0576576418678164;0.0561603662604766;0.0543690070588537;0.0524112861498674;0.0503997334532639;0.103995854069092;0.110001241802560;0.115882463138564;0.121220436514134;0.125528416210754;0.128341932897020;0.129321475391936;0.128341932897020;0.125528416210754;0.121220436514134;0.115882463138564;0.110001241802560;0.103995854069092;0.173805648558753;0.186476475894754;0.198834857896559;0.209792476703622;0.218089026619997;0.222582142825634;0.222582142825634;0.218089026619997;0.209792476703622;0.198834857896559;0.186476475894754;0.173805648558753;0.256970289056431;0.278678962080471;0.299487597046781;0.317128638346459;0.329085953506098;0.333333333333333;0.329085953506098;0.317128638346459;0.299487597046781;0.278678962080471;0.256970289056431;0.350084765549618;0.382203680198341;0.411951006533964;0.435339262304032;0.448331903668268;0.448331903668268;0.435339262304032;0.411951006533964;0.382203680198341;0.350084765549618;0.449336863239025;0.491675156854688;0.528742725475426;0.554529428756521;0.563821946760006;0.554529428756521;0.528742725475426;0.491675156854688;0.449336863239025;0.550663136760975;0.601184336282628;0.641773724956084;0.664679107085624;0.664679107085624;0.641773724956084;0.601184336282628;0.550663136760975;0.649915234450382;0.704819625809520;0.743507500235625;0.757559126971731;0.743507500235625;0.704819625809520;0.649915234450382;0.743029710943569;0.797198267897582;0.827957170600959;0.827957170600959;0.797198267897582;0.743029710943569;0.826194351441247;0.873904113219911;0.891261985882293;0.873904113219911;0.826194351441247;0.896004145930908;0.931765437358897;0.931765437358897;0.896004145930908;0.949600266546736;0.968948006912267;0.949600266546736;0.984784023135109;0.984784023135109;1];
    otherwise
        t=GaussLobattoR(n, 0, 1);
        [S, T]=TrigLobatto(t);
end

if nargout==3
    R=1-S-T;
end

%% Demo
% % Generate nodes
% n=6;
% [S, T]=TrigNodeVect(n);
% 
% % Plot unit triangle
% figure; hold on;
% plot([0,1,0,0], [0,0,1,0]);
% plot(S, T, 'ro');
% axis equal;
% triplotunit(S,T,0*S);




