% 绘制三种类型的B样条曲线，需要前面所给的所有.m文件
clear all;
%控制顶点
% % 火焰
% P{1}=[91.5,414.2;122,244;232,86;358.5,9];
% P{2}=[358.5,9;260,172;405,277;357.4,479.5];
% P{3}=[357.4,479.5;398,412;428,333;421.1,201.3];
% P{4}=[421.1,201.3;482,258;509,324;529.3,406.4];
% P{5}=[529.3,406.4;543,386;554,366;567.1,358.3];
% P{6}=[567.1,358.3;525,502;728,515;459.9,751.1];
% P{7}=[459.9,751.1;544,571.5;475,485;443.9,456.4];
% P{8}=[443.9,456.4;425,538.5;393,614.4;341.6,666.6];
% P{9}=[341.6,666.6;318,502;269,555;256.2,292.2];
% P{10}=[256.2,292.2;234,349;199.5,444;185,493.4];
% P{11}=[185,493.4;176,481;163,462;155.1,460.7];
% P{12}=[155.1,460.7;195.5,561.5;58.5,618;159,746.2];
% P{13}=[159,746.2;-121,530;96,486;46,358.2];
% P{14}=[46,358.2;69.5,380;73,384;91.5,414.2];

% % 鸟
% P{1}=[302.9,302.2;278,167.5;432,106;512,190];
% P{2}=[512,190;547,183;570,172;589,161];
% P{3}=[589,161;580,191;564,211;536,228];
% P{4}=[536,228;565,224;579,221;606,209];
% P{5}=[606,209;594,228;572,254;545,273];
% P{6}=[545,273;546,573;252,729;9,582];
% P{7}=[9,582;72,588;129,578;190,531];
% P{8}=[190,531;132,529;90,487;75,446];
% P{9}=[75,446;92,450;114,450;131,444];
% P{10}=[131,444;56,429;31,362;33,321];
% P{11}=[33,321;45,330;32,339;87,338];
% P{12}=[87,338;20,291;29,210;50,175];
% P{13}=[50,175;95,231;176,296;302.9,302.2];
% 
% % 叶子
% P{1}=[311,86;332,174;373,212;429,179];
% P{2}=[429,179;406,219;401,262;404,323];
% P{3}=[404,323;404,323;470,264;470,264];
% P{4}=[470,264;481,308;524,311;561,301];
% P{5}=[561,301;531,372;569,402;597,415];
% P{6}=[597,415;550,427;445,479;462,594];
% P{7}=[462,594;418,541;367,535;316,568];
% P{8}=[316,568;316,568;316,681;316,681];
% P{9}=[316,681;316,681;302,681;302,681];
% P{10}=[302,681;302,681;302,565;302,565];
% P{11}=[302,565;273,542;205,552;153,581];
% P{12}=[153,581;181,495;92,432;14,408];
% P{13}=[14,408;61,402;82,379;68,300];
% P{14}=[68,300;118,318;148,308;156,240];
% P{15}=[156,240;163,276;162,287;217,324];
% P{16}=[217,324;232,283;207,208;196,187];
% P{17}=[196,187;246,199;294,180;311,86];



% 
% % 狮子
% P{1}=[192,49;192,49;190,44;190,44];
% P{2}=[190,44;190,44;197,36;197,36];
% P{3}=[197,36;197,36;194,34;194,34];
% P{4}=[194,34;194,34;219,21;219,21];
% P{5}=[219,21;219,21;215,19;215,19];
% P{6}=[215,19;215,19;239,12;239,12];
% P{7}=[239,12;239,12;236,10;236,10];
% P{8}=[236,10;236,10;260,6;260,6];
% P{9}=[260,6;260,6;256,4;256,4];
% P{10}=[256,4;275,2;292,12;295,18];
% P{11}=[295,18;294,24;298,25;298,28];
% P{12}=[298,28;311,37;313,40;311,44];
% P{13}=[311,44;303,61;302,62;287,55];
% P{14}=[287,55;285,57;282,64;280,75];
% P{15}=[280,75;280,75;277,73;277,73];
% P{16}=[277,73;277,73;273,97;273,97];
% P{17}=[273,97;273,97;270,95;270,95];
% P{18}=[270,95;265,111;265,114;257,114];
% P{19}=[257,114;257,155;259,192;271,184];
% P{20}=[271,184;275,188;276,194;275,195];
% P{21}=[275,195;275,195;251,195;251,195];
% P{22}=[251,195;241,166;240,162;236,143];
% P{23}=[236,143;226,182;223,187;233,185];
% P{24}=[233,185;236,188;239,192;237,195];
% P{25}=[237,195;237,195;213,195;213,195];
% P{26}=[213,195;213,195;209,123;209,123];
% P{27}=[209,123;163,138;139,115;134,117];
% P{28}=[134,117;119,167;122,167;129,185];
% P{29}=[129,185;139,184;140,186;141,195];
% P{30}=[141,195;141,195;117,195;117,195];
% P{31}=[117,195;98,155;106,148;103,137];
% P{32}=[103,137;88,157;87,184;88,185];
% P{33}=[88,185;98,184;99,185;100,195];
% P{34}=[100,195;100,195;76,195;76,195];
% P{35}=[76,195;73,179;72,158;79,129];
% P{36}=[79,129;56,182;5,159;12,135];
% P{37}=[12,135;25,139;31,158;51,146];
% P{38}=[51,146;67,134;69,113;70,97];
% P{39}=[70,97;79,30;135,52;192,49];

% %凤凰
% P{1}=[317,219;212,169;209,408;431,367];
% P{2}=[431,367;397,366;391,363;344,336];
% P{3}=[344,336;385,329;418,341;440,356];
% P{4}=[440,356;426,340;406,327;390,301];
% P{5}=[390,301;417,313;438,338;459,369];
% P{6}=[459,369;426,325;459,285;502,278];
% P{7}=[502,278;486,282;449,319;472,347];
% P{8}=[472,347;469,340;463,309;499,297];
% P{9}=[499,297;491,301;466,329;483,354];
% P{10}=[483,354;482,338;484,328;489,320];
% P{11}=[489,320;485,365;504,364;512,364];
% P{12}=[512,364;544,363;544,338;569,345];
% P{13}=[569,345;569,341;568,337;544,336];
% P{14}=[544,336;554,332;573 332;579 345];
% P{15}=[579 345;577 324;531 331;531 287];
% P{16}=[531 287;538 318;556 316;580 332];
% P{17}=[580 332;594 344;581 357;590 369];
% P{18}=[590 369;538,346;574 455;477 409];
% P{19}=[477 409;450 395;423 427;462 444];
% P{20}=[462 444;445 444;425 432;433 402];
% P{21}=[433 402;424 408;403 464;487 462];
% P{22}=[487 462;449 480;383 447;424 391];
% P{23}=[424 391;351 463;491 538;518 449];
% P{24}=[518 449;515 533;396 509;388 449];
% P{25}=[388 449;348 520;321 530;196 470];
% P{26}=[196 470;149 443;132 452;123 460];
% P{27}=[123 460;116 469;85 478;104 455];
% P{28}=[104 455;63 463;107 488;81 493];
% P{29}=[81 493;49 491;69 430;139 443];
% P{30}=[139 443;87 414;36 426;28 480];
% P{31}=[28 480;32 373;131 420;220 457];
% P{32}=[220 457;204 448;183 426;181 400];
% P{33}=[181 400;178 368;214 367;214 380];
% P{34}=[214 380;210 402;183 385;196 416];
% P{35}=[196 416;211 387;217 428;214 443];
% P{36}=[214 443;232 468;319 512;350 457];
% P{37}=[350 457;317 488;230 443;230 426];
% P{38}=[230 426;301 482;364 433;405 381];
% P{39}=[405 381;399 383;384 384;379 391];
% P{40}=[379 391;362 416;344 400;331 426];
% P{41}=[331 426;329 405;341 399;349 391];
% P{42}=[349 391;309 391;335 432;276 419];
% P{43}=[276 419;312 416;282 398;315 387];
% P{44}=[315 387;281 380;267 421;242 412];
% P{45}=[242 412;263 399;235 368;272 370];
% P{46}=[272 370;259 356;242 361;236 361];
% P{47}=[236 361;226 361;199 348;212 322];
% P{48}=[212 322;216 338;230 343;233 336];
% P{49}=[233 336;239 326;200 278;150 289];
% P{50}=[150 289;90 306;100 348;116 352];
% P{51}=[116 352;127 317;154 361;123 361];
% P{52}=[123 361;108 360;82 328;112 298];
% P{53}=[112 298;101 299;95 285;87 298];
% P{54}=[87 298;87 265;127 286;134 284];
% P{55}=[134 284;116 264;137 270;129 247];
% P{56}=[129 247;156 249;131 268;167 279];
% P{57}=[167 279;168 276;159 257;164 251];
% P{58}=[164 251;175 256;179 277;193 285];
% P{59}=[193 285;191 285;187 273;190 258];
% P{60}=[190 258;205 267;221 287;219 298];
% P{61}=[219 298;232 302;241 318;245 339];
% P{62}=[245 339;293 391;354 384;373 380];
% P{63}=[373 380;167 367;225 152;317 219];
% 



%一段
P{1}=[358.5,9;260,172;405,277;357.4,479.5];


% 
% %两段
% P{1}=[91.5,414.2;122,244;232,86;358.5,9];
% P{2}=[358.5,9;260,172;405,277;91.5,414.2];

p = 3;
num=100;

for i=1:length(P)
cp{i,1}=redraw(p,P{i});%根据控制顶点重新绘制曲线
end


curve=cell2mat(cp);
curve(:,2)=max(curve(:,2))-curve(:,2);
for i=1:length(P)
    cp{i}(:,2)=max(curve(:,2))-cp{i}(:,2);
end
% jianju(1)=norm(curve(end,:)-curve(1,:));
% l(1)=jianju(1);
% for i=2:length(curve)
%     jianju(i)=norm(curve(i,:)-curve(i-1,:));
%     l(i)=l(i-1)+jianju(i);
% end
jianju(length(curve))=norm(curve(end,:)-curve(1,:));

for i=2:length(curve)
    jianju(i-1)=norm(curve(i,:)-curve(i-1,:));    
end
l(1)=jianju(1);
for i=2:length(curve)
     l(i)=l(i-1)+jianju(i);  
end

for i=1:length(cp)
%     j=1+(i-1)*200:i*200;
    lcp(1,i)=jianju(1+(i-1)*5000);
    for j=2:5000
        lcp(j,i)=lcp(j-1,i)+jianju(j+(i-1)*5000);
    end
%     j=1+(i-1)*200:i*200;
%     lcp(i)=sum(jianju(j));
    ratecp(i)=lcp(end,i)/l(end);
    gpcp{i,1}=cp{i}(1,:);
end
ncptemp=ratecp*num;
ncp(1)=round(ncptemp(1));
for i=2:length(ncptemp)-1
    ncp(i)=round(ncp(i-1)+ncptemp(i))-ncp(i-1);
end
ncp(end+1)=num-sum(ncp);

cpinterval=lcp(end,:)./ncp;
for i=1:size(lcp,2)
    counter=1;
    for j=1:size(lcp,1)-1
        if floor(lcp(j,i)/(counter*cpinterval(i)))
            counter=counter+1;
%             if j==size(lcp,1)
%             else
            gpcp{i}(counter,:)=cp{i}(j+1,:);    
%             end
        end
    end
end

gpoint=cell2mat(gpcp);
fenduandian=cumsum(ncp);
fenduandian=1+[0 fenduandian(1:end-1)]';
save yiduan100.mat gpoint num fenduandian
% 
% for i=1:length(gpoint)
%     pause(0.1)
%     plot(gpoint(i,1),gpoint(i,2),'.')
%     hold on
% end









% gpinterval=l(end)/num;
% for i=1:size(lcp,2)
%     counter=1;
%     for j=1:size(lcp,1)
%         if floor(lcp(j,i)/(counter*gpinterval))
%             counter=counter+1;
%             gpcp{i}(counter,:)=cp{i}(j,:);
% %             if lcp(end,i)-lcp(j,i)<gpinterval
% %                 gpcp{i}(counter,:)=[];
% %             end
%         end
%     end
% end
% if norm(gpcp{1}(1,:)-gpcp{end}(end,:))<gpinterval
%     gpcp{end}(end,:)=[];
% end
% for i=2:length(gpcp)
%     if norm(gpcp{i}(1,:)-gpcp{i-1}(end,:))<gpinterval
%         gpcp{i-1}(end,:)=[];
%     end
% end
    


















% for i=1:length(curve)
% pause(0.02)
% plot(curve(i,1),curve(i,2),'.');
% hold on
% end

function p_u=redraw(p,P)

n = length(P)-1;  

ui=jiedianxiangliang(n,p,0,1); % 准均匀B样条的节点矢量

Njp_u = zeros(1, n+1);
j=0;
for u = 0 : 0.0002 : 1-0.0002
    j=j+1;
%     p_u=zeros(length(u),2);
    for i = 1 : n+1
        Njp_u(j, i) = Njp(i, p , u, ui);
    end
    p_u(j,:) = Njp_u(j,:)*P;       
    
end
end