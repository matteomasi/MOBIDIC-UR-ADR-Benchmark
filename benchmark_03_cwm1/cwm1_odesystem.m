function dCdt = cwm1_odesystem(t,C,p)
%% Definition of the system of ODEs
% From Langergraber et al. (2009)
%   N = 17 model processes
%   M = 16 model components
%
% Last update: 28/02/2024

%% STOICHIOMETRIC MATRIX
S = zeros(17,16);
% comp#1 So
S(2,1) = 1 - 1/p(49);
S(4,1) = 1 - 1/p(49);
S(7,1) = -(4.57-p(50))/p(50);
S(15,1) = -(2-p(54))/p(54);
% comp#2 Sf
S(1,2) = 1 - p(46);
S(2,2) = -1/p(49);
S(3,2) = -1/p(49);
S(6,2) = p(47);
S(8,2) = p(47);
S(9,2) = -1/p(51);
S(10,2) = p(47);
S(12,2) = p(47);
S(14,2) = p(47);
S(17,2) = p(47);
% comp#3 Sa
S(4,3) = -1/p(49);
S(5,3) = -1/p(49);
S(9,3) = (1-p(51))/p(51);
S(11,3) = -1/p(52);
S(13,3) = -1/p(53);
% comp#4 Si
S(1,4) = p(46);
% comp#5 Snh
S(1,5) = p(57) - (1-p(46))*p(55) - p(46)*p(56);
S(2,5) = p(55)/p(49) - p(59);
S(3,5) = S(2,5);
S(4,5) = -p(59);
S(5,5) = S(4,5);    
S(6,5) = p(59) - p(47)*p(55) - (1 - p(47) - p(48))*p(57) - p(48)*p(58);
S(7,5) = -p(59) - 1/p(50);
S(8,5) = S(6,5);
S(9,5) = p(55)/p(51) - p(59);
S(10,5) = S(6,5);
S(11,5) = S(4,5);
S(12,5) = S(6,5);
S(13,5) = S(4,5);
S(14,5) = S(6,5);
S(15,5) = S(4,5);
S(16,5) = S(4,5);
S(17,5) = S(6,5);
% comp#6 Sno
S(3,6) = -(1-p(49))/(2.86*p(49));
S(5,6) = S(3,6);
S(7,6) = 1/p(50);
S(16,6) = -(1-p(54))/(0.875*p(54));
% comp#7 Sso4
S(13,7) = -(1-p(53))/(2*p(53));
S(15,7) = 1/p(54);
S(16,7) = 1/p(54);
% comp#8 Sh2s
S(13,8) = (1-p(53))/(2*p(53));
S(15,8) = -1/p(54);
S(16,8) = -1/p(54);
% comp#9 Xs
S(1,9) = -1;
S(6,9) = 1 - p(47) - p(48);
S(8,9) = S(6,9);
S(10,9) = S(6,9);
S(12,9) = S(6,9);
S(14,9) = S(6,9);
S(17,9) = S(6,9);
% comp#10 Xi
S(6,10) = p(48);
S(8,10) = p(48);
S(10,10) = p(48);
S(12,10) = p(48);
S(14,10) = p(48);
S(17,10) = p(48);
% comp#11 Xh
S(2,11) = 1;
S(3,11) = 1;
S(4,11) = 1;
S(5,11) = 1;
S(6,11) = -1;
% comp#12 Xa
S(7,12) = 1;
S(8,12) = -1;
% comp#13 Xfb
S(9,13) = 1;
S(10,13) = -1;
% comp#14 Xamb
S(11,14) = 1;
S(12,14) = -1;
% comp#15 Xasrb
S(13,15) = 1;
S(14,15) = -1;
% comp#16 Xsob
S(15,16) = 1;
S(16,16) = 1;
S(17,16) = -1;


%% PROCESS RATES
P = zeros(1,17);
% process#1 - Hydrolysis
P(1) = p(1)* C(9) / (C(11)+C(13)) / ( p(2)+( C(9)/(C(11)+C(13)) )  ) * ( C(11)+p(3)*C(13)  ); %
% process#2 - Aerobic growth of Xh on Sf
P(2) = p(4) * C(2)/(p(8)+C(2)) * C(2)/(C(2)+C(3)) * C(1)/(p(7)+C(1)) * C(5)/(p(11)+C(5)) * p(12)/(p(12)+C(8)) * C(11); %
% process#3 - Anoxic growth of Xh on Sf
P(3) = p(5)*p(4) * C(2)/(p(8)+C(2)) * C(2)/(C(2)+C(3)) * p(7)/(p(7)+C(1)) * C(6)/(p(10)+C(6)) * C(5)/(p(11)+C(5)) * p(12)/(p(12)+C(8)) * C(11); %
% process#4 - Aerobic growth of Xh on Sa
P(4) = p(4) * C(3)/(p(9)+C(3)) * C(3)/(C(2)+C(3)) * C(1)/(p(7)+C(1)) * C(5)/(p(11)+C(5)) * p(12)/(p(12)+C(8)) * C(11); %
% process#5 - Anoxic growth of Xh on Sa
P(5) = p(5)*p(4) * C(3)/(p(9)+C(3)) * C(3)/(C(2)+C(3)) * p(7)/(p(7)+C(1)) * C(6)/(p(10)+C(6)) * C(5)/(p(11)+C(5)) * p(12)/(p(12)+C(8)) * C(11); %
% process#6 - Lysis of Xh
P(6) = p(6)*C(11); %
% process#7 - Aerobic growth of Xa on Snh
P(7) = p(13) * C(5)/(p(16)+C(5)) * C(1)/(p(15)+C(1)) * p(17)/(p(17)+C(8)) * C(12); %
% process#8 - Lysis of Xa
P(8) = p(14)*C(12); %
% process#9 - Growth of Xfb
P(9) = p(18) * C(2)/(p(21)+C(2)) * p(24)/(p(24)+C(8)) * p(20)/(p(20)+C(1)) * p(22)/(p(22)+C(6)) * C(5)/(p(23)+C(5)) * C(13); %
% process#10 - Lysis of Xfb
P(10) = p(19)*C(13); %
% process#11 - Growth of Xamb
P(11) = p(25) * C(3)/(p(28)+C(3)) * p(31)/(p(31)+C(8)) * p(27)/(p(27)+C(1)) * p(29)/(p(29)+C(6)) * C(5)/(p(30)+C(5)) * C(14); %
% process#12 - Lysis of Xamb
P(12) = p(26)*C(14); %
% process#13 - Growth of Xasrb
P(13) = p(32) * C(3)/(p(35)+C(3)) * C(7)/(p(38)+C(7)) * p(39)/(p(39)+C(8)) * p(34)/(p(34)+C(1)) * p(36)/(p(36)+C(6)) * C(5)/(p(37)+C(5)) * C(15); %
% process#14 - Lysis of Xasrb
P(14) = p(33)*C(15); %
% process#15 - Aerobic growth of Xsob on Sh2s
P(15) = p(40) * C(8)/(p(45)+C(8)) * C(1)/(p(42)+C(1)) * C(5)/(p(44)+C(5)) * C(16); %
% process#16 - Anoxic growth of Xsob on Sh2s
P(16) = p(40)* p(60) * C(8)/(p(45)+C(8)) * C(6)/(p(43)+C(6)) * p(42)/(p(42)+C(1)) * C(5)/(p(44)+C(5)) * C(16); %
% process#17 - Lysis of Xsob
P(17) = p(41)*C(16); %


%% System of ODEs
dCdt = S'*P';


end