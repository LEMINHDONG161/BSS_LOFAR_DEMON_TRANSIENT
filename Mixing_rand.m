function [R, Td] = Mixing_rand(Target_1, Target_2, Fs, x, y, x_sen, y_sen, p)

% S1: Source 1, S2: Source 2, S3 Source 3
% R1: Sensor 1, R2: Sensor 2 R3 Sensor 3, R4 Sensor 4

dep = 80; % Sensor depth

%% Distance (Target - Sensor)
d_S1_R1 = sqrt(((abs(y(1)-y_sen(1)))^2)+((abs(x(1)-x_sen(1)))^2));
d_S1_R2 = sqrt(((abs(y(1)-y_sen(2)))^2)+((abs(x(1)-x_sen(2)))^2));
d_S1_R3 = sqrt(((abs(y(1)-y_sen(3)))^2)+((abs(x(1)-x_sen(3)))^2));
d_S1_R4 = sqrt(((abs(y(1)-y_sen(4)))^2)+((abs(x(1)-x_sen(4)))^2));
d_S2_R1 = sqrt(((abs(y(2)-y_sen(1)))^2)+((abs(x(2)-x_sen(1)))^2));
d_S2_R2 = sqrt(((abs(y(2)-y_sen(2)))^2)+((abs(x(2)-x_sen(2)))^2));
d_S2_R3 = sqrt(((abs(y(2)-y_sen(3)))^2)+((abs(x(2)-x_sen(3)))^2));
d_S2_R4 = sqrt(((abs(y(2)-y_sen(4)))^2)+((abs(x(2)-x_sen(4)))^2));
%% Transmission Loss
k = [d_S1_R1 d_S1_R2 d_S1_R3 d_S1_R4 d_S2_R1 d_S2_R2 d_S2_R3 d_S2_R4];

S1_R1= Target_1*(1/k(1));
S1_R2= Target_1*(1/k(2));
S1_R3= Target_1*(1/k(3));
S1_R4= Target_1*(1/k(4));
S2_R1= Target_2*(1/k(5));
S2_R2= Target_2*(1/k(6));
S2_R3= Target_2*(1/k(7));
S2_R4= Target_2*(1/k(8));
%% Time delay
c = 1500; % sound speed

TD_S1_R1 = d_S1_R1/c;
TD_S1_R2 = d_S1_R2/c;
TD_S1_R3 = d_S1_R3/c;
TD_S1_R4 = d_S1_R4/c;
TD_S2_R1 = d_S2_R1/c;
TD_S2_R2 = d_S2_R2/c;
TD_S2_R3 = d_S2_R3/c;
TD_S2_R4 = d_S2_R4/c;
% Zero pading
L_S1_R1 = length(0:1/Fs:TD_S1_R1);
L_S1_R2 = length(0:1/Fs:TD_S1_R2);
L_S1_R3 = length(0:1/Fs:TD_S1_R3);
L_S1_R4 = length(0:1/Fs:TD_S1_R4);
L_S2_R1 = length(0:1/Fs:TD_S2_R1);
L_S2_R2 = length(0:1/Fs:TD_S2_R2);
L_S2_R3 = length(0:1/Fs:TD_S2_R3);
L_S2_R4 = length(0:1/Fs:TD_S2_R4);

S1_R1_M = S1_R1;
S1_R2_M = S1_R2;
S1_R3_M = S1_R3;
S1_R4_M = S1_R4;
S2_R1_M = S2_R1;
S2_R2_M = S2_R2;
S2_R3_M = S2_R3;
S2_R4_M = S2_R4;



S1_R1(1,1:L_S1_R1) = 0;
S1_R2(1,1:L_S1_R2) = 0;
S1_R3(1,1:L_S1_R3) = 0;
S1_R4(1,1:L_S1_R4) = 0;
S2_R1(1,1:L_S2_R1) = 0;
S2_R2(1,1:L_S2_R2) = 0;
S2_R3(1,1:L_S2_R3) = 0;
S2_R4(1,1:L_S2_R4) = 0;

Td = max([L_S1_R1 L_S1_R2 L_S1_R3 L_S1_R4 L_S2_R1 L_S2_R2 L_S2_R3 L_S2_R4])

%% Multi path
S1_R1_M1 = 0.7*S1_R1_M; S1_R1_M2 = 0.5*S1_R1_M;
S2_R2_M1 = 0.7*S2_R2_M; S2_R2_M2 = 0.5*S2_R2_M;
S1_R3_M1 = 0.7*S1_R3_M; S1_R3_M2 =  S1_R3_M;
S2_R3_M1 = 0.7*S2_R3_M; S2_R3_M2 = 0.5*S2_R3_M;
S2_R1_M = 0.8*S2_R1_M;
S1_R2_M = 0.8*S1_R2_M;
S1_R3_M = 0.8*S1_R3_M;
S2_R3_M = 0.8*S2_R3_M;

LL1 = 700; LL2 = 900; LL3 = 1400; LL4 = 2000;
TT1 = LL1/c; TT2 = LL2/c; TT3 = LL2/c; TT4 = LL4/c;
LTT1 = length(0:1/Fs:TT1); LTT2 = length(0:1/Fs:TT2); LTT3 = length(0:1/Fs:TT3);
LTT4 = length(0:1/Fs:TT4);

S1_R1_M1(1,1:LTT1) = 0; S1_R3_M2(1,1:LTT4) = 0;
S1_R3_M1(1,1:LTT1) = 0; S1_R1_M2(1,1:LTT2) = 0;
S2_R2_M1(1,1:LTT1) = 0; S2_R2_M2(1,1:LTT2) = 0;
S2_R1_M(1,1:LTT3) = 0; S1_R2_M(1,1:LTT3) = 0; S2_R3_M(1,1:LTT3) = 0;

%% Surface Noise
WN = sqrt(db2mag(p-100))*randn(1,size(Target_1,2));

L = size(WN,2);
f = Fs*(0:(L/2))/L;
No_1k = round(1000/mean(diff(f)));

NL1 = 17*log10(f./1000); f2 = linspace(max(f), min(f), length(f));
NL2 = 17*log10(f2./1000);
NL = [NL1 NL2(2:end)];
NL(1:No_1k+1) = 0;
NL(length(NL)-No_1k:end) = 0;

Y = fft(WN);
y = Y./db2mag(NL);
x = real(ifft(y));

WNs = x;

%% Noise Loss (Sruface - sensor1)
TL_N1 = 20*log10(dep);
TL_N2 = 20*log10(dep+x_sen(2));
TL_N3 = 20*log10(dep+x_sen(3));
TL_N4 = 20*log10(dep+x_sen(4));

TL_N = [db2mag(-TL_N1) db2mag(-TL_N2) db2mag(-TL_N3) db2mag(-TL_N4)];

%% Mix
% R1 = S1_R1+S2_R1+(WNs(1,:)*(1/dep));
% R2 = S1_R2+S2_R2+(WNs(1,:)*(1/(dep-x_sen(2))));
% 
% R1 = S1_R1+S2_R1+S1_R1_M1+(WNs(1,:)*(1/dep));
% R2 = S1_R2+S2_R2+S2_R2_M1+(WNs(1,:)*(1/(dep-x_sen(2))));
% 
% R1 = S1_R1+S2_R1+S1_R1_M1+S1_R1_M2+(WNs(1,:)*(1/dep));
% R2 = S1_R2+S2_R2+S2_R2_M1+S2_R2_M2+(WNs(1,:)*(1/(dep-x_sen(2))));

% size(S1_R1)
% Random channel
H = rand(100,3)
%K= rand(1000)*WNs;
% Assume Target are the same
R = H*[S1_R1; S2_R1; rand(1)*WNs];

% R1 = S1_R1+S2_R1;
% R2 = S1_R2+S2_R2 ;
% R3 = S1_R3+S2_R3 ;
% R4 = S1_R4+S2_R4 ;

% R1 = S1_R1+S2_R1+(WNs(1,:));
% R2 = S1_R2+S2_R2+(WNs(1,:));

