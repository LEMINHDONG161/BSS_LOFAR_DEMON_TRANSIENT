function [R, Td] = Mixing_100(Target_1, Target_2, Fs, x, y, x_sen, y_sen, p)

% S1: Source 1, S2: Source 2, S3 Source 3
% R1: Sensor 1, R2: Sensor 2 R3 Sensor 3, R4 Sensor 4

dep = 80; % Sensor depth

%% Distance (Target1 - Sensor)
d_S1_R = zeros(100,1);
for i=1:100
    d_S1_R(i,:)= sqrt(((abs(y(1)-y_sen(i)))^2)+((abs(x(1)-x_sen(i)))^2));
end
clear i
%% Distance (Target2- sensor)
d_S2_R = zeros(100,1);
for j=1:100
    d_S2_R(j,:)= sqrt(((abs(y(2)-y_sen(j)))^2)+((abs(x(2)-x_sen(j)))^2));
end
clear j
%% Transmission Loss

for i =1:100
    k= Target_1*(1/d_S1_R(i,:));
    S1_R(i,:) = k;
end
clear i

for j=1:100
    h=Target_2*(1/d_S2_R(j,:));
    S2_R(j,:)=h;
end
clear j

%% Distance sensor 1

%% Time delay
c = 1500; % sound speed
TD_S1_R = zeros(100,1);
for i =1:100
    TD_S1_R(i,:)=d_S1_R(i,:)/c;
end
clear i
TD_S2_R = zeros(100,1);
for j =1:100
    TD_S2_R(j,:)=d_S2_R(j,:)/c;
end
clear j
% Zero pading
L_S1_R = zeros(100,1);
for i=1:100
    L_S1_R(i,:) = length(0:1/Fs:TD_S1_R(i,:));
end
clear i
L_S2_R = zeros(100,1);
for j=1:100
    L_S2_R(j,:) = length(0:1/Fs:TD_S2_R(j,:));
end
clear j
%for i=1:100
%    S1_R(i,:)(1,1:L_S1_R(i,:)) = 0;
%end
%clear i

%for j=1:100
%    S1_R(j,:)(1,1:L_S1_R(j,:)) = 0;
%end
%clear j
Td1 = max(L_S1_R);
Td2 = max(L_S2_R);
Td = max([Td1 Td2]);



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
H = rand(100,2)
% Assume Target are the same

R = H*[S1_R(1); S2_R(1)];
% R1 = S1_R1+S2_R1;
% R2 = S1_R2+S2_R2 ;
% R3 = S1_R3+S2_R3 ;
% R4 = S1_R4+S2_R4 ;

% R1 = S1_R1+S2_R1+(WNs(1,:));
% R2 = S1_R2+S2_R2+(WNs(1,:));

