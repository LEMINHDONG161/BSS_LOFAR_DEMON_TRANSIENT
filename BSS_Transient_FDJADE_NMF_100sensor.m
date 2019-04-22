clc; clear all; 
close all;


Data_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program\SaveData';
Prgm_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program';

%% Load
cd(Data_path)
disp('Load...')
load('Target_interrup1'); Target_1 = Target_sig;
load('Target_Transient140'); Target_2 = Target_sig;
%load('Target_2_3'); Target_3 = Target_sig;

Fs = 10e3;
time = 180;
t = 0:1/Fs:time; % Time (sec)

%% Set
cd(Prgm_path)
% 좌표값 (표적 위치)
x = [-50 50 ];
% x2 = [-50 50];
% y = [300 350];
y = [500 500 ];
% y2 = [500 500]; 

% 좌표값 (센서 위치)
x_sen = 10*rand(100,1);
y_sen = 5*rand(100,1);

p = 0; % Surface noise (0~4 stat: 44.5~66.5)

win = 1024; % window (1.024sec)
overlap = win*0.5; % overlap
nfft = 1024;

%% Mixing
% disp('Mixing...')
% [R1, R2,R3,R4,Td] = Mixing3(Target_1, Target_2,  Fs, x, y, x_sen, y_sen, p);
% Tar = [Target_1; Target_2];
% R = [R1; R2; R3; R4];


%% Random Mixing
disp('Random Mixing...')
[R,Td] = Mixing_100(Target_1, Target_2,  Fs, x, y, x_sen, y_sen, p);
Tar = [Target_1; Target_2];
R = R(:,Td:end);

%figure,plot(R(1,1e4:2e4))
%figure,plot(R(2,1e4:2e4))
%figure,plot(R(1,:)-R(2,:))


%% Variance of sensor = Signal Power
% According to these you can add Noise
%Sig_power = var(R')

%% LOFAR
disp('LOFAR...')

dfs =1000; % Down sampling frequency
sc = Fs/dfs; % Down sampling rate
Ds_T = Tar;


Ds_T_deci_1 = decimate(Ds_T(1,:),sc); % Down sampling 1
Ds_T_deci_2 = decimate(Ds_T(2,:),sc); % Down sampling 2
Ds_T_deci = [Ds_T_deci_1; Ds_T_deci_2]; % lofar (Resampling)

tt = 0:1/dfs:(size(Ds_T_deci,2)-1)/dfs;

[Ds_T_STFT1,F_LOFAR_T,T_LOFAR_T] = stft(Ds_T_deci(1,:),win,overlap,nfft,dfs); % STFT 1
[Ds_T_STFT2,F_LOFAR_T,T_LOFAR_T] = stft(Ds_T_deci(2,:),win,overlap,nfft,dfs); % STFT 2
Ds_T_STFT = cell(1,2); 

Ds_T_STFT{1,1} = Ds_T_STFT1; Ds_T_STFT{1,2} = Ds_T_STFT2; % LOFAR (STFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
Ds = R;
Ds_dec_array =cell(1,100);
Ds_STFT_array = cell(1,100);
for i=1:100  
M = decimate(Ds(i,:),sc); % Down sampling 1
Ds_dec_array{1,i}=M;

[K,F_LOFAR,T_LOFAR] = stft(Ds_dec(i,:),win,overlap,nfft,dfs); % STFT 
Ds_STFT_array{1,i} =K;
end
clear i



%% JADE (Frequency domain)
disp('FD-JADE...')
FD_sep1_F = []; FD_sep2_F = [];FD_sep3_F = [];  FD_sep_F = cell(1,2,3);
for i = 1:size(Ds_STFT_array,2)
    ST = [Ds_STFT(:,2)'; Ds_STFT(:,3)';Ds_STFT(:,5)'];
    [Ae, FD_jade] = jade(ST,3);
    FD_sep1_F = [FD_sep1_F FD_jade(1,:)'];
    FD_sep2_F = [FD_sep2_F FD_jade(2,:)'];
    FD_sep3_F = [FD_sep3_F FD_jade(3,:)'];
   
end
clear i

FD_sep_F{1,1} = FD_sep1_F; FD_sep_F{1,2} = FD_sep2_F;FD_sep_F{1,3} = FD_sep3_F; % STFT (Seperated)

% IFFT
[FD_sep1, T_ifft1] = istft(FD_sep1_F, win, overlap, nfft, dfs);
[FD_sep2, T_ifft2] = istft(FD_sep2_F, win, overlap, nfft, dfs);
[FD_sep3, T_ifft3] = istft(FD_sep3_F, win, overlap, nfft, dfs);

FD_sep = [FD_sep1; FD_sep2; FD_sep3]; % Seperated signals

%% NMF
% NMF with output FDJADE
N_FD_sep1_F = abs(FD_sep1_F);
N_FD_sep2_F = abs(FD_sep2_F);
N_FD_sep3_F = abs(FD_sep3_F);

% NMF with received signal
%N_FD_sep1_F = abs(Ds_STFT1);
%N_FD_sep2_F = abs(Ds_STFT2);
% NMF with target signal
%N_FD_sep1_F = abs(Ds_T_STFT1);
%N_FD_sep2_F = abs(Ds_T_STFT2);

[W1, H1] = nmf(N_FD_sep2_F, 2, 100);

[W2, H2] = nmf(N_FD_sep2_F, 2, 100);

phi1 = angle(FD_sep1_F); phi2 = angle(FD_sep2_F);

XmagHat1 = cell(1,2); XmagHat2 = cell(1,2);
for j = 1:2
XmagH1 = W1(:,j)*H1(j,:);
XmagH2 = W2(:,j)*H2(j,:);

XmagHat1{1,j} = XmagH1;
XmagHat2{1,j} = XmagH2;
end    

%% Supervised NMF
A = (W2(:,2));
B = (W2(:,1));

h_W1 = (B.*(1./A));
h_W2 = (A.*(1./B));

[Wa1, Ha1]   = nmf_s(N_FD_sep2_F, 2, [h_W1 randn(size(A))], 100, 1);
[Wa2, Ha2]   = nmf_s(N_FD_sep2_F, 2, [h_W2 randn(size(B))], 100, 1);
    Wa2(:,1) = log((1./Wa2(:,1)));
    XmagHat1 = Wa2(:,1)*(Ha1(1,:));
for i=1
   %Wa1(:,i) = ((1./Wa1(:,i)));
  
 
    
  
    XmagHat2 = Wa1(:,i)*(Ha1(i,:));
    %XmagHat1 = Wa2(:,i)*(Ha2(i,:));
    
    
end
%XmagHat1 = Wa2(:,1)*(Ha2(1,:));
%XmagHat2 = Wa2(:,1)*(Ha2(1,:));
%% NMF + FD-JADE
 %HH = XmagHat./max(max(XmagHat));

 HH2 = XmagHat1;
 AA1 = abs(FD_sep_F{1,2});
 AA2 = abs(FD_sep_F{1,2});

 KKK = HH2.*AA2;
 %KKK2 = KKK.*AA2;
 
 %NK = KKK2./max(max(KKK2));
 


%% Figure
%% Figure Target Signal
figure,
set(gcf,'numbertitle','off','name', 'Target signal');
for ii = 1:size(Ds_T_STFT,2)
    subplot(2,1,ii)
      surf(T_LOFAR_T,F_LOFAR_T,10*log10(abs(Ds_T_STFT{1,ii}./max(abs(Ds_T_STFT{1,ii}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-40 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
end
figure,
set(gcf,'numbertitle','off','name', ' Frequency Target signal');
for ii = 1:size(Ds_T_STFT,2)
    subplot(2,1,ii)
    plot(F_LOFAR_T,db(abs(Ds_T_STFT{1,ii}(:,2))/max(abs(Ds_T_STFT{1,ii}(:,2)))))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',6);
    ylabel('Normalized (dB)','fontsize',6);
    xlim([0 500])
    ylim([-40 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
end

%% Figure Received Signal
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal');

    surf(T_LOFAR,F_LOFAR,10*log10(abs(Ds_STFT{1,1}./max(abs(Ds_STFT{1,1}))))... % we can choose which sensor you want to show repplace 1 by n
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-40 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',6); ylabel('Frequency (Hz)','fontsize',6);
    set(gca,'fontsize',12)
    set(gcf,'color','w')

figure,
set(gcf,'numbertitle','off','name', ' Frequency Recieved signal');

    plot(F_LOFAR,db(abs(Ds_STFT{1,1}(:,2))/max(abs(Ds_STFT{1,1}(:,2)))))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',6);
    ylabel('Normalized(dB)','fontsize',6);
    xlim([0 500])
    ylim([-40 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on

%% Figure FD - JADE
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE');
for ii = 1:size(FD_sep_F,3)
    subplot(3,1,ii)
    surf(T_LOFAR,F_LOFAR,10*log10(abs(FD_sep_F{1,ii})./max(abs(FD_sep_F{1,ii})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-40 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',6); ylabel('Frequency (Hz)','fontsize',6);
    set(gca,'fontsize',6)
    set(gcf,'color','w')
end
figure,
set(gcf,'numbertitle','off','name', ' frequency FD-JADE');
for ii = 1:size(FD_sep_F,3)
    subplot(3,1,ii)
    plot(F_LOFAR,db(abs(FD_sep_F{1,ii}(:,100))/max(abs(FD_sep_F{1,ii}(:,100)))))
    xlabel('Frequency (Hz)','fontsize',6);
    ylabel('Normalized(dB)','fontsize',6);
    xlim([0 500])
    ylim([-40 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
end

figure,
set(gcf,'numbertitle','off','name', 'FD-JADE_transient');

    surf(T_LOFAR,F_LOFAR,10*log10(abs(FD_sep_F{1,1})./max(abs(FD_sep_F{1,1})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-40 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')


%% Figure  NNMF
%% Display  NNMF
figure,
surf(T_LOFAR,F_LOFAR,db(XmagHat1./max(max(XmagHat1))),'edgecolor','none'); axis tight;
%surf(T_LOFAR_T,F_LOFAR_T,db(XmagHat./max(max(XmagHat))),'edgecolor','none'); axis tight; % display for target input nnmf
shading interp;
view(0,90);
colorbar;
caxis([-140 0])
xlim([0 180])
ylim([0 500])
xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
set(gca,'fontsize',12)
set(gcf,'color','w')
%% Display  NNMF
figure,
surf(T_LOFAR,F_LOFAR,db(XmagHat2./max(max(XmagHat2))),'edgecolor','none'); axis tight;
%surf(T_LOFAR_T,F_LOFAR_T,db(XmagHat./max(max(XmagHat))),'edgecolor','none'); axis tight; % display for target input nnmf
shading interp;
view(0,90);
colorbar;
caxis([-140 0])
xlim([0 180])
ylim([0 500])
xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
set(gca,'fontsize',12)
set(gcf,'color','w')
%% Display combine NNMF + FDJADE
figure,
surf(T_LOFAR,F_LOFAR,db(KKK./max(max(KKK))),'edgecolor','none'); axis tight;
%surf(T_LOFAR_T,F_LOFAR_T,db(XmagHat./max(max(XmagHat))),'edgecolor','none'); axis tight; % display for target input nnmf
shading interp;
view(0,90);
colorbar;
caxis([-140 0])
xlim([0 180])
ylim([0 500])
xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
set(gca,'fontsize',12)
set(gcf,'color','w')





figure,
surf(T_LOFAR,F_LOFAR,db(XmagHat2./max(max(XmagHat2))),'edgecolor','none'); axis tight;
%surf(T_LOFAR_T,F_LOFAR_T,db(XmagHat./max(max(XmagHat))),'edgecolor','none'); axis tight; % display for target input nnmf
shading interp;
view(0,90);
colorbar;
caxis([-40 0])
xlim([0 180])
ylim([0 500])
xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
set(gca,'fontsize',12)
set(gcf,'color','w')
figure,
plot(T_LOFAR,db(XmagHat2./max(max(XmagHat2))));
%plot(T_LOFAR_T,db(XmagHat./max(max(XmagHat)))); % display for target input nnmf
xlabel('Time XmagHat2 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
ylim([ -60 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR,db(XmagHat2./max(max(XmagHat2))));
%plot(F_LOFAR_T,db(XmagHat./max(max(XmagHat)))); % display for target input nnmf
xlabel('Frequency XmagHat2 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 500])
ylim([-60 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on


%% Test component H
figure,
plot(T_LOFAR,db(Ha1./max(max(Ha1))));
%plot(T_LOFAR,Ha1);
%plot(T_LOFAR_T,Ha1);% display for target input nnmf
xlabel('Time  Ha1','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
ylim([-60 30])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(T_LOFAR,db(Ha1(1,:)./max(max(Ha1(1,:)))));
%plot(T_LOFAR,(Ha1(1,:)));
%plot(T_LOFAR_T,Ha1);% display for target input nnmf
xlabel('Time  ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 180])
ylim([-60 30])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
%plot(T_LOFAR,(Ha1./max(max(Ha1))),'r-');
plot(T_LOFAR,(Ha1(2,:)));
%plot(T_LOFAR_T,Ha1);% display for target input nnmf
xlabel('Time Ha1_2 ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 180])
%ylim([-60 30])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on

figure,
%plot(T_LOFAR,(Ha2./max(max(Ha2))));
plot(T_LOFAR,Ha2);
%plot(T_LOFAR_T,Ha2);% display for target input nnmf
xlabel('Time Ha2 ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 180])
%ylim([-1 1])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,

plot(T_LOFAR,db(Ha2(1,:)./max(max(Ha2(1,:)))));
%plot(T_LOFAR,(Ha1(1,:)));
%plot(T_LOFAR_T,Ha1);% display for target input nnmf
xlabel('Time  ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 180])
ylim([-60 30])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
%plot(T_LOFAR,(Ha2./max(max(Ha2))));
plot(T_LOFAR,(Ha2(2,:)));
%plot(T_LOFAR_T,Ha2);% display for target input nnmf
xlabel('Time Ha2_2  ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 180])
%ylim([-60 30])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on


figure,
%plot(T_LOFAR,(H2./max(max(H2))));
subplot(2,1,1);
plot(T_LOFAR,H2(1,:));
%plot(T_LOFAR_T,(H2./max(max(H2))));% display for target input nnmf
xlabel('Time H2 ','fontsize',12);
ylabel('Amplitude  ','fontsize',12);
xlim([0 180])
%ylim([-1 1])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
%figure,
%plot(T_LOFAR,(H2./max(max(H2))));
subplot(2,1,2);
plot(T_LOFAR,H2(2,:));
%plot(T_LOFAR_T,(H2./max(max(H2))));% display for target input nnmf
xlabel('Time H2 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
%ylim([-1 1])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
%% Test component W
figure,
plot(F_LOFAR,(Wa1./max(max(Wa1))));
%plot(F_LOFAR,Wa1);
%plot(F_LOFAR_T,Wa1);% display for target input nnmf
xlabel('Frequency Wa1 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 500])
ylim([0 1])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,

plot(F_LOFAR,1-(Wa1(:,1)./max(max(Wa1(:,1)))));
%plot(F_LOFAR,Wa2);
%plot(F_LOFAR_T,Wa2);% display for target input nnmf
xlabel('Frequency  ','fontsize',12);
ylabel(' Amplitude  ','fontsize',12);
xlim([0 500])
%ylim([-800 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR,(Wa1(:,2)));
%plot(F_LOFAR,Wa2);
%plot(F_LOFAR_T,Wa2);% display for target input nnmf
xlabel('Frequency Wa1_2 ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 500])
%ylim([-60 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR,(Wa2./max(max(Wa2))));
%plot(F_LOFAR,Wa2);
%plot(F_LOFAR_T,Wa2);% display for target input nnmf
xlabel('Frequency Wa2 ','fontsize',12);
ylabel('Amplitude db','fontsize',12);
xlim([0 500])
%ylim([-60 10])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR,(Wa2(:,1)));
%plot(F_LOFAR,Wa2);
%plot(F_LOFAR_T,Wa2);% display for target input nnmf
xlabel('Frequency Wa2_1 ','fontsize',12);
ylabel('Amplitude db ','fontsize',12);
xlim([0 500])
%ylim([-1100 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR,(Wa2(:,2)));
%plot(F_LOFAR,Wa2);
%plot(F_LOFAR_T,Wa2);% display for target input nnmf
xlabel('Frequency Wa2_2','fontsize',12);
ylabel('Amplitude  db ','fontsize',12);
xlim([0 500])
%ylim([-60 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
%plot(F_LOFAR,(W2./max(max(W2))),'r');
plot(F_LOFAR,W2(:,1),'r-');
%plot(F_LOFAR_T,W2,'b-');% display for target input nnmf
xlabel('Frequency W2 ','fontsize',12);
%ylabel('Amplitude ','fontsize',12);
xlim([0 500])
%ylim([0 0.3])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
hold on
plot(F_LOFAR,W2(:,2),'b-');
%plot(F_LOFAR_T,W2,'b-');% display for target input nnmf
xlabel('Frequency W2 ','fontsize',12);
%ylabel('Amplitude ','fontsize',12);
xlim([0 500])
%ylim([0 0.3])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
%hold on
%figure,


%% 


 

