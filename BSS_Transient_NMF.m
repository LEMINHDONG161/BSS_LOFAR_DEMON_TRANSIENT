clc; clear; 
% close all;

Data_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program\SaveData';
Prgm_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program';

%% Load
cd(Data_path)
disp('Load...')
load('Target_nnf9'); Target_1 = Target_sig;
load('Target_nnf9_2'); Target_2 = Target_sig;
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
x_sen = [-5 5];
y_sen = [0 0];

p = 0; % Surface noise (0~4 stat: 44.5~66.5)

win = 1024; % window (1.024sec)
overlap = win*0.5; % overlap
nfft = 1024;

%% Mixing
disp('Mixing...')
[R1, R2,Td] = Mixing(Target_1, Target_2,  Fs, x, y, x_sen, y_sen, p);
Tar = [Target_1; Target_2];
R = [R1; R2];

R = R(:,Td:end);


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = R;


Ds_deci_1 = decimate(Ds(1,:),sc); % Down sampling 1
Ds_deci_2 = decimate(Ds(2,:),sc); % Down sampling 2
Ds_deci = [Ds_deci_1; Ds_deci_2]; %  (Resampling)

tt = 0:1/dfs:(size(Ds_deci,2)-1)/dfs;

[Ds_STFT1,F_LOFAR,T_LOFAR] = stft(Ds_deci(1,:),win,overlap,nfft,dfs); % STFT 1
[Ds_STFT2,F_LOFAR,T_LOFAR] = stft(Ds_deci(2,:),win,overlap,nfft,dfs); % STFT 2

Ds_STFT = cell(1,2); 

Ds_STFT{1,1} = Ds_STFT1; Ds_STFT{1,2} = Ds_STFT2; % LOFAR (STFT)


%% JADE (Frequency domain)
disp('FD-JADE...')
FD_sep1_F = []; FD_sep2_F = []; FD_sep_F = cell(1,2);
for i = 1:size(Ds_STFT{1,1},2)
    ST = [Ds_STFT{1,1}(:,i)'; Ds_STFT{1,2}(:,i)'];
    [Ae, FD_jade] = jade(ST,2);
    FD_sep1_F = [FD_sep1_F FD_jade(1,:)'];
    FD_sep2_F = [FD_sep2_F FD_jade(2,:)'];
end
clear i

FD_sep_F{1,1} = FD_sep1_F; FD_sep_F{1,2} = FD_sep2_F; % STFT (Seperated)
FD_sep1 = []; FD_sep2 = []; FD_sep = cell(1,2);
% IFFT
[FD_sep1, T_ifft] = istft(FD_sep1_F, win, overlap, nfft, dfs);
[FD_sep2, T_ifft] = istft(FD_sep2_F, win, overlap, nfft, dfs);


FD_sep{1,1} = FD_sep1; FD_sep{1,2} = FD_sep2; % Seperated signals

%% NMF
N_FD_sep1_F = abs(Ds_T_STFT1);
N_FD_sep2_F = abs(Ds_T_STFT2);

[W1, H1] = nmf(N_FD_sep1_F, 2, 100);

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

[Wa1, Ha1]   = nmf_s(N_FD_sep1_F, 2, [h_W1 randn(size(A))], 100, 1);
[Wa2, Ha2]   = nmf_s(N_FD_sep2_F, 2, [h_W2 randn(size(B))], 100, 1);

for i=1
    XmagHat = Wa1(:,i)*(Ha1(i,:));
end

%% NMF + FD-JADE
% HH = XmagHat./max(max(XmagHat));
% 
% HH2 = XmagHat;
% AA1 = abs(FD_sep_F{1,1});
% AA2 = abs(FD_sep_F{1,2});
% 
% KKK = HH2.*AA1;
% KKK2 = KKK.*AA2;
% 
% NK = KKK2./max(max(KKK2));
% 
% figure,
% surf(T_DEMON,F_DEMON,db(NK),'edgecolor','none'); axis tight;
% shading interp;
% view(0,90);
% colorbar;
%     caxis([-60 0])
% ylim([0 300])
% xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
% set(gca,'fontsize',12)
% set(gcf,'color','w')

%% Figure
%% Figure
figure,
set(gcf,'numbertitle','off','name', 'Target signal');
for ii = 1:size(Ds_T_STFT,2)
    subplot(2,1,ii)
    surf(T_LOFAR_T,F_LOFAR_T,10*log10(abs(Ds_T_STFT{1,ii}./max(abs(Ds_T_STFT{1,ii}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-20 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
end
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal');
for ii = 1:size(Ds_STFT,2)
    subplot(2,1,ii)
    surf(T_LOFAR,F_LOFAR,10*log10(abs(Ds_STFT{1,ii}./max(abs(Ds_STFT{1,ii}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-20 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
end

figure,
set(gcf,'numbertitle','off','name', 'FD-JADE');
for ii = 1:size(FD_sep_F,2)
    subplot(2,1,ii)
    surf(T_LOFAR,F_LOFAR,10*log10(abs(FD_sep_F{1,ii})./max(abs(FD_sep_F{1,ii})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-20 0])
    xlim([0 180])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
end

% 


figure,
surf(T_LOFAR_T,F_LOFAR_T,db(XmagHat./max(max(XmagHat))),'edgecolor','none'); axis tight;
shading interp;
view(0,90);
colorbar;
caxis([-20 0])
xlim([0 180])
ylim([0 500])
xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
set(gca,'fontsize',12)
set(gcf,'color','w')
figure,
subplot(2,1,1);
%plot(T_LOFAR_T,(H1./max(max(H1))));
plot(T_LOFAR_T,H1);
xlabel('Time H1','fontsize',12);
ylabel('Amplitude','fontsize',12);
xlim([0 180])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
subplot(2,1,2);
%plot(T_LOFAR_T,(H2./max(max(H2))));
plot(T_LOFAR_T,H2);
xlabel('Time H2 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(T_LOFAR_T,(XmagH1./max(max(XmagH1))));
xlabel('Time XmagH1 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(T_LOFAR_T,(XmagH2./max(max(XmagH2))));
xlabel('Time XmagH2 ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(T_LOFAR_T,(XmagHat./max(max(XmagHat))));
xlabel('Time XmagHat ','fontsize',12);
ylabel('Amplitude ','fontsize',12);
xlim([0 180])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
subplot(2,1,1);
%plot(F_LOFAR_T,(W1./max(max(W1))));
plot(F_LOFAR_T,W1,'r-');
xlabel(' Basic Vector Frequency W1 ','fontsize',12);
ylabel('Amplitude','fontsize',12);
%xlim([0 500])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on


%plot(F_LOFAR_T,(W2./max(max(W2))));
subplot(2,1,2);
plot(F_LOFAR_T,W2,'b-');
xlabel('Basic Vector Frequency W2 ','fontsize',12);
ylabel('Amplitude','fontsize',12);
%xlim([0 500])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR_T,W2,'b-');
hold
figure,
plot(F_LOFAR_T,(XmagH1./max(max(XmagH1))));
xlabel('Frequency XmagH1 ','fontsize',12);
ylabel('Amplitude','fontsize',12);
%xlim([0 500])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR_T,(XmagH2./max(max(XmagH2))));
xlabel('Frequency XmagH2 ','fontsize',12);
ylabel('Amplitude','fontsize',12);
%xlim([0 500])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
figure,
plot(F_LOFAR_T,(XmagHat./max(max(XmagHat))));
xlabel('Frequency XmagHat ','fontsize',12);
ylabel('Amplitude','fontsize',12);
%xlim([0 500])
%ylim([-10 0])
set(gca,'fontsize',12)
set(gcf,'color','w')
grid on
 %plot(T_DEMON,(H1./max(max(H1))));
 %xlabel('Frequency (Hz)','fontsize',12);
    %ylabel('Normalized spectrum (dB)','fontsize',12);
    %xlim([0 500])
    %ylim([-10 0])
    %set(gca,'fontsize',12)
    %set(gcf,'color','w')
    %grid on

