clc; clear; 
close all;

Data_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program\SaveData';
Prgm_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program';

%% Load
cd(Data_path)
disp('Load...')
load('Target_old1'); Target_1 = Target_sig;
load('Target_old2'); Target_2 = Target_sig;

Fs = 10000;
time = 180;
t = 0:1/Fs:time; % Time (sec)

%% Set
cd(Prgm_path)
% 좌표값 (표적 위치)
x = [-5 3 ];
y = [2600 2300];

% 좌표값 (센서 위치)
x_sen = [4 1 7];
y_sen = [0 0 0];

p = 44.5; % Surface noise (0~4 stat: 44.5~66.5)

win = 1024; % window (1.024sec)
overlap = 512; % overlap
nfft = 1024;

%% Mixing
disp('Mixing...')
wind = hanning(1180000)';
wind = wind(590001:1180000);
winn = [ones(1,10001) wind];
[R1_a, R2_a,R3_a,Td] = Mixing2(Target_1, Target_2, Fs, x, y, x_sen, y_sen, p);
R1 = R1_a;
R2 = R2_a;
R3 = R3_a;
R = [R1; R2;R3];

R = R(:,Td:end);

%% DEMON
disp('LOFAR...')

dfs =1000; % Down sampling frequency
sc = Fs/dfs; % Down sampling rate

Ds = R;
% + A(:,Td:end); % DEMON (Envelope [LPF])

Ds_deci_1 = decimate(Ds(1,:),sc); % Down sampling 1
Ds_deci_2 = decimate(Ds(2,:),sc); % Down sampling 2
Ds_deci_3 = decimate(Ds(3,:),sc);
Ds_deci = [Ds_deci_1; Ds_deci_2; Ds_deci_3 ]; % DEMON (Resampling)

tt = 0:1/dfs:(size(Ds_deci,3)-1)/dfs;

[Ds_STFT1,F_DEMON,T_DEMON] = stft(Ds_deci(1,:),win,overlap,nfft,dfs); % STFT 1
[Ds_STFT2,F_DEMON,T_DEMON] = stft(Ds_deci(2,:),win,overlap,nfft,dfs); % STFT 2
[Ds_STFT3,F_DEMON,T_DEMON] = stft(Ds_deci(3,:),win,overlap,nfft,dfs); % STFT 2
%[S_T1,F_T,T_T,P_T] = spectrogram(R1,2^12,512,1024*30,Fs,'yaxis');
%[S_T2,F_T,T_T,P_T] = spectrogram(R2,2^12,512,1024*30,Fs,'yaxis');
Ds_STFT = cell(1,2,3); 

Ds_STFT{1,1} = Ds_STFT1; Ds_STFT{1,2} = Ds_STFT2; Ds_STFT{1,3} = Ds_STFT3;% LOFAR (STFT)

%% JADE (Frequency domain)
disp('FD-JADE...')
FD_sep1_F = []; FD_sep2_F = []; FD_sep3_F = []; FD_sep_F = cell(1,2,3);
for i = 1:size(Ds_STFT{1,1},2)
    ST = [Ds_STFT{1,1}(:,i)'; Ds_STFT{1,2}(:,i)'; Ds_STFT{1,3}(:,i)'];
    [Ae, FD_jade] = jade(ST,3);
    FD_sep1_F = [FD_sep1_F FD_jade(1,:)'];
    FD_sep2_F = [FD_sep2_F FD_jade(2,:)'];
    FD_sep3_F = [FD_sep3_F FD_jade(3,:)'];
end
clear i

FD_sep_F{1,1} = FD_sep1_F; FD_sep_F{1,2} = FD_sep2_F; FD_sep_F{1,3} = FD_sep3_F % STFT (Seperated)

% IFFT
[FD_sep1, T_ifft1] = istft(FD_sep1_F, win, overlap, nfft, dfs);
[FD_sep2, T_ifft2] = istft(FD_sep2_F, win, overlap, nfft, dfs);
[FD_sep3, T_ifft3] = istft(FD_sep3_F, win, overlap, nfft, dfs);

FD_sep = [FD_sep1; FD_sep2; FD_sep3]; % Seperated signals

%% Figure
figure, spectrogram(R1,Fs,'yaxis');
figure, spectrogram(R2,Fs,'yaxis');
 figure, spectrogram(R1,2^12,512,Fs,'yaxis');

 [S_T,F_T,T_T,P_T] = spectrogram(R1,2^12,512,1024*30,Fs,'yaxis');
 
figure,
plot(F_T,P_T(:,2)/max(P_T(:,2)))
 xlabel('Frequency (Hz)','fontsize',12); ylabel('Normalized spectrum','fontsize',12);
 set(gca,'fontsize',12)
 set(gcf,'color','w')
 grid on
 xlim([0 500])
figure, spectrogram(R2,2^12,512,Fs,'yaxis');


 [S_T,F_T,T_T,P_T] = spectrogram(R2,2^12,512,1024*30,Fs,'yaxis');
 
figure,
plot(F_T,P_T(:,2)/max(P_T(:,2)))
 xlabel('Frequency (Hz)','fontsize',12); ylabel('Normalized spectrum','fontsize',12);
 set(gca,'fontsize',12)
 set(gcf,'color','w')
 grid on
 xlim([0 300])
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal');
for ii = 1:size(Ds_STFT,3)
    subplot(3,1,ii)
    surf(T_DEMON,F_DEMON,10*log10(abs(Ds_STFT{1,ii}./max(abs(Ds_STFT{1,ii}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-25 5])
   ylim([0 300])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
end

figure,
set(gcf,'numbertitle','off','name', 'FD-JADE');
for ii = 1:size(FD_sep_F,3)
    subplot(3,1,ii)
    surf(T_DEMON,F_DEMON,10*log10(abs(FD_sep_F{1,ii})./max(abs(FD_sep_F{1,ii})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-25 5])
    ylim([0 300])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
end

figure,
set(gcf,'numbertitle','off','name', 'Recieved signal');
for ii = 1:size(Ds_STFT,3)
    subplot(3,1,ii)
    plot(F_DEMON,db(abs(Ds_STFT{1,ii}(:,2))/max(abs(Ds_STFT{1,ii}(:,2)))))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum (dB)','fontsize',12);
    xlim([0 300])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
end

figure,
set(gcf,'numbertitle','off','name', 'FD-JADE');
for ii = 1:size(FD_sep_F,3)
    subplot(3,1,ii)
    plot(F_DEMON,db(abs(FD_sep_F{1,ii}(:,100))/max(abs(FD_sep_F{1,ii}(:,100)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum (dB)','fontsize',12);
    xlim([0 300])
    ylim([-40 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
end
