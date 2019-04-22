clc; clear; 
close all;

Data_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program\SaveData';
Prgm_path = 'C:\Users\jys\Desktop\BSS and Antenna\Passive sonar\Program_final\Target_BSS_Program';

%% Load
cd(Data_path)
disp('Load...')
load('Target_narrow1'); Target_1 = Target_sig;
load('Target_narrow2'); Target_2 = Target_sig;

Fs = 10000;
time = 180;
t = 0:1/Fs:time; % Time (sec)

%% Set
cd(Prgm_path)
% 좌표값 (표적 위치)
x = [-50 50];
y = [4000 4000];

% 좌표값 (센서 위치)
x_sen = [-5 5];
y_sen = [0 0];

p = 66.5; % Surface noise (0~4 stat: 44.5~66.5)

win = 1024; % window (1.024sec)
overlap = 512; % overlap
nfft = 1024;

%% Mixing
disp('Mixing...')

[R1, R2,Td] = Mixing(Target_1, Target_2, Fs, x, y, x_sen, y_sen, p);
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

% IFFT
[FD_sep1, T_ifft1] = istft(FD_sep1_F, win, overlap, nfft, dfs);
[FD_sep2, T_ifft2] = istft(FD_sep2_F, win, overlap, nfft, dfs);

FD_sep = [FD_sep1; FD_sep2]; % Seperated signals

%% Figure Target 
% Spectrogram time - frequency information
figure,
set(gcf,'numbertitle','off','name', 'Target signal 1');
    surf(T_LOFAR_T,F_LOFAR_T,10*log10(abs(Ds_T_STFT{1,1}./max(abs(Ds_T_STFT{1,1}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-3 0])
   ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
figure,
set(gcf,'numbertitle','off','name', 'Target signal 2');
    surf(T_LOFAR_T,F_LOFAR_T,10*log10(abs(Ds_T_STFT{1,2}./max(abs(Ds_T_STFT{1,2}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-3 0])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
% Frequency information
figure,
    set(gcf,'numbertitle','off','name', 'Target signal 1');

    plot(F_LOFAR_T,db(abs(Ds_T_STFT{1,1}(:,2))/max(abs(Ds_T_STFT{1,1}(:,2)))))
    %     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum dB','fontsize',12);
    xlim([0 500])
    ylim([-50 0 ])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
    set(gcf,'numbertitle','off','name', 'Target signal 1');

    plot(F_LOFAR_T,abs(Ds_T_STFT{1,1}(:,2))/max(abs(Ds_T_STFT{1,1}(:,2))))
    %     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum ','fontsize',12);
    xlim([0 500])
    ylim([0 1])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
    set(gcf,'numbertitle','off','name', 'Target signal 2');

    plot(F_LOFAR_T,db(abs(Ds_T_STFT{1,2}(:,2))/max(abs(Ds_T_STFT{1,2}(:,2)))))
    %     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum dB','fontsize',12);
    xlim([0 500])
    ylim([-50 0 ])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
 figure,
    set(gcf,'numbertitle','off','name', 'Target signal 2');

    plot(F_LOFAR_T,abs(Ds_T_STFT{1,2}(:,2))/max(abs(Ds_T_STFT{1,2}(:,2))))
    %     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum ','fontsize',12);
    xlim([0 500])
    ylim([0 1])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
%% Figure Received Signal
% Spectrogram time - frequency information
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 1');
    surf(T_LOFAR,F_LOFAR,10*log10(abs(Ds_STFT{1,1}./max(abs(Ds_STFT{1,1}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-3 0])
   ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 1');
    surf(T_LOFAR,F_LOFAR,abs(Ds_STFT{1,1}./max(abs(Ds_STFT{1,1})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([0 10])
   ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 2');
    surf(T_LOFAR,F_LOFAR,10*log10(abs(Ds_STFT{1,2}./max(abs(Ds_STFT{1,2}))))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-3 0])
   ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
% Frequency information
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 1');

    plot(F_LOFAR,db(abs(Ds_STFT{1,1}(:,2))/max(abs(Ds_STFT{1,1}(:,2)))))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum (dB)','fontsize',12);
    xlim([0 500])
    ylim([-50 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 1');

    plot(F_LOFAR,abs(Ds_STFT{1,1}(:,2))/max(abs(Ds_STFT{1,1}(:,2))))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum ','fontsize',12);
    xlim([0 500])
    ylim([0 1])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
    figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 1');

    plot(F_LOFAR,abs(Ds_STFT{1,1}(:,2)))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel(' spectrum ','fontsize',12);
    xlim([0 500])
    %ylim([0 10])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
set(gcf,'numbertitle','off','name', 'Recieved signal 2');

    plot(F_LOFAR,db(abs(Ds_STFT{1,2}(:,2))/max(abs(Ds_STFT{1,2}(:,2)))))
%     plot(F_R,db(abs((Ds_STFT{1,ii}(:,2)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum (dB)','fontsize',12);
    xlim([0 500])
    ylim([-50 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
%% Figure FD-JADE     
% Spectrogram time - frequency information
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE 1');
    surf(T_LOFAR,F_LOFAR,10*log10(abs(FD_sep_F{1,1})./max(abs(FD_sep_F{1,1})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-3 0])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE 2');
    surf(T_LOFAR,F_LOFAR,10*log10(abs(FD_sep_F{1,2})./max(abs(FD_sep_F{1,2})))...
        ,'edgecolor','none'); axis tight;
    shading interp;
    view(0,90);
    colorbar;
    caxis([-3 0])
    ylim([0 500])
    xlabel('Time (Sec)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12);
    set(gca,'fontsize',12)
    set(gcf,'color','w')
 % Frequency information   
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE 1');
    plot(F_LOFAR,db(abs(FD_sep_F{1,1}(:,100))/max(abs(FD_sep_F{1,1}(:,100)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum dB ','fontsize',12);
    xlim([0 500])
    ylim([-50 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE 1');
    plot(F_LOFAR,abs(FD_sep_F{1,1}(:,100))/max(abs(FD_sep_F{1,1}(:,100))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum ','fontsize',12);
    xlim([0 500])
    ylim([0 1])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE 2');
    plot(F_LOFAR,db(abs(FD_sep_F{1,2}(:,100))/max(abs(FD_sep_F{1,2}(:,100)))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum dB ','fontsize',12);
    xlim([0 500])
    ylim([-50 0])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on
figure,
set(gcf,'numbertitle','off','name', 'FD-JADE 2');
    plot(F_LOFAR,abs(FD_sep_F{1,2}(:,100))/max(abs(FD_sep_F{1,2}(:,100))))
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Normalized spectrum ','fontsize',12);
    xlim([0 500])
    ylim([0 1])
    set(gca,'fontsize',12)
    set(gcf,'color','w')
    grid on

