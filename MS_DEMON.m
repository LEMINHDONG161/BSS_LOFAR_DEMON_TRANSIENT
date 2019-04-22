function [y square_sig] = MS_DEMON(x)

y = [];
for i = 1:size(x,1)
load BPF_500_4500_Fs_10k
f_sig = filter(Num,1,x(i,:)); clear Num % BPF
square_sig = f_sig.^2; % ^2
load LPF_500_200_Fs_10k
Lf_sig = filter(Num,1,square_sig); clear Num %LPF
% Lf_sig = envelope(square_sig); clear Num %LPF
sig = detrend(Lf_sig); % Dc removal
y = [y; sig];
end
