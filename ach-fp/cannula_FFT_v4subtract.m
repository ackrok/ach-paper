load('C:\Users\Anya\Desktop\FP_LOCAL\4_glutamate\raw_data\cannula_demod_noFilter.mat');
% Processing: digitalLIA, but commented out all filtering (lines 74, 78, 80)

%% FFT
p1_mat = [];
for x = 1:length(new)
vec = new(x).demod(9000001:end); % keep only 30min, to match scop length
Fs = new(x).rawFs;

T = 1/Fs;               % Sampling period
L = length(vec);        % Length of signal
vec(isnan(vec)) = [];
fftACh = fft(vec);      % Discrete Fourier Transform of photometry signal
P2 = abs(fftACh/L);     % Two-sided spectrum P2
P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;     % Frequency domain vector
P1 = medfilt1(P1);      % Median filter initial FFT

p1_mat(:,x) = movmean(P1,500);
end
p1_all = p1_mat;
fprintf('FFT done!\n');

%% Normalize FFT
p1_mat = p1_all;
tmp = [];
r = [find(f == 0.01):find(f == 100)]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 1):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(p1_mat,2)
    vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100]
    auc_tmp = trapz(vec_norm(r_14));
    tmp(:,x) = vec_norm;
end

%% PLOT FFT w/o subtraction
figure; hold on
plot(flog, nanmean(norm(:,[1:4]),2), 'Color', [0.05 0.75 0.45]);
plot(flog, nanmean(norm(:,[5:8]),2), 'm');
plot(flog, nanmean(norm(:,[9:12]),2), 'Color', [0.85 0.35 0.1]);
plot(flog, nanmean(scop,2), 'k');

legend({'aCSF','d1d2','glu','mAChRant'})
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});

%% PLOT SUBTRACTION
tmp = [];
for x = 1:12; tmp(:,x) = norm(:,x) - nanmean(scop,2); end
figure; hold on
plot([-2 2],[0 0],'--k');
plot(flog, tmp(:,[1:4]), 'Color', [0.05 0.75 0.45]);
plot(flog, tmp(:,[5:8]), 'Color', 'm');
plot(flog, tmp(:,[9:12]), 'Color', [0.85 0.35 0.1]);
legend({'line','aCSF','d1d2','glu'})
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting scopolamine');

%% AUC SUBTRACTION
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 1):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:12
    auc(x) = trapz(tmp(r_14,x))/length(r_14);
end
auc = reshape(auc, [4 3]);
figure; violinplot(auc); xticklabels({'aCSF','d1d2','glu'})