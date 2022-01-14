%% NOTE
% important: make sure to run adjProcessFP to cut 30Hz
% lpCut = 30; filtOrder = 8;
% rawFs = 2000; dsRate = 40; dsType = 2;
% interpType = 'linear'; fitType = 'interp'; basePrc = 5; winSize = 10; winOv = 0;

%% FFT
p1_mat = [];
for x = 1:length(raw)
    vec = []; for y = 1:size(raw(x).fp_prc1,2)
        vec = [vec; raw(x).fp_prc1(:,y)]; end
Fs = raw(x).Fs;

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

%% PLOT all traces on logarithmic axis
ds = 10;
figure;
loglog(f(1:ds:end), p1_mat((1:ds:end),1)); hold on
for x = 2:length(raw)
    loglog(f(1:ds:end), p1_mat((1:ds:end),x));
end
xlim([0.01 100]);
xlabel('Frequency'); ylabel('Power (a.u.)');
title(sprintf('%s FFT (n = %d recs)',raw(1).FPnames{1},size(p1_adj,2)));

%% normalize [0 1]
p1_adj = [];
r = [0.01 100];
for x = 1:size(p1_mat,2)
    a = p1_mat(:,x);
    b = log10(a);
    fft_min = b(f == r(2));
    fft_max = b(1);
    c = (b - fft_min) / (fft_max - fft_min);
    p1_adj(:,x) = c;
end

%% PLOT ADJUSTED FFT
figure; hold on
shadederrbar(log10(f), nanmean(p1_adj,2), SEM(p1_adj,2), 'g'); hold on
% for x = 1:length(new)
%     plot(log10(f), tmp(:,x));
% end
xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
xlabel('Frequency'); ylabel('Power (a.u.)');
title(sprintf('%s adjusted FFT (n = %d recs)',raw(1).FPnames{1},size(p1_adj,2)));