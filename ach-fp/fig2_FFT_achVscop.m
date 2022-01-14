load('fig2_ach+scop_demod_noFilter');
% Processing: digitalLIA, but commented out all filtering (lines 74, 78, 80)

%% FFT
p1_mat = [];
for x = 1:length(new)
vec = new(x).demod;
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

%% PLOT all traces on logarithmic axis
ds = 10;
figure;
loglog(f(1:ds:end), p1_mat((1:ds:end),1)); hold on
for x = 2:length(new)
    loglog(f(1:ds:end), p1_mat((1:ds:end),x));
end
xlim([0.01 100]);
xlabel('Frequency'); ylabel('Power (a.u.)');
title('ACh3.0 control + scopolamine (n = 3)');

%% normalize [0 1]
tmp = [];
r = [0.01 100];
for x = 1:length(new)
    a = p1_mat(:,x);
    b = log10(a);
    fft_min = b(f == r(2));
    fft_max = b(1);
    c = (b - fft_min) / (fft_max - fft_min);
    tmp(:,x) = c;
end

%% PLOT ADJUSTED FFT
figure; hold on
shadederrbar(log10(f), nanmean(tmp(:,[1:3]),2), SEM(tmp(:,[1:3]),2), 'g'); hold on
shadederrbar(log10(f), nanmean(tmp(:,[4:6]),2), SEM(tmp(:,[4:6]),2), 'k');
% for x = 1:length(new)
%     plot(log10(f), tmp(:,x));
% end
xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
xlabel('Frequency'); ylabel('Power (a.u.)');
title('ACh3.0 control-scopolamine (n = 3)');

%% PLOT SUBTRACTED FFT
figure; hold on
ds = 10;
tmp = []; flog = log10(f);
for x = 1:3
    tmp(:,x) = normalize(log10(p1_mat(:,x)),'range') - normalize(log10(p1_mat(:,x+3)),'range');
    plot(flog(1:ds:end), tmp(1:ds:end,x));
end
shadederrbar(flog(1:ds:end), nanmean(tmp(1:ds:end,:),2), SEM(tmp(1:ds:end,:),2), 'k');
xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
xlabel('Frequency'); ylabel('Power (a.u.)');
title('ACh3.0 control-scopolamine (n = 3)');
