%% Load data
% Data structure should include field 'fp_prc1', otherwise run adjProcessFP
% on raw data first

%% Single sided ACh amplitude spectrum
figure; 
plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
for ii = 1:length(beh)
    sp(ii) = subplot(plm,pln,ii); hold on
    for jj = 1:size(beh(ii).fp_prc1,2)
        y = beh(ii).fp_prc1(:,jj);
        Fs = 50;                % Sampling frequency
        T = 1/Fs;               % Sampling period
        L = length(y);          % Length of signal
        y(isnan(y)) = [];
        fftACh = fft(y);        % Discrete Fourier Transform of photometry signal
        P2 = abs(fftACh/L);     % Two-sided spectrum P2
        P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;     % Frequency domain vector
        P1 = medfilt1(P1);      % Median filter initial FFT
        
        plot(f,P1,'--'); hold on
    end
    title(sprintf('%s: FFT amp spectrum', beh(ii).rec))
    xlabel('f (Hz)'); ylabel('|P1(f)|'); xlim([0 6]); % set(gcf,'color','w')
end

%% 
% beh = raw_ach([1:14, 18:22]); % remove high noise recordings
% beh = modACh;
% beh = raw_ko([3,8:10,12:18]);
% beh = raw_gfp([1:2,4:5,7:12]); 
% beh = raw_ant;

% for x = 1:length(beh); fp_prc1 = [];
%     for y = 1:size(beh(x).fp_prc1,2)
%         fp_prc1 = [fp_prc1; beh(x).fp_prc1(:,y)];
%     end
%     beh(x).fp_prc1 = fp_prc1;
% end

p1_da = [];
h = waitbar(0, 'fft');
for ii = 1:length(beh)
for jj = 1:size(beh(ii).rawfp,2) %jj = 1:size(beh(ii).fp_prc1,2)
    % y = beh(ii).fp_prc1(:,jj);
    y = beh(ii).rawfp(:,jj); Fs = 2000;
    
    % Fs = 50;                % Sampling frequency
    T = 1/Fs;               % Sampling period
    L = length(y);          % Length of signal
    y(isnan(y)) = [];
    fftACh = fft(y);        % Discrete Fourier Transform of photometry signal
    P2 = abs(fftACh/L);     % Two-sided spectrum P2
    P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;     % Frequency domain vector
    P1 = medfilt1(P1);      % Median filter initial FFT

    p1_da = [p1_da, P1];
end
waitbar(ii/length(beh),h);
end; close(h);

%%
ds = 100;

figure;
y = movmean(nanmean(p1_ach,2),50);
loglog(f(1:ds:length(y)-1), y(1:ds:length(y)-1),'--k');

y = movmean(nanmean([p1_ach],2),50); z = movmean(SEM([p1_ach],2),50);
shadederrbar(f(1:ds:length(y)-1), y(1:ds:length(y)-1), z(1:ds:length(y)-1), 'g'); hold on

%y = movmean(nanmean(p1_gfp,2),50); z = movmean(SEM(p1_gfp,2),50);
%shadederrbar(f(1:ds:length(y)-1), y(1:ds:length(y)-1), z(1:ds:length(y)-1), 'k'); hold on

y = movmean(nanmean(p1_ant,2),50); z = movmean(SEM(p1_ant,2),50);
shadederrbar(f(1:ds:length(y)-1), y(1:ds:length(y)-1), z(1:ds:length(y)-1), 'm'); hold on

y = movmean(nanmean(p1_ko,2),50); z = movmean(SEM(p1_ko,2),50);
shadederrbar(f(1:ds:length(y)-1), y(1:ds:length(y)-1), z(1:ds:length(y)-1), 'r'); hold on

% y = movmean(nanmean(p1_mod,2),50); z = movmean(SEM(p1_mod,2),50);
% shadederrbar(f(1:ds:length(y)-1), y(1:ds:length(y)-1), z(1:ds:length(y)-1), 'c'); hold on

% title('Average FFT: ACh (17mice), GFP(3 mice), KO (3mice), ANT (5mice)');
xlim([0.1 30]); xlabel('Frequency (Hz)'); ylabel('Power'); set(gcf,'color','w')
legend({'data1','ACh CW','ACh CW','ACh ANT','ACh ANT','ACh KO','ACh KO'});

%% PLOT ACh - ANT
% figure; ds = 100;
% y = movmean(nanmean(p1_ant,2),50);
% loglog(f(1:ds:45000), y(1:ds:45000),'--k');
% y = movmean(nanmean([p1_mod],2),50); z = movmean(SEM([p1_mod],2),50);
% shadederrbar(f(1:ds:45000), y(1:ds:45000), z(1:ds:45000), 'g'); hold on
% y = movmean(nanmean(p1_ant,2),50); z = movmean(SEM(p1_ant,2),50);
% shadederrbar(f(1:ds:45000), y(1:ds:45000), z(1:ds:45000), 'k'); hold on
% xlim([0.1 8.5]); xlabel('Frequency (Hz)'); ylabel('Power');
% title('Average FFT: ACh vs ANT (5mice)'); legend({'data1','ACh','ACh','ANT','ANT'});

sm = 100; ds = 100;
b = movmean(nanmean(p1_ach,2) - nanmean(p1_ach,2),sm);
figure; plot(f(1:ds:45000), b(1:ds:45000), 'g');
xlim([0.1 10]); ylim([0 0.04]); yticks([0:0.01:0.04]);
xlabel('Frequency (Hz)'); ylabel('Power (ACh - ANT)'); grid on
title('Average FFT: ACh - ANT (n = 5 mice). sm = 100, ds = 100');

%%
figure; 
% shadederrbar(f, nanmean(mat_ach,2), SEM(mat_ach,2), 'b'); hold on; 
loglog(f, nanmean(mat_ko,2), 'b'); xlim([0.1 8.5]);
% ylim([0 0.15]);
xlabel('f (Hz)'); ylabel('|P1(f)|'); % xlim([0 6]); set(gcf,'color','w')
title('Average FFT: ACh');

%%
range = [181:14401];
mdl = fitlm(f(range), (nanmean(mat_gfp(range,:),2)))