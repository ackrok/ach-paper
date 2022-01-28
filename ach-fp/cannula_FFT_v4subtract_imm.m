%% new_ach new_da new_scop s raw

%% Extract IMMOBILITY times
adjFs = new(1).rawFs/50; % Adjust imm onset times to match rawFs
for x = 1:length(new)
    if x >= 1 && x <= 4
        onRest = s(1).s(x).onRest*adjFs;
        offRest = s(1).s(x).offRest*adjFs; % Adjust imm offset times to match rawFs
    elseif x >= 5 && x <= 8
        onRest = s(2).s(x-4).onRest*adjFs;
        offRest = s(2).s(x-4).offRest*adjFs; % Adjust imm offset times to match rawFs
    elseif x >= 9 && x <= 12
        onRest = s(3).s(x-8).onRest*adjFs;
        offRest = s(3).s(x-8).offRest*adjFs; % Adjust imm offset times to match rawFs
    end
    new(x).onRest = onRest;
    new(x).offRest = offRest;
    new(x).L_imm = sum(offRest - onRest);
end

% load('R:\homes\ack466\ACh paper\2_ACh_fluctuations\GZ_mAChRant_wheel2by2.mat')
% new_scop(4).onRest = beh(5).onRest*adjFs; new_scop(4).offRest = beh(5).offRest*adjFs;
% new_scop(5).onRest = beh(6).onRest*adjFs; new_scop(5).offRest = beh(6).offRest*adjFs;
% new_scop(6).onRest = beh(8).onRest*adjFs; new_scop(6).offRest = beh(8).offRest*adjFs;

%% DEMOD extract during IMMOBILITY + during INFUSION window
%new AK190_210909: new(x).onRest(z) < 6000000
for x = 13:15
    tmp = [];
    for z = 1:length(new(x).onRest) % Iterate over each immobility period
        % if new(x).onRest(z) < 6000000; continue; end % Extract only from t = +30min post infusion start time
        % if length(tmp) > 5000000; continue; end
        tmp = [tmp; new(x).demod(new(x).onRest(z):new(x).offRest(z))]; % Extract and concatenate immobility periods
    end
    new(x).demodImmSub = tmp;
end


%% FFT
%new = new_ach([1:8 13 14 11 12]);
p1_mat = [];
for x = 1:length(new)
    vec = new(x).demodImmSub(1:2500000); % keep only 30min, to match scop length
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
fprintf('FFT done! \n');

%% Normalize FFT
p1_mat = p1_all;
tmp = [];
r = [find(f == 0.01):find(f == 100)]; % Restrict to [0.01 100]
% r = [find(f == 0.5):find(f == 100)]; % Restrict to [0.5 50]Hz
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 1):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:size(p1_mat,2)
    vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    % vec_norm = normalize(log10(p1_mat(r,x)),'zscore'); % Normalize z-score
    % vec_norm = normalize(log10(p1_mat(r,x)),'range');
    % vec_norm = normalize(vec_norm,'zscore');
    tmp(:,x) = vec_norm;
end
norm = tmp;
norm = tmp(:,[1:12]);
scop = tmp(:,[13:15]);

%% PLOT FFT w/o subtraction
figure; hold on
plot(flog, nanmean(norm(:,[1:4]),2), 'Color', [0.05 0.75 0.45]);
plot(flog, nanmean(norm(:,[5:8]),2), 'm');
plot(flog, nanmean(norm(:,[9:12]),2), 'Color', [0.85 0.35 0.1]);
plot(flog, nanmean(scop,2), 'k');
legend({'aCSF','d1d2','glu','mAChRant'})
xlabel('Frequency'); % xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});

%% PLOT SUBTRACTION
tmp = []; for x = 1:12; tmp(:,x) = norm(:,x) - nanmean(scop,2); end % Subtract avg FFT for mAChR antagonist
figure; hold on
plot([-2 2],[0 0],'--k');
plot([0 0],[-0.2 0.3],'--k'); plot([0.6021 0.6021],[-0.2 0.3],'--k'); 
plot(flog, nanmean(tmp(:,[1:4]),2), 'Color', [0.05 0.75 0.45]);
plot(flog, nanmean(tmp(:,[5:8]),2), 'Color', 'm');
plot(flog, nanmean(tmp(:,[9:12]),2), 'Color', [0.85 0.35 0.1]);
legend({'line','1 Hz','4 Hz','aCSF','d1d2','glu'})
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting mAChRantag');

%% PLOT SUBTRACTION + SHADEDERRBAR
tmp = [];
for x = 1:12; tmp(:,x) = norm(:,x) - nanmean(scop,2); end
figure; hold on
plot([-2 2],[0 0],'--k');
plot([0 0],[-0.2 0.3],'--k'); plot([0.6021 0.6021],[-0.2 0.3],'--k'); 
shadederrbar(flog, nanmean(tmp(:,[1:4]),2), SEM(tmp(:,[1:4]),2), [0.05 0.75 0.45]); hold on
shadederrbar(flog, nanmean(tmp(:,[5:8]),2), SEM(tmp(:,[5:8]),2), 'm'); hold on
shadederrbar(flog, nanmean(tmp(:,[9:12]),2), SEM(tmp(:,[9:12]),2), [0.85 0.35 0.1]); hold on
legend({'line','aCSF','d1d2','glu','1 Hz','4 Hz'})
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
title('FFT subtracting mAChRantag');

%% AUC post-subtraction
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 1):find(f_sub == 4)]; % AUC from [1 4] Hz
for x = 1:12
    auc(x) = trapz(tmp(r_14,x))/length(r_14);
end
auc = reshape(auc, [4 3]);

figure; violinplot(auc);  xticklabels({'aCSF','d1d2','glu'})
hold on; plot([0.5 3.5],[0 0],'--k');
title('AUC post-subtraction');

%% AUC PLOT
a = auc + 0.15;
figure; hold on
plot([0.5 3.5],[0 0],'--k');
errorbar(nanmean(a)',SEM(a,1)','.k','MarkerSize',20);
plot([1.15; 1.85; 2.85].*ones(3,4),a','.m','MarkerSize',20);
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'aCSF','D1R/D2R','NMDA/AMPA'});