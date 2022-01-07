%Cross-correlation between ACh and DA photometry signals acquired during
%dual recordings
%
%   Description: This script is for running cross-correlation analysis on
%   two continuous photometry signals (ACh, DA) after isolating time
%   periods when animal is immobile (no movement, no reward)
%
%   Author: Anya Krok, December 2021

%% 
% good = [1:length(modAChDA)]; rmv = [3 5 18 22:24 27:28 30 33:34 36 41:42]; good(rmv) = []; % acceleration
% good = [10:16,20,22:29,32:35,40:41,44:47]; % reward
good = [1:length(modAChDA)]; rmv = 42; good(rmv) = [];
beh = modAChDA(good);

%%
corr_mat = nan(1001,length(beh));
corr_5 = nan(1001,length(beh)); corr_50 = corr_5; corr_95 = corr_5; 

h = waitbar(0, 'cross corr');
for x = 1:length(beh); y = [1 2]; %CHANGE - which FP signal to run CCG or ACG on
    
    %% acetylcholine
    fp_1 = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mu = nanmean(fp_1); % mean of entire photometry signal
    fp_1 = fp_1 - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_1)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_1)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_imm = extractEventST([1:length(fp_1)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    fp_1_imm = fp_1(idx_imm_nonRew); 
    
    fp_2 = beh(x).FP{y(2)}; % dopamine
    fp_2 = fp_2 - nanmean(fp_2);
    fp_2_imm = fp_2(idx_imm_nonRew);
    %%
    if isempty(fp_1_imm) || isempty(fp_2_imm); continue; end
    %% cross-correlation, whole signal
    % [xcf, lags, bounds] = crosscorr(fp_ach_imm, fp_da_imm,'NumLags',100,'NumSTD',3);
    % corr_achda(:,x) = xcf; corr_achda_bounds(:,x) = bounds; % STORE
    % [corr_achda_shuff(:,x),~,~] = crosscorr(fp_ach_imm(randperm(length(fp_ach_imm))), fp_da_imm(randperm(length(fp_ach_imm))),'NumLags',100,'NumSTD',3);
    [corr_mat(:,x),lags] = xcorr(fp_1_imm, fp_2_imm, 10*Fs, 'coeff');
    tmp_shuff = [];
    for s = 1:50
        tmp_shuff(:,s) = xcorr(fp_1_imm, circshift(fp_2_imm, 30*Fs), 10*Fs, 'coeff');
        % tmp_shuff(:,s) = xcorr(fp_1_imm(randperm(length(fp_1_imm))), fp_2_imm(randperm(length(fp_2_imm))), 60*Fs, 'coeff');
    end
    corr_5(:,x) = prctile(tmp_shuff, [5], 2);
    corr_50(:,x) = prctile(tmp_shuff, [50], 2);
    corr_95(:,x) = prctile(tmp_shuff, [95], 2);
    
%%
    waitbar(x/length(beh),h);
end
close(h);

%% PLOT: average full cross-correlation between ACh, DA signals
corr_adj = corr_mat - nanmean(corr_mat([1:50],:));
% corr_adj = corr_mat;

figure; hold on
plot([0 0],[-0.7 0.2],'--k');
shadederrbar(lags/Fs, nanmean(corr_shuff,2), SEM(corr_shuff,2), 'k'); hold on
plot(lags/Fs, corr_adj, 'Color', [0 0 0 0.1]);
shadederrbar(lags/Fs, nanmean(corr_adj,2), SEM(corr_adj,2), 'b');
xlabel('Lag (s)'); ylabel('Sample Cross Correlation');
title(sprintf('ACh/DA cross-correlation (n = %d recs)',size(corr_mat,2)));
% Interpretation:
% peak corr has POSITIVE lag, so DA signal is equal to ACh signal shifted
% by 10 samples (200ms) to the left


%% N = X mice
% corr_adj = corr_mat - nanmean(corr_mat([1:50],:));
corr_adj = corr_mat;

tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);
corr_an = []; corr_an_shuff = [];
corr_an_5 = []; corr_an_50 = []; corr_an_95 = [];
for x = 1:nAn
    ii = find(strcmp(tmp,uni{x}));
    corr_an(:,x) = nanmean([corr_adj(:,ii)],2);
    % corr_an_shuff(:,x) = nanmean([corr_shuff(:,ii)],2);
    % corr_an_shuff(:,x) = corr_shuff(:,ii(randi(length(ii),1)));
    corr_an_5(:,x) = nanmean([corr_5(:,ii)],2); % corr_5(:,ii(randi(length(ii),1))); %
    corr_an_50(:,x) = nanmean([corr_50(:,ii)],2); % corr_50(:,ii(randi(length(ii),1)));
    corr_an_95(:,x) = nanmean([corr_95(:,ii)],2); % corr_95(:,ii(randi(length(ii),1)));
end
[min_val, ii] = max(corr_an);
min_time = lags(ii)/50;

%% PLOT SHADEDERRBAR
figure; hold on
plot([0 0],[-0.6 0.6],'--k');
% shadederrbar(lags/Fs, nanmean(corr_an_shuff,2), SEM(corr_an_shuff,2), 'k'); hold on
shadederrbar(lags/Fs, nanmean(corr_an_50,2), nanmean(corr_an_5,2), 'k'); hold on
plot(lags/Fs, corr_an, 'Color', [0 0 0 0.1]);
% shadederrbar(lags/Fs, nanmean(corr_an,2), SEM(corr_an,2), 'b');
xlabel('Lag (s)'); xlim([-1 1]);
ylabel('Sample Cross Correlation'); ylim([-0.2 1.2]); yticks([-0.2:0.2:1.2]);
title(sprintf('ACh/DA cross-correlation (n = %d mice)',nAn));
% Interpretation:
% peak corr has POSITIVE lag, so DA signal is equal to ACh signal shifted
% by 10 samples (200ms) to the left

%%
figure; 
subplot(1,2,1); violinplot(min_val); 
ylabel('max corr coeff'); % ylim([-1 0]); yticks([-1:0.2:0]);
title(sprintf('ACh/DA coeff: %1.2f',nanmean(min_val)));
subplot(1,2,2); violinplot(min_time*1000); 
ylabel('max corr lag (ms)'); % ylim([0 250]); yticks([0:50:250]);
title(sprintf('ACh/DA lag: %1.1f ms',nanmean(min_time*1000)));

%% PLOT proportion of CCG > 95% CI
above95 = []; below5 = [];
for x = 1:length(uni)
    a = corr_an(:,x);
    above95 = [above95, a > corr_an_95(:,x)]; %binary vector where CCG passed 95% confidence interval
    below5 = [below5, a < corr_an_5(:,x)]; %binary vector where CCG below 5% confidence interval
end

figure; hold on
bar(lags/Fs, 100*sum(above95,2)/size(above95,2),'FaceColor','b','FaceAlpha',0.5);
bar(lags/Fs, -100*sum(below5,2)/size(below5,2),'FaceColor','b','FaceAlpha',0.5);
xlabel('Lag (s)'); ylabel('Prop of CCG > 95% CI (%)'); ylim([-100 100]);
title(sprintf('Proportion of CCG > 95p CI (n = %d mice)',size(above95,2)))