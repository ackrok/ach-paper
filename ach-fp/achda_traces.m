%% ACh leads DA
x = 27;

fp_ach = beh(x).FP{1}; Fs = beh(x).Fs; % extract photometry signal from structure
fp_mu = nanmean(fp_ach); % mean of entire photometry signal
fp_ach = fp_ach - fp_mu; % SUBTRACT baseline% from fp%, now centered on BASELINE
idx_rew = extractEventST([1:length(fp_ach)]', floor(beh(x).reward)-100, floor(beh(x).reward)+100, 1); % identify sample during reward
idx_imm = extractEventST([1:length(fp_ach)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
fp_ach_imm = fp_ach(idx_imm_nonRew); 
fp_da = beh(x).FP{2}; % dopamine
fp_da = fp_da - nanmean(fp_da);
fp_da_imm = fp_da(idx_imm_nonRew);

figure; hold on; 
plot(beh(x).time(idx_imm_nonRew), fp_ach_imm, 'g'); 
plot(beh(x).time(idx_imm_nonRew), fp_da_imm, 'm');
ylabel('Photometry (%dF/F)'); ylim([-20 20])
yyaxis right; 
acc = getAcc(beh(x).vel); 
plot(beh(x).time(idx_imm_nonRew), acc(idx_imm_nonRew), 'k');
ylabel('Acceleration (cm/s^2)'); ylim([-2 6]);
title(sprintf('%s',beh(x).rec))

% x = 41; xlim([676 681.5]); xticks([676:681]); % JM007 ACh leads DA
% x = 30; xlim([1581.5 1587.5]); xticks([1582:1587]); % AK190 ACh leads DA