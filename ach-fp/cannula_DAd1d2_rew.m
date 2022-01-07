%% create structure with infusion photometry data
for y = 1:3
    [s(y).a, t, s(y).ev] = plot_fp2event(s(y).s, [-6 2], 0);
end
 
%% mu, std of rest ACh distribution during infusions
mu = []; sigma = []; % initialize/clear matrix
fp_cell = {}; fp_max = {};
% fp_1 = {}; fp_2 = {};
 
z = 2; % ACh or DA
for y = 1:length(s) % iterate over infusion
    beh = s(y).s;
    for x = 1:length(beh) % iterate over animal
        ii = find(strcmp({cannula.rec},beh(x).rec)); ii = ii(1); % identify matching baseline recording from cannula larger structure
        fp_base = cannula(ii).FP{z}; % extract full baseline photometry signal
        idx = extractEventST([1:length(fp_base)]', cannula(ii).onRest, cannula(ii).offRest, 1); % rest baseline
        fp_base_mu = nanmean(fp_base(idx)); % mean of baseline photometry recording
        
        a_rew = s(y).a{x,z}; % extract reward-aligned signal
        a_rew = a_rew - fp_base_mu; % subtract baseline%
        ev = s(y).ev{x}; % extract event times
        inf_win = s(y).win(x,:).*60;
        idx_inf = ev > inf_win(1) & ev < inf_win(2); % exclude events outside of infusion window
        a_rew = a_rew(:,idx_inf); % exclude events outside of infusion window
        % a_rew = a_rew - nanmean(a_rew(find(t == -6):find(t == -1),:)); % adjust for own baseline
        m = max(a_rew(find(t == 0):find(t == 0.5),:));
        
        fp_cell{x,y} = a_rew;
        fp_max{x,y} = m(:);
    end
end

%%
fig = figure; fig.Position(3) = 1375;
y = 1; sp(y) = subplot(1,3,y); hold on; axis('square');
plot([0 0],[-2 6],'--k');
for x = 1:4
    plot(t, nanmean(fp_cell{x,y},2), 'Color', [1 0 1]);
end
xlabel('Latency to Reward (s)'); 
ylabel('relative change rDA1m (%dF/F'); xlim([-1 2]);
title(sprintf('rDA1m - %s',s(y).inf));

y = 2; sp(y) = subplot(1,3,y); hold on; axis('square');
plot([0 0],[-2 6],'--k');
for x = 1:4
    plot(t, nanmean(fp_cell{x,y},2), 'Color', [0.5 0.2 0.55]);
end
xlabel('Latency to Reward (s)'); ylabel('relative change rDA1m (%dF/F');
title(sprintf('rDA1m - %s',s(y).inf));
xlim([-1 2]);
 
sp(3) = subplot(1,3,3); hold on; axis('square');
tmp_mu = cellfun(@nanmean, fp_max);
tmp_std = cellfun(@nanstd, fp_max);
plot([1.15; 1.85].*ones(2,4),tmp_mu(:,[1:2])','.:k');
errorbar(nanmean(tmp_mu(:,[1:2]),1),SEM(tmp_mu(:,[1:2]),1),'.k')
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','D1R/D2Rant'});
ylabel('rDA1m peak amplitude (%dF/F)');
title(sprintf('rDA1m - rew peak amplitude'));
linkaxes(sp,'y');
movegui(gcf,'center');

%%
fig = figure; fig.Position(3) = 1375;

sp(1) = subplot(1,3,1); hold on; axis('square');
plot([0 0],[-2 6],'--k');
y = 1; tmp = []; for x = 1:4; tmp(:,x) = nanmean(fp_cell{x,y},2); end
shadederrbar(t, nanmean(tmp,2), SEM(tmp,2), 'm');
y = 2; tmp = []; for x = 1:4; tmp(:,x) = nanmean(fp_cell{x,y},2); end
shadederrbar(t, nanmean(tmp,2), SEM(tmp,2), 'r');
xlabel('Latency to Reward (s)'); ylabel('relative change rDA1m (%dF/F');
title(sprintf('rDA1m - %s',s(y).inf));
xlim([-1 2]);
 
sp(2) = subplot(1,3,2); hold on; axis('square');
tmp_mu = cellfun(@nanmean, fp_max);
tmp_std = cellfun(@nanstd, fp_max);
plot([1.15; 1.85].*ones(2,4),tmp_mu(:,[1:2])','.:k');
errorbar(nanmean(tmp_mu(:,[1:2]),1),SEM(tmp_mu(:,[1:2]),1),'.k')
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'aCSF','D1R/D2Rant'});
ylabel('rDA1m peak magnitude (%dF/F)');
title(sprintf('rDA1m - rew peak magnitude'));
linkaxes(sp,'y');
movegui(gcf,'center');
