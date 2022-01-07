%%
h = waitbar(0, 'pauses');
for y = 1:length(s) % iterate over infusion
    %%
    beh = s(y).s;
    [align_rew, t_rew, ev_rew] = plot_fp2event(beh,[-6 2],0);
    [align_acc, t_acc, ev_acc] = plot_fp2event(beh,[-6 2],0);

    a_rew_avg = []; a_acc_avg = [];
    for x = 1:length(beh) % iterate over animal
        fp = beh(x).FP{1}; Fs = 50; % extract full photometry signal
        idx_inf = [s(y).win(x,:)*60]; % infusion window
        
        ii = find(ev_rew{x} > idx_inf(1) & ev_rew{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
        a_rew = align_rew{x,1}(:,ii); % extract aligned traces that are in infusion window
        a_rew = a_rew - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        a_rew = a_rew - nanmean(a_rew([find(t_rew == -6):find(t_rew == -1)],:),1);
        
        ii = find(ev_acc{x} > idx_inf(1) & ev_acc{x} < idx_inf(end)); % index, wrt ev vector, where events are in infusion window
        a_acc = align_acc{x,1}(:,ii); % extract aligned traces that are in infusion window
        a_acc = a_acc - s(y).base_mu_rest(x); % subtract mean of baseline recording to center around 0%
        a_acc = a_acc - nanmean(a_acc([find(t_acc == -6):find(t_acc == -1)],:),1);
        
        s(y).a_rew{x} = a_rew;
        s(y).a_acc{x} = a_acc;
        a_rew_avg(:,x) = nanmean(a_rew,2); 
        a_acc_avg(:,x) = nanmean(a_acc,2);
    end
    s(y).a_rew_avg = a_rew_avg;
    s(y).a_acc_avg = a_acc_avg;
end

%% PLOT AVERAGE
fig = figure; fig.Position(3) = 1375;
clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1];
subplot(1,3,1); hold on
plot([0 0],[-4 4], '--k');
for y = 1:length(s)
    shadederrbar(t_rew, nanmean(s(y).a_rew_avg,2), SEM(s(y).a_rew_avg,2), clr(y,:)); hold on
end
xlabel('Latency to Reward (s)'); xlim([-1 2]); xticks([-1:2]);
ylabel('ACh3.0 (%dF/F)'); ylim([-4 4]); yticks([-4:2:4]);

title(sprintf('ACh to reward'));
%
subplot(1,3,2); hold on
plot([0 0],[-2 7], '--k');
for y = 1:length(s)
    shadederrbar(t_rew, nanmean(s(y).a_acc_avg,2), SEM(s(y).a_acc_avg,2), clr(y,:)); hold on
end
xlabel('Latency to Acc Peak (s)'); xlim([-1.5 1.5]); xticks([-1:1]);
ylabel('ACh3.0 (%dF/F)'); ylim([-2 7]); yticks([-2:2:6])
title(sprintf('ACh to acceleration'));

movegui(gcf,'center');

%% 
amp_rew = []; lat_rew = []; amp_acc = [];
win_rew = [find(t_rew == 0):find(t_rew == 1.5)];
win_acc = [find(t_rew == 0):find(t_rew == 0.5)];
for y = 1:3; for x = 1:4
    [mm, ii] = min(s(y).a_rew{x}(win_rew,:));
    ii = ii + win_rew(1) - 1;
    amp_rew(x,y) = nanmean(mm);
    lat_rew(x,y) = nanmean(t_rew(ii));
    [mm, ii] = max(s(y).a_acc{x}(win_acc,:));
    amp_acc(x,y) = nanmean(mm);
    end; end

fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
violinplot(-1.*amp_rew); xticklabels({s.inf});
ylabel('reward pause amplitude'); ylim([0 5]); yticks([0:2:5]);
title(sprintf('(a/d = %1.3f) (a/n = %1.3f)',signrank(amp_rew(:,1),amp_rew(:,2)),signrank(amp_rew(:,1),amp_rew(:,3))))
subplot(1,3,2); hold on
violinplot(lat_rew); xticklabels({s.inf});
ylabel('reward pause latency'); ylim([0 1]); yticks([0:0.5:1]);
title(sprintf('(a/d = %1.3f) (a/n = %1.3f)',signrank(lat_rew(:,1),lat_rew(:,2)),signrank(lat_rew(:,1),lat_rew(:,3))))
subplot(1,3,3); hold on
violinplot(amp_acc); xticklabels({s.inf});
ylabel('acceleration amplitude'); ylim([0 10]); yticks([0:2:10]);
title(sprintf('(a/d = %1.3f) (a/n = %1.3f)',signrank(amp_acc(:,1),amp_acc(:,2)),signrank(amp_acc(:,1),amp_acc(:,3))))