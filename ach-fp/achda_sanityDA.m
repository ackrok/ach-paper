%% Plot DA to reward, only JM004 signal
beh = d2ant([3,4]); % extract only JM004 signal
[align_1, t] = plot_fp2event(beh,[-6 2]);
[align_2, ~] = plot_fp2event(behLes(7),[-6 2]);
align_full = {};
align_full{1} = align_1{1,2} - nanmean(beh(1).FP{2}); 
align_full{2} = align_1{2,2} - nanmean(beh(2).FP{2}); 
align_full{3} = align_2{1,2} - nanmean(behLes(7).FP{2}); 

%%
fig = figure; fig.Position([3 4]) = [1000 650];

x = 1; sp1(x) = subplot(2,2,x); hold on
plot([0 0],[0 15],'--k');
shadederrbar(t, nanmean(align_full{x},2), SEM(align_full{x},2), 'm');
xlabel('Latency to Reward (s)'); xlim([-0.5 2]); xticks([0:2]);
ylabel('rDA1m (%dF/F)');
axis('square')
title(sprintf('%s rDA1m',strtok(beh(x).rec,'_')));
 
subplot(2,2,2); hold on
plot([0 0],[-5 30],'--k');
plot(t, align_full{x}, 'Color', [1 0 1 0.5]);
xlabel('Latency to Reward (s)'); xlim([-0.5 2]); xticks([0:2]);
axis('square')
 
x = 2; sp1(x) = subplot(2,2,x+1); 
plot([0 0],[0 15],'--k');
shadederrbar(t, nanmean(align_full{x}(:,[1:100]),2), SEM(align_full{x}(:,[1:100]),2), 'r');
xlabel('Latency to Reward (s)'); xlim([-0.5 2]); xticks([0:2]);
axis('square')
title('rDA1m + D2R antagonist');
 
x = 3; sp1(x) = subplot(2,2,x+1); 
plot([0 0],[0 15],'--k');
shadederrbar(t, nanmean(align_full{x},2), SEM(align_full{x},2), 'k');
xlabel('Latency to Reward (s)'); xlim([-0.5 2]); xticks([0:2]);
axis('square')
title('rDA1m + 6ohda lesion');
 
linkaxes(sp1,'x'); linkaxes(sp1,'y')
movegui(gcf,'center');