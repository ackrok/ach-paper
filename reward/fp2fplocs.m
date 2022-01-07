[align, t] = plot_fp2fp;

%%
figure;
plm = floor(sqrt(size(align,1))); pln = ceil(size(align,1)/plm); % Subplot size depending on number of recordings
clr = {'g','r','b'}; 
lbl = 'Reward Delivery'; lbl = 'Movement Onset'; lbl = 'Acceleration Peak'; lbl = 'ACh peak'; lbl = 'ACh trough';
for x = 1:length(align)
    sp(x) = subplot(plm,pln,x); y = 1;
    %x = 13; y = 1;
    %plot(time, align{x,y}, ':'); hold on % Plot average STA for each photometry signal
    %shadederrbar(t, nanmean(align{x,y},2), SEM(align{x,y},2), clr{y}); hold on % Plot average STA for each photometry signal
    shadederrbar(t, nanmean(align_mov{x,y},2), SEM(align_mov{x,y},2), clr{1}); hold on % Plot average STA for each photometry signal
    shadederrbar(t, nanmean(align_rest{x,y},2), SEM(align_rest{x,y},2), clr{2}); hold on % Plot average STA for each photometry signal
    xlabel(sprintf('Latency to %s (s)',lbl)); 
    ylabel('DA (z-score)'); grid on; xlim([t(1) t(end)]);
    title(sprintf('%s - %s',beh(x).rec,beh(x).site)); 
    xlim([-1 1])
end; linkaxes(sp,'y'); 

%%
y = 1; % CHANGE: 1 = align to ACh peak, 2 = align to DA peak

full_mov = []; full_rest = []; full_rew = []; % full_lick = [];
for x = 1:length(align_mov)
    full_mov = [full_mov, nanmean(align_mov{x,y},2)];
    full_rest = [full_rest, nanmean(align_rest{x,y},2)];
end
for x = 1:length(align_rew)
    full_rew = [full_rew, nanmean(align_rew{x,y},2)];
    % full_lick = [full_lick, align_lick{x,y}];
end
%
figure; hold on
clr = {'g','r','b','m'};
shadederrbar(t, nanmean(full_mov,2), SEM(full_mov,2), clr{1}); hold on % Plot average STA for each photometry signal
shadederrbar(t, nanmean(full_rest,2), SEM(full_rest,2), clr{2}); hold on % Plot average STA for each photometry signal
% shadederrbar(time, nanmean(full_lick,2), SEM(full_lick,2), clr{4}); hold on % Plot average STA for each photometry signal
% yyaxis right
shadederrbar(t, nanmean(full_rew,2), SEM(full_rew,2), clr{3}); hold on % Plot average STA for each photometry signal
switch y
    case 2; xlabel('Latency to DA peak (s)'); ylabel('ACh (% dF/F)'); title('ACh fp to DA peak (n = 5 mice) - DMS');
    case 1; xlabel('Latency to ACh peak (s)'); ylabel('DA (% dF/F)'); title('DA fp to ACh peak (n = 5 mice) - DMS');
end
xlim([-1 1]);
    
%%
align_GZ3 = cell(1,3); y = 2;
align_GZ3{1} = [align_mov{7,y}, align_mov{8,y}, align_mov{9,y}];
align_GZ3{2} = [align_rest{7,y}, align_rest{8,y}, align_rest{9,y}];
align_GZ3{3} = [align_rew{2,y}, align_rew{5,y}, align_rew{8,y}, align_rew{11,y}, align_rew{14,y}];

% align_GZ2 = cell(1,3); y = 2;
% align_GZ2{1} = [align_mov{4,y}, align_mov{5,y}, align_mov{6,y}];
% align_GZ2{2} = [align_rest{4,y}, align_rest{5,y}, align_rest{6,y}];
% align_GZ2{3} = [align_rew{1,y}, align_rew{4,y}, align_rew{7,y}, align_rew{10,y}, align_rew{13,y}];

figure; clr = {'g','r','b'}; 
for x = 1:3
    shadederrbar(time, nanmean(align_GZ3{x},2), SEM(align_GZ3{x},2), clr{x}); hold on % Plot average STA for each photometry signal
end
% ylabel('DA (z-score)');
% yyaxis right; x = 3; shadederrbar(time, nanmean(align_GZ2{x},2), SEM(align_GZ2{x},2), clr{x}); 
lbl = 'ACh peak'; xlabel(sprintf('Latency to %s (s)',lbl)); 
ylabel('DA (z-score)'); grid on; xlim([time(1) time(end)]);
title('GZ003 - DLS DA to ACh peaks');
xlim([-1 1])

%%
[max_mov, i_mov] = max(full_mov([51:76],:),[],1); i_mov = t(i_mov + 50);
[max_rest, i_rest] = max(full_rest([51:76],:),[],1); i_rest = t(i_rest + 50);
[max_rew, i_rew] = max(full_rew([51:76],:),[],1); i_rew = t(i_rew + 50);

[min_mov, i2_mov] = min(full_mov([1:51],:),[],1); i2_mov = t(i2_mov);
[min_rest, i2_rest] = min(full_rest([1:51],:),[],1); i2_rest = t(i2_rest);
[min_rew, i2_rew] = min(full_rew([1:51],:),[],1); i2_rew = t(i2_rew);
%%
figure; hold on
histogram(i2_mov,'BinWidth',0.02,'Normalization','Probability','FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.1,'DisplayName','MOV');
histogram(i2_rest,'BinWidth',0.02,'Normalization','Probability','FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.1,'DisplayName','REST');
histogram(i2_rew,'BinWidth',0.02,'Normalization','Probability','FaceColor','b','FaceAlpha',0.2,'EdgeAlpha',0.1,'DisplayName','REW');
switch y
    case 2; xlabel('Latency to DA peak (s)'); ylabel('probability'); title('ACh fp to DA peak (n = 8 mice)');
    case 1; xlabel('Latency to ACh peak (s)'); ylabel('probability'); title('DA fp to ACh peak (n = 8 mice)');
end
%%
figure; hold on
violinplot({i_rest, i_mov, i_rew, i2_rest, i2_mov, i2_rew});
