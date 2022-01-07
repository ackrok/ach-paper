%% characterizing recordings

%% UNIT 3D PLOT -- Peak Latency (ms), Prop (ISI > 2s), Firing Rate (Hz)
roi = WT;
%
peakL = [roi.peakL]; 
fr    = [roi.fr]; 
pISI2 = [roi.pISI2];
%
clr = {'k','b','g'};
figure; hold on
for x = unique([roi.label])
    idx = find([roi.label] == x);
    %scatter3(pISI2(idx), peakL(idx), fr(idx), 'DisplayName', sprintf('Type %d',x));
    plot3(pISI2(idx),peakL(idx),fr(idx),'.','MarkerSize',10,'Color',clr{x},'DisplayName', sprintf('Type %d',x));
end
xlabel('prop (ISI) > 2s'); ylabel('Peak Latency (ms)'); zlabel('Firing Rate (Hz)');
title(sprintf('Control Units (n = %d), t1-%d, t2-%d, t3-%d',length(roi),length(find([roi.label] == 1)),length(find([roi.label] == 2)),length(find([roi.label] == 3))));
legend; grid on; xlim([0 1]); ylim([0 3]); zlim([0 50]); view(45,45); 
xticks([0:0.25:1]);

%% Waveform Comparison
figure;
subplot(1,2,1);
wf_cin = [WT(find([WT.label] == 2)).wf];

%% Firing Rate Comparison
sub = WT(find([WT.label] == 2));
spn = WT(find([WT.label] == 1));

plotcell = cell(1,2); plotcell{1} = log([spn.fr]); plotcell{2} = log([sub.fr]);
figure; violinplot(plotcell); 
xticklabels({'pSPN','pCIN'}); 
ylabel('Firing Rate (Hz)'); yticks([-3:3]); yticklabels({'0.001','0.01','0.1','0','1','10','100'}); grid on
title(sprintf('Firing Rate of pSPNs (n=%d), pCINs (n=%d)',length(spn),length(sub)));

%% Plot Average Waveforms
fig = figure; fig.Position(3) = 1800;
sp(1) = subplot(1,3,1);
sub = WT(find([WT.label] == 1));
shadederrbar([-4:1/30:4],nanmean([sub.wf],2),nanstd([sub.wf],[],2),'k');
xlabel('Time (ms)'); ylabel('Voltage (uV)'); title('SPN');
sp(2) = subplot(1,3,2);
sub = WT(find([WT.label] == 2));
wf = [sub.wf];
a = [];
for x = 1:length(sub)
    [m,ii] = max(abs(sub(x).wf));
    if sub(x).wf(ii) < 0; a = [a,x]; end
end
shadederrbar([-4:1/30:4],nanmean(wf(:,a),2),nanstd(wf(:,a),[],2),'b');
xlabel('Time (ms)'); ylabel('Voltage (uV)'); title('CIN');
sp(3) = subplot(1,3,3);
sub = WT(find([WT.label] == 3));
shadederrbar([-4:1/30:4],nanmean([sub.wf],2),nanstd([sub.wf],[],2),'k');
xlabel('Time (ms)'); ylabel('Voltage (uV)'); title('FSI');
linkaxes(sp,'y'); linkaxes(sp,'x'); xlim([-2 3])

%% Plot ACGs
%load acg_cin, acg_spn, xLin from file R:tritsn01Lab/In Vivo/DATAcomb/acg_cin+spn_210227.mat
%also need WT to be loaded
acg = cell(1,2); acg{1} = acg_spn; acg{2} = acg_cin;
fig = figure; fig.Position(3) = 1800;
clr = {'k','b'}; ttl = {'SPN','CIN'};
for z = 1:2
    sp(z) = subplot(1,3,z); 
    sub = WT(find([WT.label] == z));
    plot_mat = [];
    for x = 1:size(acg{z},2)
        plot_mat(:,x) = movmean(acg{z}(:,x),20)./length(sub(x).st);
    end
    x = xLin'; y = nanmean(plot_mat,2)'; sem = nanstd(plot_mat,[],2)';
    color = char2rgb(clr{z}); patchcolor = color+(1-color)*.8;
    yerru = y+sem; yerrl = y-sem;
    xpatch=[x,fliplr(x)]; ypatch=[yerru,fliplr(y)]; ypatch2=[y,zeros(1,length(y))];
    hold on
    fill(xpatch,ypatch,patchcolor,'FaceAlpha',0.5,'EdgeAlpha',0);
    fill(xpatch,ypatch2,patchcolor,'FaceAlpha',0.85,'EdgeAlpha',0);
    main = plot(x,y,'-','Color',color); hold off
    %shadederrbar(xLin, nanmean(plot_mat,2), nanstd(plot_mat,[],2), 'k');
    xlim([-500 500]); ylim([0 30]); xlabel('Lag (ms)'); title(ttl{z});
end
xticks([-500:200:500]);


%%load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullWT_2-10-21.mat')

figure; hold on
[~, edges] = histcounts(log10([0.01:0.05:100]));
clr = {'k','b'};
x = 2;
roi = cinWT; sub = roi(find([roi.label] == x)); fr = [sub.fr];
%[~, edges] = histcounts(log10(fr));
histogram(fr,10.^edges,'Normalization','probability','FaceColor',clr{x},'FaceAlpha',0.2,'EdgeAlpha',0.5); %histogram of ISI distribution
set(gca,'xscale','log') %set log-scale for x-axis
xlim([0.01 100]); xticks([10^-2, 10^-1, 10^0, 10^1, 10^2]);
xticklabels({'0.01','0.1','1','10','100'});
xlabel('spikes/s'); ylabel('proportion of units');
title(sprintf('pCIN firing rate: mu = %1.2f Hz (n = %d)',nanmean(fr),length(sub)));
axis('square')
 
%% Plotting Waveforms
wf_wt = []; group = 2; roi = cinWT;
idx = find([roi.label] == group);
for x = 1:length(idx); n = idx(x) ;
    thiswf = roi(n).wf; 
    switch length(thiswf)
        case 181
            wf_wt(:,x) = zscore(thiswf(:)); 
        case 241
            thiswf(212:end) = []; thiswf(1:30) = []; 
            wf_wt(:,x) = zscore(thiswf(:)); 
    end
end
wf_t = [-3:1/30:3]';
%%
figure; 
shadederrbar(wf_t,nanmean(wf_wt,2),nanstd(wf_wt,[],2),'b'); 
xlabel('Time (ms)'); ylabel('(z-score)'); ylim([-5 1]);
title(sprintf('pCIN waveform (n = %d units)',size(wf_wt,2)));
axis('square')

