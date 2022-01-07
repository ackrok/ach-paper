load('R:\homes\ack466\ACh paper\2_ACh_fluctuations\allACh_mod.mat')
load('R:\homes\ack466\ACh paper\2_ACh_fluctuations\allACh_beh.mat')
load('C:\Users\Anya\Desktop\IV_LOCAL\data_comb\fullcin_6-2021_WT.mat')

%% MOD ACh
beh = modACh;   % CHANGE
fp_mu = [];     % initialize
for x = 1:length(beh)
    fp = beh(x).FP{1};
    idx = {};
    idx{1} = extractEventST([1:length(fp)]',beh(x).onRest,beh(x).offRest,1);
    idx{2} = extractEventST([1:length(fp)]',beh(x).on,beh(x).off,1);
    for ii = 1:2
        fp_mu(ii,x) = nanmean(fp(idx{ii})); end
end
fp_delta = (fp_mu(2,:) - fp_mu(1,:))./fp_mu(1,:);

tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp);
fp_mu_an = []; fp_delta_an = [];
for x = 1:length(uni)
    ii = find(strcmp(tmp, uni{x}));
    fp_mu_an(:,x) = nanmean(fp_mu(:,ii),2); 
    fp_delta_an(x) = nanmean(fp_delta(:,ii));
end

%% CW ACh
beh = wt_ach;   % CHANGE
fp_mu = [];     % initialize
for x = 1:22
    fp_sub = cell(length(beh(x).on),1);
    for y = 1:length(beh(x).on)
        fp = beh(x).fp(:,y); fp = fp - prctile(fp, 1);
        idx = {};
        idx{1} = extractEventST([1:length(fp)]',beh(x).onRest{y},beh(x).offRest{y},1);
        idx{2} = extractEventST([1:length(fp)]',beh(x).on{y},beh(x).off{y},1);
        for ii = 1:2; fp_sub{ii} = [fp_sub{ii}; fp(idx{ii})]; end % concatenate photometry from multiple short recordings
    end
    for ii = 1:2
        fp_mu(ii,x) = nanmean(fp_sub{ii}); end
end
fp_delta = (fp_mu(2,:) - fp_mu(1,:))./fp_mu(1,:);

tmp = {}; for x = 1:22; tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp);
% fp_mu_an = []; fp_delta_an = [];
for x = 1:length(uni)
    ii = find(strcmp(tmp, uni{x}));
    % fp_mu_an(:,x) = nanmean(fp_mu(:,ii),2); fp_delta_an(x) = nanmean(fp_delta(:,ii));
    fp_mu_an = [fp_mu_an, nanmean(fp_mu(:,ii),2)]; % concatenate with MOD
    fp_delta_an = [fp_delta_an, nanmean(fp_delta(:,ii))]; % concatenate with MOD
end

%% CIN firing rate
sub = cinWT; beh = behWT;
fr_mu = nan(2,length(sub));
for x = 1:length(sub)
    b = find(strcmp({beh.rec},sub(x).rec)); 
    if isempty(beh(b).on); continue; end
    st = sub(x).st;
    fr_mu(1,x) = 1/mean(diff(extractEventST(st,beh(b).onRest,beh(b).offRest,0)));
    fr_mu(2,x) = 1/mean(diff(extractEventST(st,beh(b).on,beh(b).off,0)));
end
fr_delta = (fr_mu(2,:) - fr_mu(1,:))./fr_mu(1,:);

tmp = {}; for x = 1:length(sub); tmp{x} = strtok(sub(x).rec,'_'); end
uni = unique(tmp);
for x = 1:length(uni)
    ii = find(strcmp(tmp, uni{x}));
    fr_mu_an(:,x) = nanmean(fr_mu(:,ii),2); 
    fr_delta_an(x) = nanmean(fr_delta(:,ii));
end

%% PLOT
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
plot([1.15; 1.85].*ones(2,length(fp_delta_an)),fp_mu_an, ':k');
errorbar(nanmean(fp_mu_an,2),nanstd(fp_mu_an,[],2),'k');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'imm','mov'});
ylabel('ACh3.0 (% dF/F)');
axis('square');
title(sprintf('ACh3.0 (n = %d mice) (p = %1.3f)',size(fp_mu_an,2),signrank(fp_mu_an(1,:),fp_mu_an(2,:))));

subplot(1,3,2); hold on
% plot([1.15; 1.85].*ones(2,length(fr_delta)),fr_mu, ':k');
% errorbar(nanmean(fr_mu,2),SEM(fr_mu,2),'k');
plot([1.15; 1.85].*ones(2,length(fr_delta_an)),fr_mu_an, ':k');
errorbar(nanmean(fr_mu_an,2),nanstd(fr_mu_an,[],2),'k');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'imm','mov'});
ylabel('CIN Firing Rate (Hz)'); ylim([2 10]); yticks([2:2:10]);
axis('square');
title(sprintf('pCIN (n = %d mice) (p = %1.3f)',size(fr_mu_an,2),signrank(fr_mu_an(1,:),fr_mu_an(2,:))));

subplot(1,3,3); hold on
plot([1 2],[0 0],'--k');
violinplot({fp_delta_an,fr_delta_an}); xlim([0.5 2.5]); xticklabels({'ACh3.0','pCIN'});
ylabel('delta increase mvmt - rest'); yticks([-0.5:0.5:1]);
axis('square');
title('ACh3.0 vs pCIN FR delta increase');
movegui(gcf,'center');