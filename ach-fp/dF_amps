load('allACh_beh.mat'); load('alltdT_beh.mat'); load('ChATKO-ACh_beh.mat');
beh = wt_ach(1:22);
beh = tdt(1:4);
beh = ach_chatko(2:end);

%% Proper Rest Threshold
velThres = 2.5;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;
for x = 1:length(beh)
    beh(x).onRest = cell(1,length(beh(x).on)); beh(x).offRest = beh(x).onRest;
    for y = 1:size(beh(x).vel,2)
        vel = beh(x).vel(:,y); vel = abs(vel);
        [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
        [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
        onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
        beh(x).onRest{y} = onsetInd; beh(x).offRest{y} = offsetInd;
    end
end; fprintf('Done %s: new rest threshold = 0.25.\n',beh(1).FPnames{1});

%% Extract ACh dF/F values for REST vs MOV
mat = struct;
for x = 1:length(beh)
    mat(x).rec = beh(x).rec;
    mat(x).FPnames = beh(x).FPnames;
    mat(x).fp = raw(x).fp(:);
    mat(x).fp_mvmt = []; mat(x).fp_rest = [];
    for y = 1:length(beh(x).on)
%%         
%         fp = beh(x).fp(:,y); % Extract photometry signal for this sweep
%         fp = (fp - min(fp))/(max(fp) - min(fp)); % Min-Max Normalization
          fp = raw(x).fp(:,y); % Extract photometry signal for this sweep
%         fp = fpfinal{x,y};

%%        
        t_mr = cell(2,2); t_mr{1,1} = beh(x).on{y}; t_mr{2,1} = beh(x).off{y}; 
        t_mr{1,2} = beh(x).onRest{y}; t_mr{2,2} = beh(x).offRest{y}; 
        tmp_sig = cell(1,2);
        for a = 1:2 % Repeat over movement, rest
            tmp_vec = []; % Clear/initialize vector
            for z = 1:length(t_mr{1,a}) % Iterate over all mov/rest periods
                range = [t_mr{1,a}(z), t_mr{2,a}(z)];  % Sample range
                tmp_vec = [tmp_vec; fp(range(1):range(2))]; % Extract values within range
            end
            tmp_sig{a} = tmp_vec;
        end
        mat(x).fp_mvmt = [mat(x).fp_mvmt; tmp_sig{1}]; % Concatenate to output matrix
        mat(x).fp_rest = [mat(x).fp_rest; tmp_sig{2}]; % Concatenate to output matrix
    end
end
fprintf('Done %s: extracted REST + MOV photometry signal.\n',beh(1).FPnames{1});

fp_rest_all = []; fp_mvmt_all = [];
for x = 1:length(mat)
    fp_rest_all = [fp_rest_all; mat(x).fp_rest]; 
    fp_mvmt_all = [fp_mvmt_all; mat(x).fp_mvmt];
end

%%
figure; plm = floor(sqrt(length(mat))); pln = ceil(length(mat)/plm);
for x = 1:length(mat)
    if isempty(mat(x).fp_mvmt) || isempty(mat(x).fp_rest); continue; end
    sp(x) = subplot(plm,pln,x); hold on
    histogram(mat(x).fp_rest,'BinWidth',0.1,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability');
    histogram(mat(x).fp_mvmt,'BinWidth',0.1,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability');
    title(sprintf('%s',mat(x).rec));
end; linkaxes(sp,'x');

%% HISTOGRAM: REST vs MOV
figure; hold on
yyaxis left
histogram(fp_rest_all,'BinWidth',0.1,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','REST');
yyaxis right
histogram(fp_mvmt_all,'BinWidth',0.1,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','MOV');
xlabel('% dF/F'); ylabel('Probability'); grid on; legend; %xlim([0 1]) %xlim([-5 20])
title(sprintf('%s',beh(1).FPnames{1}))
%title(sprintf('%s: muREST = %1.3f, muMOV = %1.3f (p = %1.3f)',beh(1).FPnames{1},nanmean(fp_rest_all),nanmean(fp_mvmt_all),ranksum(fp_rest_all,fp_mvmt_all)));

%% HISTOGRAM: REST WT vs REST KO
figure; y1 = wt_fprest; y2 = ko_fprest;
yyaxis left; ylabel('Probability');
histogram(y1,'BinWidth',0.1,'FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','WT');
yyaxis right; ylabel('Probability');
histogram(y2,'BinWidth',0.1,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','ChATcKO');
xlabel('% dF/F'); grid on; legend;
title(sprintf('dF/F - 25p : muWT = %1.3f, muKO = %1.3f (p = %1.3f)',nanmean(y1),nanmean(y2),ranksum(y1,y2)));


%% COMPARE PLOT: tdTomato vs ACh
% tdt_fp = cell(1,2); tdt_fp{1} = fp_mvmt_all; tdt_fp{2} = fp_rest_all;
% wt_fp = cell(1,2); wt_fp{1} = fp_mvmt_all; wt_fp{2} = fp_rest_all;
figure; hold on
plot(nanmean(tdt_fp{2})*ones(5,1),[0:0.01:0.04],'--r','DisplayName','tdT-muREST');
plot(nanmean(tdt_fp{1})*ones(5,1),[0:0.01:0.04],'--g','DisplayName','tdT-muMOV');
plot(nanmean(wt_fp{2})*ones(5,1),[0:0.01:0.04],'r','DisplayName','ACh-muREST');
plot(nanmean(wt_fp{1})*ones(5,1),[0:0.01:0.04],'g','DisplayName','ACh-muMOV');
histogram(fp_rest_all,'BinWidth',0.1,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','ACh-REST');
histogram(fp_mvmt_all,'BinWidth',0.1,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','ACh-MOV');
xlabel('% dF/F'); ylabel('Probability'); grid on; legend; xlim([-5 20])
title(sprintf('%s: muREST = %1.3f, muMOV = %1.3f || %s: muREST = %1.3f, muMOV = %1.3f',...
    'ACh',nanmean(wt_fp{2}),nanmean(wt_fp{1}),'tdT',nanmean(tdt_fp{2}),nanmean(tdt_fp{1})));

% figure; hold on
% x = 1; y = tdt_fp_mu; bar(x,nanmean(y),'r','FaceAlpha',0.2,'DisplayName','tdT-muREST'); 
% x = 2; y = tdt_fp{1}; bar(x,nanmean(y),'g','FaceAlpha',0.2,'DisplayName','tdT-muMOV'); 
% x = 3; y = wt_fp{2}; bar(x,nanmean(y),'r','FaceAlpha',0.5,'DisplayName','ACh-muREST'); 
% x = 4; y = wt_fp{1}; bar(x,nanmean(y),'g','FaceAlpha',0.5,'DisplayName','ACh-muMOV');
% xticks([1:4]); xticklabels({'tdT-muREST','tdT-muMOV','ACh-muREST','ACh-muMOV'});

%% BAR PLOT: ACh vs tdT -- mean of distribution
wt_fp_mu = nan(length(mat_ach),2);
for x = 1:length(mat_ach) 
    wt_fp_mu(x,1) = nanmean(mat_ach(x).fp_mvmt); wt_fp_mu(x,2) = nanmean(mat_ach(x).fp_rest); end
tdt_fp_mu = nan(length(mat_tdt),2);
for x = 1:length(mat_tdt) 
    tdt_fp_mu(x,1) = nanmean(mat_tdt(x).fp_mvmt); tdt_fp_mu(x,2) = nanmean(mat_tdt(x).fp_rest); end

figure; hold on
x = 1; y = tdt_fp_mu(:,2); bar(x,nanmean(y),'r','FaceAlpha',0.1,'DisplayName','tdT-muREST'); errorbar(x,nanmean(y),SEM(y,1),'k');
    %scatter(x*ones(length(y),1),y,'r');
x = 2; y = tdt_fp_mu(:,1); bar(x,nanmean(y),'g','FaceAlpha',0.1,'DisplayName','tdT-muMOV'); errorbar(x,nanmean(y),SEM(y,1),'k');
    %scatter(x*ones(length(y),1),y,'g'); 
    %plot([1*ones(1,length(y));2*ones(1,length(y))],[tdt_fp_mu(:,2)';tdt_fp_mu(:,1)'],':k');
x = 3; y = wt_fp_mu(:,2); bar(x,nanmean(y),'r','FaceAlpha',0.4,'DisplayName','ACh-muREST'); errorbar(x,nanmean(y),SEM(y,1),'k');
    %scatter(x*ones(length(y),1),y,'r');
x = 4; y = wt_fp_mu(:,1); bar(x,nanmean(y),'g','FaceAlpha',0.4,'DisplayName','ACh-muMOV'); errorbar(x,nanmean(y),SEM(y,1),'k');
    %scatter(x*ones(length(y),1),y,'g'); 
    %plot([3*ones(1,length(y));4*ones(1,length(y))],[wt_fp_mu(:,2)';wt_fp_mu(:,1)'],':k');
xticks([1:4]); xticklabels({'tdT-muREST','tdT-muMOV','ACh-muREST','ACh-muMOV'});
ylabel('% dF/F'); grid on; 

%% BAR PLOT: ACh vs ChATcKO -- mean of distribution
wt_fp_mu = nan(length(mat_ach),2);
for x = 1:length(mat_ach) 
    wt_fp_mu(x,1) = nanmean(mat_ach(x).fp_mvmt); wt_fp_mu(x,2) = nanmean(mat_ach(x).fp_rest); end
chat_fp_mu = nan(length(mat_chat),2);
for x = 1:length(mat_chat) 
    chat_fp_mu(x,1) = nanmean(mat_chat(x).fp_mvmt); chat_fp_mu(x,2) = nanmean(mat_chat(x).fp_rest); end

figure; hold on
x = 1; y = chat_fp_mu(:,2); bar(x,nanmean(y),'r','FaceAlpha',0.1,'DisplayName','ChATcKO-muREST'); errorbar(x,nanmean(y),SEM(y,1),'k');
x = 2; y = chat_fp_mu(:,1); bar(x,nanmean(y),'g','FaceAlpha',0.1,'DisplayName','ChATcKO-muMOV'); errorbar(x,nanmean(y),SEM(y,1),'k');
x = 3; y = wt_fp_mu(:,2); bar(x,nanmean(y),'r','FaceAlpha',0.4,'DisplayName','WT-ACh-muREST'); errorbar(x,nanmean(y),SEM(y,1),'k');
x = 4; y = wt_fp_mu(:,1); bar(x,nanmean(y),'g','FaceAlpha',0.4,'DisplayName','WT-ACh-muMOV'); errorbar(x,nanmean(y),SEM(y,1),'k');
xticks([1:4]); xticklabels({'ChATcKO','ChATcKO','WT-ACh','WT-ACh'});
ylabel('% dF/F'); grid on; 

%% HISTOGRAM: 95% confidence intervals
%mat = mat_ach;
%bin = 0.1; edges = [-5:bin:20]; mid = edges(1:end-1)+(bin/2); %for %dF/F
bin = 0.01; edges = [0:bin:1]; mid = edges(1:end-1)+(bin/2); %for %dF/F min max normalization

y_mov = []; y_rest = [];
for x = 1:length(mat)
    y_mov(x,:) = histcounts(mat(x).fp_mvmt,edges,'Normalization','probability');
    y_rest(x,:) = histcounts(mat(x).fp_rest,edges,'Normalization','probability');
end
y = y_rest; % CHANGE
N = size(y,1); yMean = nanmean(y); ySEM = nanstd(y)/sqrt(N);
CI95 = tinv([0.025 0.975], N-1);                        % Calculate 95% Probability Intervals of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));                  % Calculate 95% Confidence Intervals of all experiments at each value of 'x'
figure;
shadederrbar(mid, yMean, yCI95(2,:), 'r');

y = y_mov; % CHANGE
N = size(y,1); yMean = nanmean(y); ySEM = nanstd(y)/sqrt(N);
CI95 = tinv([0.025 0.975], N-1);                        % Calculate 95% Probability Intervals of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));                  % Calculate 95% Confidence Intervals of all experiments at each value of 'x'
hold on
shadederrbar(mid, yMean, yCI95(2,:), 'g');

xlabel('dF/F (%)'); ylabel('probability'); legend({'REST 95%CI','REST mean','MOV 95%CI','MOV mean'})
title(sprintf('%s: mean +/- 95p CI',mat(1).FPnames{1}));
%ylim([0 0.035]); yticks([0:0.01:0.03]); grid on;

%% BOX PLOT: 95% confidence intervals
plotvec = [tdt_fp{2};tdt_fp{1};wt_fp{2};wt_fp{1}]; % Combine data into a single column vector
groups = [repmat({'tdT-rest'},length(tdt_fp{2}),1); repmat({'tdT-mov'},length(tdt_fp{1}),1);...
    repmat({'ACh-rest'},length(wt_fp{2}),1); repmat({'ACh-mov'},length(wt_fp{1}),1)];
figure; boxplot(plotvec,groups,'Notch','on','symbol','');
ylabel('dF/F (%)'); title('Summary Statistics');
