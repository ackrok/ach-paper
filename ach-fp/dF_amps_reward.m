beh = behPav;
%%
mat = struct;
for x = 1:length(beh)
    mat(x).rec = beh(x).rec;
    mat(x).FPnames = beh(x).FPnames;
    fp = beh(x).FP{1}; Fs = beh(x).Fs; % Extract photometry from data
    mat(x).fp = fp; % All photometry dF/F values
    staTmp = getSTA(fp, beh(x).reward/Fs, Fs, [0 1]); % Reward alignment [0,1] from reward delivery
    mat(x).fp_rew = staTmp(:); % Reward photometry dF/F values
    
    acc = [beh(x).vel(1); diff(movmean(beh(x).vel,10))];
    [pks,locs] = findpeaks(acc,'MinPeakProminence',50,'MinPeakDistance',1); % Acceleration peaks
    locs = beh(x).time(locs); % Acceleration peaks, now in seconds
    locs(locs < 1 + beh(x).reward(end)/Fs) = []; % Remove accelerations during reward period
    staTmp = getSTA(fp, locs, Fs, [-0.4 1]); % Acceleration alignment [-0.4 1] from acc peak
    mat(x).fp_acc = staTmp(:); % Acceleration photometry dF/F values
    
    mat(x).fp_mvmt = []; mat(x).fp_rest = []; 
    t_mr = cell(2,2); t_mr{1,1} = beh(x).on; t_mr{2,1} = beh(x).off; 
    t_mr{1,2} = beh(x).onRest; t_mr{2,2} = beh(x).offRest; 
    tmp_sig = cell(1,2);
    for a = 1:2 % Repeat over movement, rest
        tmp_vec = []; % Clear/initialize vector
        for z = 1:length(t_mr{1,a}) % Iterate over all mov/rest periods
            if t_mr{1,a}(z) < beh(x).reward(end)+Fs; continue; end
            range = [t_mr{1,a}(z), t_mr{2,a}(z)];  % Sample range
            tmp_vec = [tmp_vec; fp(range(1):range(2))]; % Extract values within range
        end
        tmp_sig{a} = tmp_vec;
    end
    mat(x).fp_mvmt = [mat(x).fp_mvmt; tmp_sig{1}]; % Concatenate to output matrix
    mat(x).fp_rest = [mat(x).fp_rest; tmp_sig{2}]; % Concatenate to output matrix
end
% Extract and concatenate across all recordings:
fp_rew = [mat.fp_rew]; fp_mvmt = []; fp_rest = []; fp_acc = [];
for x = 1:length(mat)
    fp_mvmt = [fp_mvmt; mat(x).fp_mvmt]; 
    fp_rest = [fp_rest; mat(x).fp_rest]; 
    fp_acc = [fp_acc; mat(x).fp_acc]; 
end

%% HISTOGRAM
figure; hold on
histogram(fp_rew(:),'BinWidth',0.2,'FaceColor','c','FaceAlpha',0.5,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','REWARD');
histogram(fp_acc,'BinWidth',0.2,'FaceColor','k','FaceAlpha',0.5,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','ACC');
%histogram(fp_mvmt,'BinWidth',0.2,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','MOV');
histogram(fp_rest,'BinWidth',0.2,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability','DisplayName','REST');
xlabel('% dF/F'); ylabel('Probability'); grid on; legend; %xlim([0 1]) %xlim([-5 20])

%% HISTOGRAM SUB PLOTS
figure; plm = floor(sqrt(length(mat))); pln = ceil(length(mat)/plm);
for x = 1:length(mat)
    sp(x) = subplot(plm,pln,x); hold on
    histogram(mat(x).fp_rew,'BinWidth',0.2,'FaceColor','c','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability');
    title(sprintf('%s',mat(x).rec));
    if isempty(mat(x).fp_mvmt) || isempty(mat(x).fp_rest); continue; end
    histogram(mat(x).fp_rest,'BinWidth',0.2,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability');
    histogram(mat(x).fp_mvmt,'BinWidth',0.2,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.5,'Normalization','probability');
end; linkaxes(sp,'x');
%% PLOT 95% confidence intervals
bin = 0.2; edges = [-5:bin:20]; mid = edges(1:end-1)+(bin/2); %for %dF/F

y_val = cell(1,3);
for x = 1:length(mat)
    y_val{1}(x,:) = histcounts(mat(x).fp_rest,edges,'Normalization','probability');
    y_val{2}(x,:) = histcounts(mat(x).fp_rew,edges,'Normalization','probability');
    y_val{3}(x,:) = histcounts(mat(x).fp_acc,edges,'Normalization','probability');
    %y_val{4}(x,:) = histcounts(mat(x).fp_mvmt,edges,'Normalization','probability');
end

figure; clr = {'r','c','k','g'};
for z = 1:length(y_val)
    y = y_val{z};
    N = size(y,1); yMean = nanmean(y); ySEM = nanstd(y)/sqrt(N);
    CI95 = tinv([0.025 0.975], N-1);                        % Calculate 95% Probability Intervals of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:));                  % Calculate 95% Confidence Intervals of all experiments at each value of 'x'
    shadederrbar(mid, yMean, yCI95(2,:), clr{z}); hold on
end
xlabel('dF/F (%)'); ylabel('probability'); 
legend({'REST 95%CI','REST mean','REW 95%CI','REW mean','ACC 95%CI','ACC mean','MOV 95%CI','MOV mean'})
title(sprintf('%s: mean +/- 95p CI',mat(1).FPnames{1}));
%ylim([0 0.035]); yticks([0:0.01:0.03]); grid on;
