
%% CCG DELTA - AVERAGE
figure; hold on; sm = 5;
shadederrbar(time, movmean(nanmean(ccg50,2),sm), movmean(nanmean(ccg95,2),sm), 'k'); hold on
shadederrbar(time, movmean(nanmean(ccgDelta,2),sm), movmean(SEM(ccgDelta,2),sm), 'b'); 
xlabel('Lag (s)'); ylabel('CCG (deltaFR)'); title(sprintf('GATcKO: CCG (n = %d pairs)',size(ccgDelta,2)));
xlim([-1 1]);

%% CCG - UNIT PAIRS
% n = [mat.n]; m = [mat.m];
% figure; plm = floor(sqrt(size(ccg_mvmt,2))); pln = ceil(size(ccg_mvmt,2)/plm);
% for x = 1:size(ccg_mvmt,2)
%     sp(x) = subplot(plm,pln,x); hold on
%     plot(time, ccg_full(:,x), 'b'); %plot(time, ccg_mvmt(:,x), 'g');
%     xlim([-1 1]); title(sprintf('m #%d | ref #%d',m(x),n(x)));
% end
% linkaxes(sp,'y');

for x = 1:length(mat)
    figure; plm = floor(sqrt(length(mat(x).n))); pln = ceil(length(mat(x).n)/plm);
    for y = 1:length(mat(x).n)
        subplot(plm,pln,y); 
        shadederrbar(time, mat(x).shuffPrc{y}(:,2), mat(x).shuffPrc{y}(:,3)-mat(x).shuffPrc{y}(:,2), 'k'); hold on
        plot(time, mat(x).ccg(:,y), 'b');
        xlim([-1 1]); title(sprintf('%s: m#%d|ref#%d',mat(x).rec,mat(x).m(y),mat(x).n(y)));
    end
end

%% FIRING RATE - MOV vs REST
fr = [];
for x = 1:length(sub)
    idx_b = find(strcmp({beh.rec},sub(x).rec));
    fr(x,1) = 1/mean(diff(extractEventST(sub(x).st,beh(idx_b).on/50,beh(idx_b).off/50,1)));
    fr(x,2) = 1/mean(diff(extractEventST(sub(x).st,beh(idx_b).onRest/50,beh(idx_b).offRest/50,1)));
end

figure; hold on
plot([ones(1,length(fr));2*ones(1,length(fr))],fr',':k');
clr = {'g','r'};
for x = 1:2
    y = fr(:,x);
    bar(x, nanmean(y), clr{x}, 'FaceAlpha', 0.2);
    scatter(x*ones(length(y),1), y, clr{x});
    errorbar(x, nanmean(y), SEM(y,1), 'k');
end
xticks([1 2]); xticklabels({'MOV','REST'}); xlim([0.5 2.5])
ylabel('Firing Rate (Hz)');
title(sprintf('GAT1KO - MOV vs REST (n = %d units)',length(y)));

%% PLOT versus DISTANCE
x_var = [mat.dist]'; %x_var = [mat.fr]';
max_ccg = max(ccgDelta,[],1)';

tbl = table(x_var,max_ccg,'VariableNames',{'dist','maxCCG'}); %Create a table
mdl = fitlm(tbl,'maxCCG ~ dist'); %Fit a linear regression model with max(CCG) as response variable, unitDistance as predictor variable
ci = coefCI(mdl); %Find confidence intervals for the coefficients of the model
x_hat = [1:max(x_var)]; 
y_hat = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}*x_hat; %Model fit y
yerrl = ci(1,1) + ci(2,1)*x_hat; yerru = ci(1,2) + ci(2,2)*x_hat; %Lower and upper confidence intervals
xpatch = [x_hat,fliplr(x_hat)]; ypatch = [yerru,fliplr(yerrl)]; %Generate patch for plotting shaded confidence interval

figure; hold on
plot(x_hat,nanmean(ccg95(:))*ones(length(x_hat),1),'--k'); %Dashed line for 95% confidence interval
plot(x_var,max_ccg,'.b','MarkerSize',10);
fill(xpatch,ypatch,[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0); %Plotting shaded confidence interval
plot(x_hat,y_hat,'-r'); %Plotting model fit line
xlabel('Unit Distance (um)'); ylabel('Max CCG rest spikes (deltaFR)');
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)',size(ccgDelta,2)));
title(sprintf('CCG peak ~ Unit Distance (n = %d pairs)\nR-squared = %1.3f',size(ccgDelta,2),mdl.Rsquared.Ordinary));
xlim([-50 1150]); 
