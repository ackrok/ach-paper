%% Load data
load('allACh_beh.mat');
load('allACh_rawfp.mat')
raw = raw_ach; raw_fs = 2000;
beh = wt_ach(1:22); Fs = 50;

%% Load data
load('ChATKO-ACh_beh.mat')
load('ChATKO-ACh_rawfp.mat')
raw = raw_chatko;
beh = ach_chatko;

%% Peak to peak noise in raw data
noise_10ms = []; noise_max = []; noise_min = []; % noise_2s = [];
for x = 1:length(beh)
    for y = 1:size(beh(x).fp,2)
        raw_fp = raw(x).rawfp(:,y);
        raw_fp_rest = [];
        if isempty(beh(x).onRest{y})
            noise_10ms(x,y) = nan; noise_2s(x,y) = nan;
            continue; end
        for z = 1:length(beh(x).onRest{y})
            raw_fp_rest = [raw_fp_rest; raw_fp([beh(x).onRest{y}(z)*dsRate : beh(x).offRest{y}(z)*dsRate])];
        end
        raw_fp = raw_fp_rest; % overwrite raw_fp to be only rest periods
        
        bin = 0.01; % in seconds
        tmp = [];
        for aa = 1:[(length(raw_fp)/raw_fs)/bin]
            raw_fp_seg = raw_fp([aa*(bin*raw_fs)-(bin*raw_fs-1)]:[aa*(bin*raw_fs)]);
            tmp(aa,1) = max(raw_fp_seg);
            tmp(aa,2) = min(raw_fp_seg);
            tmp(aa,3) = mean(raw_fp_seg);
        end
        noise_10ms(x,y) = mean([tmp(:,1) - tmp(:,2)]./tmp(:,3));
        noise_max(x,y) = mean([tmp(:,1)]);
        noise_min(x,y) = mean([tmp(:,2)]);
        
        
%         bin = 2; % in seconds
%         tmp = [];
%         for aa = 1:[(length(raw_fp)/raw_fs)/bin]
%             raw_fp_seg = raw_fp([aa*(bin*raw_fs)-(bin*raw_fs-1)]:[aa*(bin*raw_fs)]);
%             tmp(aa,1) = max(raw_fp_seg);
%             tmp(aa,2) = min(raw_fp_seg);
%             tmp(aa,3) = mean(raw_fp_seg);
%         end
%         noise_2s(x,y) = mean([tmp(:,1) - tmp(:,2)]./tmp(:,3));
    end
end

%fprintf('%s-swp%d: max-min noise in %1.2fs bins = %1.3f\n',beh(x).rec,y,bin,mean([test(:,1) - test(:,2)]));

% vec = noise_2s - noise_10ms; vec = vec(:);      % Create Data
% SEM = nanstd(vec)/sqrt(length(vec));            % Standard Error
% ts = tinv([0.025  0.975],length(vec)-1);        % T-Score
% CI = nanmean(vec) + ts*SEM;                     % Confidence Intervals
% 
% fprintf('%s 95p CI: [%1.4f %1.4f](V)\n',beh(x).FPnames{1},CI(1),CI(2));
% WT 95p CI: [0.1164 0.1473](V)
% ChATcKO ACh 95p CI: [0.1451 0.2972](V)
