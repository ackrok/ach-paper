%% Load data
% load('allACh_beh.mat');
% load('allACh_rawfp.mat')
raw = raw_ach; raw_fs = 2000;
beh = wt_ach(1:22); Fs = 50;

%% Load data
% load('ChATKO-ACh_beh.mat')
% load('ChATKO-ACh_rawfp.mat')
raw = raw_chatko;
beh = ach_chatko;

%% AUC during rest 
vec = []; 
time = [1/Fs:1/Fs:length(raw(1).fp(:,1))/Fs];
for x = 1:length(raw)
    for y = 1:size(raw(x).fp,2)
        % fp = fpfinal{x,y};
        fp = raw(x).fp(:,y);
        % prcval = prctile(fp, 1); fp = fp - prcval; % Subtract lower 1% of dF/F values
        if isempty(raw(x).onRest{y}); continue; end
        vec_tmp = [0];
        for z = 1:length(raw(x).onRest{y})
            range = [raw(x).onRest{y}(z):raw(x).offRest{y}(z)]; %Range of time, in samples, for this rest bout
            seg = fp(range); %Extract segment of photometry
            vec_tmp = [vec_tmp + trapz(seg)/(raw(x).offRest{y}(z)-raw(x).onRest{y}(z))]; % Adjust each rest period trapz by length of rest period
        end
%         vec = [vec, vec_tmp];
        vec = [vec, vec_tmp/length(raw(x).onRest{y})]; %Adjust for number of rest periods per sweep
    end
end
%
figure; violinplot(vec); ylabel('AUC per second'); title('AUC REST');
%figure; violinplot({vec_wt, vec_ko}); xticklabels({'WT','KO'});

%% AUC during MOVEMENT 
vec = []; 
for x = 1:length(raw)
    for y = 1:size(raw(x).fp,2)
        if isempty(beh(x).on{y}); continue; end
        vec_tmp = 0;
        for z = 1:length(beh(x).on{y})
            range = [beh(x).on{y}(z) : beh(x).off{y}(z)]; %Range of time, in samples, for this rest bout
%             seg = raw(x).fp(range,y);   %Extract segment of photometry
            seg = fpfinal{x,y}(range); %Extract segment of photometry and time vectors
            vec_tmp = [vec_tmp + trapz(seg)/(beh(x).off{y}(z) - beh(x).on{y}(z))]; %Adjust each rest period trapz by length of rest period
        end
%         vec = [vec, vec_tmp];
        vec = [vec, vec_tmp/length(beh(x).on{y})]; %Adjust for number of rest periods per sweep
    end
end

figure; violinplot(vec); ylabel('AUC (dF/F) per second'); title('AUC MOV'); grid on

%% REST vs MOV
figure; 
violinplot({vec_wt_rest, vec_ko_rest, vec_wt_mov, vec_ko_mov}); 
xticklabels({'WT REST','KO REST','WT MOV','KO MOV'});
ylabel('AUC (per second)'); grid on
title('AUC: WT vs KO, subtract bottom 1%');
