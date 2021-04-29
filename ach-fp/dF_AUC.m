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
    b = find(strcmp({beh.rec},raw(x).rec));
    for y = 1:length(beh(b).onRest)
        % fp = fpfinal{x,y};
        fp = raw(x).fp_prc1(:,y);
        %
        if isempty(beh(b).onRest{y}); continue; end
        vec_tmp = [0];
        for z = 1:length(beh(b).onRest{y})
            range = [beh(b).onRest{y}(z):beh(b).offRest{y}(z)]; %Range of time, in samples, for this rest bout
            seg = fp(range); %Extract segment of photometry
            vec_tmp = [vec_tmp + trapz(seg)/(beh(b).offRest{y}(z)-beh(b).onRest{y}(z))]; % Adjust each rest period trapz by length of rest period
        end
%         vec = [vec, vec_tmp];
        vec = [vec, vec_tmp/length(beh(b).onRest{y})]; %Adjust for number of rest periods per sweep
    end
end
%
 figure; violinplot(vec); ylabel('AUC per second'); title('AUC REST');
% figure; violinplot({vec_wt, vec_gfp, vec_ko}); xticklabels({'WT','GFP','KO'}); ylabel('AUC per second'); title('AUC REST');

%% AUC during MOVEMENT 
vec = []; 
time = [1/Fs:1/Fs:length(raw(1).fp(:,1))/Fs];
for x = 1:length(raw)
    b = find(strcmp({beh.rec},raw(x).rec));
    for y = 1:length(beh(b).on)
        % fp = fpfinal{x,y};
        fp = raw(x).fp_prc1(:,y);
        %
        if isempty(beh(b).on{y}); continue; end
        vec_tmp = [0];
        for z = 1:length(beh(b).on{y})
            range = [beh(b).on{y}(z):beh(b).off{y}(z)]; %Range of time, in samples, for this rest bout
            seg = fp(range); %Extract segment of photometry
            vec_tmp = [vec_tmp + trapz(seg)/(beh(b).off{y}(z)-beh(b).on{y}(z))]; % Adjust each rest period trapz by length of rest period
        end
%         vec = [vec, vec_tmp];
        vec = [vec, vec_tmp/length(beh(b).on{y})]; %Adjust for number of rest periods per sweep
    end
end
%
% figure; violinplot(vec); ylabel('AUC (dF/F) per second'); title('AUC MOV');
% figure; violinplot({vec_wtMov, vec_gfpMov, vec_koMov}); xticklabels({'WT','GFP','KO'}); ylabel('AUC (dF/F) per second'); title('AUC MOV');

%% REST vs MOV
figure; 
violinplot({vec_wt_rest, vec_ko_rest, vec_wt_mov, vec_ko_mov}); 
xticklabels({'WT REST','KO REST','WT MOV','KO MOV'});
ylabel('AUC (per second)'); grid on
title('AUC: WT vs KO, subtract bottom 1%');
