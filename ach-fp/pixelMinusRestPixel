%% subtracting rest values for pixel values of FP
% load full FP data structure
% cd '/Users/akrok/Desktop/FP_local/FP_comb'
% load('allACh_accPk_fpAlign_12-20.mat')
%% extract pixel values, subtracting mean of rest pixel values
mat = cell(length(beh),2); 
diffFs = 2000/50;
for x = 1:length(beh)
    mat{x,1} = []; mat{x,2} = [];
    for y = 1:length(beh(x).onRest) %iterate over all sweeps for recording
        fp = beh(x).fp(:,y); %extract photometry signal for this sweep
        restSub = [];
        for z = 1:length(beh(x).onRest{y})
            fp_seg = fp(beh(x).onRest{y}(z)*diffFs : beh(x).offRest{y}(z)*diffFs);
            restSub = [restSub; fp_seg/nanmean(fp_seg)];
        end
        movSub = [];
        for z = 1:length(beh(x).on{y})
            fp_seg = fp(beh(x).on{y}(z)*diffFs:beh(x).off{y}(z)*diffFs);
            movSub = [movSub; fp_seg/nanmean(fp_seg)];
        end
        mat{x,1} = [mat{x,1}; movSub - nanmean(restSub)]; %subtract mean of all rest pixel values
        mat{x,2} = [mat{x,2}; restSub - nanmean(restSub)]; %subtract mean of all rest pixel values
    end
end
%% extract pixel values
mat = cell(length(beh),2); diffFs = 1;
for x = 1:length(beh)
    mat{x,1} = []; mat{x,2} = [];
    for y = 1:length(beh(x).onRest) %iterate over all sweeps for recording
        fp = beh(x).fp(:,y); %extract photometry signal for this sweep
        fp = normalize(fp);
        restSub = [];
        for z = 1:length(beh(x).onRest{y})
            fp_seg = fp(beh(x).onRest{y}(z)*diffFs : beh(x).offRest{y}(z)*diffFs);
            restSub = [restSub; fp_seg];
        end
        movSub = [];
        for z = 1:length(beh(x).on{y})
            fp_seg = fp(beh(x).on{y}(z)*diffFs:beh(x).off{y}(z)*diffFs);
            movSub = [movSub; fp_seg];
        end
        mat{x,1} = [mat{x,1}; movSub - nanmean(restSub)]; %subtract mean of all rest pixel values
        mat{x,2} = [mat{x,2}; restSub - nanmean(restSub)]; %subtract mean of all rest pixel values
    end
end
%% plot all into subplots
figure; plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
for x = 1:length(beh)
sp(x) = subplot(plm,pln,x); hold on
histogram(mat{x,2},'BinWidth',0.25,'Normalization','probability','FaceColor','r','DisplayName','rest');
histogram(mat{x,1},'BinWidth',0.25,'Normalization','probability','FaceColor','g','DisplayName','mov');
title(sprintf('%s',beh(x).rec));
end
%% plot all together
vec_mov = []; vec_rest = [];
for x = 1:length(mat)
vec_mov = [vec_mov; mat{x,1}]; vec_rest = [vec_rest; mat{x,2}];
end
figure; hold on
histogram(vec_rest,'BinWidth',0.25,'Normalization','probability','FaceColor','r','DisplayName','rest');
histogram(vec_mov,'BinWidth',0.25,'Normalization','probability','FaceColor','g','DisplayName','mov');
title(sprintf('Pixel Values - mean(Rest Pixel Values) (n = %d recs)',length(mat)));
xlabel('dF/F'); ylabel('probability');

%%
figure; hold on
histogram(mat{x,1},'BinWidth',0.25,'Normalization','probability','FaceColor','r','DisplayName','rest');
histogram(mat{x,1},'BinWidth',0.25,'Normalization','probability','FaceColor','g','DisplayName','mov');
title(sprintf('Pixel Values - mean(Rest Pixel Values) (n = %d recs)',length(mat)));
xlabel('dF/F'); ylabel('probability');
