%% GLM: ACh photometry ~ CIN spike times

%%
% fPath = '/Users/akrok/Desktop/IV_local';
% fName = uigetfile([fPath,'/*.mat']); fName = {fName};
% load(fullfile(fPath,fName{1})); 
out = runglm_AChCIN(data);

%%
figure; hold on
for x = 1:length(out)
    scatter(x.*ones(length(out(x).X_var)-1,1), out(x).p_avg(2:end), 'filled');
end
line([0.5, length(out)+0.5], [0.05 0.05]);
xticks([1:length(out)]); xticklabels({out.date});
ylabel('p avg');
title('IV043 GLM: ACh ~ CINst')

%%
mat = struct;
for x = 1:length(out)
    for y = 1:length(out(x).X_var)-1
        z = 1+length(mat);
        mat(z).mouse = out(x).mouse;
        mat(z).date = out(x).date;
        mat(z).n = out(x).X_var{y+1};
        mat(z).p_all = out(x).p(y+1,:);
        mat(z).p_avg = mean([out(x).p(y+1,:)],2); 
        mat(z).p_sem = SEM([out(x).p(y+1,:)],2);
        mat(z).rmse = out(x).rmse;
        mat(z).rmse_rmv = out(x).rmv_rmse(y);
        mat(z).rmse_cons = out(x).rmse_constant;
        mat(z).r2 = out(x).r2;
        mat(z).r2_rmv = out(x).rmv_r2(y);
        mat(z).r2_cons = out(x).r2_constant;
        
    end
end
mat(1) = []; 
        
        
