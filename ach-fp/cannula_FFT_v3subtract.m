load('C:\Users\Anya\Desktop\FP_LOCAL\4_glutamate\raw_data\cannula_demod_noFilter.mat');
% Processing: digitalLIA, but commented out all filtering (lines 74, 78, 80)

%% FFT
p1_mat = [];
for x = 1:length(new)
vec = new(x).demod(9000001:end);
Fs = new(x).rawFs;

T = 1/Fs;               % Sampling period
L = length(vec);        % Length of signal
vec(isnan(vec)) = [];
fftACh = fft(vec);      % Discrete Fourier Transform of photometry signal
P2 = abs(fftACh/L);     % Two-sided spectrum P2
P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;     % Frequency domain vector
P1 = medfilt1(P1);      % Median filter initial FFT

p1_mat(:,x) = movmean(P1,500);
end
p1_all = p1_mat;
fprintf('FFT done!\n');

%% PLOT all traces on logarithmic axis
ds = 100;
figure;
loglog(f(1:ds:end), p1_mat((1:ds:end),1)); hold on
for x = 2:length(new)
    loglog(f(1:ds:end), p1_mat((1:ds:end),x));
end
xlim([0.01 100]);
xlabel('Frequency'); ylabel('Power (a.u.)');
title('rDA1m cannula (n = 4)'); % title('ACh3.0 cannula (n = 4)');

%% PLOT SUBTRACTED FFT
p1_mat = p1_all(:,[1:8]);

fig = figure; fig.Position([3 4]) = [442 362]; hold on
ds = 50;
tmp = []; flog = log10(f); 
r = [find(flog == -2):find(flog == 2)];
flog = flog(r);
for x = 1:4
    tmp(:,x) = normalize(log10(p1_mat(r,x)),'range') - normalize(log10(p1_mat(r,x+4)),'range');
    % plot(flog(1:ds:end), tmp(1:ds:end,x), 'Color', [0 0 0 0.2]);
end
plot([-2 2],[0 0],'--k');
shadederrbar(flog(1:ds:end), nanmean(tmp(1:ds:end,:),2), SEM(tmp(1:ds:end,:),2), 'm'); %[0.05 0.75 0.45]);
xlabel('Frequency'); xlim([-2 2]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
ylabel('Power (a.u.)'); ylim([-0.2 0.2]); yticks([-0.2:0.1:0.2]);
% title(sprintf('ACh3.0 %s - %s (n = 4), ds %d',new(1).rx,new(5).rx,ds)); 
title(sprintf('rDA1m %s - %s (n = 4), ds %d',new(1).rx,new(5).rx,ds));
axis('square')

%% PLOT AUC [1 4] Hz
p1_mat = p1_all;

tmp = [];
r = [find(f == 0.01):find(f == 100)]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_14 = [find(f_sub == 1):find(f_sub == 4)]; % AUC from [1 4] Hz

for x = 1:size(p1_mat,2)
    vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100]
    auc_tmp = trapz(vec_norm(r_14));
    
    tmp(:,x) = vec_norm;
    % tmp(x) = auc_tmp/length(r_14);
end
% auc = reshape(tmp, [4 3]);

%%
a = auc_ach(:,[1:3]) - nanmean(auc_ach(:,4));
figure; violinplot(a); 
hold on; plot([0.5 3.5],[0 0],'--k');
ylabel('AUC [1 4]Hz subtract scopolamine avg')
xticklabels({'aCSF','d1d2','glu'});
p = []; for x = 1:3; [~,p(x)] = ttest(a(:,x)); end
title(sprintf('ACh AUC [1 4]Hz subtract scopolamine \n t-test aCSF (%1.3f), d1d2 (%1.3f), glu (%1.3f)',...
    p(1),p(2),p(3)));
