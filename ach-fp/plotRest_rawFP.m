%% SUBPLOTS: rawFP during each rest period (1x sweep, all rests)
x = 1; y = 1; Fs = 50; raw_fs = 2000;

figure; plm = ceil(length(beh(x).onRest{y})/3); 
for z = 1:length(beh(x).onRest{y})
    sp(z) = subplot(plm,3,z);
    raw_y = raw(x).rawfp(:,y); raw_t = [1/raw_fs : 1/raw_fs : length(raw_y)/raw_fs]'; plot(raw_t, raw_y)
%     final_y = raw(x).fp(:,y); final_t = beh(x).time; plot(final_t,final_y);
    xlim([beh(x).onRest{y}(z)/Fs beh(x).offRest{y}(z)/Fs]); 
    title(sprintf('rest #%d',z)); grid on
end; linkaxes(sp,'y');

%% OVERLAY: rawFP for all rest periods (1x sweep, all rests)
x = 5; y = 1; dsRate = 40;
signal = raw(x).rawfp(:,y); %[beh(x).vel(1,y); diff(movmean(beh(x).vel(:,y),10))]; 
figure; hold on
for z = 1:length(raw(x).onRest{y})
    %plot(raw(x).rawfp([beh(x).onRest{y}(z)*dsRate:beh(x).offRest{y}(z)*dsRate],y));
    plot(signal([raw(x).onRest{y}(z)*dsRate:raw(x).offRest{y}(z)*dsRate])); 
end
title(sprintf('%s-swp#%d - rawFP',raw(x).rec,y)); ylabel('Voltage (V)');


%% Plot: velocity + rest-bouts + fp
x = 1; 
figure;
for y = 1:size(beh(x).fp,2)
    subplot(size(beh(x).fp,2),1,y); hold on
    bouts = zeros(1,length(beh(x).time));
    for ii = 1:length(beh(x).onRest{y})
        bouts(beh(x).onRest{y}(ii):beh(x).offRest{y}(ii)) = 1;
    end
    area(beh(x).time, max(beh(x).vel(:,y)).*bouts,'FaceColor','r','FaceAlpha',0.2,'EdgeAlpha',0,'ShowBaseLine','off','DisplayName','bouts')
    plot(beh(x).time, beh(x).vel(:,y), 'k');
    yyaxis right; plot(beh(x).time, beh(x).fp(:,y), 'g'); 
end

%% Extract FP during rest
x = 1;
fp = {}; t = {}; 
for y = 1:size(beh(x).fp,2)
    for z = 1:length(beh(x).onRest{y})
        fp{z,y} = beh(x).fp([beh(x).onRest{y}(z):beh(x).offRest{y}(z)],y);
        t{z,y}  = beh(x).time(beh(x).onRest{y}(z):beh(x).offRest{y}(z));
    end
end

%% Plot: subplots for each rest-period
x = 1; y = 2;

figure; hold on
plm = floor(sqrt(length(beh(x).onRest))); pln = ceil(length(beh(x).onRest)/plm);
for z = 1:length(beh(x).onRest)
    subplot(plm,pln,z);
    plot(beh(x).time(beh(x).onRest(z):beh(x).offRest(z)), ...
        beh(x).FP{1}([beh(x).onRest(z):beh(x).offRest(z)]), 'k');
   % plot(raw(x).rawfp([beh(x).onRest{y}(z)*dsRate:beh(x).offRest{y}(z)*dsRate],y))
end

%% PLOT: single rest period
x = 1; y = 3; z = 2;

lpCut = 10; filtOrder = 8;
dsRate = 40; dsType = 2;
interpType = 'linear'; fitType = 'line'; basePrc = 10; winSize = 10; winOv = 0;

figure; 
sp(1) = subplot(4,1,1);
    final_t = beh(x).time; final_y = beh(x).fp(:,y);
    plot(final_t, final_y, 'k');
    title(sprintf('%s-acq#%d: rest #%d',beh(x).rec,y,z),'Interpreter','none'); grid on
sp(2) = subplot(4,1,2);
    raw_y = raw(x).rawfp(:,y); raw_t = [1/raw_fs : 1/raw_fs : length(raw_y)/raw_fs]';
    plot(raw_t, raw_y)
    title('raw voltage'); grid on
sp(3) = subplot(4,1,3);     
    nbFP = filterFP(raw_y,raw_fs,lpCut,filtOrder,'lowpass');
    plot(raw_t, nbFP)
    title('low-pass filter 10Hz'); grid on
sp(4) = subplot(4,1,4);
    plot_fp = raw(x).rawfp(:,y)./noise_10ms(x,y);
    plot(raw_t, plot_fp)
    title('divide rawV by average 10ms bin noise'); grid on
% sp(5) = subplot(5,1,5);
%     nbFP = downsampleTLab(nbFP,dsRate,dsType);
%     [dF] = baselineFP(nbFP,interpType,fitType,basePrc,winSize,winOv,Fs);
%     plot(final_t(1:length(dF)), dF)
%     title('baseline (dF/F) after downsampling'); grid on

linkaxes(sp,'x'); xlim([beh(x).onRest{y}(z)/Fs beh(x).offRest{y}(z)/Fs]); 

%% OVERLAY: photometry + spike times
x = 39; y = 1; z = 16; Fs = 50;

figure; hold on
final_t = beh(x).time; final_y = beh(x).fp(:,y);
plot(final_t, final_y, 'k');
sub = cinwt(find(strcmp({cinwt.rec},beh(x).rec)));
for aa = 1:length(sub)
    plot(sub(aa).st, aa/10*ones(length(sub(aa).st),1), '.b');
end
title(sprintf('%s: rest #%d',beh(x).rec,z),'Interpreter','none'); grid on
xlim([beh(x).onRest{y}(z)/Fs beh(x).offRest{y}(z)/Fs]); 
