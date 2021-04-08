%%
velThres = 0.25;
Fs = 50;
minRestTime = 4*Fs; minRunTime = 1*Fs;
timeThres = 4*Fs; timeShift = 0.5*Fs;

for x = 1:length(beh)
    vel = beh(b).vel; vel = abs(vel);
    [onsetInd,offsetInd] = getOnsetOffset(-vel,-velThres,minRunTime,minRestTime,1);
    [onsetInd,offsetInd] = adjOnsetOffset(onsetInd,offsetInd,timeThres,vel);
    onsetInd = onsetInd+timeShift; offsetInd = offsetInd-timeShift;
    beh(b).onRest = onsetInd; beh(b).offRest = offsetInd;
end

%% Instantaneous Firing Rate of REST spikes
%mat = struct; 
bin = 1; Fs = 50; acqFs = 30000;
h = waitbar(0,'Inst-FR of Rest Spikes');
for x = 1:length(sub)
    b = find(strcmp({beh.rec},sub(x).rec)); %Find matching behavior data
    st = sub(x).st; 
    instfr_rest = [];
    for z = 1:length(beh(b).onRest)
        bout_length = (beh(b).offRest(z)-beh(b).onRest(z))/Fs; %Length of bout
        logicalIndexes = st < beh(b).offRest(z)/Fs & st > beh(b).onRest(z)/Fs;
        tmp = get_inst_fr(st(logicalIndexes)-beh(b).onRest(z)/Fs, bout_length, bin, acqFs);
        instfr_rest = [instfr_rest; tmp(1:bout_length/bin)];
    end
    %instfr = get_inst_fr(st_rest, beh(b).time(end), bin, acqFs);
    %mat(x).rec = sub(x).rec; mat(x).n = sub(x).n;
    %mat(x).instfr_rest = instfr(1:beh(b).time(end)/bin);
    mat(x).instfr_rest = instfr_rest;
    full = movmean(mat(x).instfr_raw,10);
    mat(x).instfr_rest_norm = (movmean(instfr_rest,10) - min(full)) / (max(full) - min(full)); %Min-max normalization
    waitbar(x/length(sub),h);
end; close(h); fprintf('Done! \n');

%% IV069-rec01: inst fr for each unit
figure;
for x = 1:9
    sp(x) = subplot(9,1,x); plot(movmean(mat(x+85).instfr_rest,10),'b');
    ylabel('Inst Freq (Hz)'); grid on;
    title(sprintf('IV069-rec01 unit #%d',sub(x+85).n));
end
linkaxes(sp,'y'); ylim([0 10])
%% IV069-rec01: overlay all units + velocity trace
tmp = [];
for x = 1:9; tmp(:,x) = movmean(mat(x+85).instfr_rest,10); end
figure; 
sp(1) = subplot(4,1,1); %Overlay Inst Frequency (Hz)
plot(tmp); ylabel('Inst Freq (Hz)'); title('IV069-rec01'); grid on;

tmp = [mat(86:94).instfr_rest_norm];
sp(2) = subplot(4,1,2); %Overlay Normalized Inst Frequency
plot(tmp); ylabel('Normalized Inst Freq'); grid on; ylim([0 1]);

sp(3) = subplot(4,1,3); %Average Normalized Inst Frequency
plot(sum(tmp,2)/size(tmp,2),'k'); ylabel('Normalized Inst Freq'); grid on; ylim([0 1]);

vel_rest = []; b = 38; %Velocity trace
for z = 1:length(beh(b).onRest)
    bout_length = (beh(b).offRest(z)-beh(b).onRest(z))/Fs; %Length of bout
    vel_tmp = beh(b).vel(beh(b).onRest(z):beh(b).offRest(z)); 
    vel_tmp = vel_tmp(1:50:end);
    vel_rest = [vel_rest; vel_tmp(1:bout_length/bin)]; 
end
sp(4) = subplot(4,1,4); plot(vel_rest); ylabel('Velocity (cm/s)'); 
ylim([-1 5]); grid on; linkaxes(sp,'x');


%% IV069-rec01: Correlation
tmp = [mat(86:94).instfr_rest_norm];
%tmp = [mat2(86:94).instfr];
figure; 
h = heatmap(corrcoef(tmp),'CellLabelColor','none'); colormap parula; h.ColorLimits = [-1 1];
title('IV069-rec01: corr-coef InstFR rest spikes')

%% ALL CORRELATION
r = [];
uni = unique({mat.rec});
for x = 1:length(uni)
    idx = find(strcmp({mat.rec},uni{x}));
    if length(idx) < 2; continue; end %If there are no possible pairs to be made, continue to next unique recording ID
    a = corrcoef([mat(idx).instfr_rest_norm]); %Correlation coefficient for inst freq of unit pairs in this unique recording
    a = triu(a,1); %Extract only elements above the main diagonal (auto correlation)
    a(a == 0) = nan; %Remove elements from below the main diagonal
    a = reshape(a, [size(a,1)*size(a,2), 1]); %Reshape matrix into vector
    a = reshape(a(~isnan(a)), [nchoosek(length(idx),2),1]); %Reshape non-NaN elements into a vector of length = nPairs, serves as check
    r = [r; a]; %Save correlation coefficients for unit pairs from this unique recording
end
r_square = reshape(r,[15,11]); %Reshape into square shape

figure;
h = heatmap(r_square,'CellLabelColor','none'); colormap parula; h.ColorLimits = [-1 1];
title(sprintf('corr-coef InstFR rest spikes (n = %d pairs)',length(r)));

figure; violinplot(r);
ylabel('Correlation Coefficient'); ylim([-1 1]); xlim([0.6 1.4]); grid on
title(sprintf('corr-coef InstFR rest spikes (n = %d pairs)',length(r)));
