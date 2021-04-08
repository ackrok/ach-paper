%% Align FP to REWARDED TRIALS
mat = struct; % Initialize output structure
gen = struct; gen.bin = 0.25; gen.window = [0 0.25]; nShuff = 10; %CHANGE: window for PETH 
h = waitbar(0, 'PETH: CIN spikes to acceleration peaks');
for x = 1:length(beh)
lick = beh(x).lick(:)/beh(x).Fs; % Lick times, in seconds
stim = beh(x).reward/beh(x).Fs; % Cue or Reward times, in second

peth = getClusterPETH(lick,stim,gen); % PETH: licks to stimulus
lick_yes = find(peth.cts{1} > 0); % Reward trials where mouse licks
if isempty(lick_yes); continue; end

fp = [beh(x).vel(1); diff(movmean(beh(x).vel,10))]; tmpAlign = cell(1);

[sta, time, staZ] = getSTA(fp, stim(lick_yes), beh(x).Fs, [-2 2]);
tmpAlign{1} = staZ;


mat(x).rec = beh(x).rec;
mat(x).stimLick = stim(lick_yes);
mat(x).align = tmpAlign;
waitbar(x/length(beh),h);
end; fprintf('Done! \n'); close(h);


figure; clr = {'k'};
for x = 1:length(mat)
sp(x) = subplot(3,5,x);
if isempty(mat(x).stimLick); continue; end

shadederrbar(time, nanmean(mat(x).align{y},2), SEM(mat(x).align{y},2), clr{y}); hold on

xlabel('Latency from Reward (s)'); ylabel('FP (z-score)'); grid on
title(sprintf('%s-DLS (%d trials)',mat(x).rec,length(mat(x).stimLick)));
end; linkaxes(sp,'x');

align2 = cell(1);
for x = 1:15
if isempty(mat(x).stimLick); align2{y} = [align2{y},nan(201,1)]; continue; end
align2{1} = [align2{1}, nanmean(mat(x).align{y},2)]; end
figure;
shadederrbar(time, nanmean(align2{1},2), SEM(align2{1},2), 'k'); hold on
xlabel('Latency from Reward (s)'); ylabel('FP (z-score)');
title('ACh, DA aligned to Reward (lick within 250ms)'); xlim([0 1]);
