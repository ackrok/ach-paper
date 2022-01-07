fPath = '/Users/akrok/Desktop/IV_local/data_temp/';
[fName,fPath] = uigetfile([fPath,'*.mat'],'Multiselect','On');

%% Extract Behavioral Data
roi = struct;
for x = 1:length(fName)
    load(fullfile(fPath,fName{x})); 
    a = length(roi) + 1;
    roi(a).rec = data.mouserec; %data.mouserec;
    roi(a).vel = data.final.vel;
    roi(a).time = data.final.time;
    if isfield(data.final,'FP')
        roi(a).FP = data.final.FP;
        roi(a).FPnames = data.final.FPnames;
    end
    if isfield(data.final,'mov'); subst = data.final.mov;
    elseif isfield(data.final,'beh'); subst = data.final.beh; end
    try
        roi(a).on = subst.onsets;
        roi(a).off = subst.offsets;
        roi(a).onRest = subst.onsetsRest;
        roi(a).offRest = subst.offsetsRest;  
    bouts = zeros(1,length(data.final.time)); vec = [];
    for ii = 1:subst.numBouts
        bouts(subst.onsets(ii):subst.offsets(ii)) = 1;
        vec = [vec; mean(data.final.vel(subst.onsets(ii):subst.offsets(ii)))];
    end
    roi(a).boutProp   = length(find(bouts == 1))/length(bouts); %proportion of time spent in a bout
    roi(a).numBouts   = subst.numBouts; %number of bouts in recording
    roi(a).avgBoutDur = subst.avgBoutDuration; %average bout duration
    roi(a).stdBoutDur = subst.stdBoutDuration; %standard dev of bout durations
    roi(a).maxVel     = max(data.final.vel); %maximum velocity achieved
    roi(a).avgVel     = mean(vec); %average of average velocity during bouts
    end
end
if isempty(roi(1).time); roi(1) = []; end
beh = roi;

%% Extract CIN data
roi = struct;
for x = 1:length(fName)
    load(fullfile(fPath,fName{x})); 
    if ~isfield(data,'lbl'); continue; end
    if isfield(data.lbl,'cin')
    for y = 1:length(data.lbl.cin)
        a = length(roi) + 1;
        roi(a).rec = data.mouserec;
        roi(a).n    = data.lbl.cin(y); roi(a).label = 2;
        roi(a).st   = data.clusters(roi(a).n).spikeTimes;
        roi(a).wf   = data.WF(roi(a).n).topMu;
        roi(a).fr   = data.spk(roi(a).n).fr;
        roi(a).pISI2 = data.spk(roi(a).n).pISI2;
        roi(a).CV   = data.spk(roi(a).n).CV;
        roi(a).halfW = data.WF(roi(a).n).halfWidth;
        roi(a).peakL = data.WF(roi(a).n).peakLatency;
        roi(a).coor  = data.clusters(roi(a).n).coordinates;
    end; end
    %%
    if isfield(data.lbl,'spn')
    for y = 1:length(data.lbl.spn)
        a = length(roi) + 1;
        roi(a).rec = data.mouserec;
        roi(a).n    = data.lbl.spn(y); roi(a).label = 1;
        roi(a).st   = data.clusters(roi(a).n).spikeTimes;
        roi(a).wf   = data.WF(roi(a).n).topMu;
        roi(a).fr   = data.spk(roi(a).n).fr;
        roi(a).pISI2 = data.spk(roi(a).n).pISI2;
        roi(a).CV   = data.spk(roi(a).n).CV;
        roi(a).halfW = data.WF(roi(a).n).halfWidth;
        roi(a).peakL = data.WF(roi(a).n).peakLatency;
        roi(a).burst = data.spk(roi(a).n).burst;
        roi(a).coor  = data.clusters(roi(a).n).coordinates;
    end; end
end
if isempty(roi(1).n); roi(1) = []; end
units = roi;

%% Burst/Pause detection for CINs
h = waitbar(0,'burst/pause detection CINs');
for x = 1:length(roi)
  try
    [Bursts, Pauses] = RGSDetect(roi(x).st, 2, [-3.0:0.005:1.5], 0.05, 0.05);
    try roi(x).burst = Bursts; catch roi(x).burst = []; end
    try roi(x).pause = Pauses; catch roi(x).pause = []; end
    close all
  end
  waitbar(x/length(roi),h);
end; close(h);

for y = 1:length(roi)
    if isempty(roi(y).burst); continue; end
    IBF = []; a = 0; 
    for x = 1:length(roi(y).burst.NumSpikes)
        a = [a(end)+1 : a(end)+roi(y).burst.NumSpikes(x)];
        burstSpikes = roi(y).burst.BurstingSpikes(a);
        IBF(x) = 1/mean(diff(burstSpikes));
    end
    roi(y).burst.bHz = IBF(:);
end
