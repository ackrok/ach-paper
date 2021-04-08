fPath = '/Users/akrok/Desktop/IV_local/data_temp/';
fName = uigetfile([fPath,'*.mat'],'Multiselect','On');
%%
roi = struct;
for x = 1:length(fName)
    load(fullfile(fPath,fName{x})); 
    roi(x).rec = data.mouserec;
    roi(x).vel = data.final.vel;
    roi(x).time = data.final.time;
    
    if isfield(data.final,'mov'); subst = data.final.mov;
    elseif isfield(data.final,'beh'); subst = data.final.beh; end
    
    try
        roi(x).on = subst.onsets;
        roi(x).off = subst.offsets;
        roi(x).onRest = subst.onsetsRest;
        roi(x).offRest = subst.offsetsRest;
        
    bouts = zeros(1,length(data.final.time)); vec = [];
    for ii = 1:subst.numBouts
        bouts(subst.onsets(ii):subst.offsets(ii)) = 1;
        vec = [vec; mean(data.final.vel(subst.onsets(ii):subst.offsets(ii)))];
    end

    roi(x).boutProp   = length(find(bouts == 1))/length(bouts); %proportion of time spent in a bout
    roi(x).numBouts   = subst.numBouts; %number of bouts in recording
    roi(x).avgBoutDur = subst.avgBoutDuration; %average bout duration
    roi(x).stdBoutDur = subst.stdBoutDuration; %standard dev of bout durations
    roi(x).maxVel     = max(data.final.vel); %maximum velocity achieved
    roi(x).avgVel     = mean(vec); %average of average velocity during bouts
    end
end


%% Subplots of multiple comparisons
full = struct([]); full(1).cohort = beh; full(2).cohort = beh_gat; %full(3).cohort = beh_ldopa; full(4).cohort = beh_chat; 
vec_avg = cell(4,length(full)); 
for y = 1:length(full)
    roi = full(y).cohort;
    for x = 1:length(roi)
        try
              z = 1; vec_avg{z,y}(x) = mean(roi(x).boutProp).*100;  name{z} = 'Percent of Total Time in Bouts (%)';
              z = 2; vec_avg{z,y}(x) = mean(roi(x).avgBoutDur);     name{z} = 'Average Bout Duration (s)';
              z = 3; vec_avg{z,y}(x) = mean(roi(x).maxVel);         name{z} = 'Maximal Velocity Achieved (cm/s)';
              z = 4; vec_avg{z,y}(x) = mean(roi(x).avgVel);         name{z} = 'Average Velocity within Bout (cm/s)';
        end
    end
end;  fprintf('Done! \n');
%% Plotting comparison of multiple cohorts
figure;
for z = 1:4
    subplot(2,2,z);
    m = []; group = [];
    for x = 1:length(full)
        m = [m, vec_avg{z,x}]; 
        group = [group, x.*ones(1,length(vec_avg{z,x}))];
    end
    [~,~,stats] = kruskalwallis(m,group,'off');
    c = multcompare(stats,'display','off');

    violinplot(vec_avg(z,:), {'WT','GAT'}); title(name{z}); ylabel(name{z});
    %title(sprintf('%s',name{z}));
    title(sprintf('%s || wt/gat (%1.3f)',name{z},c(6)));
    %title(sprintf('%s || wt/gat (%1.3f)',name{z},signrank(vec_avg{1},vec_avg{2}));
    %title(sprintf('%s \n WT/6ohda (%1.3f) || WT/chat-ko (%1.3f) || 6ohda/chat-ko (%1.3f)',name{z},c(1,6),c(2,6),c(3,6))); 
    %title(sprintf('%s \n WT/6ohda (%1.3f) || WT/Ldopa (%1.3f) || WT/chatKO (%1.3f) \n 6ohda/Ldopa (%1.3f) || 6ohda/chatKO (%1.3f) || Ldopa/chatKO (%1.3f)',...
        %name{z},c(1,6),c(2,6),c(3,6),c(4,6),c(5,6),c(6,6))); 
end


%%
vec = {3}; roi1 = beh_wt; roi2 = beh_6ohda; roi3 = beh_ldopa; 
vec{:,1} = [roi1.boutProp];
vec{:,2} = [roi2.boutProp];
vec{:,3} = [roi3.boutProp];
 
figure; 
violinplot(vec, {'WT','6ohda','Ldopa'}); 
ylabel('Velocity (cm/s)'); 
%%
m = [vec{:,1},vec{:,2},vec{:,3}]; 
group = [ones(1,length(vec{:,1})),2*ones(1,length(vec{:,2})),3*ones(1,length(vec{:,3}))];
[p,t,stats] = kruskalwallis(m,group,'off');
c = multcompare(stats,'display','off');
title(sprintf('Percent of Total Time in Bouts \n wt/6ohda: p=%1.3f || wt/ldopa: p=%1.3f || 6ohda/ldopa: p =%1.3f',c(1,6),c(2,6),c(3,6))); 


%% Plot velocity
figure; roi = beh_gat; plotme = [1:length(roi)]; 
for x = plotme
    subplot(length(plotme),1,x-plotme(1)+1); 
    plot(roi(x).time,roi(x).vel); title(roi(x).rec,'Interpreter','none')
end
    
