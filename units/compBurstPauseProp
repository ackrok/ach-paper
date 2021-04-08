%% Identifying bursts during MOV vs. REST
lbl = {'CIN pause'};
time = [-1:1/50:1]; 

mat_b = cell(length(sub),2); 
%mat_c_mov = []; mat_c_rest = []; mat_c_shuff = [];
for x = 1:length(sub)
    if lbl{1} == 'CIN burst'; events = sub(x).burst.Windows(:,1);
    elseif lbl{1} == 'CIN pause'; events = sub(x).pause.Windows(:,1); 
    else; continue;
    end
    evMov_idx = []; evRest_idx = [];
    f = find(strcmp({beh.rec},sub(x).rec));
    for a = 1:length(beh(f).on)
        logicalIndexes = events <= beh(f).off(a) & events >= beh(f).on(a);
        evMov_idx = [evMov_idx; find(logicalIndexes)];
    end
    for a = 1:length(beh(f).onRest)
        logicalIndexes = events <= beh(f).offRest(a) & events >= beh(f).onRest(a);
        evRest_idx = [evRest_idx; find(logicalIndexes)];
    end
    mat_b{x,1} = evMov_idx; mat_b{x,2} = evRest_idx;
end 
%% comparing BURST/PAUSE properties
vec_full = cell(1,2); vec_avg = cell(1,2); 
for x = 1:length(sub)
    n = {2}; n{1} = mat_a_mov{x}; n{2} = mat_a_rest{x}; 
    for y = 1:length(n)
        %a = [sub(x).burst.Windows(n{y},2) - sub(x).burst.Windows(n{y},1)]; name = 'burst length (s)';
        %a = [sub(x).burst.NumSpikes(n{y})]; name = '# spikes in each burst';
        
        vec_full{:,y} = [vec_full{:,y}; a];
        vec_avg{:,y}(x) = mean(a);
    end
end

m = [vec_avg{:,1},vec_avg{:,2}]; 
group = [ones(1,length(vec_avg{:,1})),2*ones(1,length(vec_avg{:,2}))];
[~,~,stats] = kruskalwallis(m,group,'off');
c = multcompare(stats,'display','off');
 
figure; violinplot(vec_avg, {'MOV','REST'}); title(name); ylabel(name);
title(sprintf('%s || mov/rest (p = %1.3f)',name,c(6))); 


%% comparing BURST/PAUSE properties
full = struct([]); full(1).cohort = cinwt; full(2).cohort = lesion2; full(3).cohort = ldopa; 
vec_avg = {length(full)}; vec_all = {length(full)};
for y = 1:length(full)
    roi = full(y).cohort; vec_avg{:,y} = []; 
    for x = 1:length(roi)
        try
%             vec_full{:,y} = [vec_full{:,y}; roi(x).pause.AllLengths]; name = 'pause length (s)';
%             vec_full{:,y} = [vec_full{:,y}; [roi(x).burst.Windows(:,2) - roi(x).burst.Windows(:,1)]]; name = 'burst length (s)';
%             vec_full{:,y} = [vec_full{:,y}; roi(x).burst.IBF]; name = 'inter-burst frequency (Hz)';
%             vec_avg{:,y}(x) = mean(roi(x).pause.IPF); vec_full{:,y} = [vec_full{:,y}; roi(x).pause.IPF]; name = 'inter-pause frequency (Hz)';
%             vec_full{:,y} = [vec_full{:,y}; roi(x).burst.NumSpikes]; name = '# spikes in each burst';
%             vec_full{:,y} = [vec_full{:,y}; roi(x).pause.NumSpikes]; name = '# spikes in each pause';
%             vec_avg{:,y}(x) = mean(roi(x).pause.AllLengths); name = 'pause length (s)';
%             vec_avg{:,y}(x) = mean([roi(x).burst.Windows(:,2) - roi(x).burst.Windows(:,1)]); name = 'burst length (s)';
%             vec_avg{:,y}(x) = length(roi(x).pause.PausingSpikes)/length(roi(x).st)*100; name = 'proportion of spikes in pauses (%)';
%             vec_avg{:,y}(x) = length(roi(x).burst.BurstingSpikes)/length(roi(x).st)*100; name = 'proportion of spikes in bursts (%)';
        end
    end
end
%figure; violinplot(vec); xticklabels({'WT','6ohda','Ldopa'}); title(name); ylabel(name);
 
m = [vec_avg{:,1},vec_avg{:,2},vec_avg{:,3}]; 
group = [ones(1,length(vec_avg{:,1})),2*ones(1,length(vec_avg{:,2})),3*ones(1,length(vec_avg{:,3}))];
[~,~,stats] = kruskalwallis(m,group,'off');
c = multcompare(stats,'display','off');
 
figure; violinplot(vec_avg, {'WT','6ohda','Ldopa'}); title(name); ylabel(name);
title(sprintf('%s \n WT/6ohda (%1.3f) || WT/Ldopa (%1.3f) || 6ohda/Ldopa (%1.3f)',name,c(1,6),c(2,6),c(3,6))); 

