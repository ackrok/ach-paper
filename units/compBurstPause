
%%
full = struct([]); full(1).cohort = cinwt; full(2).cohort = lesion2; full(3).cohort = chat; 
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

%%
full = struct([]); full(1).cohort = cinwt; full(2).cohort = lesion2; %full(3).cohort = ldopa; full(4).cohort = chat;
vec_avg = cell(6,length(full)); 
for y = 1:length(full)
    roi = full(y).cohort;
    for x = 1:length(roi)
        try
            z = 1; vec_avg{z,y}(x) = roi(x).fr; name{z} = 'firing rate (Hz)';
            z = 2; vec_avg{z,y}(x) = roi(x).CV; name{z} = 'CV';
              z = 2; vec_avg{z,y}(x) = mean(roi(x).pause.AllLengths);                                   
                    name{z} = 'pause length (s)';
              z = 3; vec_avg{z,y}(x) = roi(x).pprop*100; %(sum(roi(x).pause.AllLengths)/roi(x).st(end))*100;              
                    name{z} = 'prop of time in pauses (%)';
              z = 4; vec_avg{z,y}(x) = mean([roi(x).burst.Windows(:,2) - roi(x).burst.Windows(:,1)]);   
                    name{z} = 'burst length (s)';
              z = 4; vec_avg{z,y}(x) = mean(roi(x).burst.bHz);                                          
                    name{z} = 'inter-burst frequency (Hz)';
              z = 6; vec_avg{z,y}(x) = roi(x).bprop*100; %(sum(roi(x).burst.Windows(:,2)-roi(x).burst.Windows(:,1))/roi(x).st(end))*100;   
                    name{z} = 'prop of time in bursts (%)';
              z = 5; vec_avg{z,y}(x) = mean(roi(x).burst.NumSpikes);                                    
                    name{z} = '# spikes in each burst';
        end
    end
    vec_avg{3,y}(find(vec_avg{3,y} == 0)) = NaN; vec_avg{6,y}(find(vec_avg{6,y} == 0)) = NaN;
end; fprintf('Done! \n');
%%
figure;
for z = 1:6
    subplot(2,3,z);
%     m = []; group = [];
%     for x = 1:length(full)
%         m = [m, vec_avg{z,x}]; 
%         group = [group, x.*ones(1,length(vec_avg{z,x}))];
%     end
%     [~,~,stats] = kruskalwallis(m,group,'off');
%     c = multcompare(stats,'display','off');
    
    violinplot(vec_avg(z,:), {'WT','6ohda'}); title(name{z}); ylabel(name{z});
    title(sprintf('%s || (%1.3f)',name{z},ranksum(vec_avg{z,1},vec_avg{z,2})));
    %title(sprintf('%s \n WT/6ohda (%1.3f) || WT/Ldopa (%1.3f) || 6ohda/Ldopa (%1.3f)',name{z},c(1,6),c(2,6),c(3,6))); 
    %title(sprintf('%s \n WT/6ohda (%1.3f) || WT/Ldopa (%1.3f) || WT/chatKO (%1.3f) \n 6ohda/Ldopa (%1.3f) || 6ohda/chatKO (%1.3f) || Ldopa/chatKO (%1.3f)',...
    %    name{z},c(1,6),c(2,6),c(3,6),c(4,6),c(5,6),c(6,6))); 
end

