gen = struct; gen.bin = 0.01; gen.window = [-0.015 0.015]; 

x = 2;
idx = find(strcmp({sub.rec},sub(x).rec)); idx(idx == x) = [];
st = [sub(x).st]; % Extract spike times of this unit
    st_other = {sub(idx).st}; % Extract spike times of all other units from this recording
    fr_other = [sub(idx).fr];
Fs = 50;

%%
for aa = [2:5]
    
gen = struct; gen.bin = aa/100; gen.window = [-gen.bin/2 gen.bin/2]; 
st_rest = extractEventST(st, beh(1).onRest/Fs, beh(1).offRest/Fs, 0); % Event times during movement
peth_rest = getClusterPETH(st_other, st_rest, gen); % PETH: spike times aligned to spike times
t_0 = 1; cts0_rest = []; ll_pair = [];
for y = 1:length(idx)
    tmp = peth_rest.cts{y}(t_0,:);
    cts0_rest(y,:) = tmp;
    ll_pair(y) = length(find(tmp))/length(tmp);
end  
b = sum(cts0_rest,1);
ll = length(find(b >= 1))/length(b); % Probability of reference unit and at least 1 other unit active during same bin

nShuff = 10; ll_p = []; ll_pair_p = []; % CHANGE - number of times to generate poisson spike trains
h = waitbar(0, 'synchSpikes - poisson Spike trains');
for z = 1:nShuff
    %fr_rest = 1/mean(diff(extractEventST(st, beh(1).onRest/Fs, beh(1).offRest/Fs, 1)));
    %tmp = poissonSpikeGen(fr_rest, beh(1).time(end), 1); 
    tmp = shuffleST(st, 1); 
    st_rest_p = tmp{1};
    st_rest_other_p = {};
    for y = 1:length(idx)
        %fr_rest = 1/mean(diff(extractEventST(st_other{y}, beh(1).onRest/Fs, beh(1).offRest/Fs, 1)));
        %tmp = poissonSpikeGen(fr_rest, beh(1).time(end), 1); 
        tmp = shuffleST(st_other{y}, 1); 
        st_rest_other_p{y} = tmp{1};
    end
    peth_rest_p = getClusterPETH(st_rest_other_p, st_rest_p, gen); % PETH: spike times aligned to spike times

    t_0 = 1; cts0_rest_p = [];
    for y = 1:length(idx)
        tmp = peth_rest_p.cts{y}(t_0,:);
        cts0_rest_p(y,:) = tmp;
        ll_pair_p(z,y) = length(find(tmp))/length(tmp);
    end  
    % sta_rest_p = getSTA(fp, st_rest_p, Fs, [-1, 1]);
    
    b = sum(cts0_rest_p,1);
    ll_p(z) = length(find(b >= 1))/length(b); % Probability of poisson unit and at least 1 other poisson unit active during same bin
    waitbar(z/nShuff, h);
end
close(h);

test(aa).rec = mat(2).rec; test(aa).n = mat(2).n; test(aa).m = mat(2).m; 
test(aa).bin = gen.bin;
test(aa).ll = ll; test(aa).ll_p = ll_p; 
test(aa).ll_pair = ll_pair; test(aa).ll_pair_p = ll_pair_p; 
test(aa).cts0_rest = cts0_rest;
test(aa).cts0_rest_p = cts0_rest_p;
end
    
%% LIKELIHOOD of 2+ UNITS CO-ACTIVE
figure;
for aa = 1:length(test)
    sp(aa) = subplot(1,3,aa); hold on
    scatter([1:4],test(aa).ll_pair,'b','filled','DisplayName','data');
    for y = 1:nShuff
        scatter([1:4],test(aa).ll_pair_p(y,:),'k','DisplayName','poisson');
    end
    title(sprintf('%d ms bins',aa*10)); xticks([1:4]); xlim([0.5 4.5])
    ylabel('Proportion of Spikes where 2+ Units are Spiking');
end
linkaxes(sp,'y');
