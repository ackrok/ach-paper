function stSub = extractEventST_2(st,on,off)
% Extract spike times within a specified event time range, artificially
% removing gaps between spikes generated by separation in event times
%
% stSub = extractEventST_2(st,on,off)
%
% Description: This function will take inputs of spike times and the start
% and end times for an event, and return a vector of spike times that occur
% within the event time range.
%
% INPUT
%   'st' - vector of spike times, in units that match that of event times
%   'on' - vector of start times for event(s)
%   'off' - vector of end times for event(s)
%
% OUTPUT
%   'stSub' - vector of subset of spike times occuring within event range,
%               artificially remove gaps between spikes
%
% Anya Krok, January 2022
% adapted from extractEventST: adjust spike times for length of time between
% events so that there are no large gaps between spikes, as is generated
% when using extractEventST with flag 1 for concatenating spikes
        

stSub = []; % Initialize output vector
t_adj = []; % Initialize time adjustment value
for z = 1:length(on)
    if z == 1
        t_adj = on(1); % Start t_adj value at the start time of 1st event
    else
        t_adj = t_adj + on(z) - off(z-1);
    end
    st_original = extractEventST(st, on(z), off(z), 1); % Extract spikes during this event times, concatenating spikes
    st_idx = ismember(st,st_original); % Find index of spikes that are within this event window
    st_adj = st(st_idx) - t_adj; % Adjust spike times
    stSub = [stSub; st_adj]; % Concatenate into output vector
end
end