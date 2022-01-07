function [cross_cell, cross_pts] = consecutive_vec2cell(cross)
% Purpose: convert vector with regions of consecutive numbers into a cell
% array where each cell contains of vector of single sequence of
% consecutive values.
%
% INPUTS
%   'cross' - vector column of numerical values
% OUPUTS
%   'cross_cell' - cell array, where each cell contains a single sequence
%       of consecutive values
%   'cross_pts' - vector with 2 columns, where 1st column contains start
%   value of a single sequence and 2nd column is the end value
%
% Anya Krok, December 2021

cross = cross(:); % check that vector is a column
idlFwd = [diff(cross) == 1;0]; % elements are consecutive if they the next element is one larger than they are
idlBwd = [0;flip(diff(flip(cross)) == -1)]; % elements are also consecutive if the one before them is one less than they are
idl = idlFwd | idlBwd; % elements are consecutive if either criteria is met
iStart = diff([0;idlFwd])==1; % mark where beginning of each group of consecutive elements occurs, these are jumps from 0 to 1
groupCount = cumsum(iStart); % count up jumps to form group numbers
group = groupCount.*idl ; % assign group number, but only to ones that belong to a group (the consecutive values)
    % (for non-members groupCount.*idl is zero so these are assigned to group zero)
numGroups = max(group); 
cross_cell = cell(numGroups,1); % loop through to assign to cell array by group
cross_pts = nan(numGroups,2); % start and stop point of each crossing
for k = 1:numGroups 
    cross_cell{k} = cross(k == group);
    cross_pts(k,1) = cross_cell{k}(1);
    cross_pts(k,2) = cross_cell{k}(end);
end
    
end