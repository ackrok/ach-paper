function [idxMax, idxGood, valMag, crossGroups, crossStart] = findIdxPause(fpInputVec, flag, thres, width)
% Photometry pause/peak identification
%   Description: Identify peak/pause events in continuous photometry signal
%   using defined threshold and width
%
% [idxMax, idxStart, valMag, crossGroups, crossStart] = findIdxPause(fpInputVec, 'pause', thres, width)
% [idxMax, idxStart, valMag, crossGroups, crossStart] = findIdxPause(fpInputVec, 'peak', thres, width)
%
% INPUT
%   'fpInputVec' - vector with photometry signal
%   'flag' - character string, either 'pause' or 'peak'
%   'thres' - numerical value, in %dF/F
%       examples: pause -2, peak +4
%   'width'  -  numerical values, in samples
%       examples: pause 7, peak 5
%
% OUTPUT
%   'idxMax' - vector of indices for point of maximum deflection, in samples
%   'idxGood' - vector of indices for crossGroups that satisfy criteria
%   'valMag' - vector of values for magnitude of maximum deflection, in %dF/F
%   'crossGroups' - cell array with vectors of consecutive indices for points
%   crossing threshold
%   'crossStart' - two-column matrix with first and last index of
%   crossGroups
%
% Created by: Anya Krok, December 2021
%

    % thresholds = [1:20]; % thresholds: defined as deviation from new baseline
    % peak_test = cell(1,length(thresholds));
    % mu = []; sigma = [];
    % for t = 1:length(thresholds) % Iterate over possible thresholds
    %     thres = thresholds(t);
    
    switch flag
        case 'pause'; cross = find(fpInputVec < thres); % Find all samples below threshold
        case 'peak'; cross = find(fpInputVec > thres); % Find all samples below threshold
    end
    [crossGroups, crossStart] = consecutive_vec2cell(cross); % Extract grouped crossings
    valMag = nan(length(crossGroups),1); % Initialize vector for amplitudes
    idxGood = []; idxMax = [];
    for c = 1:length(crossGroups) % Iterate over discrete pause/peaks
        fp_vec = [fpInputVec(crossGroups{c})]; % Vector for photometry values that cross during this pause/peak
        if length(crossGroups{c}) > width % Only find amplitudes for peaks that are X samples long (1 sample = 20ms if Fs = 50Hz)
            switch flag
                case 'pause'; [valMag(c),ii] = min(fp_vec); % Pause magnitude
                case 'peak'; [valMag(c),ii] = max(fp_vec); % Peaks magnitude
            end
            idxGood = [idxGood; c]; % Index of crosses that satisfy parameters (width)
            idxMax = [idxMax; crossGroups{c}(ii)]; % Index of max point in cross
        end
    end
    %         mu(t) = nanmean(peak_test{t});
    %         sigma(t) = nanstd(peak_test{t});
    %     end
    %     figure; errorbar(mu, sigma); hold on; errorbar(0,nanmean(m),nanstd(m))
end