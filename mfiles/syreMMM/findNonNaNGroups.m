function [groups, gmaxIdxs] = findNonNaNGroups(data)
    % Find indices of non-NaN values
    nonNanIdx = find(~isnan(data));
    
     if isempty(nonNanIdx)
        groups = {};
        gmaxIdxs = [];
        return;
    end
    
    % Find boundaries between groups
    gaps = diff(nonNanIdx) > 1;
    
    % Group start indices
    groupStarts = [nonNanIdx(1) nonNanIdx(find(gaps==1) + 1)];
   
    % Group end indices  
    groupEnds = [nonNanIdx(find(gaps==1))  nonNanIdx(end)];
    
    % Create cell array of groups
    numGroups = length(groupStarts);
    groups = cell(1, numGroups);
    
    for i = 1:numGroups
        groups{i} = groupStarts(i):groupEnds(i);
    end
    
    [~, maxIdx] = max(cellfun(@length, groups));
    gmaxIdxs = groups{maxIdx};
end