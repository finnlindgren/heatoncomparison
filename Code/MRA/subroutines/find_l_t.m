function [ l,t ] = find_l_t( i, nRegions )
% find_l_t Finds l (level) and t (tile number)
%   This function finds the level and tile number given the continuous
%   index i as an input
switch  i
    case 1
        l=1; t=1; % Special case for zeroth region
    otherwise
        cumRegions=cumsum(nRegions);
        ind_smaller=find(cumRegions<i,1,'last');
        l=ind_smaller+1;
        t=i-cumRegions(ind_smaller);
end
end

