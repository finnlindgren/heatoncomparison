function [ l_parent,t_parent,i_parent ] = find_parent( i, nRegions,J )
% find_parent Finds descriptors of parent
%   This function finds the level, tile number and continuous index of the
%   parent of the region with index i

switch  i
    case 1
        % Special case for zeroth region
        t_parent=NaN;
        l_parent=NaN;
        i_parent=NaN;
    otherwise
        cumRegions=cumsum(nRegions);
        ind_smaller=find(cumRegions<i,1,'last');
        l=ind_smaller+1;
        t=i-cumRegions(ind_smaller);
         % figure out ID of parent
        t_parent=ceil(t/J);
        l_parent=l-1;
        i_parent=sum(nRegions(1:l-2))+t_parent;
end


end