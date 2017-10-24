function [ i_ancestry ] = find_ancestry( i, nRegions,J )
% find_parent Finds descriptors of parent
%   This function finds the level, tile number and continuous index of the
%   parent of the region with index i

switch  i
    case 1
        % Special case for zeroth region
        %t_ancestry=NaN;
        %l_ancestry=NaN;
        i_ancestry=[];
    otherwise
        cumRegions=cumsum(nRegions);
        ind_smaller=find(cumRegions<i,1,'last');
        %l=ind_smaller+1;
        %t=i-cumRegions(ind_smaller);
        % pre-allocate vector
        i_ancestry=zeros(ind_smaller,1, 'int64'); % number of ancestry members is one smaller than level
         % fill ancestry by looping through parents parents
         %while i_parent > 1
         for k=1:ind_smaller
             [ ~,~,i_parent ] = find_parent( i, nRegions,J );
             i_ancestry(k,1)=i_parent;
             i=i_parent;
         end
         i_ancestry=flip(i_ancestry); % From coarsest to finest
      
end


end
