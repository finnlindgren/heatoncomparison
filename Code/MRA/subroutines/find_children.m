function [ i_children ] = find_children( i,nRegions,J  )
% find_children Find immediate children
%   Find children at resolution right below current

cumRegions=cumsum(nRegions);

if length(cumRegions)==1  % Check that there is more than one region total
    error('myfuns:find_children:OnlyOneRegion', ...
        'no chidren if only one region');
end

if i>cumRegions(end-1)
    i_children=[];  % No children on finest resolution level
else
    %L=length(cumRegions);% Find out total number of Levels
    switch  i
        case 1
            % Special case for zeroth region
            i_children=2:cumRegions(2); % 
        otherwise% Special case for second finest level
            i_children=find_child( i,nRegions,J  ); % Only immediate children
%         otherwise
%             [ l,~ ] = find_l_t( i, nRegions );
%             % Create a cell structure where elements are children on each level
%             temp_i_children=cell(L-l,1);% There will be children on L-l levels
%             
%             
%             temp_i_children{1}=find_child( i,nRegions,J  ); % immediate children
%             
%             % Loop through lower levels
%             for k=2:L-l
%                 for j=1:length(temp_i_children{k-1})
%                     temp=find_child( temp_i_children{k-1}(j),nRegions,J  );
%                     i_children{k}=[temp_i_children{k};temp];
%                 end
%             end
%             i_children=[]
%             for k=1:L-l
%             i_children=[i_children;temp_i_children{k}];
%             end
    end
end
end


%% NOTES:
% I originally misunderstood what this function did in Matthias' code.
% find_offspring finds the entire set of connected higher resolution
% regions.


