function [ i_children ] = find_child( i, nRegions,J )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
        [ l,t ] = find_l_t( i, nRegions );
        l_children=l+1;
        t_children=(J*(t-1)+1):(J*t);
        i_children=zeros(length(t_children),1, 'int64');
        for k=1:length(t_children)
            i_children(k) = find_i( l_children,t_children(k), nRegions );
        end

end

