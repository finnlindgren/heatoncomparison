function [ i ] = find_i( l,t, nRegions )
% find_i Finds continuous index i
%   This function finds the continuous index i given the l (level) and t
%   (tile number) as inputs

        i=sum(nRegions(1:l-1))+t;

end