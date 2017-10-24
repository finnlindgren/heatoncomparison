function [ M,nRegiongsFinest, M_integer, nRegions, totalRegions ] = find_mrequired( n,r,J )
% find_mrequired Finds number of levels
%   Finds the number of levels for a given number of observations, number
%   of knots per region (r) and number of partitions per region (J). The
%   idea is to work backwards from the finest level by making the average
%   number of observations per region similar to the number of knots.

% Find number of regions required at finest level

nRegiongsFinest=n/r;

% Solve for M which is m at the finest level using the J^M is the number of
% regions at the finest level

M=log(nRegiongsFinest)/log(J);
M_integer=ceil(M);

% Comment: Matlab doesn't have a direct command to solve a logarithm for
% any base, so we use log(x)/log(b) in place of log_base_b(x)

mLevels=0:M_integer; % vector of levels
nRegions=J.^mLevels; % Regions (partitions) by level
totalRegions=sum(nRegions);
end

