function  R = co( locs1,locs2,theta )
% CO Covariance function
%   Generic function for testing


  h=pdist2(locs1,locs2)/theta(2);
  %R= theta(1)*(1 + sqrt(3).*h).*exp(-sqrt(3).*h); % nu=3/2
  R= theta(1)*exp(-h); % nu=1/2 exponential

end

