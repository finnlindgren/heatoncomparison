function [ X,Y ] = create_knots( xmin,xmax,nx,ymin,ymax,ny, offsetpercentage )
% create_knots for each partition 
%   Determines knot locations for each partition

offsetx=(xmax-xmin)*offsetpercentage/100;
offsety=(ymax-ymin)*offsetpercentage/100;
[X,Y]=meshgrid(linspace(xmin+offsetx,xmax-offsetx,nx),linspace(ymin+offsety,ymax-offsety,ny));
end

