function [ knots, partitions,nRegions, outputdata, predlocs ] = build_structure( M,J,r,domainBoundaries, offsetpercentage, varargin )
% build_structure Build nested structure
%   Build partitioning structure and knots as a function of levels

% M= Total number of levels (l)
% J= Number of partitions for each region
% r= number of knots per partition

% Check number of optional input arguments
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:build_structure:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
% Optional argument that can be passed is data
% [data]=varargin{1};

optargs(1:numvarargs) = varargin;
[data, predvec]=optargs{:};
%% Calculate quantities of interest

mLevels=0:M-1; % vector of levels
nRegions=J.^mLevels; % Regions (partitions) by level
totalRegions=sum(nRegions);

% Calculate number of knots in each direction
if isinteger(sqrt(r))  % Simple way to assign knots for J =30 Generalize!
nx0=sqrt(r);nx=sqrt(r); % Could be done differently
ny0=sqrt(r);ny=sqrt(r); % Could be done differently
else
nx0=ceil(sqrt(r));nx=ceil(sqrt(r)); % Could be done differently
ny0=r/nx0;ny=r/nx; % Could be done differently    
    
end

%% Create knots for partitions

%Pre-allocate 
knots=cell(totalRegions,1);
outputdata=cell(totalRegions,1);
partitions=cell(totalRegions,1);
%indices=cell(totalRegions,1); % legacy from Julia setup
%indres=cell(length(mLevels),1);% legacy from Julia setup

% Construct zeroth level
xmin0=domainBoundaries(1);xmax0=domainBoundaries(2);
ymin0=domainBoundaries(3);ymax0=domainBoundaries(4);
[knots_x,knots_y]=create_knots(xmin0,xmax0, nx0, ymin0,ymax0,ny0,offsetpercentage);
knots{1,1}=[knots_x(:),knots_y(:)];

[ xmin, xmax,ymin, ymax ] = create_partition( xmin0,xmax0, ymin0, ymax0,J );
partitions{1,1}=[ xmin, xmax,ymin, ymax ];
% tempx=linspace(xmin0,xmax0, J/2+1);
% xmin=repmat(temp(1:end-1), J/2,1);xmin=xmin(:);
% xmax=repmat(temp(2:end),J/2,1);xmax=xmax(:);
% 
% tempy=linspace(ymin0,ymax0, J/2+1);
% ymin=repmat(temp(1:end-1),1, J/2);ymin=ymin';
% ymax=repmat(temp(2:end),1, J/2);ymax=ymax';
%i_Regions=2;
if numvarargs==2
    finestknotlevel=M-1;
else
    finestknotlevel=M;
end

for l= 2:finestknotlevel
    disp(['Building Level ',num2str(l),' starting']);
    for t=1:J:nRegions(l)
        
        % figure out ID of parent
         [ i ] = find_i( l,t, nRegions );
         [ ~,~,i_parent ] = find_parent( i, nRegions,J );
        % get partition coordinates of parent
        xmin=partitions{i_parent,1}(:,1);
        xmax=partitions{i_parent,1}(:,2);
        ymin=partitions{i_parent,1}(:,3);
        ymax=partitions{i_parent,1}(:,4);
        for j=1:J
            i_current=sum(nRegions(1:l-1))+t+j-1;
            [knots_x,knots_y]=create_knots(xmin(j),xmax(j), nx, ymin(j),ymax(j),ny,offsetpercentage);
            [ xmin_temp, xmax_temp,ymin_temp, ymax_temp ] = create_partition( xmin(j),xmax(j), ymin(j), ymax(j),J );
            knots{i_current,1}=[knots_x(:),knots_y(:)];
            partitions{i_current,1}=[ xmin_temp, xmax_temp,ymin_temp, ymax_temp ];
            %i_Regions=i_Regions+1;
        end
    end
end
% Special construct to find knots for finest resolution region
if numvarargs==2  % The last level is built using data as knots
    l=M;
    disp(['Building finest resolution Level ',num2str(l),' starting']);
    n_tilesfinestlevel=nRegions(end)/J;
    for t=1:J:nRegions(l)
        
        % figure out ID of parent
        [ i ] = find_i( l,t, nRegions );
        [ ~,~,i_parent ] = find_parent( i, nRegions,J );
        % get partition coordinates of parent
        xmin=partitions{i_parent,1}(:,1);
        xmax=partitions{i_parent,1}(:,2);
        ymin=partitions{i_parent,1}(:,3);
        ymax=partitions{i_parent,1}(:,4);
        for j=1:J
            i_current=sum(nRegions(1:l-1))+t+j-1;
            ind=find(data(:,1)>xmin(j)& data(:,1)<xmax(j) & data(:,2)>ymin(j)& data(:,2)<ymax(j));
            knots_x=data(ind,1);knots_y=data(ind,2);
            %[knots_x,knots_y]=create_knots(xmin(j),xmax(j), nx, ymin(j),ymax(j),ny,offsetpercentage);
            %[ xmin_temp, xmax_temp,ymin_temp, ymax_temp ] = create_partition( xmin(j),xmax(j), ymin(j), ymax(j),J );
            knots{i_current,1}=[knots_x(:),knots_y(:)];
            outputdata{i_current,1}=data(ind,3);% This is to only pass the data, not the location to MRA
            %partitions{i_current,1}=[ xmin_temp, xmax_temp,ymin_temp, ymax_temp ];
            %i_Regions=i_Regions+1;
            data(ind,:)=[];% Eliminate the data that has already been assigned to a region, speeds up subsequent searching
            
            %% Vinay's addition to partition the prediction locations
            if ~isnan(predvec)  % If predicting
            predInd = find(predvec(:,1)>=xmin(j) & predvec(:,1)<=xmax(j) & predvec(:,2)>=ymin(j)& predvec(:,2)<=ymax(j));
            predlocs{i_current,1} = predvec(predInd,:);
            else
                predlocs = NaN;
            end
            
        end
        disp(['Tile ',num2str(floor(t/4)+1),' of ', num2str(n_tilesfinestlevel), ' completed']);
    end
else
end
end

