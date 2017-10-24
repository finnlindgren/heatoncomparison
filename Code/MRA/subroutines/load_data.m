function [ data, regressionModel, domainBoundaries, predVec, theta, varEps ] = load_data(dataType)
%load_data Loads data for Heaton comparison project.
%   Data files are loaded as a function of a variable called dataType

switch dataType
    
    case 'satellite'
        load('./satelliteData.mat')
        % Values of parameters of covariance function
        theta=[5.57,0.12]; varEps=0.01;
        
    case 'simulated'
        load('./simulatedData.mat')
        % Values of parameters of covariance function
        theta=[8.13,0.72]; varEps=0.1;
        
    otherwise
        error('Error. Specified dataType is not a valid data set.');
end
% Determine the boundaries of the domain spanded by the data.
xmin0 = min(lon);
xmax0 = max(lon);
ymin0 = min(lat);
ymax0 = max(lat);
domainBoundaries=[xmin0, xmax0, ymin0, ymax0];

% Find observation locations.
logicalInd = ~isnan(obs);

% Declare all locations as prediction locations.
xPred = lon;
yPred = lat;
predVec = [xPred(:),yPred(:)];

% Assign lon, lat and observations to data matrix.
data(:,1) = lon(logicalInd);
data(:,2) = lat(logicalInd);
data(:,4) = obs(logicalInd);

% Detrend data.
regressionModel=fitlm(data(:,1:2),data(:,4),'linear');
residuals=table2array(regressionModel.Residuals(:,1));
data(:,3)=residuals;

end

