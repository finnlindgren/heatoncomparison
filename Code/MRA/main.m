%% MRA main script for Heaton project
% Executes the MRA for the simulated and satellite data. Calculation
% options include prediction, parameter optimization and likelihood only
% evaluation.

%% Preliminary settings
clear all; addpath('subroutines');

%% Options to specify [Options to change here.]
% Choose dataType 'satellite' or 'simulated'
dataType='simulated'; 
% Choose calculationType 'prediction', 'optimize', 'likelihood'
calculationType='optimize'; 
plotting=1; % Only use '1' for testing.
%% Setup options [Do NOT change.]
M=9; % Total number of levels
J=2; % Number of partitions by level
r=64; % number of knots per partition
offsetPercentage = 0.01; % offset percentage from partition boundaries

%% Load data
[data, regressionModel, domainBoundaries, predVec, theta, varEps] = load_data(dataType);

%% Build hierarchical grid structure
[knots, ~,nRegions, outputdata, predlocs] = build_structure( M,J,r,domainBoundaries, offsetPercentage,data(:,1:3),predVec );

%% Potential optimization
switch calculationType
    case {'optimize', 'prediction'}
        switch calculationType
            case 'optimize'
                fun=@(thetaOpt)MRA( [thetaOpt(1) thetaOpt(2)], outputdata, knots, M, J, nRegions,thetaOpt(3) );
                % Dummy values required by optimization routine
                A = [];b = [];Aeq = [];    beq = [];
                % Limits and initial values for parameter search
                lb = [0,0,0]; ub = [10,1,5]; x0 = [5,0.3,0.1];
                
                tic; x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
                elapsedTimeOptimization=toc;
                
                % Assign values from optimization to theta and varEps
                theta=[x(1) x(2)];
                varEps=x(3);
        end
        % Prediction
        tic
        [sum_loglik,preds] = MRA( theta, outputdata, knots, M, J, nRegions,varEps, predlocs);
        elapsedTimePrediction=toc;
        % Reformat data for plotting
        pred = cell2mat(preds);predloc=cell2mat(predlocs);predVariance=pred(:,2);
        % Add the prediction from the regression
        predRegression=predict(regressionModel,predloc);
        predMean=pred(:,1)+predRegression;
        save(['ResultsMra',dataType],'predloc','predMean','predVariance');
    case 'likelihood'
        % Prediction
        tic
        [sum_loglik] = MRA( theta, outputdata, knots, M, J, nRegions,varEps);
        elapsedTimeLikelihood=toc;
    otherwise
        warning('Unexpected calculationType. Code is not executed.')
end

%% plots
if plotting==1
figure
scatter(data(:,1),data(:,2),5,data(:,4));
colormap(parula)
colorbar
[cmin,cmax] = caxis;
caxis([cmin,cmax])
title('Ground truth')
xlabel('x locations')
ylabel('y locations')
figure
scatter(predloc(:,1),predloc(:,2),5,predMean);
colormap(parula)
colorbar
caxis([cmin,cmax])
title('Predicted values')
xlabel('x locations')
ylabel('y locations')
figure
scatter(predloc(:,1),predloc(:,2),5,predVariance);
colormap(flip(autumn))
colorbar
title('Prediction variance')
xlabel('x locations')
ylabel('y locations')
end

%% Notes:

% calculationType specifies what is calculated. 
% Option 'prediction' uses a default values for the parameters and just conducts 
% the prediction. 
% Option 'optimize' optimizes over the range, variance and measurement error 
% and then predicts using the values obtained from the optimization.
% Option 'likelihood' only calculates the likelihood.