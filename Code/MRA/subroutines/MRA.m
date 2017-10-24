function [ sum_loglik,preds ] = MRA( theta, data, knots, M, J, nRegions,varargin )
% MRA Main MRA function
%   Detailed explanation goes here

% Check number of optional input arguments
numvarargsin = length(varargin);
if numvarargsin > 2
    error('myfuns:createPrior:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {0 NaN};
% overwrite the ones specified in varargin.
optargs(1:numvarargsin) = varargin;
 [varEps, predlocsj]=optargs{:};

% Calculate key quantities
totalRegions=sum(nRegions);
cumRegions=cumsum(nRegions);

% Pre-allocate quantities
RPriorChol=cell(totalRegions,1);
KcB=cell(totalRegions,1);
Atildeprev=cell(totalRegions,1);
wtildeprev=cell(totalRegions,1);
loglik=NaN(totalRegions,1);
if iscell(predlocsj)==1  % If predicting
    postmean=cell(totalRegions,1);
    postvar=cell(totalRegions,1);
    Btilde=cell(totalRegions,1);
    preds=cell(totalRegions,1);
    Rpostchol=cell(totalRegions,1);
    KcA=cell(totalRegions,1); % Could be matrices? Entries are vectors
    Kcw=cell(totalRegions,1);  % Could be matrices? Entries are vectors   
else  % New DH: 06/24/2017 workaround to pass correct object to createPrior
    predlocsj = num2cell(nan(totalRegions,1)); 
    
end
%data=cell(totalRegions,1); % NEEDS TO MOVE UP

%% loop from coarsest to finest level

for l=1:M
  
% loop through all i for a given level and find ancestry
for i=(cumRegions(l)-nRegions(l)+1):cumRegions(l)
   i_ancestry = find_ancestry( i, nRegions,J );
   [temp1, temp2, temp3, temp4, temp5]=createPrior(theta, M, knots([i_ancestry;i],1),RPriorChol(i_ancestry,1),KcB(i_ancestry,1), data{i,1}, varEps, predlocsj{i,1});
   % , temp3, temp4, temp5
   %[ RpriorCholj,KcBc,varargout ] = createPrior( theta, M,knotsb,RpriorCholb,KcBb,dataj,varargin )
   RPriorChol{i}=temp1;
   KcB{i}=temp2;
   if l==M  % If highest resolution level
      Atildeprev{i}=temp3; 
      wtildeprev{i}=temp4;
      if isnan(predlocsj{end})==0  % If predicting
          postmean{i}=temp5{1};
          postvar{i}=temp5{2};
          Btilde{i}=temp5{3};
      else % If not predicting
       loglik(i,1)=temp5;   
      end
   end
end
disp(['Prior Level ',num2str(l), ' completed'])
end
clear KcB; % To save memory

%% Posterior inference 
% From second finest to coarsest resolution

for l=M-1:-1:1
    %fprintf(['Level Posterior ',num2str(l)]);
    disp(['Posterior Level ',num2str(l),' starting']);
    Atildecur=cell(totalRegions,1);
    wtildecur=cell(totalRegions,1);
    % loop through  all regions i for the level 
    for i=(cumRegions(l)-nRegions(l)+1):cumRegions(l)
        
        %for i=cumRegions(end-1):-1:1 % DELETE
        [ i_children ] = find_children( i,nRegions,J  );% Find indices of children
        % Calculate posterior quantities
        Rpriorcholj=RPriorChol{i};
        wtildechildren=wtildeprev(i_children);
        Atildechildren=Atildeprev(i_children);
        
        [ wtildecurj, Atildecurj,loglikj, Rpostcholj, Kcwj, KcAj ] = postInf( Rpriorcholj, wtildechildren, Atildechildren );
% Rpostcholj(1,1)
% Kcwj(1,1)
% KcAj{1,1}(1,1)

        wtildecur{i}=wtildecurj;
        Atildecur{i}=Atildecurj;
        
        if isnan(predlocsj{end})<1
            Rpostchol{i} = Rpostcholj;
            Kcw{i} = Kcwj;
            KcA{i} = KcAj;
            
        else
        loglik(i,1)=loglikj;
        end
        %i
    end
    wtildeprev=wtildecur;
    Atildeprev=Atildecur;
end
sum_loglik=sum(loglik);

%% Spatial prediction
if isnan(predlocsj{end})==0    
    for i=(cumRegions(M)-nRegions(M)+1):cumRegions(M)
        if (M > 0)
            i_ancestry = find_ancestry( i, nRegions,J );
            preds{i,1} = predfct(postmean{i,1},postvar{i,1},Btilde{i,1},Rpostchol(i_ancestry,1),KcA(i_ancestry,1),Kcw(i_ancestry,1));
        else
            preds{i,:}=[postmean(i,1),postvar(i,1)];
        end
    end
% end
end


%% Notes:

% DH: 06/24/2017: logical check if prediction locations are passed: check
% if predlocsj is a cell array. If not, it is NaN. Could potentially cause
% problems if we switch away from using cell arrays.