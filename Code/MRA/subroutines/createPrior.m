function [ RpriorCholj,KcBc,Atj, wtj, retlikpred ] = createPrior( theta, M,knotsb,RpriorCholb,KcBb,dataj,varargin )
%createPrior Creates prior values
%   for current region and ancestry; this function contains optional input
%   and output arguments depending on wheter the level is the last level
%   and/or if predictions are made, optional outputs are Atj,wtj,Btildej

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
 
 
 % Preliminaries
 L=length(knotsb); % Find current level
 mlM=(L<M); % Indicator if current level is less than last level
 
 % Create prior quantities
 RpriorCholj=[];
 KcBc=cell(L-1,1);
 V=cell(L,1);
 for l=1:L
     V{l,1}=co(knotsb{L,1}, knotsb{l,1},theta);
     switch l  % Switch construct to avoid going in the loop for highest level
         case 1
             %KcBc{l}=V{l}'; %%% CHECK TEMPORARY FIX
         otherwise
             for k=1:(l-1)
                 if l<L
                     V{l}=V{l}-KcBc{k}'*KcBb{l,1}{k,1};% ????
                 else
                     V{l}=V{l}-KcBc{k}'*KcBc{k};
                 end
             end
     end
     if l<L
         KcBc{l}=RpriorCholb{l}\V{l}';
     else
         RpriorCholj=chol(V{l}, 'lower'); % MISSING Varepsilon!
         % New 6/14/2017 DMH
         RpriorCholj=RpriorCholj+diag(linspace(varEps,varEps,size(RpriorCholj,1)));
         
     end
     
 end
 %varargout{1}=V; Just for testing. REMOVE


if mlM<1  % Check if region is at lowest level
     % Begin inference at lowest level
        % precompute solves
        Sicy=RpriorCholj\dataj;
        SicB=cell(L-1,1);
        for l=1:L-1
            SicB{l}=RpriorCholj\V{l};
        end
        % inference quantities
        wtj=cell(L-1,1);
        Atj=cell(L-1,L-1);
        for l=1:L-1
            wtj{l}=SicB{l}'*Sicy;
            
            for k =l:L-1
                Atj{l,k}=SicB{l}'*SicB{k};
            end
        end
        %varargout{1}=Atj;
        %varargout{2}=wtj;
        % Check if predicting
        if isnan(predlocsj)  % If NOT predicting
            loglikj=2*sum(log(diag(RpriorCholj)))+Sicy'*Sicy;
            retlikpred=loglikj;
            %varargout{3}=retlikpred;
        else  % If predicting
            RpriorChol=[RpriorCholb;{RpriorCholj}];
            KcB=[KcBb;{KcBc}];
            % calculate Bp and L
            KcBp=cell(L,1);
            Vp=cell(L,1);
            for l=1:L
                Vp{l}=co(predlocsj,knotsb{l},theta);
                for k=1:l-1
                    Vp{l}=Vp{l}-KcBp{k}'*KcB{l,1}{k,1};
                end
                KcBp{l}=RpriorChol{l}\Vp{l}';
            end
            
            Vpp=co(predlocsj,predlocsj,theta); % Covariance matrix of prediction locations
            for l=1:(L-1) % CHANGED to L-1 from L 10 July 2017 DH
                Vpp=Vpp-KcBp{l}'*KcBp{l};
            end
            % Initialize prediction inference
            postmeanj=NaN(size(predlocsj,1),L); % prediction mean matrix for all levels
            postvarj=NaN(size(predlocsj,1),L); % prediction variance matrix for all levels
            postmeanj(:,L)=KcBp{L}'*Sicy;
            postvarj(:,L)=diag(Vpp - KcBp{L}'*KcBp{L});
            Btildej=cell(L,1);
            for k=1:L-1
                Btildej{L+1}{k}=Vp{k}-KcBp{L}'*SicB{k};
            end
            retlikpred={postmeanj, postvarj,Btildej};
        end
else  % if NOT lowest level
    wtj=NaN;Atj=NaN;retlikpred=NaN; % Assign NaNs to outputs
end

end


%% Notes

% Important change in line 96: Largest value of loop changed to L-1 from L