function [ wtildecurj, Atildecurj,loglikj, Rpostcholj, Kcwj, KcAj ] = postInf( Rpriorcholj, wtildechildren, Atildechildren )
% postInf Posterior inference
%   From second finest to coarsest resolution
% Prediction option is NOT implemented yet, that requires one additional
% input argument
l=length(wtildechildren{1}); % Gives resolution; help construct to get resolution
w=cell(l,1);
A=cell(l,l); % This is a square cell structure to hold tiles

for j=1:l  % j is meaningless counter, in Matthias' code this counter is l
    temp=[];% helper construct to duplicate Julia functionality
    for c=1:length(wtildechildren)
        temp=[temp wtildechildren{c}{j}];
    end
    w{j}=sum(temp,2); clear temp;
    
    for k=j:l
        temp=NaN(size(Atildechildren{1, 1}{1, 1}));% helper construct to duplicate Julia functionality
        temp(:,:,length(Atildechildren))=NaN(size(Atildechildren{1, 1}{1, 1}));% helper construct to duplicate Julia functionality
        for c=1:length(Atildechildren)
            temp(:,:,c)=Atildechildren{c}{j,k};
        end
        A{j,k}=sum(temp,3);clear temp;
        
    end
end

% Calculate Cholesky of K.inv
Rpost=Rpriorcholj*Rpriorcholj'+A{l,l};
Rpostcholj=chol(Rpost,'lower');

% pre-compute the solves required later
Kcwj=Rpostcholj\w{l};
KcAj=cell(l-1,1);
for j=1:l-1
    KcAj{j}=Rpostcholj\A{j,l}';
end

% calculate w.tilde and A.tilde

if l==1
    wtildecurj=NaN;
    Atildecurj=NaN;
else
    wtildecurj=cell(l-1,1); % Each cell holds vectors
    Atildecurj=cell(l-1,l); % This is a square cell structure to hold tiles
    for j=1:l-1
        wtildecurj{j}=w{j}-KcAj{j}'*Kcwj;
        for k=j:l-1
            Atildecurj{j,k}=A{j,k}-KcAj{j}'*KcAj{k};    
        end
    end   
end

loglikj=2*sum(log(diag(Rpostcholj)))-2*sum(log(diag(Rpriorcholj)))-Kcwj'*Kcwj;
% Values to return: wtildecurj, Atildecurj,loglikj

end
%% Notes: 
% This function contains complicated constructs to duplicate Julia
% functionality to sum up over children.

% The prediction option for this function is not yet programmed.