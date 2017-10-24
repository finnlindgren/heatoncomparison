function [ xmin, xmax,ymin, ymax ] = create_partition( xmin0,xmax0, ymin0, ymax0,J )
% create_partition specifies child partitions
%   Partitions are specified by their mininal and maximal x and y values
switch  J
    case 2
        if (xmax0-xmin0)>=(ymax0-ymin0) % X dimension is larger or both equal
            tempx=linspace(xmin0,xmax0, J+1);
            xmin=repmat(tempx(1:end-1), J/2,1);xmin=xmin(:);
            xmax=repmat(tempx(2:end),J/2,1);xmax=xmax(:);
            
            tempy=linspace(ymin0,ymax0, J/2+1);
            ymin=repmat(tempy(1:end-1),1, J/2+1);ymin=ymin';
            ymax=repmat(tempy(2:end),1, J/2+1);ymax=ymax';
            
        else
            tempx=linspace(xmin0,xmax0, J/2+1);
            xmin=repmat(tempx(1:end-1),1, J/2+1);xmin=xmin(:);
            xmax=repmat(tempx(2:end),1,J/2+1);xmax=xmax(:);
            
            tempy=linspace(ymin0,ymax0, J+1);
            ymin=repmat(tempy(1:end-1),1, J/2);ymin=ymin';
            ymax=repmat(tempy(2:end),1, J/2);ymax=ymax';        
        end
    case {4,16}
        tempx=linspace(xmin0,xmax0, J/2+1);
        xmin=repmat(tempx(1:end-1), J/2,1);xmin=xmin(:);
        xmax=repmat(tempx(2:end),J/2,1);xmax=xmax(:);
        
        tempy=linspace(ymin0,ymax0, J/2+1);
        ymin=repmat(tempy(1:end-1),1, J/2);ymin=ymin';
        ymax=repmat(tempy(2:end),1, J/2);ymax=ymax';
        
        
end


end


%% Notes: this function should work for all even Js. Idea is to separate between J=2,4,8
