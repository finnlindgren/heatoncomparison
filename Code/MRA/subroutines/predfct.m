function [ predsj ] = predfct( postmeanj,postvarj, Btildej, Rpostcholb, KcAb, Kcwb )
%Prediction at finest resolution at the end
M = size(postmeanj,2);
for k=M-1:-1:1
    KcBtilde = Rpostcholb{k}\Btildej{k+2}{k}';
    Btildej{k+1} = cell(k,1);
    for l=(k-1):-1:1
        Btildej{k+1}{l} = Btildej{k+2}{l} - KcBtilde'*KcAb{k}{l};
    end
end
for k=1:M-1
    KcBtildecur = Rpostcholb{k}\Btildej{k+2}{k}';
    postmeanj(:,k) = KcBtildecur'*Kcwb{k};
    postvarj(:,k) = sum(KcBtildecur.^2,1);
end

% predsj = [sum(postmeanj(1:end),2),sum(postvarj(1:end),2)];
predsj= [sum(postmeanj,2),sum(postvarj,2)];
end

