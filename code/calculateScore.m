function [P,R,rmse1,rmse2,F] = calculateScore(GT, TT, dataLength, margin)
% For segmentation:
% 
% A true positive is where you estimate one boundary at a true boundary position
% 
% A false positive is where you estimate an additional boundary (number of estimates - 1) at a true boundary position
% 
% A false negative is where you don't estimate a boundary at a true boundary position

n = length(GT);				
k = length(TT);
if k>0

diff1=[];diff2=[];

ind(1:n) = -1; positives(1:n)=0;
selected(1:k)=false;
minV(1:n) = inf;
TP=0;
for(j = 1:1:n)
    for(i = 1:1:k)
        if(abs(TT(i)-GT(j)) < abs(minV(j)))
            minV(j) = abs(TT(i) - GT(j));
            ind(j) = i;
        end
    end
    if(minV(j) <= margin)
        TP = TP+1;
        selected(ind(j))=true;
    end
end

dminV2=minV;
for(j = 1:n)
    if(minV(j) >  margin)
        nearestEstimate=max(GT);
        for t=1:k
            if (~selected(t) && abs(TT(t)-GT(j)) < nearestEstimate )
                nearestEstimate =abs(TT(t)-GT(j));
                ind(j) = t;
            end
        end
        dminV2(j)=nearestEstimate;
        selected(ind(j))=true;
    end
end

P=0; F=0; R=0;
FP=0; FN=0;

%RMSE1
rmse1 = sqrt(sum(minV.^2)/length(GT))/dataLength;

%RMSE2
rmse2 = sqrt(sum(dminV2.^2)/length(GT))/dataLength;

%F-score

P=TP/length(GT);

detected(1:k)=0;
for(j=1:n)
    for(i=1:k)
        if(abs(TT(i)-GT(j)) <= margin)
            positives(j)= positives(j)+1;
            detected(i) = detected(i)+1;
        end
    end
end

FN = sum(positives==0);
FP = sum(positives(positives>=2)-1)+ sum(detected==0);
TP = sum(positives>0);

R = TP / (TP + FN);
P = TP / (TP + FP);
if (P==0 && R==0)F=0;else
    F = 2 *((P * R)/(P + R));
end
else
    P=0;F=0;R=0;rmse1=1;rmse2=1;
end

