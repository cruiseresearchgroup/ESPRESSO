function [score , score2] = calcScore(groundTruth, detectedSegLoc, dataLength)
[m n] = size(groundTruth);				%SHOHREH-- m and l always 1 . yeah?
[l k] = size(detectedSegLoc);
if(l ~=1)
    detectedSegLoc= detectedSegLoc';
    [l k] = size(detectedSegLoc);
end

sumDiff = 0;
ind(1:n) = -1;
minV(1:n) = inf;
indNearest = 1;
%while(indNearest ~= -1)
for(j = 1:1:n)
    for(i = 1:1:k)
        if(abs(detectedSegLoc(i)-groundTruth(j)) < abs(minV(j)))
            minV(j) = abs(detectedSegLoc(i) - groundTruth(j));
            ind(j) = i;
        end
    end
end
sumOfDiff = sum(minV);
score = sumOfDiff/(dataLength*length(minV));

score2 = sqrt(sumOfDiff.^2/dataLength);
end
