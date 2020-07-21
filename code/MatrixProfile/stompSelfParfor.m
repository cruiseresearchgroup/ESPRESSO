% Compute the self-similarity join of a given time series A
% Yan Zhu 03/09/2016 modified by Chin-Chia Michael Yeh 03/10/2016
%
% [matrixProfile, matrixProfileIndex] = stompSelfParfor(A, subLen, workerNum)
% Output:
%     matrixProfile: matrix porfile of the self-join (vector)
%     matrixProfileIndex: matrix porfile index of the self-join (vector)
% Input:
%     A: input time series (vector)
%     subLen: interested subsequence length (scalar)
%     workerNum: the number of worker for MATLAB's Parallel Computing
%                Toolbox (scalar)
%
% Yan Zhu, Zachary Zimmerman, Nader Shakibay Senobari, Chin-Chia Michael Yeh, Gareth Funning,
% Abdullah Mueen, Philip Brisk and Eamonn Keogh, "Matrix Profile II: Exploiting a Novel 
% Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and 
% Joins," ICDM 2016, http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [matrixProfile, profileIndex] = stompSelfParfor(data, subLen, workerNum)
%% setup parpool
if isempty(which('parpool'))
    if matlabpool('size') <= 0 %#ok<*DPOOL>
        matlabpool(workerNum);
    elseif matlabpool('size')~= workerNum
        matlabpool('close');
        matlabpool(workerNum);
    end
else
    parProfile = gcp('nocreate');
    if isempty(gcp('nocreate'))
        parpool(workerNum);
    elseif parProfile.NumWorkers ~= workerNum
        delete(gcp('nocreate'));
        parpool(workerNum);
    end
end

%% set trivial match exclusion zone
excZone = round(subLen/2);

%% check input
if subLen > length(data)/2
    error('Error: Time series is too short relative to desired subsequence length');
end
if subLen < 4
    error('Error: Subsequence length must be at least 4');
end
if length(data) == size(data, 2)
    data = data';
end

%% check skip position
proLen = length(data) - subLen + 1;
skipLoc = false(proLen, 1);
for i = 1:proLen
    if any(isnan(data(i:i+subLen-1))) || any(isinf(data(i:i+subLen-1)))
        skipLoc(i) = true;
    end
end
data(isnan(data)) = 0;
data(isinf(data)) = 0;

%% initialization
proLen = length(data) - subLen + 1;
matrixProfiles = cell(1, workerNum);
profileIndices = cell(1, workerNum);
[dataFreq, dataLen, dataMean, dataSig] = fastFindNNPre(data, subLen);

%% seperate index into different job
jobIdx = zeros(workerNum, 2);
jobIdx(1, 1) = 1;
jobIdx(1, 2) = round(proLen/workerNum);
for i = 2:workerNum
    jobIdx(i, 1) = jobIdx(i-1, 2)+1;
    jobIdx(i, 2) = jobIdx(i, 1) + round(proLen/workerNum) - 1;
end
jobIdx(workerNum, 2) = proLen;

%% compute the matrix profile
query = data(1:1+subLen-1);
[~, firstProduct, ~, ~] = ...
    fastFindNN(dataFreq, query, dataLen, subLen, dataMean, dataSig);

parfor jobi = 1:workerNum
    querySum = 0;
    query2Sum = 0;
    dropVal = 0;
    distProfile = zeros(proLen, 1);
    lastProduct = zeros(proLen, 1);
    pickedIdx = jobIdx(jobi, 1):jobIdx(jobi, 2);
    for i = 1:length(pickedIdx)
        % compute the distance profile
        idx = pickedIdx(i);
        query = data(idx:idx+subLen-1);
        if i == 1
            [distProfile(:), lastProduct(:), querySum, query2Sum, querySig] = ...
                fastFindNN(dataFreq, query, dataLen, subLen, dataMean, dataSig);
            distProfile(:) = real(distProfile);
        else
            querySum = querySum-dropVal+query(subLen);
            query2Sum = query2Sum-dropVal^2+query(subLen)^2;
            queryMean=querySum/subLen;
            querySig2 = query2Sum/subLen-queryMean^2;
            querySig = sqrt(querySig2);
            lastProduct(2:dataLen-subLen+1) = lastProduct(1:dataLen-subLen) - ...
                data(1:dataLen-subLen)*dropVal + data(subLen+1:dataLen)*query(subLen);% change here
            lastProduct(1)=firstProduct(idx);
            distProfile(:) = 2*(subLen - ...
                (lastProduct-subLen*dataMean*queryMean)./(dataSig*querySig));
            distProfile(:) = real(distProfile);
        end
        dropVal=query(1);
        
        % apply exclusion zone
        excZoneStart = max(1, idx-excZone);
        excZoneEnd = min(proLen, idx+excZone);
        distProfile(excZoneStart:excZoneEnd) = inf;
        distProfile(dataSig<eps) = inf;
        if skipLoc(i) || (querySig < eps)
            distProfile = inf(size(distProfile));
        end
        distProfile(skipLoc) = inf;
        
        % figure out and store the neareest neighbor
        if i == 1
            matrixProfiles{jobi} = inf(proLen, 1);
            profileIndices{jobi} = inf(proLen, 1);
        end
        updatePos = distProfile < matrixProfiles{jobi};
        profileIndices{jobi}(updatePos) = idx;
        matrixProfiles{jobi}(updatePos) = distProfile(updatePos);
    end
end

%% merge output of parfor
matrixProfile = inf(proLen, 1);
profileIndex = zeros(proLen, 1);
for jobi = 1:workerNum
    updatePos = matrixProfiles{jobi} < matrixProfile;
    matrixProfile(updatePos) = matrixProfiles{jobi}(updatePos);
    profileIndex(updatePos) = profileIndices{jobi}(updatePos);
end
matrixProfile = sqrt(matrixProfile);
matrixProfile(skipLoc) = inf;
profileIndex(skipLoc) = inf;

%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [dataFreq, dataLen, dataMean, dataSig] = fastFindNNPre(data, subLen)
% compute stats about data
dataLen = length(data);
data(dataLen+1:2*dataLen) = 0;
dataFreq = fft(data);
dataCum = cumsum(data);
data2Cum =  cumsum(data.^2);
data2Sum = data2Cum(subLen:dataLen)-[0;data2Cum(1:dataLen-subLen)];
dataSum = dataCum(subLen:dataLen)-[0;dataCum(1:dataLen-subLen)];
dataMean = dataSum./subLen;
dataSig2 = (data2Sum./subLen)-(dataMean.^2);
dataSig = sqrt(dataSig2);

function [distProfile, lastProduct, querySum, query2Sum, querySig] = ...
    fastFindNN(dataFreq, query, dataLen, subLen, dataMean, dataSig)
% proprocess query for fft
query = query(end:-1:1);
query(subLen+1:2*dataLen) = 0;

% compute the product
queryFreq = fft(query);
productFreq = dataFreq.*queryFreq;
product = ifft(productFreq);

% compute the stats about query
querySum = sum(query);
query2Sum = sum(query.^2);
queryMean=querySum/subLen;
querySig2 = query2Sum/subLen-queryMean^2;
querySig = sqrt(querySig2);

% compute the distance profile
distProfile = 2*(subLen-(product(subLen:dataLen)-subLen*dataMean*queryMean)./...
    (dataSig*querySig));
lastProduct=real(product(subLen:dataLen));