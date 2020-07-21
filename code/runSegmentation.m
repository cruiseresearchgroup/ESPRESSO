% Script for running time series segmentation. Computes the matrix profile 
% and the matrix profile index and segments the time series. Returns all
% computed values for reference
% 
% function [crosscount, splitLoc, MatrixProfile, MPindex] = 
%                                      RunSegmentation(ts, slWindow)
% Input parameters:
% ts - time series to segment
% slWindow - sliding window size
%
% Output parameters:
% crosscount - number of crossings at each point
% splitLoc - split locations
% MatrixProfile - matrix profile
% MPindex - matrix profile index
function [crosscount,wcc,semcc,semArcSet,newWCC, newCC] = ...
        runSegmentation(MatrixProfile, MPindex, slWindow )    
    
    %[MatrixProfile, MPindex] = timeseriesSelfJoinFast(ts, slWindow);
    %[MatrixProfile, MPindex] =  stompSelfParfor(ts, slWindow,4);
    [crosscount1] = segmentTimeSeries(MPindex);
    [crosscount] = normCrossCountAll(crosscount1, slWindow);
    
    [wcc1, newWCC1] = calculateMotifDensityMatrix(MatrixProfile,MPindex);
    
    [wcc] = normCrossCountAll(wcc1, slWindow);
    
    [newWCC] = normCrossCountAll(newWCC1, slWindow);
%      subplot(4,2,3);  
%      plot(crosscount,'b');
%  	title('cross counts');
%     subplot(4,2,4);  
%      plot(crosscount1,'b');
%  	title('cross counts1');
%      
%     subplot(4,2,5);    
%     plot(wcc,'b');
%  	title('Weighted crosscounts');
%     subplot(4,2,6);    
%     plot(wcc1,'b');
%  	title('Weighted crosscounts1');
    %%CALCUALTE SEMANTIC DENSITY MATRIX PROFILE
    [semArcSet, newCC] = calculateSemanticDensityMatrix(MatrixProfile,MPindex, 5,slWindow);
    [semcc] = normCrossCountAll(semArcSet, slWindow);
    
%     subplot(4,2,7);    
%     plot(semcc,'b');
%  	title('semantci cross counts');
%     subplot(4,2,8);    
%     plot(semArcSet,'b');
%  	title('semantci cross counts1');
%     figure('name','semanticc');
%     plot(semArcSet);
    
    %workOnfindingMinimum(crosscount1, crosscount, wcc1,wcc, semArcSet,semcc);

%     figure('name','wcc');
%     plot(wcc1);
    