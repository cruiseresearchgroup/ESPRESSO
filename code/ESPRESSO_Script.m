% run ESPRESSO
function ESPRESSO(K, data, subsequence, chain_len)
%%% Input 
% data : input original data ,
% K : number of segments, 
% subsequence : window size
% chain_len : maximum length of chain for each subsequence
%%% returns:
% extracted segment boundaries, F-score, MAE



    [numTS,lenTS] = size(data)    % number of time-series and length of them

    [MP, MPI] = computMP(data, subsequence);

    for i = [1:size(MP,1)]
        [~, wcac(i,:)] = calculateSemanticDensityMatrix(MP,MPI, chain_length, subsequence);
    end
    [espTT,~] = separateGreedyIG(data, K, wcac, 0.01);
    return espTT

function [MP, MPI] = computMP(Integ_TS, subsequence)
    for i=1:1:size(Integ_TS,1)
        [MP(i,:), MPI(i,:)] = timeseriesSelfJoinFast(Integ_TS(i,:),subsequence);
    end
end


%% To calculate the score use the following 
% 
%    [Precision,Recal,rmse1,rmse2,Fscore] = calculateScore(GT_TT, espTT, lenTS, margin);
%    [mae,~] = calcScoreMAE(GT_TT, espTT, lenTS);
