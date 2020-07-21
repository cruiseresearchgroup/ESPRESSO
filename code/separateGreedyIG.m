function [bestTT,bestIG] = separateGreedyIG(ts, numSegms, CC,pdist)
[numberofTS, dataLength] = size(ts);
Integ_TS = Clean_TS(ts,1);

pdist = floor(dataLength * pdist);
CC = -1 * smoothdata(CC,'gaussian',5);
bestTT = [];
bestIG = 0;
for d = 1 : numberofTS
    [pks, locs, width, proms] = findpeaks(CC(d,:),'SortStr','descend', 'minpeakdistance',pdist);
    [p, idx] = sort(proms,'descend');
    TT = [];
    maxIG = zeros(1,length(locs));
    if length(locs) > 0
        remain_locs = locs;
        for i = 1 : length(locs)
            tempTT = [];
            for j = 1 : length(remain_locs)
                c = [TT, remain_locs(j), dataLength];
                ig = IG_Cal(Integ_TS, c,i);
                if(ig > maxIG(i))
                    maxIG(i) = ig;
                    tempTT = remain_locs(j);
                end
            end
            TT = [TT, tempTT];
            remain_locs(remain_locs == tempTT) = [];
        end
        
        if(numSegms-1 > length(maxIG))
            t = length(maxIG); 
        else
            t = numSegms-1;
        end
        if(maxIG(t) > bestIG)
            bestIG = maxIG(t);
            bestTT = TT(1:t);
        end
    else
        bestIG = [];
        bestTT = [];
    end
end

end
