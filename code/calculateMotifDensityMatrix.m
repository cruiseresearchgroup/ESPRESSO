function [crosscount, newc] = calculateMotifDensityMatrix(MatrixProfile,MPindex)
profile_len = length(MPindex);
nnmark = zeros(1,profile_len);
newmark = zeros(1,profile_len);

%%FIND THE MAXIMUM LENGTH OF AN ARC
max_len =0;
min_len = inf;
for i=1:profile_len
    len = abs(i-MPindex(i));
    if(len > max_len)
        max_len = len;
    end
    if(len<min_len)
        min_len=len;
    end
end
diff=max_len-min_len;
miu = mean(MatrixProfile);
sim_diff = max(abs(min(MatrixProfile)-miu), abs(max(MatrixProfile)-miu));
%%AGGREGAT WEIGHTED-ARCS (GIVE WEIGHT BASED ON LENGTH)
%%AIM : WEIGHT ~ 1/LENGTH
for i=1:profile_len
    
        small=min(i,MPindex(i));
        large=max(i,MPindex(i));
        arc_length = abs(i-MPindex(i));
        %arc_weight = 1- ((arc_length-min_len)/diff);
        %if(arc_length == min_len) arc_weight =1;else
        arc_weight =1 - ((arc_length-min_len)/diff).^2 ; 
        nnmark(small:large)=nnmark(small:large)+arc_weight; 
        newmark(small:large) = newmark(small:large) + (1+1*((MatrixProfile(i)-miu)/sim_diff) * arc_weight); 
end
crosscount = nnmark;
newc = newmark;


