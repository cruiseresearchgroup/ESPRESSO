function[TT,maxVAR]=DP_IG(Integ_TS,k,step,double)
Integ_TS=Clean_TS(Integ_TS,double);
[Nu_TS,Le_TS]=size(Integ_TS);n=Le_TS;
%step=1;
cost(n,n,1)=0;
for i=1:step:Le_TS
    for j=i+1:step:Le_TS
        if(i==j)
            cost(i,j,1)=0;
            continue;
        end
        for l=1:Nu_TS
             TS_dist(l)=Integ_TS(l,j)-Integ_TS(l,i);
        end
        cost(i:i+step-1,j:j+step-1,1)=((j-i)/Le_TS)*Sh_Entropy(TS_dist);
        
    end
    
end


 for b=2:k+1
     for i=2:n
         cost(1,i,b)=cost(1,i,b-1);
         pos(i,b)=1;
         for j=1+step:step:i-1
             if cost(1,j,b-1)+cost(j+1,i,1)<= cost(1,i,b)
                 cost(1,i,b)=cost(1,j,b-1)+cost(j+1,i,1);
                 pos(i,b)=j;
             end
         end
     end
 end
 maxVAR=cost(1,n,k+1);
 
 id=n;
 for b=k+1:-1:2
     TT(b-1)=pos(id,b);
     id=TT(b-1);
 end
 