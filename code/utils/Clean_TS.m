%-----------------------------------------------
%This function normalizes the time series and 
%doubles the number of the time series to 
%address the hetergenousity in the time series
%-----------------------------------------------
function[N_Integ_TS] =Clean_TS(O_Integ_TS,double)
%double=1 Normalize and double 
%double=0 Normalize  
%double=2 without change 

Integ_TS=O_Integ_TS;
[m, n]=size(Integ_TS);
for i=1:m
  minVal=min(Integ_TS(i,:));
  if (double==2) minVal=0;end
  Integ_TS(i,:)=Integ_TS(i,:)-minVal;
  if double~=2
  sumVal=sum(Integ_TS(i,:))/1000;   %SHOHREH - 1000??
  Integ_TS(i,:)=Integ_TS(i,:)/sumVal;
  end
  if double==1
    maxVal=max(Integ_TS(i,:));
    Integ_TS(i+m,:)=maxVal-Integ_TS(i,:);
    sumVal=sum(Integ_TS(i+m,:))/1000;
    Integ_TS(i+m,:)=Integ_TS(i+m,:)/sumVal;
  end
end
N_Integ_TS=Integ_TS;
N_Integ_TS=cumsum(Integ_TS,2);
