function y=SH_Entropy(x)
x(x==0)=[];       
ss=sum(x);
sum1=sum((x/ss).*log(x/ss));
y=-1*sum1;
       
 