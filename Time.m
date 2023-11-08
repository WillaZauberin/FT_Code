function numt=Time(t,h)

 numt=zeros(h,1);%时间变化的数集
for i=1:h
    tt=strsplit(t{i,1},':');%在“：”处把时间字符串拆分
   miao1(i)=str2num(tt{4});%读取第4列的秒数据
    miao2(i)=str2num(tt{5});
   fen(i)=str2num(tt{3});
%    shi(i)=str2num(tt{1});
miao(i)=miao1(i)+0.001*miao2(i);
   numt(i)=(fen(i)-fen(1))*60+(miao(i)-miao(1));
end

