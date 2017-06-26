close all;clear all; 
 
disp('********************************************************'); 
disp('***********波导宽边纵向并联缝隙驻波天线设计***********'); 
disp('********************************************************'); 
 
f=input('工作中心频率(GHz):f='); 
er=input('波导内介质相对介电常数：er='); 
lambda=300/f/sqrt(er);%波长，mm 
k=2*pi/lambda; 
a=input('波导宽边尺寸(mm)：a='); 
b=input('波导窄边尺寸(mm)：b='); 
lambda_g=lambda/sqrt(1-(lambda/(2*a))^2);%波导波长 
nn=input('请输入缝隙单元数目：'); 
 
disp('请输入单元激励功率分布类型：'); 
disp('1.功率等幅激励；2.其他功率分布激励.'); 
amp_type=input('');%激励功率分布类型 
 
% aa 就代表了功率 
switch amp_type  %功率分布类型 
    case 1 %等幅激励 
        aa=ones(1,nn); 
    case 2 %其他分布 
        for n=1:nn 
            fprintf('请输入第 %d 个单元分配的功率：',n); 
            aa(n)=input(''); 
        end 
    otherwise 
        disp('WRONG SELECTION!!!'); 
        %break;  
end 
 
K=1/sum(aa); 
g=K.*aa;%每个缝隙的电导 
g0=2.09*a*lambda_g/b/lambda*(cos(lambda*pi/2/lambda_g))^2; 
%Dn为每个缝隙距宽边中心的距离 
D=asin(sqrt(g./g0))*a/pi; 
 
ls=lambda_g/2; 
l1=lambda_g/4;%并联缝隙，短路板距末端缝中心为lambda_g/4，裂缝总是位于驻波电压的波峰点 
l0=lambda/2; 
fprintf('缝隙长度(介质波导填充时，此值较小，应靠近半个自由空间波长 %4.2fmm)：l0= %4.2fmm\n',300/f/2,l0); 
fprintf('缝隙间距：ls= %4.2fmm\n',ls); 
fprintf('最后一个缝隙中心距终端短路间距：l1= %4.2fmm\n',l1); 
 
for n=1:nn 
    fprintf('第 %d 个缝隙距中心间距：D(%d)=%4.3fmm\n',n,n,(-1)^(n+1)*D(n)); 
end 
 
fprintf('缝隙宽度取为lambda/40--lambda/20：%4.3fmm ---%4.3fmm\n',lambda/40,lambda/20); 
disp('存在波导壁厚时，需调整lambda_g，缝长和缝隙距中心距离'); 
 
bandwidth=0.5/nn;%对波导谐振阵来说，驻波小于2的工作带宽只有：50%/nn 
fprintf('该波导谐振阵驻波小于2的相对带宽为：%4.2f %%\n',bandwidth*100); 