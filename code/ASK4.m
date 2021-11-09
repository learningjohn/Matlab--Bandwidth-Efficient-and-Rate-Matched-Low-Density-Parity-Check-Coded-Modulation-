N=10^5;           %要传输的符号数
s_i=[-3,-1,1,3];  %要传输的四个符号
Es_N0_dB=0:0.5:25;  %设置Es_N0的范围
err=zeros(1,length(Es_N0_dB));
sim_ser=zeros(1,length(Es_N0_dB));
s_r=zeros(1,N);   %接收到的符号
dff=zeros(1,N);   %符号之差

for i=1:length(Es_N0_dB)
s=randsrc(1,N,s_i); %令各个符号等概率出现
s_d=(1/sqrt(5))*s;  %归一化各个符号的距离,令Es=1
SER=zeros(1,length(Es_N0_dB)); %先预设好SER
n=randn(1,N);       %产生N个噪声，其均值为0，方差为1
s_n=sqrt(10^(-Es_N0_dB(i)/10))*1/sqrt(2)*n;   %给噪声加上合适的方差
y=s_d+s_n;          %给每一个符号都加上噪声
%%%*******解调部分*******%距离判定***这个是仿真部分，也就是用概率去算的
s_r(find(y<(-2/sqrt(5))))=-3;
s_r(find(y>=(2/sqrt(5))))=3;
s_r(find(y>=(-2/sqrt(5))& y<0))=-1;
s_r(find(y<(2/sqrt(5))&y>=0))=1; 
%%%*****将s和s_r相减，返回非0的数组然后用find函数找出它们的位置********%
dff=abs(s-s_r);
err(i)=size(find(dff),2);
end
%%%用概率法求出SER和Es_N0_dB的关系
sim_ser=err/N;
semilogy(Es_N0_dB,sim_ser,'mx-');
%theoryBer = 0.75*erfc(sqrt(0.2*(10.^(Es_N0_dB/10))));
theoryBer = (7/8)*erfc(sqrt((10.^(Es_N0_dB/10))/21));
axis([0 25 10^-5 1]);
hold on 
semilogy(Es_N0_dB,theoryBer,'b.-');
hold off
grid on
legend('simulation','theory');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('SER for 4-ASK modulation')
addpath('mat_data\');
load('best_P_8ASK.mat');
M = 8;
A = PX*(([M-1:2:1,1:2:M-1]'/2).^2);
%A = PX*[9/4;1/4;1/4;9/4];
%theoryBer_best = PX*[1/2;1;1;1/2].*erfc(sqrt((10.^(Es_N0_dB/10)')./(4*A)));
theoryBer_best = PX*[1/2;ones(M-2,1);1/2].*erfc(sqrt((10.^(Es_N0_dB/10)')./(4*A)));