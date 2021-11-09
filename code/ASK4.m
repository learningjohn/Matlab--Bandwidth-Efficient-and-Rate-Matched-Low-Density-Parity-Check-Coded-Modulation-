N=10^5;           %Ҫ����ķ�����
s_i=[-3,-1,1,3];  %Ҫ������ĸ�����
Es_N0_dB=0:0.5:25;  %����Es_N0�ķ�Χ
err=zeros(1,length(Es_N0_dB));
sim_ser=zeros(1,length(Es_N0_dB));
s_r=zeros(1,N);   %���յ��ķ���
dff=zeros(1,N);   %����֮��

for i=1:length(Es_N0_dB)
s=randsrc(1,N,s_i); %��������ŵȸ��ʳ���
s_d=(1/sqrt(5))*s;  %��һ���������ŵľ���,��Es=1
SER=zeros(1,length(Es_N0_dB)); %��Ԥ���SER
n=randn(1,N);       %����N�����������ֵΪ0������Ϊ1
s_n=sqrt(10^(-Es_N0_dB(i)/10))*1/sqrt(2)*n;   %���������Ϻ��ʵķ���
y=s_d+s_n;          %��ÿһ�����Ŷ���������
%%%*******�������*******%�����ж�***����Ƿ��沿�֣�Ҳ�����ø���ȥ���
s_r(find(y<(-2/sqrt(5))))=-3;
s_r(find(y>=(2/sqrt(5))))=3;
s_r(find(y>=(-2/sqrt(5))& y<0))=-1;
s_r(find(y<(2/sqrt(5))&y>=0))=1; 
%%%*****��s��s_r��������ط�0������Ȼ����find�����ҳ����ǵ�λ��********%
dff=abs(s-s_r);
err(i)=size(find(dff),2);
end
%%%�ø��ʷ����SER��Es_N0_dB�Ĺ�ϵ
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