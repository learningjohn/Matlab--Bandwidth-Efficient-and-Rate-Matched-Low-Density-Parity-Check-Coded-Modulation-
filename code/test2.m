%通过ASK的最优分布求得QAM的分布，并比较平均分布在AWGN信道下的性能
clc;clear;


load('best_P_8ASK.mat');

EsN0dB =0:1:25;
M=64;
QAM_symbol = qammod(0:M-1,M,'bin','UnitAveragePower',true).';
nSyms = 100000;
EsN0=10.^(EsN0dB/10);
Fer1 = zeros(1,length(EsN0));
Fer2= zeros(1,length(EsN0));

p_test = PX(13,:)'*PX(13,:);
QAM_dmat = (repmat([-sqrt(M)+1:2:sqrt(M)-1],sqrt(M),1)./2).^2 + (repmat([sqrt(M)-1:-2:-sqrt(M)+1]',1,sqrt(M))./2).^2;

A_mean = (M-1)/6;
A_test = sum(sum(p_test.*QAM_dmat));


for i =1:length(EsN0dB)
    
    %均匀分布理论公式计算式
    Px_theoryfun = @(SNR,A,M) 4*(1-1/sqrt(M))*(qfunc(sqrt(SNR/(2*A)))-qfunc(sqrt(SNR/(2*A)))^2);
    %均匀分布理论误符号率
    SER_mean(i) = Px_theoryfun(EsN0(i),A_mean,M);
    %带入均匀分布公式的MB分布理论误符号率
    SER_meanMB(i) = Px_theoryfun(EsN0(i),A_test,M);
    %带入推导公式的MB分布理论误符号率
    SER_best(i) = QAM_theorySer(p_test,EsN0dB(i),M);
    %均匀分布带入推导公式计算
    SER_mean2(i) = QAM_theorySer(ones(sqrt(M),sqrt(M))/M,EsN0dB(i),M);
end
%maltab自带的理论误码率公式
[BER,SER] = berawgn(EsN0dB-10*log10(log2(M)),'qam',M);
%绘制分布
subplot(1,2,1)
bar3(ones(8,8)/64);xlabel('x');ylabel('y');zlabel('p');
title('64QAM平均分布,A=10.5');
subplot(1,2,2)
bar3(PX(13,:)'*PX(13,:));xlabel('x');ylabel('y');zlabel('p');
title('64QAM在13dB下的M-B分布，A=4.6644')
%绘制不同分布理论误码率
figure()
semilogy(EsN0dB,SER,'-*k','LineWidth',2,'MarkerSize',10);hold on;
semilogy(EsN0dB,SER_meanMB,'-ob','LineWidth',2,'MarkerSize',10);
semilogy(EsN0dB,SER_best,'-+r','LineWidth',2,'MarkerSize',10);hold off
ylabel('Ser');xlabel('信噪比（EsN0）')
legend('均匀分布理论误符号率','带入均匀分布公式的MB分布理论误符号率','MB分布理论误符号率','FontSize',10);
%绘制平均分布理论误码率
figure()
semilogy(EsN0dB,SER,'-*k','LineWidth',2,'MarkerSize',10);hold on;
semilogy(EsN0dB,SER_mean,'-ob','LineWidth',2,'MarkerSize',10);
semilogy(EsN0dB,SER_mean2,'-+r','LineWidth',2,'MarkerSize',10);hold off
ylabel('Ser');xlabel('信噪比（EsN0）')
legend('matlab自带误符号率','当前公式误符号率','前面报告误符号率','FontSize',10);
title('64QAM均匀分布');

