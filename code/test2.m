%通过ASK的最优分布求得QAM的分布，并比较平均分布在AWGN信道下的性能
clc;clear;
addpath("Compute_fun\");
addpath("mat_data\");
load('best_P_8ASK.mat');
addpath('ccdm\');
EsN0dB =8:0.1:18;
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
    SER_mean(i) = QAM_theorySer(ones(sqrt(M),sqrt(M))/M,EsN0dB(i),M);
end
%maltab自带的理论误码率公式
[BER,SER] = berawgn(EsN0dB-10*log10(log2(M)),'qam',M);

figure()
semilogy(EsN0dB,SER,'-k');hold on;
semilogy(EsN0dB,SER_meanMB,'-b');
semilogy(EsN0dB,SER_best,'-r');hold off
ylabel('Ser');xlabel('信噪比（EsN0）')
legend('均匀分布理论误符号率','带入均匀分布公式的MB分布理论误符号率','MB分布理论误符号率','FontSize',10);