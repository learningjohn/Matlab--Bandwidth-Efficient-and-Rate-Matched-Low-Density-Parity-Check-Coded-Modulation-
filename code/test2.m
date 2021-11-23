%ͨ��ASK�����ŷֲ����QAM�ķֲ������Ƚ�ƽ���ֲ���AWGN�ŵ��µ�����
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

    %���ȷֲ����۹�ʽ����ʽ
    Px_theoryfun = @(SNR,A,M) 4*(1-1/sqrt(M))*(qfunc(sqrt(SNR/(2*A)))-qfunc(sqrt(SNR/(2*A)))^2);
    %���ȷֲ������������
    SER_mean(i) = Px_theoryfun(EsN0(i),A_mean,M);
    %������ȷֲ���ʽ��MB�ֲ������������
    SER_meanMB(i) = Px_theoryfun(EsN0(i),A_test,M);
    %�����Ƶ���ʽ��MB�ֲ������������
    SER_best(i) = QAM_theorySer(p_test,EsN0dB(i),M);
    %���ȷֲ������Ƶ���ʽ����
    SER_mean(i) = QAM_theorySer(ones(sqrt(M),sqrt(M))/M,EsN0dB(i),M);
end
%maltab�Դ������������ʹ�ʽ
[BER,SER] = berawgn(EsN0dB-10*log10(log2(M)),'qam',M);

figure()
semilogy(EsN0dB,SER,'-k');hold on;
semilogy(EsN0dB,SER_meanMB,'-b');
semilogy(EsN0dB,SER_best,'-r');hold off
ylabel('Ser');xlabel('����ȣ�EsN0��')
legend('���ȷֲ������������','������ȷֲ���ʽ��MB�ֲ������������','MB�ֲ������������','FontSize',10);