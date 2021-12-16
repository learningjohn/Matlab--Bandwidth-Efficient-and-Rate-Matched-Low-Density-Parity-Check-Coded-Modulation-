%ͨ��ASK�����ŷֲ����QAM�ķֲ������Ƚ�ƽ���ֲ���AWGN�ŵ��µ�����
clc;clear;
addpath("Compute_fun\");
addpath("mat_data\");
load('best_P_8ASK.mat');
addpath('ccdm\');
EsN0dB =0:1:25;
M=64;
QAM_symbol = qammod(0:M-1,M,'bin','UnitAveragePower',true).';
nSyms = 100000;
EsN0=10.^(EsN0dB/10);
Es=1;
N0=Es./EsN0;
variance=N0;
Standard_variance=sqrt(variance);
Fer1 = zeros(1,length(EsN0));
Fer2= zeros(1,length(EsN0));
QAM_dmat = (repmat([-sqrt(M)+1:2:sqrt(M)-1],sqrt(M),1)./2).^2 + (repmat([sqrt(M)-1:-2:-sqrt(M)+1]',1,sqrt(M))./2).^2;
p_test = PX(13,:)'*PX(13,:);
for i =1:length(EsN0dB)
%     pOpt = PX(i,:)'*PX(i,:);
    pOpt = p_test;
    pOpt = pOpt(:);
    [p_quant,nBitsInfo,n_i] = ccdm.initialize(pOpt,nSyms);
    
    meanConstPower = sum(abs(QAM_symbol).^2.*p_quant);
    txBits = randi(2,1,nBitsInfo)-1;
    i_TX=ccdm.encode(txBits(1,:),n_i).'+1;
    txSyms = i_TX-1;
    IQ = QAM_symbol(i_TX);
    Es = meanConstPower;
    N0=Es./EsN0(i);
    variance=N0;
    Standard_variance1=sqrt(variance);
    n=(randn(1,nSyms).'+1i*randn(1,nSyms).')./sqrt(2);
    n_w=Standard_variance1*n;
    y = IQ+n_w;
    R_C = qamdemod(y,M,'bin','UnitAveragePower',true);
    Fer1(i) = sum(R_C ~= txSyms)/nSyms;
    
    src_symbols_hat = ccdm.decode(R_C',n_i,nBitsInfo);
    Ber1(i) = sum(src_symbols_hat ~= txBits)/length(src_symbols_hat);
    
    SER_best(i) = QAM_theorySer(p_test,EsN0dB(i),M);
    
end
[BER,SER] = berawgn(EsN0dB-10*log10(log2(M)),'qam',M);

figure()
semilogy(EsN0dB,Fer1,'-*k');hold on;
semilogy(EsN0dB,SER,'b');
semilogy(EsN0dB,SER_best,'-+r');hold off
ylabel('��֡��');xlabel('����ȣ�EsN0��')
legend('MB�ֲ�����������','���ȷֲ�����������','MB�ֲ�����������')

%
figure()
semilogy(EsN0dB,BER,'-*k');hold on;
semilogy(EsN0dB,Ber1,'b');hold off
ylabel('��֡��');xlabel('����ȣ�EsN0��')
legend('MB�ֲ�����������','���ȷֲ�����������')

acit_ASK = 1:2:sqrt(M);
ASK_symbol = [-acit_ASK(end:-1:1),acit_ASK];%ASK������
max_pow = mean(abs(ASK_symbol).^2);
for j =1:length(EsN0dB)
    P=10.^(EsN0dB(j)/10);  %�źŹ��ʣ�Ĭ����������Ϊ1
    delta = sqrt(P/(ASK_symbol.^2*PX(j,:)'));
    I(j) = mutualinfo(PX(j,:),delta,ASK_symbol);
    I_mean(j) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol);
    C(j) = log2(1+P)/2;
end
figure()
hold on;
plot(EsN0dB,I*2,'-r','LineWidth',1);
plot(EsN0dB,I_mean*2,'k','LineWidth',1);
plot(EsN0dB,C*2,'--k','LineWidth',2);hold off;
xlabel('EsN0');ylabel('Rate');title('64QAM')
legend('���ŷֲ�����Ϣ��','ƽ���ֲ�����Ϣ��','�ŵ�����');