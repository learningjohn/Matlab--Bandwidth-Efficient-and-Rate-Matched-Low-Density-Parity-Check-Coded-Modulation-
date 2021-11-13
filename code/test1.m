%通过ASK的最优分布求得QAM的分布，并比较平均分布在AWGN信道下的性能
clc;clear;
addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
load('best_P_8ASK.mat');
addpath('ccdm\');
EsN0dB =0:0.5:25; 
M=64;
C = qammod(0:M-1,M,'bin','UnitAveragePower',true).';
nSyms = 100000;
EsN0=10.^(EsN0dB/10);
Es=1;
N0=Es./EsN0;
variance=N0;
Standard_variance=sqrt(variance);
Fer1 = zeros(1,length(EsN0));
Fer2= zeros(1,length(EsN0));
QAM_dmat = (repmat([-sqrt(M)+1:2:sqrt(M)-1],sqrt(M),1)./2).^2 + (repmat([sqrt(M)-1:-2:-sqrt(M)+1]',1,sqrt(M))./2).^2;
for i =1:length(EsN0dB)
    pOpt = PX(i,:)'*PX(i,:);
    pOpt = pOpt(:);
    [p_quant,nBitsInfo,n_i] = ccdm.initialize(pOpt,nSyms);
    
    meanConstPower = sum(abs(C).^2.*p_quant);
    txBits = randi(2,1,nBitsInfo)-1;
    i_TX=ccdm.encode(txBits(1,:),n_i).'+1;
    txSyms = i_TX-1;
    IQ = C(i_TX);
    Es = meanConstPower;
    N0=Es./EsN0(i);
    variance=N0;
    Standard_variance1=sqrt(variance);
    n=(randn(1,nSyms).'+1i*randn(1,nSyms).')./sqrt(2);
    n_w=Standard_variance1*n;
    y = IQ+n_w;
    R_C = qamdemod(y,M,'bin','UnitAveragePower',true);
    Fer1(i) = sum(R_C ~= txSyms)/nSyms;

    src_symbols_hat = ccdm.decode(txSyms',n_i,nBitsInfo);
    Ber1(i) = sum(src_symbols_hat ~= txBits)/length(src_symbols_hat);
    
    SER_best(i) = QAM_theorySer(PX(i,:)'*PX(i,:),EsN0dB(i),M);
    
end
[BER,SER] = berawgn(EsN0dB-10*log10(log2(M)),'qam',M);

figure()
semilogy(EsN0dB,Fer1,'-*k');hold on;
semilogy(EsN0dB,SER,'b');
semilogy(EsN0dB,SER_best,'-+r');hold off
ylabel('误帧率');xlabel('信噪比（EsN0）')
legend('MB分布仿真误码率','均匀分布理论误码率','MB分布理论误码率')

% 
% figure()
% semilogy(EsN0dB,BER,'-*k');hold on;
% semilogy(EsN0dB,Ber1,'b');hold off
% ylabel('误帧率');xlabel('信噪比（EsN0）')
% legend('MB分布仿真误码率','均匀分布理论误码率')