%LDPC译码仿真

clc;clear;
%仿真参数设置
% EbN0_dB = 0:0.5:10;
%QAM参数
M = 16;
Es = 1;

% EsN0_dB = EbN0_dB + 10*log10(log2(M));
EsN0_dB = 0:0.5:10;
EsN0=10.^(EsN0_dB/10);

N0=Es./EsN0;
variance=N0;
sigma=sqrt(variance);
%仿真次数
Sim_times = 1e3;
%码长，可在648,1296,1944中选择
N =648;
N_qam = N/log2(M);
%编码码率，可在1/2,2/3,3/4中选择
R = 1/2;
%LDPC的H矩阵，使用802.11n方法构造，来自https://github.com/simgunz/802.11n-ldpc
H =sparse(buildH( N, R));
%初始化LDPC编解码函数，系统自带
ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H);

for i = 1:length(EsN0_dB)
    errbit_count(i) = 0;
    for j = 1:Sim_times
        %产生随机比特
        txBits = randi(2,1,N*R)-1;
        %LDPC编码
        enc_bits = ldpcEncoder(txBits');
        %QAM调制
        X_mod=qammod(enc_bits,M,'Gray','UnitAveragePower',true,'InputType','bit'); 
        %AWGN噪声
        y = X_mod + sigma(i).*(randn(N_qam,1)+1i*randn(N_qam,1))./sqrt(2);
        %QAM解调
        y_demod=qamdemod(y,M,'Gray','UnitAveragePower',true,'OutputType','approxllr','NoiseVariance',variance); 
        %LDPC译码
        dec_bits = ldpcDecoder(y_demod);
        errbit_count(i) = errbit_count(i) +sum(dec_bits~=txBits');
    end
    
end
biterr_rate = errbit_count./(N*R*Sim_times);
semilogy(EsN0_dB,biterr_rate);