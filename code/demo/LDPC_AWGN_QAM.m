%LDPC�������

clc;clear;
%�����������
% EbN0_dB = 0:0.5:10;
%QAM����
M = 16;
Es = 1;

% EsN0_dB = EbN0_dB + 10*log10(log2(M));
EsN0_dB = 0:0.5:10;
EsN0=10.^(EsN0_dB/10);

N0=Es./EsN0;
variance=N0;
sigma=sqrt(variance);
%�������
Sim_times = 1e3;
%�볤������648,1296,1944��ѡ��
N =648;
N_qam = N/log2(M);
%�������ʣ�����1/2,2/3,3/4��ѡ��
R = 1/2;
%LDPC��H����ʹ��802.11n�������죬����https://github.com/simgunz/802.11n-ldpc
H =sparse(buildH( N, R));
%��ʼ��LDPC����뺯����ϵͳ�Դ�
ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H);

for i = 1:length(EsN0_dB)
    errbit_count(i) = 0;
    for j = 1:Sim_times
        %�����������
        txBits = randi(2,1,N*R)-1;
        %LDPC����
        enc_bits = ldpcEncoder(txBits');
        %QAM����
        X_mod=qammod(enc_bits,M,'Gray','UnitAveragePower',true,'InputType','bit'); 
        %AWGN����
        y = X_mod + sigma(i).*(randn(N_qam,1)+1i*randn(N_qam,1))./sqrt(2);
        %QAM���
        y_demod=qamdemod(y,M,'Gray','UnitAveragePower',true,'OutputType','approxllr','NoiseVariance',variance); 
        %LDPC����
        dec_bits = ldpcDecoder(y_demod);
        errbit_count(i) = errbit_count(i) +sum(dec_bits~=txBits');
    end
    
end
biterr_rate = errbit_count./(N*R*Sim_times);
semilogy(EsN0_dB,biterr_rate);