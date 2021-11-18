addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
m_ASK = 1:4;              %ASK����
SNR_dB =0:1:30;     
figure()
hold on;
for i=1:length(m_ASK)
    acit_ASK = 1:2:2^(m_ASK(i)); 
    ASK_symbol = [-acit_ASK(end:-1:1),acit_ASK];%ASK������
    filename_bestP = ['best_P_',num2str(2^m_ASK(i)),'ASK.mat'];
    max_pow = mean(abs(ASK_symbol).^2);
    load(filename_bestP);%��ȡ��ǰ�����µ����Ÿ��ʷֲ�
    for j =1:length(SNR_dB)
        P=10.^(SNR_dB(j)/10);  %�źŹ��ʣ�Ĭ����������Ϊ1
        I(j) = mutualinfo(PX(i,:),delta,ASK_symbol);
        I_mean(j) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol);
    end
    clear(PX)
end