%����2��3��8��16ASM�����Ż���Ϣ��������ƽ������Ϣ����Ƚ�

addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
m_ASK = 3;              %ASK����
SNR_dB =0:1:35;     
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
        delta = sqrt(P/(ASK_symbol.^2*PX(j,:)'));
        I(j) = mutualinfo(PX(j,:),delta,ASK_symbol);
        I_mean(j) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol);
        C(j) = log2(1+P)/2;
    end
    plot(SNR_dB,I,'-r','LineWidth',1);plot(SNR_dB,I_mean,'k','LineWidth',1);
    clear PX 
end
plot(SNR_dB,C,'--k','LineWidth',2);hold off;
    parfor j =1:length(SNR_dB)
        I_QAM_mean(j) = mutualinfo_QAM(ones(8,8)/(8^2),SNR_dB(j),8^2);
    end