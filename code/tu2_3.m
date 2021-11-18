addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
m_ASK = 1:4;              %ASK长度
SNR_dB =0:1:30;     
figure()
hold on;
for i=1:length(m_ASK)
    acit_ASK = 1:2:2^(m_ASK(i)); 
    ASK_symbol = [-acit_ASK(end:-1:1),acit_ASK];%ASK星座点
    filename_bestP = ['best_P_',num2str(2^m_ASK(i)),'ASK.mat'];
    max_pow = mean(abs(ASK_symbol).^2);
    load(filename_bestP);%读取当前调制下的最优概率分布
    for j =1:length(SNR_dB)
        P=10.^(SNR_dB(j)/10);  %信号功率，默认噪声功率为1
        I(j) = mutualinfo(PX(i,:),delta,ASK_symbol);
        I_mean(j) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol);
    end
    clear(PX)
end