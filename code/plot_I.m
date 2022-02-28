%该函数用来绘制信道容量，均匀分布互信息熵和最优分布互信息熵
clc;clear;
M_ASK = [4,8,16,32,64];   %ASK长度
for i =1:length(M_ASK)
    ASK_half = 1:2:(M_ASK(i));
    ASK_symbol = [-ASK_half(end:-1:1),ASK_half];%ASK星座点
    mean_pow = mean(abs(ASK_symbol).^2);
    snr_dB = 0:1:40;
    best_P = [];
    for j =1:length(snr_dB)
        best_P(j,:) = initialize_PX(snr_dB(j),M_ASK(i));
        P=10.^(snr_dB(j)/10);  %信号功率，默认噪声功率为1
        delta = sqrt(P/(ASK_symbol.^2*best_P(j,:)'));
        I(i,j) = mutualinfo(best_P(j,:),delta,ASK_symbol);
        I_mean(i,j) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/mean_pow),ASK_symbol);
        C(j) = log2(1+P)/2;
    end
    
end
figure()
hold on;
plot(snr_dB,C,'--k','LineWidth',2);
plot(snr_dB,I_mean(1,:),'k','LineWidth',1);
plot(snr_dB,I_mean(2,:),'k','LineWidth',1);
plot(snr_dB,I_mean(3,:),'k','LineWidth',1);
plot(snr_dB,I_mean(4,:),'k','LineWidth',1);
plot(snr_dB,I_mean(5,:),'k','LineWidth',1);
xlabel('SNR');ylabel('Rate');
legend('信道容量','平均分布互信息量4-ASK','8-ASK','16-ASK','32-ASK','64-ASK');grid on;