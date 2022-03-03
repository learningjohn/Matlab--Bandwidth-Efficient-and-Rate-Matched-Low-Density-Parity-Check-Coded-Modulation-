addpath("Compute_fun\");
addpath("mat_data\");
addpath('ccdm\');
%该函数用来绘制信道容量，均匀分布互信息熵和最优分布互信息熵
clc;clear;
M_ASK = [4,8,16,32,64];   %ASK长度
% snr_dB = 10:30;
snr_dB = 6:16;
for i =1:length(M_ASK)
    ASK_half = 1:2:(M_ASK(i));
    ASK_symbol = [-ASK_half(end:-1:1),ASK_half];%ASK星座点
    mean_pow = mean(abs(ASK_symbol).^2);
    
    best_P = [];
    for j =1:length(snr_dB)
        best_P(j,:) = initialize_PX(snr_dB(j),M_ASK(i));
        P=10.^(snr_dB(j)/10);  %信号功率，默认噪声功率为1
        delta = sqrt(P/(ASK_symbol.^2*best_P(j,:)'));
        I(i,j) = mutualinfo(best_P(j,:),delta,ASK_symbol);
        I_mean(i,j) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/mean_pow),ASK_symbol);
        C(j) = log2(1+P)/2;
        H_MB(i,j) = -(log2(best_P(j,:)))*best_P(j,:)';
        
        H_PAS(i,j) =-(log2(2*best_P(j,M_ASK(i)/2+1:end)))*2*best_P(j,M_ASK(i)/2+1:end)';
        H(i,j) = log2(M_ASK(i));
    end
    
end

figure()
hold on;
plot(snr_dB,C,'--k','LineWidth',2);
plot(snr_dB,I_mean(1,:),'k','LineWidth',1);
plot(snr_dB,I_mean(2,:),'k','LineWidth',1);
plot(snr_dB,I_mean(3,:),'k','LineWidth',1);
plot(snr_dB,I_mean(4,:),'k','LineWidth',1);
plot(snr_dB,I_mean(5,:),'r','LineWidth',1);

xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');
legend('信道容量','平均分布互信息量');grid on;
hold off;

figure()
hold on;
plot(snr_dB,C,'--k','LineWidth',2);

plot(snr_dB,I_mean(3,:),'k','LineWidth',1);
plot(snr_dB,I(3,:),'r','LineWidth',1);

plot(snr_dB,H(3,:),'k-.','LineWidth',1);
plot(snr_dB,H_MB(3,:),'r-.','LineWidth',1);
xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');
legend('信道容量','平均分布互信息量','MB最优分布互信息','均匀分布自信息熵','MB分布自信息熵');grid on;
title('16-ASK');
hold off;

figure()
hold on;
plot(snr_dB,C,'--k','LineWidth',2);

% plot(snr_dB,H_MB(2,:),'r-.','LineWidth',1);
% plot(snr_dB,H_PAS(2,:),'b-*','LineWidth',1);

plot(snr_dB,I(3,:),'k','LineWidth',1);
plot(snr_dB,H_MB(3,:),'r-.','LineWidth',1);
plot(snr_dB,H_PAS(3,:),'b-.*','LineWidth',1);
plot(snr_dB,I(2,:),'y','LineWidth',1);
xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');grid on;
legend('信道容量','16-ASK M-B最优分布I(X;Y)','16-ASK M-B最优分布H(X)','16-ASK PAS速率R','8-ASK M-B最优分布I(X;Y)')
hold off

figure()
hold on;
plot(snr_dB,C,'k','LineWidth',2);
plot(snr_dB,I(2,:),'--r','LineWidth',1);
% plot(snr_dB,H_MB(3,:),'r-.','LineWidth',1);
plot(snr_dB,H_PAS(2,:)+1/4,'g-.','LineWidth',1);
plot(snr_dB,H_PAS(2,:),'b-.','LineWidth',1);


xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');grid on;
legend('信道容量','I(X;Y)','PAS速率H(A)+1/4','PAS速率H(A)')
hold off