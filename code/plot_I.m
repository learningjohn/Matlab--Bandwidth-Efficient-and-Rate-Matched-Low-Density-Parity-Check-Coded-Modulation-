addpath("Compute_fun\");
addpath("mat_data\");
addpath('ccdm\');
%�ú������������ŵ����������ȷֲ�����Ϣ�غ����ŷֲ�����Ϣ��
clc;clear;
M_ASK = [4,8,16,32,64];   %ASK����
% snr_dB = 10:30;
snr_dB = 6:16;
for i =1:length(M_ASK)
    ASK_half = 1:2:(M_ASK(i));
    ASK_symbol = [-ASK_half(end:-1:1),ASK_half];%ASK������
    mean_pow = mean(abs(ASK_symbol).^2);
    
    best_P = [];
    for j =1:length(snr_dB)
        best_P(j,:) = initialize_PX(snr_dB(j),M_ASK(i));
        P=10.^(snr_dB(j)/10);  %�źŹ��ʣ�Ĭ����������Ϊ1
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
legend('�ŵ�����','ƽ���ֲ�����Ϣ��');grid on;
hold off;

figure()
hold on;
plot(snr_dB,C,'--k','LineWidth',2);

plot(snr_dB,I_mean(3,:),'k','LineWidth',1);
plot(snr_dB,I(3,:),'r','LineWidth',1);

plot(snr_dB,H(3,:),'k-.','LineWidth',1);
plot(snr_dB,H_MB(3,:),'r-.','LineWidth',1);
xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');
legend('�ŵ�����','ƽ���ֲ�����Ϣ��','MB���ŷֲ�����Ϣ','���ȷֲ�����Ϣ��','MB�ֲ�����Ϣ��');grid on;
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
legend('�ŵ�����','16-ASK M-B���ŷֲ�I(X;Y)','16-ASK M-B���ŷֲ�H(X)','16-ASK PAS����R','8-ASK M-B���ŷֲ�I(X;Y)')
hold off

figure()
hold on;
plot(snr_dB,C,'k','LineWidth',2);
plot(snr_dB,I(2,:),'--r','LineWidth',1);
% plot(snr_dB,H_MB(3,:),'r-.','LineWidth',1);
plot(snr_dB,H_PAS(2,:)+1/4,'g-.','LineWidth',1);
plot(snr_dB,H_PAS(2,:),'b-.','LineWidth',1);


xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');grid on;
legend('�ŵ�����','I(X;Y)','PAS����H(A)+1/4','PAS����H(A)')
hold off