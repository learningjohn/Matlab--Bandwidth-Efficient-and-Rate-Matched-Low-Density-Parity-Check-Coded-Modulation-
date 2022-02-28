addpath("Compute_fun\");
addpath("mat_data\");
addpath('ccdm\');
%�ú������������ŵ����������ȷֲ�����Ϣ�غ����ŷֲ�����Ϣ��
clc;clear;
M_ASK = [4,8,16,32,64];   %ASK����
for i =1:length(M_ASK)
    ASK_half = 1:2:(M_ASK(i));
    ASK_symbol = [-ASK_half(end:-1:1),ASK_half];%ASK������
    mean_pow = mean(abs(ASK_symbol).^2);
    snr_dB = 0:1:40;
    best_P = [];
    for j =1:length(snr_dB)
        best_P(j,:) = initialize_PX(snr_dB(j),M_ASK(i));
        P=10.^(snr_dB(j)/10);  %�źŹ��ʣ�Ĭ����������Ϊ1
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
plot(snr_dB,I(1,:),'r','LineWidth',1);
plot(snr_dB,I_mean(3,:),'k','LineWidth',1);
plot(snr_dB,I(3,:),'r','LineWidth',1);
plot(snr_dB,I_mean(5,:),'k','LineWidth',1);
plot(snr_dB,I(5,:),'r','LineWidth',1);
xlabel('SNR(dB)');ylabel('Rate(bit/s/hz)');
legend('�ŵ�����','ƽ���ֲ�����Ϣ��','MB���ŷֲ�����Ϣ');grid on;
hold off;