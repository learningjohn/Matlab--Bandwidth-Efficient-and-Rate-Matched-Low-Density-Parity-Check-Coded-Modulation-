clc;clear;
snr_dB = 18:0.25:26;
P=10.^(snr_dB/10);
C = log2(1+P)./2;
%PAS-32ASK-5/6
bitrate = [3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9];
err_snr = [18.82,19.37,20.13,20.81,21.4,21.97,22.6,23.25,23.85,24.5];
%16ASK
bitrate1 = [repmat(3,1,length(20.3:0.01:(21.35-0.01))),...
    repmat(3.2,1,length(21.35:0.01:(22.2-0.01))),...
    repmat(4*5/6,1,length(22.2:0.01:(23.75-0.01))),...
    repmat(4*8/9,1,length(23.75:0.01:(24.08-0.01))),...
    repmat(4*9/10,1,length(24.08:0.01:26))];

err_snr1 = [20.3:0.01:(21.35-0.01),21.35:0.01:(22.2-0.01),22.2:0.01:(23.75-0.01),23.75:0.01:(24.08-0.01),24.08:0.01:26];
figure()
hold on;
plot(snr_dB,C,'-.k','LineWidth',2);grid on;
plot(err_snr,bitrate,'-ob','LineWidth',1);
plot(err_snr1,bitrate1,'--k','LineWidth',1);
hold off;
xlabel('SNR');ylabel('Bit/Channel use (Bit/S/Hz)');
legend('信道容量','32-ASK,PAS整形信号，Fer=10^-3','16-ASK LDPC均匀分布信号，Fer=10^-3');