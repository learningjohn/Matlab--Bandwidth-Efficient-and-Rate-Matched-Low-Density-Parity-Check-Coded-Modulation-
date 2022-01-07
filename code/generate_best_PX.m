%��ͬ������£����㲻ͬ���ƽ���ASK��I(x,y)���ŵ�MB�ֲ�
%���н����󱣴�ֲ�PX
addpath("Compute_fun\");
addpath("mat_data\");
addpath("ccdm\");
clc;clear;
m_ASK = 32;              %ASK����
ASK_half = 1:2:(m_ASK); 
ASK_symbol = [-ASK_half(end:-1:1),ASK_half];%ASK������
max_pow = mean(abs(ASK_symbol).^2);
min_pow = 1;
% delta_range = [0,10]; %�ƽ�ָ��Χ
err_min = 0.000000001;    %��С���
a = 0.618;b=1-a;    
err = 1;
SNR_dB =15:1:28;      
for i =1:length(SNR_dB)
         
P=10.^(SNR_dB(i)/10);  %�źŹ��ʣ�Ĭ����������Ϊ1
count = 0;
delta_range = [sqrt(P/max_pow),sqrt(P/min_pow)];
err = 1;
while abs(err)>err_min
    delta_left =delta_range(1)+ (delta_range(2)-delta_range(1))*b;
    delta_right = delta_range(1) + (delta_range(2)-delta_range(1))*a;
    v_1 = Mid_way(delta_left,P,ASK_symbol);
    v_2 = Mid_way(delta_right,P,ASK_symbol);
    if any([v_1,v_2] < 0)
        flag = 1;
    end
    PX_left = PXv(ASK_symbol,v_1,ASK_symbol);
    PX_right = PXv(ASK_symbol,v_2,ASK_symbol); 
    I_left = mutualinfo(PX_left,delta_left,ASK_symbol);
    I_right = mutualinfo(PX_right,delta_right,ASK_symbol);
    err = I_left - I_right;
    if err>0
        delta_range(2) = delta_right;
   
    else
         delta_range(1) = delta_left;
    end
    count = count+1

end
delta = (delta_left+delta_right)/2;
v= Mid_way(delta,P,ASK_symbol);
PX(i,:) = PXv(ASK_symbol,v,ASK_symbol);
I(i) = mutualinfo(PX(i,:),delta,ASK_symbol);
H(i) = -(log2(PX(i,:))*PX(i,:)');
H1(i) = -(log2(2*PX(i,1:length(PX(i,:))/2))*2*PX(i,1:length(PX(i,:))/2)');
I_mean(i) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol);
C(i) = log2(1+P)/2;
[p_quant,nBitsInfo,n_i] = ccdm.initialize(PX(i,:),10000);
R(i) = nBitsInfo/(10000);
end
figure()
hold on;grid on;
plot(SNR_dB,C);
plot(SNR_dB,I);
% plot(SNR_dB,I_mean);
plot(SNR_dB,H1);
plot(SNR_dB,H1+1/6);
% plot(SNR_dB,R*2/3);
legend('�ŵ�����','���ŷֲ�����Ϣ��','2/3������Ϣ��','3/4������Ϣ��')
hold off;
%  filename_bestP = ['best_P_',num2str(m_ASK),'ASK.mat'];
%  save(filename_bestP,'PX')
%snr_Px = [SNR_dB',PX];
% filename_snrbestP = ['snr_bestP_',num2str(m_ASK),'ASK.mat'];
%save(filename_snrbestP,'snr_Px')
%������ϵõ�������,32ASK,��������Ϊ5/6ʱ��ʵ�����ʶ�Ӧ�������
R_SNR = [3,15.3;3.1,16.4;3.2,17.6;3.3,18.8;3.4,20;3.5,21.1;3.6,22.3;3.7,23.5;3.8,24.7;3.9,26.1;];