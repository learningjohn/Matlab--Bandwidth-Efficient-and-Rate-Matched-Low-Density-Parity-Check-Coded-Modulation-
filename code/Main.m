%��ͬ������£�����ASK��I(x,y)����MB�ֲ�
%
addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
m_ASK = 3;              %ASK����
acit_ASK = 1:2:2^(m_ASK); 
ASK_symbol = [-acit_ASK(end:-1:1),acit_ASK];%ASK������
max_pow = mean(abs(ASK_symbol).^2);
min_pow = 1;
% delta_range = [0,10]; %�ƽ�ָ��Χ
err_min = 0.000000001;    %��С���
a = 0.618;b=1-a;    
err = 1;
SNR_dB =0:1:30;      
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
I_mean(i) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol);
C(i) = log2(1+P)/2;
end
figure()
plot(SNR_dB,C);hold on;plot(SNR_dB,I);plot(SNR_dB,I_mean);
legend('�ŵ�����','���ŷֲ�����Ϣ��','ƽ���ֲ�����Ϣ��')
filename_bestP = ['best_P_',num2str(2^m_ASK),'ASK.mat'];
save(filename_bestP,'PX')
