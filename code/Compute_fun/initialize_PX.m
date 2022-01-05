function [Px] = initialize_PX(snr,M_ASK)
%��ʼ����ǰ���������SNR��Ӧ�����ŷֲ�������ļ��洢�ĸ���û�ж�Ӧ������ȷֲ�
%�����¼��㲢���´洢�ļ�
%   �˴���ʾ��ϸ˵��
filename_bestP = ['snr_bestP_',num2str(M_ASK),'ASK.mat'];
load(['mat_data/',filename_bestP]);
snr_exit = snr_Px(:,1);
Px_exit = snr_Px(:,2:end);
ASK_symbol = (-M_ASK+1):2:(M_ASK-1);%ASK������
max_pow = mean(abs(ASK_symbol).^2);
min_pow = 1;
% delta_range = [0,10]; %�ƽ�ָ��Χ
err_min = 0.000000001;    %��С���
a = 0.618;b=1-a;
for i = 1:length(snr)
    if any(snr(i)==snr_exit)
        Px(i,:) = Px_exit(snr(i)==snr_exit,:);
    else

        SNR_dB =snr(i);
        P=10.^(SNR_dB/10);  %�źŹ��ʣ�Ĭ����������Ϊ1
        count = 0;
        delta_range = [sqrt(P/max_pow),sqrt(P/min_pow)];
        err = 1;
        while abs(err)>err_min
            delta_left =delta_range(1)+ (delta_range(2)-delta_range(1))*b;
            delta_right = delta_range(1) + (delta_range(2)-delta_range(1))*a;
            v_1 = Mid_way(delta_left,P,ASK_symbol);
            v_2 = Mid_way(delta_right,P,ASK_symbol);
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
        end
        delta = (delta_left+delta_right)/2;
        v= Mid_way(delta,P,ASK_symbol);
        Px(i,:) = PXv(ASK_symbol,v,ASK_symbol);
        snr_Px = [snr_Px;[snr(i),Px(i,:)]];

    end
end
save(['mat_data/',filename_bestP],'snr_Px');

end

