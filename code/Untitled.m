%不同信噪比下，计算ASK的I(x,y)最大的MB分布
addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
m_ASK = 4;              %ASK长度
acit_ASK = 1:2:2^(m_ASK-1);
ASK_symbol = [-acit_ASK(end:-1:1),acit_ASK];%ASK星座点
max_pow = mean(abs(ASK_symbol).^2);
min_pow = 1;
% delta_range = [0,10]; %黄金分割法范围
err_min = 0.000000001;    %最小误差
a = 0.618;b=1-a;
err = 1;
SNR_dB =5;
P=10.^(SNR_dB/10);  %信号功率，默认噪声功率为0
count = 0;
delta_range = [sqrt(P/max_pow),sqrt(P/min_pow)];
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
    count = count+1
    
end
delta = (delta_left+delta_right)/2;
v= Mid_way(delta,P,ASK_symbol);
PX = PXv(ASK_symbol,v,ASK_symbol);
I = mutualinfo(PX,delta_left,ASK_symbol);
C = log2(1+P)/2
delta  = sqrt(P/max_pow):0.001:sqrt(P/min_pow)
for i  = 1:length(delta)
    v(i)= Mid_way(delta(i),P,ASK_symbol);
    PX = PXv(ASK_symbol,v(i),ASK_symbol);
    I(i) = mutualinfo(PX,delta(i),ASK_symbol);
    i
end