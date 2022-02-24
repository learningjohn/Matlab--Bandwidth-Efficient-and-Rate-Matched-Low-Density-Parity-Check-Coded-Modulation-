function [Li] = calculate_ibit_llr(Rx_I,ASK_of_biti0,ASK_of_biti1,Standard_variance1,P_ASK_of_biti0,P_ASK_of_biti1)
%计算ASK接收信号,使用文中公式（60）（61）（62）编写
%   此处显示详细说明
num_of_llr = length(Rx_I);
%先验信息部分
%公式（61）
PBi0 = sum(P_ASK_of_biti0);
PBi1 = sum(P_ASK_of_biti1);
pri_info = log(PBi0/PBi1);
%信道似然比部分
%公式（62）
 for n = 1:num_of_llr

     Pybi0 = sum(exp(-(Rx_I(n)-ASK_of_biti0).^2./(2*Standard_variance1^2)).*P_ASK_of_biti0)/PBi0;
     Pybi1 = sum(exp(-(Rx_I(n)-ASK_of_biti1).^2./(2*Standard_variance1^2)).*P_ASK_of_biti1)/PBi1;
     Li(n,1) = pri_info + log(Pybi0/Pybi1);
 end 

end

