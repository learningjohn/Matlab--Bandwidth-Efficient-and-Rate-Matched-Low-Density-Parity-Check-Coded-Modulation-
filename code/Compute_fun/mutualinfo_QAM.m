function [I] = mutualinfo_QAM(px,SNR_dB,M)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
I=0;
px = px(:);
qam_symbol = qammod(0:M-1,M,'bin');
Es = 10^(SNR_dB/10);
delta = sqrt(Es/(abs(qam_symbol).^2*px));
for i = 1:M
    
    p_y_x = @(x_r,x_i,u) (1/sqrt(2*pi)*exp(-(x_r-real(u)).^2)/2).*(1/sqrt(2*pi)*exp(-(x_i-imag(u)).^2)/2);
    %向量化，实现累加部分
    p_y = @(x_r,x_i) p_y_x(x_r,x_i,delta*qam_symbol)*(px);
    %二重积分上下限，由于-inf inf 的区间粒度问题，需要手动设置区间
    %由于二重积分函数integral2不支持向量化操作，手动两次一维积分
    inter_realmin = delta*real(qam_symbol(i))-5;
    inter_realmax = delta*real(qam_symbol(i))+5;
    inter_imagmin = delta*imag(qam_symbol(i))-5;
    inter_imagmax = delta*imag(qam_symbol(i))+5;
    %min() + max() 的原因是当值过小时会产生NaN,该方法可以去掉NaN的同时保留正负值
    p_yx_imag = @(x_i) integral(@(x_r) min(p_y_x(x_r,x_i,delta*qam_symbol(i)) * px(i) .* log2(p_y_x(x_r,x_i,delta*qam_symbol(i))./p_y(x_r,x_i)),0) ...
              + max(p_y_x(x_r,x_i,delta*qam_symbol(i)) * px(i) .* log2(p_y_x(x_r,x_i,delta*qam_symbol(i))./p_y(x_r,x_i)),0),...
              inter_realmin,inter_realmax,'ArrayValued',true);
    I(i) = integral(p_yx_imag,inter_imagmin,inter_imagmax,'ArrayValued',true);
end

end


