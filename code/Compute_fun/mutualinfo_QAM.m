function [I] = mutualinfo_QAM(px,SNR_dB,M)
%UNTITLED6 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I=0;
px = px(:);
qam_symbol = qammod(0:M-1,M,'bin');
Es = 10^(SNR_dB/10);
delta = sqrt(Es/(abs(qam_symbol).^2*px));
for i = 1:M
    
    p_y_x = @(x_r,x_i,u) (1/sqrt(2*pi)*exp(-(x_r-real(u)).^2)/2).*(1/sqrt(2*pi)*exp(-(x_i-imag(u)).^2)/2);
    %��������ʵ���ۼӲ���
    p_y = @(x_r,x_i) p_y_x(x_r,x_i,delta*qam_symbol)*(px);
    %���ػ��������ޣ�����-inf inf �������������⣬��Ҫ�ֶ���������
    %���ڶ��ػ��ֺ���integral2��֧���������������ֶ�����һά����
    inter_realmin = delta*real(qam_symbol(i))-5;
    inter_realmax = delta*real(qam_symbol(i))+5;
    inter_imagmin = delta*imag(qam_symbol(i))-5;
    inter_imagmax = delta*imag(qam_symbol(i))+5;
    %min() + max() ��ԭ���ǵ�ֵ��Сʱ�����NaN,�÷�������ȥ��NaN��ͬʱ��������ֵ
    p_yx_imag = @(x_i) integral(@(x_r) min(p_y_x(x_r,x_i,delta*qam_symbol(i)) * px(i) .* log2(p_y_x(x_r,x_i,delta*qam_symbol(i))./p_y(x_r,x_i)),0) ...
              + max(p_y_x(x_r,x_i,delta*qam_symbol(i)) * px(i) .* log2(p_y_x(x_r,x_i,delta*qam_symbol(i))./p_y(x_r,x_i)),0),...
              inter_realmin,inter_realmax,'ArrayValued',true);
    I(i) = integral(p_yx_imag,inter_imagmin,inter_imagmax,'ArrayValued',true);
end

end


