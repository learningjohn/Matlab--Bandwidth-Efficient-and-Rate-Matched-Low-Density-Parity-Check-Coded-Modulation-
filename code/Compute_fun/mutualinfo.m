function [I] = mutualinfo(px,delta,X)
%UNTITLED6 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I=0;

for i = 1:length(X)
    
    p_y_x = @(x,u) (1/sqrt(2*pi)*exp(-(x'-u).^2./2));
    p_y = @(x) p_y_x(x,delta*X)*(px');

    p_yxi =@(x) min(p_y_x(x,delta*X(i)) * px(i) .* log2(p_y_x(x,delta*X(i))./p_y(x)),0) ...
              + max(p_y_x(x,delta*X(i)) * px(i) .* log2(p_y_x(x,delta*X(i))./p_y(x)),0);
%min() + max() ��ԭ���ǵ�ֵ��Сʱ�����NaN,�÷�������ȥ��NaN��ͬʱ��������ֵ

%     I = I + integral(p_yxi,-inf,inf,'ArrayValued',true);
    I = I + integral(p_yxi,delta*X(i)-abs(10*delta*X(i)),delta*X(i)+abs(10*delta*X(i)),'ArrayValued',true);
end

end


