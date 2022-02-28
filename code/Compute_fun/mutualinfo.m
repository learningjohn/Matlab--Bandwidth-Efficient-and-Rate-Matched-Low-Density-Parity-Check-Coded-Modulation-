function [I] = mutualinfo(px,delta,X)
%ASK���ŵ��ŵ�����Ϣ�ؼ��㡣
%   �˴���ʾ��ϸ˵��
I=0;

    p_y_x = @(x,u) (1/sqrt(2*pi)*exp(-(x-u).^2./2));
    p_y = @(x) p_y_x(x,delta*X)*(px');
    
for i = 1:length(X)
    


    p_yxi =@(x) min(p_y_x(x,delta*X(i)) * px(i) .* (log2(p_y_x(x,delta*X(i)))-log2(p_y(x))),0) ...
              + max(p_y_x(x,delta*X(i)) * px(i) .* (log2(p_y_x(x,delta*X(i)))-log2(p_y(x))),0);
%min() + max() ��ԭ���ǵ�ֵ��Сʱ�����NaN,�÷�������ȥ��NaN��ͬʱ��������ֵ
%�����𿪼�����Ա���һЩinfֵ�Ĳ���

%     I = I + integral(p_yxi,-inf,inf,'ArrayValued',true);
    I = I + integral(p_yxi,delta*X(i)-abs(10*delta*X(i)),delta*X(i)+abs(10*delta*X(i)),'ArrayValued',true);
    %����integral����������м�����ȡһ���ʵ��Ļ������������߾���
end

end


