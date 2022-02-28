function [I] = mutualinfo(px,delta,X)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
I=0;

    p_y_x = @(x,u) (1/sqrt(2*pi)*exp(-(x-u).^2./2));
    p_y = @(x) p_y_x(x,delta*X)*(px');
    
for i = 1:length(X)
    


    p_yxi =@(x) min(p_y_x(x,delta*X(i)) * px(i) .* (log2(p_y_x(x,delta*X(i)))-log2(p_y(x))),0) ...
              + max(p_y_x(x,delta*X(i)) * px(i) .* (log2(p_y_x(x,delta*X(i)))-log2(p_y(x))),0);
%min() + max() 的原因是当值过小时会产生NaN,该方法可以去掉NaN的同时保留正负值
%对数拆开计算可以避免一些inf值的产生

%     I = I + integral(p_yxi,-inf,inf,'ArrayValued',true);
    I = I + integral(p_yxi,delta*X(i)-abs(10*delta*X(i)),delta*X(i)+abs(10*delta*X(i)),'ArrayValued',true);
end

end


