function [I] = mutualinfo(px,delta,X)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
I=0;

for i = 1:length(X)
    
    p_y_x = @(x,u) (1/sqrt(2*pi)*exp(-(x'-u).^2./2));
    p_y = @(x) p_y_x(x,delta*X)*(px');

    p_yxi =@(x) min(p_y_x(x,delta*X(i)) * px(i) .* log2(p_y_x(x,delta*X(i))./p_y(x)),0) + max(p_y_x(x,delta*X(i)) * px(i) .* log2(p_y_x(x,delta*X(i))./p_y(x)),0);
%min() + max() 的原因是去掉NaN，同时保留正负值

%     I = I + integral(p_yxi,-inf,inf,'ArrayValued',true);
    I = I + integral(p_yxi,-abs(delta*X(i)*10),abs(delta*X(i)*10),'ArrayValued',true);
end
end


