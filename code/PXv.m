function [Px] = PXv(x,v,X)
%MB分布概率计算
%   此处显示详细说明
Px = exp(-v.*(x.^2))./sum(exp(-v.*(X.^2)));
end