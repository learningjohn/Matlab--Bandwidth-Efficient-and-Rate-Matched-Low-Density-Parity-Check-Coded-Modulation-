function [Px] = PXv(x,v,X)
%MB�ֲ����ʼ���
%   �˴���ʾ��ϸ˵��
Px = exp(-v.*(x.^2))./sum(exp(-v.*(X.^2)));
end