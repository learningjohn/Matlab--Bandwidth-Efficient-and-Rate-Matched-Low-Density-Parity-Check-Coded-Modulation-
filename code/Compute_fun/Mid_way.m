function [v_out] = Mid_way(delta,P,X)
%���ַ�����vֵ
%   �˴���ʾ��ϸ˵��
v_range = [0,10]; %���ַ���Χ
err_min = 0.0001;    %��С���
% Pv = exp(-v*(x^2))/sum(exp(-v.*(X^2)));
err=1;
while abs(err)>err_min
    v = (v_range(2)-v_range(1))/2  + v_range(1);
%     v = sum(v_range)/2;
    err = P/(delta^2) - sum(PXv(X,v,X).*(X.^2));
    if err>0
    v_range(2) = v;
    elseif err ==0
        v_out = v;
        return;
    else
        v_range(1) = v;
    end
end
v_out = v;
end