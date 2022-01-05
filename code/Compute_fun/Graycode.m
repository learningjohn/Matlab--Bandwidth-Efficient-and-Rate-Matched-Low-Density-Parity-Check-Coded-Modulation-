function [Gray_X,Gray_X_bit] = Graycode(X,n)
%���ױ��룬���ظ��ױ����ı��غͶ�Ӧ����
%Graycode(X,n)  XΪ���������nΪ���λ������Ŀ����ָ��ʱĬ��Ϊlog2(max(x))
X = X-1;
n_max = floor(log2(max(X)))+1;
if nargin<2
    n = n_max;
end
bi_x = logical(de2bi(X,n,'left-msb'));
temp_bi_x = circshift(bi_x',1)';

temp_bi_x(:,1) = 1; %�����������ӳ�� 2��3��1��0
% temp_bi_x(:,1) = 0;%�����Ƹ���ӳ�� 0��1��3��2

Gray_X_bit = double(xor(bi_x,temp_bi_x));
Gray_X = bi2de(Gray_X_bit,'left-msb');
end

