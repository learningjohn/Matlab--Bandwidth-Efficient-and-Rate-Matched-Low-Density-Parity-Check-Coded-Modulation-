function [Gray_X,Gray_X_bit] = Graycode(X,n)
%格雷编码，返回格雷编码后的比特和对应整数
%Graycode(X,n)  X为输入变量，n为最高位比特数目，不指定时默认为log2(max(x))
X = X-1;
n_max = floor(log2(max(X)))+1;
if nargin<2
    n = n_max;
end
bi_x = logical(de2bi(X,n,'left-msb'));
temp_bi_x = circshift(bi_x',1)';

temp_bi_x(:,1) = 1; %二进制逆格雷映射 2，3，1，0
% temp_bi_x(:,1) = 0;%二进制格雷映射 0，1，3，2

Gray_X_bit = double(xor(bi_x,temp_bi_x));
Gray_X = bi2de(Gray_X_bit,'left-msb');
end

