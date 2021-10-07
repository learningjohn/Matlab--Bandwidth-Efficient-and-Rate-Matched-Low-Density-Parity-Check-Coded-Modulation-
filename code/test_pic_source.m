%使用图像信源与随机信源比较信道互信息熵
addpath("Compute_fun\");
addpath("mat_data\");
clc;clear;
% test_source = imread('lena512.bmp');
test_source = imread('lena512color.tiff');
test_source_1_line = test_source(:,1);
test_source_bit = double(de2bi(test_source_1_line));
test_source_bit_1_line = reshape(test_source_bit,[],1);
M = 64;
QamSymbol_pic= qammod(test_source_bit_1_line(1:682*log2(M)),M,'InputType','Bit');
QamSymbol_pic_I = real(QamSymbol_pic);
QamSymbol_pic_Q = imag(QamSymbol_pic);
ASK_symbol=unique(QamSymbol_pic_I);
[T_I,~]  = hist(QamSymbol_pic_I,ASK_symbol);
P_I = T_I/length(QamSymbol_pic_I);
meanConstPower_I = sum(abs(ASK_symbol').^2.*P_I);

[T_Q,~]  = hist(QamSymbol_pic_Q,ASK_symbol);
P_Q = T_Q/length(QamSymbol_pic_Q);
meanConstPower_Q = sum(abs(ASK_symbol').^2.*P_Q);

SNR_dB =0:0.5:25;      
max_pow = mean(abs(ASK_symbol).^2);

for i =1:length(SNR_dB)
         
P=10.^(SNR_dB(i)/10);  %信号功率，默认噪声功率为1
count = 0;

I_mean(i) = mutualinfo(ones(1,length(ASK_symbol))/length(ASK_symbol),sqrt(P/max_pow),ASK_symbol')*2;
I_I(i) = mutualinfo(P_I,sqrt(P/meanConstPower_I),ASK_symbol');
I_Q(i) = mutualinfo(P_Q,sqrt(P/meanConstPower_Q),ASK_symbol');
I(i) = I_I(i)+I_Q(i);
C(i) = log2(1+P);

end
figure()
plot(SNR_dB,C);hold on;plot(SNR_dB,I);plot(SNR_dB,I_mean);
legend('信道容量','图像分布互信息量（按位取）','平均分布互信息量')
% legend('信道容量','图像分布互信息量（按像素点取）','平均分布互信息量')