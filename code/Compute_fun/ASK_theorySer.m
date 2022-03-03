function [Ser] = ASK_theorySer(P,EsN0_dB,M)
%计算理论QAM误符号率，适用条件M>=4
%输入变量解释
%P: QAM信号各个星座点的概率，矩阵大小应与QAM星座图对应
%EsN0_dB： Es/N0的dB形式
%M: ASK点数
%输出变量 Ser:当前的理论误符号率


    EsN0=10.^(EsN0_dB/10);
    %M-QAM信号各个点在坐标轴上的位置，以d为坐标点单位
    ASK_dmat = [-(M-1):2:(M-1)]./2;
    %计算M-ASK在当前概率P下的A
    A =ASK_dmat.^2*P';
    
    Ser= [1/2,ones(1,length(P)-2),1/2]*P'* erfc(sqrt(EsN0/(4*A)));
    
    %Ser = (1-1/M)*erfc(sqrt(3/(M^2-1)*EsN0));
end

