function [Ser] = QAM_theorySer(P,EsN0_dB,M)
%计算QAM理论误符号率
%输入变量解释
%P: QAM信号各个星座点的概率，矩阵大小应与QAM星座图对应
%EsN0_dB： Es/N0的dB形式
%M: QAM点数
%输出变量 Ser:当前的理论误符号率


    EsN0=10.^(EsN0_dB/10);
    %M-QAM信号各个点在坐标轴上的位置，以d为坐标点单位
    QAM_dmat = (repmat([-sqrt(M)+1:2:sqrt(M)-1],sqrt(M),1)./2).^2 + (repmat([sqrt(M)-1:-2:-sqrt(M)+1]',1,sqrt(M))./2).^2;
    %计算M-QAM在当前概率P下的取值
    A =sum(sum( QAM_dmat.*P));
    ErrA = erfc(sqrt(EsN0(i)/(4*A)))-(0.5*erfc(sqrt(EsN0(i)/(4*A))))^2;
    ErrB = (3/2)*erfc(sqrt(EsN0(i)/(4*A)))-2*(0.5*erfc(sqrt(EsN0(i)/(4*A))))^2;
    ErrC = 2*erfc(sqrt(EsN0(i)/(4*A)))-4*(0.5*erfc(sqrt(EsN0(i)/(4*A))))^2;
    
    matrix_err = [[ErrA;ones(sqrt(M)-2,1)*ErrB;ErrA],repmat([ErrB;ones(sqrt(M)-2,1)*ErrC;ErrB],1,sqrt(M)-2),[ErrA;ones(sqrt(M)-2,1)*ErrB;ErrA]];
    SER= sum(sum(matrix_err.*P));
end

