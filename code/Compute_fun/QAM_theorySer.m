function [Ser] = QAM_theorySer(P,EsN0_dB,M)
%��������QAM������ʣ���������M>=16
%�����������
%P: QAM�źŸ���������ĸ��ʣ������СӦ��QAM����ͼ��Ӧ
%EsN0_dB�� Es/N0��dB��ʽ
%M: QAM����
%������� Ser:��ǰ�������������


    EsN0=10.^(EsN0_dB/10);
    %M-QAM�źŸ��������������ϵ�λ�ã���dΪ����㵥λ
    QAM_dmat = (repmat([-sqrt(M)+1:2:sqrt(M)-1],sqrt(M),1)./2).^2 + (repmat([sqrt(M)-1:-2:-sqrt(M)+1]',1,sqrt(M))./2).^2;
    %����M-QAM�ڵ�ǰ����P�µ�ȡֵ
    A =sum(sum( QAM_dmat.*P));
    %�������
    ErrA = erfc(sqrt(EsN0/(4*A)))-(0.5*erfc(sqrt(EsN0/(4*A))))^2;
    ErrB = (3/2)*erfc(sqrt(EsN0/(4*A)))-2*(0.5*erfc(sqrt(EsN0/(4*A))))^2;
    ErrC = 2*erfc(sqrt(EsN0/(4*A)))-4*(0.5*erfc(sqrt(EsN0/(4*A))))^2;
    
    matrix_err = [[ErrA;ones(sqrt(M)-2,1)*ErrB;ErrA],repmat([ErrB;ones(sqrt(M)-2,1)*ErrC;ErrB],1,sqrt(M)-2),[ErrA;ones(sqrt(M)-2,1)*ErrB;ErrA]];
    Ser= sum(sum(matrix_err.*P));
end

