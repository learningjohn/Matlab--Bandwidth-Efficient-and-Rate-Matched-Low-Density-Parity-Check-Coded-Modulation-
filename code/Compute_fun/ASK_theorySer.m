function [Ser] = ASK_theorySer(P,EsN0_dB,M)
%��������QAM������ʣ���������M>=4
%�����������
%P: QAM�źŸ���������ĸ��ʣ������СӦ��QAM����ͼ��Ӧ
%EsN0_dB�� Es/N0��dB��ʽ
%M: ASK����
%������� Ser:��ǰ�������������


    EsN0=10.^(EsN0_dB/10);
    %M-QAM�źŸ��������������ϵ�λ�ã���dΪ����㵥λ
    ASK_dmat = [-(M-1):2:(M-1)]./2;
    %����M-ASK�ڵ�ǰ����P�µ�A
    A =ASK_dmat.^2*P';
    
    Ser= [1/2,ones(1,length(P)-2),1/2]*P'* erfc(sqrt(EsN0/(4*A)));
    
    %Ser = (1-1/M)*erfc(sqrt(3/(M^2-1)*EsN0));
end

