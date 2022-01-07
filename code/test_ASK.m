clc;clear;
addpath("ccdm\");addpath("Compute_fun\");addpath("mat_data\");
%��ʼ��һЩ����
M = 32;
n = log2(M);
M_ASK = 2^(n);
R = 5/6;
nSim = 10000;
%������ز���
snr = 18:0.2:19;     %�����
load('R_SNR.mat');
EsN0=10.^(snr/10);
LDPC_bitlength = 64800; %�볤����648��1296   1944��64800��ѡ��

n_ASK = LDPC_bitlength/n;  %ASK�źŵĳ���
if LDPC_bitlength == 64800
    H = dvbs2ldpc(R);            %LDPC��H����
else
    H =sparse(buildH( LDPC_bitlength, R));            %LDPC��H����
end

ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H);

%ASK�����ź�ӳ���
ASK_map = 1:2:2^(n)-1;
ASK_Graymap(Graycode(1:length(ASK_map))+1) = ASK_map;
Bin2Gray_map(Graycode(1:length(ASK_map))+1) = 1:length(ASK_map);
%LDPC��������ı��ؽ�֯
switch M
    case 4
        bit_map_pi = [1];
        bit_map_pi_Uniform = [2,1];
    case 8
        bit_map_pi = [2,1];
        bit_map_pi_Uniform = [3,2,1];
    case 16
        bit_map_pi = [2,3,1];
        bit_map_pi_Uniform = [3,4,2,1];
    case 32
        bit_map_pi = [2,1,4,3];
        bit_map_pi_Uniform = [3,2,5,4,1];    
end
%����ӳ��������Ƿ���λ��
[~,Gray_table] = Graycode([length(ASK_map):-1:1,1:length(ASK_map)]);
ASK = [-ASK_map(end:-1:1),ASK_map];
Gray_table = [zeros(1,length(ASK_map)),ones(1,length(ASK_map));Gray_table'];
Gray_table_map_bin = bi2de(Gray_table','left-msb')';
xx = 1:length(Gray_table_map_bin);
Gray_table_map(Gray_table_map_bin+1) = xx;

best_P = initialize_PX(R_SNR(1,2),M_ASK);
P = @(y,a,d,px) exp(-(y-a).^2./(2*d^2)) *px';

for ii = 1:length(snr) 
    
    best_P1 = 2*best_P(length(best_P)/2+1:end);
    [p_quant,nBitsInfo,n_i] = ccdm.initialize(best_P1,n_ASK);
    
    %ʵ�����ʼ���LDPC������

%     R_real = nBitsInfo/(n_ASK*log2(M));
%     H_real =sparse(buildH( LDPC_bitlength, R_real));
%     ldpcEncoder_real = comm.LDPCEncoder(H_real);
%     ldpcDecoder_real = comm.LDPCDecoder(H_real);
    nerr_bit = 0;  
    nerr = 0;
for k = 1:nSim
        %��Դ����
        txBits_I = randi(2,1,nBitsInfo)-1;
        txBits_other = randi(2,1,n_ASK/6)-1;
        
        i_TX=ccdm.encode(txBits_I,n_i).'+1;
        %����ӳ��
        [~,Graymap_Tx_I_Bit] = Graycode(i_TX);
        %��֯ӳ��
        pimap_Tx_I_Bit = Graymap_Tx_I_Bit(:,bit_map_pi);
        %���ϱ���
        
        %LDPC����
        data_I_bit = reshape(pimap_Tx_I_Bit',[],1);
        data_I_bit_all = [data_I_bit;txBits_other'];
        encData_I_bit = ldpcEncoder(data_I_bit_all);
        sign_symbol_I = encData_I_bit(length(data_I_bit)+1:end)*2-1;
        out_Tx_I = ASK_map(i_TX)'.*sign_symbol_I;
        
         Es = mean(abs(out_Tx_I).^2);
        N0=Es./EsN0(ii);
        variance=N0/2;
        Standard_variance1=sqrt(variance);
        noise=randn(1,n_ASK).';
        %AWGN
        Rx_I = out_Tx_I + Standard_variance1.*noise;
        llr_I = zeros(n_ASK,n);
        for j = 1:n
            %I·������λLLR
            %PS
            llr_I(:,j) = log(P(Rx_I,ASK(Gray_table(j,:)==0),Standard_variance1/sqrt(2),best_P(Gray_table(j,:)==0))./...
                             P(Rx_I,ASK(Gray_table(j,:)==1),Standard_variance1/sqrt(2),best_P(Gray_table(j,:)==1)));
            
        end
        
        %PS
        %I·����
        llr_I_other = llr_I(:,1);   %����λLLR
        llr_I_info = llr_I(:,2:end);%����λLLR
        %��֯ӳ��
        llr_I_info_pidemap = llr_I_info(:,bit_map_pi)';
        llr_I_ldpc = [llr_I_info_pidemap(:);llr_I_other];
        %����
        Rx_decodData_I_bit = ldpcDecoder(llr_I_ldpc);
        Rx_decodData_I_bit_I = reshape(Rx_decodData_I_bit(1:length(data_I_bit),1),n-1,[])';
        %�⽻֯
         Rx_pimap_I_Bit = zeros(n_ASK,n-1);
        Rx_pimap_I_Bit(:,bit_map_pi) = Rx_decodData_I_bit_I;
        %��ӳ��
        i_RX = Bin2Gray_map(bi2de(Rx_pimap_I_Bit,'left-msb')+1)';
        %CCDM����
        rxBits_I = ccdm.decode(i_RX-1,n_i,nBitsInfo);
        rxBits_other = Rx_decodData_I_bit(length(data_I_bit)+1:end)';
        
        if any(Rx_decodData_I_bit~=data_I_bit_all) 
            nerr =  nerr +1;
        end
        nerr_bit =  nerr_bit + sum(rxBits_I~=txBits_I) + sum(rxBits_other~=txBits_other);

end
    biterr_rate(ii) = nerr_bit/((nBitsInfo+n_ASK/6)*nSim)
    framerr_rate(ii) = nerr/(nSim);

end