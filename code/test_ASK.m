clc;clear;
addpath("ccdm\");addpath("Compute_fun\");addpath("mat_data\");
%初始化一些参数
M = 32;
n = log2(M);
M_ASK = 2^(n);
R = 5/6;
nSim = 2000;
%编码相关参数
snr =22.5;     %信噪比

load('R_SNR.mat');
best_P = initialize_PX(R_SNR(7,2),M_ASK);

EsN0=10.^(snr/10);
LDPC_bitlength = 64800; %码长可在648，1296   1944，64800中选择

n_ASK = LDPC_bitlength/n;  %ASK信号的长度
if LDPC_bitlength == 64800
    H = dvbs2ldpc(R);            %LDPC的H矩阵
else
    H =sparse(buildH( LDPC_bitlength, R));            %LDPC的H矩阵
end

ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H);

%ASK正的信号映射表
ASK_map = 1:2:2^(n)-1;
ASK_Graymap(Graycode(1:length(ASK_map))+1) = ASK_map;
Bin2Gray_map(Graycode(1:length(ASK_map))+1) = 1:length(ASK_map);
%LDPC符号里面的比特交织
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
%格雷映射表（不考虑符号位）
[~,Gray_table] = Graycode([length(ASK_map):-1:1,1:length(ASK_map)]);
ASK = [-ASK_map(end:-1:1),ASK_map];
Gray_table = [zeros(1,length(ASK_map)),ones(1,length(ASK_map));Gray_table'];
Gray_table_map_bin = bi2de(Gray_table','left-msb')';
xx = 1:length(Gray_table_map_bin);
Gray_table_map(Gray_table_map_bin+1) = xx;


P = @(y,a,d,px) exp(-(y-a).^2./(2*d^2)) *px';

for ii = 1:length(snr) 
    
    best_P1 = 2*best_P(length(best_P)/2+1:end);
    [p_quant,nBitsInfo,n_i] = ccdm.initialize(best_P1,n_ASK);
    
    %实际码率计算LDPC的性能

%     R_real = nBitsInfo/(n_ASK*log2(M));
%     H_real =sparse(buildH( LDPC_bitlength, R_real));
%     ldpcEncoder_real = comm.LDPCEncoder(H_real);
%     ldpcDecoder_real = comm.LDPCDecoder(H_real);
    nerr_bit = 0;  
    nerr = 0;
for k = 1:nSim
        %信源比特
        txBits_I = randi(2,1,nBitsInfo)-1;
        txBits_other = randi(2,1,n_ASK/6)-1;
        
        i_TX=ccdm.encode(txBits_I,n_i).'+1;
        %格雷映射
        [~,Graymap_Tx_I_Bit] = Graycode(i_TX);
        %交织映射
        pimap_Tx_I_Bit = Graymap_Tx_I_Bit(:,bit_map_pi);
        %加上比特
        
        %LDPC编码
        data_I_bit = reshape(pimap_Tx_I_Bit',[],1);
        data_I_bit_all = [data_I_bit;txBits_other'];
        encData_I_bit = ldpcEncoder(data_I_bit_all);
        sign_symbol_I = encData_I_bit(length(data_I_bit)+1:end)*2-1;
        out_Tx_I = ASK_map(i_TX)'.*sign_symbol_I;
        
        
        Es = mean(abs(out_Tx_I).^2);
        N=Es./EsN0(ii);
        variance=N;
        Standard_variance1=sqrt(variance);
        noise=randn(1,n_ASK).';
        %AWGN
        Rx_I = out_Tx_I + Standard_variance1.*noise;
        llr_I = zeros(n_ASK,n);
        for j = 1:n
            %I路各比特位LLR
            %PS
            llr_I(:,j) = log(P(Rx_I,ASK(Gray_table(j,:)==0),Standard_variance1,best_P(Gray_table(j,:)==0))./...
                             P(Rx_I,ASK(Gray_table(j,:)==1),Standard_variance1,best_P(Gray_table(j,:)==1)));
            
        end
        
        %PS
        %I路处理
        llr_I_other = llr_I(:,1);   %符号位LLR
        llr_I_info = llr_I(:,2:end);%编码位LLR
        %交织映射
        llr_I_info_pidemap = llr_I_info(:,bit_map_pi)';
        llr_I_ldpc = [llr_I_info_pidemap(:);llr_I_other];
        %译码
        Rx_decodData_I_bit = ldpcDecoder(llr_I_ldpc);
        Rx_decodData_I_bit_I = reshape(Rx_decodData_I_bit(1:length(data_I_bit),1),n-1,[])';
        %解交织
         Rx_pimap_I_Bit = zeros(n_ASK,n-1);
        Rx_pimap_I_Bit(:,bit_map_pi) = Rx_decodData_I_bit_I;
        %逆映射
        i_RX = Bin2Gray_map(bi2de(Rx_pimap_I_Bit,'left-msb')+1)';
        %CCDM解码
        rxBits_I = ccdm.decode(i_RX-1,n_i,nBitsInfo);
        rxBits_other = Rx_decodData_I_bit(length(data_I_bit)+1:end)';
        
        if any(Rx_decodData_I_bit~=data_I_bit_all) 
            nerr(k) =  1;
        end
        nerr_bit(k) =   sum(rxBits_I~=txBits_I) + sum(rxBits_other~=txBits_other);

end
    biterr_rate(ii) = sum(nerr_bit)/((nBitsInfo+n_ASK/6)*nSim);
    framerr_rate(ii) = sum(nerr)/(nSim)

end

figure()
semilogy(snr,framerr_rate);

%3 bit/s/hz时  误码率为10^-3信噪比在18.82dB
%3.1 bit/s/hz时  误码率为10^-3信噪比在19.37dB
%3.2 bit/s/hz时  误码率为10^-3信噪比在20.13dB
%3.3 bit/s/hz时  误码率为10^-3信噪比在20.81
%3.4 bit/s/hz时  误码率为10^-3信噪比在21.4
%3.5 bit/s/hz时  误码率为10^-3信噪比在21.97 
%3.6 bit/s/hz时  误码率为10^-3信噪比在22.6
%3.7 bit/s/hz时  误码率为10^-3信噪比在23.25
%3.8 bit/s/hz时  误码率为10^-3信噪比在23.85
%3.9 bit/s/hz时  误码率为10^-3信噪比在24.5