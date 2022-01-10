%ASK误码率测试
clc;clear;
addpath("ccdm\");addpath("Compute_fun\");addpath("mat_data\");
load('R_SNR.mat');
%调制参数
M = 16;
n = log2(M);
R = 3/4;

nSim = 10000;
%信道和编码相关的参数
snr = 20.3:0.05:20.5;     %信噪比
PX = initialize_PX(snr,M);
EsN0=10.^(snr/10);
LDPC_bitlength = 64800; %码长可在648，1296   1944，64800中选择
n_ASK = LDPC_bitlength/n;  %ASK信号的长度

if LDPC_bitlength == 64800 %如果LDPC编码长度为64800则使用DVB-S2标准
    H = dvbs2ldpc(R);            %LDPC的H矩阵
   
else                       %如果长度为其他则使用ieee标准
    H =sparse(buildH( LDPC_bitlength, R));            %LDPC的H矩阵
   
end
%初始化固定码率的LDPC矩阵
ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H);


%ASK正值信号映射表
ASK_map = 1:2:2^(n)-1;
%格雷映射表
ASK_Graymap(Graycode(1:length(ASK_map))+1) = ASK_map;
%二进制比特映射表
Bin2Gray_map(Graycode(1:length(ASK_map))+1) = 1:length(ASK_map);
%LDPC前的ASK符号比特交织映射表
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
ASK = [-ASK_map(end:-1:1),ASK_map]
Gray_table = [zeros(1,length(ASK_map)),ones(1,length(ASK_map));Gray_table']
Gray_table_map_bin = bi2de(Gray_table','left-msb')'
xx = 1:length(Gray_table_map_bin);
Gray_table_map(Gray_table_map_bin+1) = xx

%对数似然比计算公式
P = @(y,a,d,px) exp(-(y-a).^2./(2*d^2)) *px';

for ii = 1:length(snr)

    nerr_bit = 0;
    nerr = 0;
    for k = 1:nSim
        %信源比特
        txBits_I_Uniform = randi(2,1,R*LDPC_bitlength)'-1;
        %不进行概率成型
        encData_I_bit_Uniform = reshape(ldpcEncoder(txBits_I_Uniform),n,[])';
        %交织
        pimap_Tx_I_Bit_Uniform = encData_I_bit_Uniform(:,bit_map_pi_Uniform);
        %Gray映射
        out_Tx_I_Uniform = ASK(Gray_table_map(bi2de(pimap_Tx_I_Bit_Uniform,'left-msb')+1));
        %均匀分布的噪声
        Es_Uniform = mean(abs(out_Tx_I_Uniform).^2);
        N0_Uniform=Es_Uniform./EsN0(ii);
        variance_Uniform=N0_Uniform;
        Standard_variance1_Uniform=sqrt(variance_Uniform);
        noise=randn(1,n_ASK).';
        %AWGN
        Rx_I_Uniform = out_Tx_I_Uniform.' + Standard_variance1_Uniform.*noise;
        
        %计算格雷映射下每一个位的对数似然比
        llr_I_Uniform = zeros(n_ASK,n);
        for j = 1:n
            %I路各比特位LLR
            %随机
            llr_I_Uniform(:,j) = log(P(Rx_I_Uniform,ASK(Gray_table(j,:)==0),Standard_variance1_Uniform,ones(1,M/2)/M)./...
                P(Rx_I_Uniform,ASK(Gray_table(j,:)==1),Standard_variance1_Uniform,ones(1,M/2)/M));
            
        end
        %随机分布
        %解交织映射
        llr_I_info_pidemap_Uniform(bit_map_pi_Uniform,:) = llr_I_Uniform';
        llr_I_ldpc_Uniform = llr_I_info_pidemap_Uniform(:);
        %译码
        rxBits_I_Uniform = ldpcDecoder(llr_I_ldpc_Uniform);
        
        if any(rxBits_I_Uniform~=txBits_I_Uniform)
            nerr(k) =  1;
        end
        nerr_bit(k) = sum(rxBits_I_Uniform~=txBits_I_Uniform)
        
    end
    biterr_rate(ii) = sum(nerr_bit)/(R*LDPC_bitlength*nSim);
    framerr_rate(ii) = sum(nerr)/(nSim)
end
