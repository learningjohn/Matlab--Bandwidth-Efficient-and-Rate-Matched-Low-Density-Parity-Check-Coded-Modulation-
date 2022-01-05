clc;clear;
addpath("ccdm\");addpath("Compute_fun\");addpath("mat_data\");
% 初始化一些参数
% 设计一个概率成型的编码调制方案，使用LDPC编码，当调制阶数为n = log2(M)时，编码码率为R = n-2/n。
M = 16
n = log2(M)
M_ASK = 2^(n/2)
R = (n-2)/n
nSim = 10000;

% QAM由两路ASK信号构成,每路ASK信号独立同分布，根据文件名读取之前保存的分布
filename_bestP = ['best_P_',num2str(2^(n/2)),'ASK.mat'];
load(filename_bestP);

%
% 信道和编码相关的参数
snr = 3:0.5:8;     %信噪比
PX = initialize_PX(snr,M_ASK)
EsN0=10.^(snr/10);
LDPC_bitlength = 648; %码长可在648，1296   1944，64800中选择

n_ASK = LDPC_bitlength*2/n;  %ASK信号的长度
if LDPC_bitlength == 64800
    H = dvbs2ldpc(R);            %LDPC的H矩阵
else
    H =sparse(buildH( LDPC_bitlength, R));            %LDPC的H矩阵
end
% 初始化LDPC编解码函数
ldpcEncoder = comm.LDPCEncoder(H);
ldpcDecoder = comm.LDPCDecoder(H);
% 调制编码相关的一些参数
%ASK正的信号映射表
ASK_map = 1:2:2^(n/2)-1
ASK_Graymap(Graycode(1:length(ASK_map))+1) = ASK_map
Bin2Gray_map(Graycode(1:length(ASK_map))+1) = 1:length(ASK_map);
%LDPC符号里面的比特交织
switch M
    case 16
        bit_map_pi = [1];
        bit_map_pi_Uniform = [2,1]
    case 64
        bit_map_pi = [2,1];
        bit_map_pi_Uniform = [3,2,1];
    case 256
        bit_map_pi = [2,3,1];
        bit_map_pi_Uniform = [3,4,2,1];
end
%I路格雷映射表
[~,Gray_table_I] = Graycode([length(ASK_map):-1:1,1:length(ASK_map)]);
ASK = [-ASK_map(end:-1:1),ASK_map]
Gray_table_I = [zeros(1,length(ASK_map)),ones(1,length(ASK_map));Gray_table_I']
Gray_table_map_bin = bi2de(Gray_table_I','left-msb')'
xx = 1:length(Gray_table_map_bin);
Gray_table_map_I(Gray_table_map_bin+1) = xx
%Q路格雷映射表
Gray_table_Q = Gray_table_I(:,end:-1:1)
Gray_table_map_Q(Gray_table_map_bin+1) = xx(end:-1:1)

% 对数似然比的计算公式
%
% 《Bandwidth Efficient and Rate-Matched Low-Density Parity-Check Coded Modulation
% Georg》
% 化简后
%
% 《Probabilistically shaped coded modulation for IM/DD system》
P = @(y,a,d) sum(exp(-(y-a).^2./(2*d^2)).');

% 仿真环节
for ii = 1:length(snr)
    
    % 初始化CCDM, best_P为原始分布，p_quant为CCDM量化后的输出分布，n_ASK为输出符号个数，nBitsInfo为输入比特个数
    best_P = PX(ii,:);
    best_P1 = 2*best_P(length(best_P)/2+1:end);
    [p_quant,nBitsInfo,n_i] = ccdm.initialize(best_P1,n_ASK);
    
    R_CCDM = nBitsInfo/(n_ASK*log2(length(best_P1)));
    
    nerr_bit = 0;
    nerr_bit_Uniform = 0;
    
    % 发送端部分
    % 对随机比特分布匹配后，进行格雷映射与交织送入LDPC编码产生冗余位，与正的ASK信号相乘后得到需要的ASK分布
    % 然后讲两路独立的ASK信号分别调制到实部虚部，相加得到方形QAM信号
    for k = 1:nSim
        
        %不进行概率成型的信源
        txBits_I_Uniform = randi(2,1,R*LDPC_bitlength)'-1;
        txBits_Q_Uniform = randi(2,1,R*LDPC_bitlength)'-1;
        
        
        %I路ASK
        
        
        %不进行概率成型
        encData_I_bit_Uniform = reshape(ldpcEncoder(txBits_I_Uniform),n/2,[])';
        %交织
        pimap_Tx_I_Bit_Uniform = encData_I_bit_Uniform(:,bit_map_pi_Uniform);

        %Gray映射
        out_Tx_I_Uniform = ASK(Gray_table_map_I(bi2de(pimap_Tx_I_Bit_Uniform,'left-msb')+1));
        
        %Q路ASK
        
        
        %不进行概率成型
        encData_Q_bit_Uniform = reshape(ldpcEncoder(txBits_Q_Uniform),n/2,[])';
        %交织
        pimap_Tx_Q_Bit_Uniform = encData_Q_bit_Uniform(:,bit_map_pi_Uniform);

        %Gray映射
        out_Tx_Q_Uniform = ASK(Gray_table_map_Q(bi2de(pimap_Tx_Q_Bit_Uniform,'left-msb')+1));
        %I路Q路实部虚部叠加就得到了我们需要的QAM信号
        
        Tx_QAM_Uniform =  out_Tx_I_Uniform + 1i.*out_Tx_Q_Uniform;
        % 信道部分·
        %生成噪声
        %均匀分布的噪声
        Es_Uniform = mean(abs(Tx_QAM_Uniform).^2);
        N0_Uniform=Es_Uniform./EsN0(ii);
        variance_Uniform=N0_Uniform;
        Standard_variance1_Uniform=sqrt(variance_Uniform);
        noise=(randn(1,n_ASK).'+1i*randn(1,n_ASK).')./sqrt(2);
        %AWGN
        Rx_QAM_Uniform = Tx_QAM_Uniform.' + Standard_variance1_Uniform.*noise;
        
        
        % 接收端部分
        %I路Q路分开处理
        
        Rx_I_Uniform = real(Rx_QAM_Uniform);Rx_Q_Uniform = imag(Rx_QAM_Uniform);%均匀分布的
        
        %计算格雷映射下每一个位的对数似然比
        
        for j = 1:n/2
            %I路各比特位LLR
            
            
            %随机
            llr_I_Uniform(:,j) = log(P(Rx_I_Uniform,ASK(Gray_table_I(j,:)==0),Standard_variance1_Uniform/sqrt(2))./...
                P(Rx_I_Uniform,ASK(Gray_table_I(j,:)==1),Standard_variance1_Uniform/sqrt(2)));
            %Q路各比特位LLR
            
            %随机
            llr_Q_Uniform(:,j) = log(P(Rx_Q_Uniform,ASK(Gray_table_Q(j,:)==0),Standard_variance1_Uniform/sqrt(2))./...
                P(Rx_Q_Uniform,ASK(Gray_table_Q(j,:)==1),Standard_variance1_Uniform/sqrt(2)));
        end
        
        
        %PS
        %I路处理
        
        %随机分布
        %解交织映射
        llr_I_info_pidemap_Uniform = llr_I_Uniform(:,bit_map_pi_Uniform)';
        llr_I_ldpc_Uniform = llr_I_info_pidemap_Uniform(:);
        %译码
        rxBits_I_Uniform = ldpcDecoder(llr_I_ldpc_Uniform);
        
        %Q路处理
        %随机分布
        %解交织映射
        llr_Q_info_pidemap_Uniform = llr_Q_Uniform(:,bit_map_pi_Uniform)';
        llr_Q_ldpc_Uniform = llr_Q_info_pidemap_Uniform(:);
        %译码
        rxBits_Q_Uniform = ldpcDecoder(llr_Q_ldpc_Uniform);
        
        %统计误比特率
        %PS
        
        %随机
        nerr_bit_Uniform =  nerr_bit_Uniform + sum(rxBits_I_Uniform~=txBits_I_Uniform) + sum(rxBits_Q_Uniform~=txBits_Q_Uniform);
    end
    
    biterr_rate_Uniform(ii) = nerr_bit_Uniform/(R*LDPC_bitlength*nSim*2)
    
end

figure()

semilogy(snr,biterr_rate_Uniform);
legend('随机分布');grid on

