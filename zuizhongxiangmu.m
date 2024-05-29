%% 多普勒检测和补偿
clc;
clear all;
close all;
%全局图像绘图设置
set(0, 'defaultAxesXGrid','on', 'defaultAxesYGrid', 'on') %打开网格
%% 模式选择
MOD = 'mo_qpsk';                                                   %可以选择不同的映射方式bpsk、qpsk、8psk可选
switch (MOD)
    case 'mo_bpsk', Mod = 2; bitnum_per = 1;
    case 'mo_qpsk', Mod = 4; bitnum_per = 2;
    case 'mo_8psk', Mod = 8; bitnum_per = 3; 
    otherwise, disp('Unknown signal constellation!');
end
%% 基本参数
fs = 48e3;                                                            % 采样频率
fl = 8e3;                                                                 % 下限频率
fh = 12e3;                                                              % 上限频率
f0 = (fl + fh) / 2;                                                        % 中心频率（单频带传输）
Rb = 2000;                                                             % 符号率
N_up = fs / Rb;                                                        % 升采样点数
N_BS = 2000;                                                           % 发送的符号数
length_BS = N_BS * N_up;
N_bit = N_BS * bitnum_per;                                             % 需要生成的比特数
alpha = 1;                                                              % 滚降系数
N_filter = 256;                                                        % 滤波器阶数
PulseShape = rcosdesign(alpha, N_filter/N_up, N_up, 'sqrt');  % 根升余弦脉冲成型滤波器
b1 = fir1(N_filter, 2 * [fl fh] / fs);                               % 带通滤波器
%% --------------------发射机部分------------------------
%% 数据信号产生及编码
load information.mat
bit_send = information(1 : N_bit);  
bit_send_1 = scramble(bit_send);
figure
subplot(3,1,1)
plot(bit_send)
title('原始信号')
subplot(3,1,2)
plot(bit_send_1)
title('加扰后的信号')
subplot(3,1,3)
plot(bit_send-bit_send_1)
title('加扰后信号与原始信号之差')

%% 信息映射
[SymbolIn, Table] = Mapping(bit_send, Mod);   
%% IQ调制脉冲成型上采样上载波
signal_IQ_1 = IQmodulate(PulseShape, N_up, SymbolIn, f0, fs);
signal_IQ = signal_IQ_1 ./ max(abs(signal_IQ_1));
%% 多普勒测量信号 （LFM）
T_syn = 0.05; B = 4e3; K = B / T_syn;                              %LFM信号参数，B带宽，T脉宽，K调频斜率
t = 0 : 1/fs : T_syn-1/fs;
signal_measure = cos(2*pi*fl*t + pi*K*t.^2);                  
length_measure = T_syn * fs;
length_GI = 0.05 * fs;                                                %保护间隔
signal_GI = zeros(1, length_GI);
%% 发送信号构成
signal_send = [signal_measure signal_GI signal_IQ signal_GI signal_measure signal_GI];    %信号结构[测量信号 保护间隔 调制信号 保护间隔 测量信号 保护间隔]

%% 画出脉冲成型前的I路信号，脉冲成型前的q路信号
figure;
subplot(2,2,1);
plot(real(SymbolIn));
title('脉冲成型前的I路信号');
subplot(2,2,2);
plot(imag(SymbolIn));
title('脉冲成型前的q路信号');
%% 画出脉冲成型后的I路信号，脉冲成型后的q路信号
subplot(2,2,3);
plot(real(signal_IQ));
title('脉冲成型后的I路信号');
subplot(2,2,4);
plot(imag(signal_IQ));
title('脉冲成型后的q路信号');

%% 画出脉冲成型前和脉冲成型后的频率幅度图
figure;
subplot(2,1,1);
plot(abs(fft(SymbolIn)));
title('脉冲成型前的频率幅度图');
subplot(2,1,2);
plot(abs(fft(signal_IQ)));
title('脉冲成型后的频率幅度图');


%% 画出载波调制前的时域波形和频域波形
figure;
subplot(2,2,1);
plot(real(SymbolIn));
title('载波调制前的时域波形');
subplot(2,2,2);
plot(abs(fft(SymbolIn)));
title('载波调制前的频域波形');
% 画出载波调制后的时域波形和频域波形
subplot(2,2,3);
plot(real(signal_IQ));
title('载波调制后的时域波形');
subplot(2,2,4);
plot(abs(fft(signal_IQ)));
title('载波调制后的频域波形');


%% 发射端流程
% 数字信源
load information.mat
bit_send = information(1 : N_bit); 

%% --------------------信道部分------------------------
%% 加噪声
SNR=15;                                         %信噪比（dB）
signal_add_noise = BandNoiseAdd(signal_send, SNR, b1 ,length_measure+length_GI, length_measure+length_GI+length_BS);
% signal_add_noise=signal_send;

%% 加多普勒
% 多普勒估计精度
% 受限于重采样函数的多普勒估计精度，多普勒设计值不满足这个倍数时，resample大概率报错
dup_precision1 = factor_resample(fs) / fs; 
% 受限于FFT的物理分辨率的多普勒精度，不满足这个倍数时，估计会存在误差
dup_precision2 = 1 / (length_measure + 2*length_GI + length_BS);  
% 求最小公倍数，同时满足上述两种精度设置                                           
dup_precision = dup_precision1;

% 多普勒添加
m = 5;                      % 估计精度的整倍数
dup_ori= m * dup_precision;  
fs1 = fs * (1 - dup_ori); 
signal_add_dopper = resample(signal_add_noise, fs1, fs);

%% 水下噪声和传播损失
% 深海噪声谱举例和海水的吸收系数
% 这里只是一个示例，您需要根据实际情况进行修改
depth = 1000; % 深度，单位：米
distance = 10000; % 距离，单位：米
n = 20; % 传播损失参数，这是一个假设的值，您需要根据实际情况进行修改
alpha = 0.1; % 海水的吸收系数，这是一个假设的值，您需要根据实际情况进行修改

% 计算传播损失
TL = n * 10 * log10(distance) + distance * alpha;

% 添加深海噪声
noise = wgn(length(signal_add_dopper), 1, 0); % 生成白噪声
noise = filter(1, [1 -0.99], noise); % 通过一个IIR滤波器使噪声具有一定的色度
noise = noise / std(noise) * 10^(-TL/20); % 根据传播损失调整噪声功率
noise = noise';
%signal_with_noise = resample(signal_add_dopper ,noise);
signal_with_noise = signal_add_dopper+noise;


%% --------------------接收机部分------------------------
signal_receive = signal_with_noise;
%% 带通滤波
signal_bandpass = filter(b1, 1, [signal_add_dopper zeros(1,fix(length(b1)/2))]);
signal_rec_pass = signal_bandpass(fix(length(b1)/2)+1:end);

%% 多普勒测量
Res_xcorr = corr_fun(signal_rec_pass, signal_measure);

% 相关法计算信号首部尾部LFM信号后沿位置
[~, pos1] = max(Res_xcorr(1 : length_measure+signal_GI));             
[~, pos2] = max(Res_xcorr(length_GI+length_BS+1 : end));      
pos2 = pos2 + length_GI + length_BS;                                      

% 计算接收信号首尾LFM间隔，与发射间隔做对比
del_rec = pos2 - pos1;
del_send = length_measure + 2*length_GI + length_BS ;

% 利用间隔变化做多普勒测量
dup_det = (del_send - del_rec) / del_send;

% 利用重采样进行多普勒补偿
fs2 = fs*(1-dup_det);
fs2 = round(fs2 / factor_resample(fs)) * factor_resample(fs);    %使其满足resample精度要求防止报错
signal_rec_dc = resample(signal_rec_pass, fs, fs2);   %dc：Doppler compensation

%% 提取信息符号
signal_rec_nodc_information = signal_rec_pass(length_measure+length_GI+1 : length_measure+length_GI+length_BS);
signal_rec_dc_information = signal_rec_dc(length_measure+length_GI+1 : length_measure+length_GI+length_BS);
% 解扰
signa1_rec_dc_information = descramble(signal_rec_dc_information);
%% IQ解调+下载波
[symbol_demodulate_nodc] = IQdemodulate(signal_rec_nodc_information, fs, length_BS, f0, PulseShape, N_up);
[symbol_demodulate_dc] = IQdemodulate(signal_rec_dc_information, fs, length_BS, f0, PulseShape, N_up);
scatterplot(symbol_demodulate_nodc);
title('采用多普勒补偿前')
scatterplot(symbol_demodulate_dc);
title('采用多普勒补偿后')

%% LMS均衡器
mu = 0.01; % 学习率
delay = 5; % 延迟
M = 15; % 滤波器长度
eq1 = dsp.LMSFilter(M,'StepSize',mu,'Method','Normalized LMS'); % 创建一个LMS均衡器
[symbol_demodulate_dc_eq1, err1] = eq1(symbol_demodulate_dc(delay+1:end).',symbol_demodulate_dc(1:end-delay).'); % 使用均衡器

%% RLS均衡器
forgetFactor = 0.98; % 遗忘因子
eq2 = dsp.RLSFilter(M,'ForgettingFactor',forgetFactor); % 创建一个RLS均衡器
[symbol_demodulate_dc_eq2, err2] = eq2(symbol_demodulate_dc(delay+1:end).',symbol_demodulate_dc(1:end-delay).'); % 使用均衡器

%% MMSE均衡器
nWeights = 15;
nFwdWeights = 10;
nSampPerSym = 1;
algType = 'MMSE';
stepSize = 0.01;
leakageFactor = 1;
eq3 = comm.LinearEqualizer('Algorithm','RLS','NumTaps',nWeights,'ForgettingFactor',leakageFactor); % 创建一个MMSE均衡器
eq3.ReferenceTap = nFwdWeights;
[symbol_demodulate_dc_eq3, err3] = eq3(symbol_demodulate_dc.',symbol_demodulate_dc.'); % 使用均衡器

%% DFE均衡器
nFwdWeights = 10;
nFbkWeights = 5;
eq4 = comm.DecisionFeedbackEqualizer('Algorithm','RLS','NumForwardTaps',nFwdWeights,'NumFeedbackTaps',nFbkWeights); % 创建一个DFE均衡器
[symbol_demodulate_dc_eq4, err4] = eq4(symbol_demodulate_dc.',symbol_demodulate_dc.'); % 使用均衡器

%% 判决解映射计算误码
for j = 1 : length(symbol_demodulate_dc_eq1)
    Distance_all = abs(symbol_demodulate_dc_eq1(j) - Table);
    Tablemin=find(Distance_all == min(Distance_all));
    symbol_decision_dc_eq1(j) = Table(Tablemin(1));
end
bit_dc_eq1  = Demapping(symbol_decision_dc_eq1 , Table , Mod);
bit_send_eq1 = bit_send(1:length(bit_dc_eq1)); % 截取bit_send以使其与bit_dc_eq1的长度相同
BER_dc_eq1 = length(find(bit_dc_eq1 ~= bit_send_eq1)) ./ length(bit_dc_eq1);

for j = 1 : length(symbol_demodulate_dc_eq2)
    Distance_all = abs(symbol_demodulate_dc_eq2(j) - Table);
    Tablemin=find(Distance_all == min(Distance_all));
    symbol_decision_dc_eq2(j) = Table(Tablemin(1));
end
bit_dc_eq2  = Demapping(symbol_decision_dc_eq2 , Table , Mod);
bit_send_eq2 = bit_send(1:length(bit_dc_eq2)); % 截取bit_send以使其与bit_dc_eq2的长度相同
BER_dc_eq2 = length(find(bit_dc_eq2 ~= bit_send_eq2)) ./ length(bit_dc_eq2);

for j = 1 : length(symbol_demodulate_dc_eq3)
    Distance_all = abs(symbol_demodulate_dc_eq3(j) - Table);
    Tablemin=find(Distance_all == min(Distance_all));
    symbol_decision_dc_eq3(j) = Table(Tablemin(1));
end
bit_dc_eq3  = Demapping(symbol_decision_dc_eq3 , Table , Mod);
bit_send_eq3 = bit_send(1:length(bit_dc_eq3)); % 截取bit_send以使其与bit_dc_eq3的长度相同
BER_dc_eq3 = length(find(bit_dc_eq3 ~= bit_send_eq3)) ./ length(bit_dc_eq3);

for j = 1 : length(symbol_demodulate_dc_eq4)
    Distance_all = abs(symbol_demodulate_dc_eq4(j) - Table);
    Tablemin=find(Distance_all == min(Distance_all));
    symbol_decision_dc_eq4(j) = Table(Tablemin(1));
end
bit_dc_eq4  = Demapping(symbol_decision_dc_eq4 , Table , Mod);
bit_send_eq4 = bit_send(1:length(bit_dc_eq4)); % 截取bit_send以使其与bit_dc_eq4的长度相同
BER_dc_eq4 = length(find(bit_dc_eq4 ~= bit_send_eq4)) ./ length(bit_dc_eq4);

fprintf(['LMS均衡后误码率：'  num2str(BER_dc_eq1) '\n'] );
fprintf(['RLS均衡后误码率：' num2str(BER_dc_eq2) '\n']);
fprintf(['MMSE均衡后误码率：' num2str(BER_dc_eq3) '\n']);
fprintf(['DFE均衡后误码率：' num2str(BER_dc_eq4) '\n']);
fprintf(['多普勒因子实际值：'  num2str(dup_ori) '\n']);
fprintf(['多普勒因子测量值：' num2str(dup_det) '\n']);

figure('Name', '时域信号');
subplot(2,1,1)
plot(signal_send);
title('发送信号')
subplot(2,1,2)
plot(signal_receive);
title('接收信号')

% % 
% %% 判决解映射计算误码
% for j = 1 : length(symbol_demodulate_nodc)
%     Distance_all = abs(symbol_demodulate_nodc(j) - Table);
%     Tablemin=find(Distance_all == min(Distance_all));
%     symbol_decision_nodc(j) = Table(Tablemin(1));
% end
% bit_nodc  = Demapping(symbol_decision_nodc , Table , Mod);
% BER_nodc = length(find(bit_nodc ~= bit_send)) ./ N_bit;
% 
% for j = 1 : length(symbol_demodulate_dc)
%     Distance_all = abs(symbol_demodulate_dc(j) - Table);
%     Tablemin=find(Distance_all == min(Distance_all));
%     symbol_decision_dc(j) = Table(Tablemin(1));
% end
% bit_dc  = Demapping(symbol_decision_dc , Table , Mod);
% BER_dc = length(find(bit_dc ~= bit_send)) ./ N_bit;
% 
% %% 绘图及输出
% 
% 
% fprintf(['多普勒因子实际值：'  num2str(dup_ori) '\n']);
% fprintf(['多普勒因子测量值：' num2str(dup_det) '\n']);
% 
% fprintf(['多普勒不补偿误码率：'  num2str(BER_nodc) '\n'] );
% fprintf(['多普勒补偿误码率：' num2str(BER_dc) '\n']);
% 
% scatterplot(symbol_demodulate_nodc);
% title('采用多普勒补偿前')
% scatterplot(symbol_demodulate_dc);
% title('采用多普勒补偿后')

