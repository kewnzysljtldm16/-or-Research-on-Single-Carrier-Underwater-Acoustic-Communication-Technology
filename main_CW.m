%% 多普勒检测和补偿
clc;
clear all;
close all;
%全局图像绘图设置
% set(0, 'defaultfigurecolor', 'w') %绘图背景为纯白
% set(0, 'defaultAxesFontSize', 12) %绘图坐标轴字体大小修改为12
% set(0, 'defaultAxesTitleFontSizeMultiplier', 1.2) %绘图标题字体大小修改为16
% set(0, 'defaultAxesFontName','宋体加粗') %修改默认字体
% set(0, 'defaultAxesLineWidth', 1) %修改默认框架线粗细
% set(0, 'defaultLineLineWidth', 1.2);%修改默认线粗细
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
fs = 100e3;                                                            % 采样频率
fl = 8e3;                                                                 % 下限频率
fh = 12e3;                                                              % 上限频率
f0 = (fl + fh) / 2;                                                        % 中心频率（单频带传输）
Rb = 2000;                                                             % 符号率
N_up = fs / Rb;                                                        % 升采样点数
N_BS = 3000;                                                           % 发送的符号数
length_BS = N_BS * N_up;
N_bit = N_BS * bitnum_per;                         % 需要生成的比特数
alpha = 1;                                                              % 滚降系数
N_filter = 512;                                                        % 滤波器阶数
PulseShape = rcosfir(alpha, [ ], N_up, 1, 'sqrt');  % 脉冲成型滤波器（低通滤波器）
b1 = fir1(N_filter, 2 * [fl fh] / fs);                               % 带通滤波器

%% --------------------发射机部分------------------------
%% 数据信号产生及编码
load information.mat
bit_send = information(1 : N_bit);            

%% 信息映射
[SymbolIn, Table] = Mapping(bit_send, Mod);   

%% IQ调制脉冲成型上采样上载波
signal_IQ = IQmodulate(PulseShape, N_up, SymbolIn, f0, fs);
signal_IQ = signal_IQ ./ max(abs(signal_IQ));
%% 多普勒测量信号 （CW）
T_syn = 0.2;                                                           %CW测量信号参数
t = 0 : 1/fs : T_syn-1/fs;
signal_measure = cos(2*pi*f0*t);                  
length_measure = T_syn * fs;
length_GI=0.1 * fs;                                                %保护间隔
signal_GI=zeros(1, length_GI);

%% 发送信号构成
signal_send = [signal_measure signal_GI signal_IQ signal_GI];    % 信号结构[测量信号 保护间隔 调制信号 保护间隔]

%   figure
%     Data_temp = signal_send./max((signal_send));
%     subplot(211);
%     % plot(DataFilter(1:end));
%     tt=0:1/fs:(length(Data_temp)-1)/fs;
%     plot(tt, Data_temp(1:end));
%     xlabel('时间/s');ylabel('归一化幅度');
% 
%     [y,f,t,p] = spectrogram(Data_temp,2048,2040,2048,fs,'yaxis');
%     subplot(212);
%     ax2 = imagesc(t,f,10*log10(abs(p)));
%     xlabel('时间/s');ylabel('频率/Hz');colorbar;
%     set(gca,'YDir','normal');axis([0,max(tt),6000,14000]);caxis([-90 max(max(10*log10(abs(p))))]);
%     shading interp;colorbar;colormap('jet');
%% -------------------------------------------------信道部分-------------------------------------------------
%% 加噪声
SNR=15;                                         %信噪比（dB）
signal_add_noise = BandNoiseAdd(signal_send, SNR, b1, length_measure+length_GI, length_measure+length_GI+length_BS);

%% 加多普勒
% 多普勒估计精度
% 受限于重采样函数的多普勒估计精度，多普勒设计值不满足这个倍数时，resample大概率报错
dup_precision1 = factor_resample(fs) / fs; 
% 受限于FFT的物理分辨率的多普勒精度，不满足这个倍数时，估计会存在误差
dup_precision2 = (1/T_syn) / f0;  
% 求最小公倍数，同时满足上述两种精度设置
icc = 1e7;                                %Integer conversion constant：整数转换常数                                             
dup_precision = lcm(dup_precision1*icc, dup_precision2*icc) / icc;

% 多普勒添加
m = 10;
dup_ori= m * dup_precision;  
fs1 = fs * (1-dup_ori); 
signal_add_dopper = resample(signal_add_noise, fs1, fs);

%% -------------------------------------------------接收部分-------------------------------------------------
signal_receive = signal_add_noise;
%% 带通滤波
signal_bandpass = filter(b1, 1, [signal_add_dopper zeros(1,fix(length(b1)/2))]);
signal_rec_pass = signal_bandpass(fix(length(b1)/2)+1 : end);

%% 多普勒测量
signal_rec_measure = signal_rec_pass(1 : length_measure);
N_fft=length(signal_rec_measure);                      % FFT计算点数
Sf = fft(signal_rec_measure, N_fft);                      % 信息符号做fft
Sff = fftshift(Sf);                                                   % 将零频分量移动至中心   matlab的fft计算结果是0到fs范围内，我们实际需要-fs/2~fs/2
Sf0 = Sff .* conj(Sff) / N_fft;                                 % 计算频域每个点的能量（功率）
SF = Sf0(N_fft/2+1 : N_fft);                                  % 取右边一半（正频率范围）
sr = fs / N_fft;                                                      % spectral resolution 物理分辨率
ff0 = (0 : N_fft/2-1) * sr;                                       % 频率和点数的对应关系
[~, f_pos] = max(SF);                                           % 获取能量（功率）最大的位置
f_det = ff0(f_pos);                                                % 获取能量（功率）最大位置对应的频率
dup_det = (f_det-f0) / f0;                                    % 通过频偏计算多普勒系数

figure('Name', '频谱分析');
plot(ff0,SF)
xlabel('频率/Hz'); ylabel('幅度'); title('频谱分析结果')
%% 重采样法进行多普勒补偿
fs2 = fs * (1-dup_det);
fs2 = round(fs2/factor_resample(fs)) * factor_resample(fs);    %使其满足resample精度要求防止报错
signal_rec_dc = resample(signal_rec_pass, fs, fs2);                         %dc：Doppler compensation

%% 提取信息符号
signal_rec_nodc_information = signal_rec_pass(length_measure+length_GI+1 : length_measure+length_GI+length_BS);
signal_rec_dc_information = signal_rec_dc(length_measure+length_GI+1 : length_measure+length_GI+length_BS);

%% IQ解调+下载波
[symbol_demodulate_nodc] = IQdemodulate(signal_rec_nodc_information, fs, length_BS, f0, PulseShape, N_up);
[symbol_demodulate_dc] = IQdemodulate(signal_rec_dc_information, fs, length_BS, f0, PulseShape, N_up);

%% 解码计算误码
for j = 1 : length(symbol_demodulate_nodc)
    Distance_all = abs(symbol_demodulate_nodc(j) - Table);
    Tablemin=find(Distance_all == min(Distance_all));
    symbol_decision_nodc(j) = Table(Tablemin(1));
end
bit_nodc  = Demapping(symbol_decision_nodc , Table , Mod);
BER_nodc = length(find(bit_nodc ~= bit_send)) ./ N_bit;

for j = 1 : length(symbol_demodulate_dc)
    Distance_all = abs(symbol_demodulate_dc(j) - Table);
    Tablemin=find(Distance_all == min(Distance_all));
    symbol_decision_dc(j) = Table(Tablemin(1));
end
bit_dc  = Demapping(symbol_decision_dc , Table , Mod);
BER_dc = length(find(bit_dc ~= bit_send)) ./ N_bit;

%% 绘图及输出
figure('Name', '时域信号');
subplot(2,1,1)
plot(signal_send);
title('发送信号')
subplot(2,1,2)
plot(signal_receive);
title('接收信号')

fprintf(['多普勒因子实际值：'  num2str(dup_ori) '\n']);
fprintf(['多普勒因子测量值：' num2str(dup_det) '\n']);

fprintf(['多普勒不补偿误码率：'  num2str(BER_nodc) '\n'] );
fprintf(['多普勒补偿误码率：' num2str(BER_dc) '\n']);

scatterplot(symbol_demodulate_nodc);
title('采用多普勒补偿前')
scatterplot(symbol_demodulate_dc);
title('采用多普勒补偿后')
