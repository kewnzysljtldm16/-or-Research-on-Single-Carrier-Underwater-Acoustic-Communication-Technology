%% �����ռ��Ͳ���
clc;
clear all;
close all;
%ȫ��ͼ���ͼ����
% set(0, 'defaultfigurecolor', 'w') %��ͼ����Ϊ����
% set(0, 'defaultAxesFontSize', 12) %��ͼ�����������С�޸�Ϊ12
% set(0, 'defaultAxesTitleFontSizeMultiplier', 1.2) %��ͼ���������С�޸�Ϊ16
% set(0, 'defaultAxesFontName','����Ӵ�') %�޸�Ĭ������
% set(0, 'defaultAxesLineWidth', 1) %�޸�Ĭ�Ͽ���ߴ�ϸ
% set(0, 'defaultLineLineWidth', 1.2);%�޸�Ĭ���ߴ�ϸ
set(0, 'defaultAxesXGrid','on', 'defaultAxesYGrid', 'on') %������
%% ģʽѡ��
MOD = 'mo_qpsk';                                                   %����ѡ��ͬ��ӳ�䷽ʽbpsk��qpsk��8psk��ѡ
switch (MOD)
    case 'mo_bpsk', Mod = 2; bitnum_per = 1;
    case 'mo_qpsk', Mod = 4; bitnum_per = 2;
    case 'mo_8psk', Mod = 8; bitnum_per = 3; 
    otherwise, disp('Unknown signal constellation!');
end

%% ��������
fs = 100e3;                                                            % ����Ƶ��
fl = 8e3;                                                                 % ����Ƶ��
fh = 12e3;                                                              % ����Ƶ��
f0 = (fl + fh) / 2;                                                        % ����Ƶ�ʣ���Ƶ�����䣩
Rb = 2000;                                                             % ������
N_up = fs / Rb;                                                        % ����������
N_BS = 3000;                                                           % ���͵ķ�����
length_BS = N_BS * N_up;
N_bit = N_BS * bitnum_per;                                             % ��Ҫ���ɵı�����
alpha = 1;                                                              % ����ϵ��
N_filter = 512;                                                        % �˲�������
PulseShape = rcosfir(alpha, [ ], N_up, 1, 'sqrt');  % ��������˲�������ͨ�˲�����
b1 = fir1(N_filter, 2 * [fl fh] / fs);                               % ��ͨ�˲���

%% --------------------���������------------------------
%% �����źŲ���������
load information.mat
bit_send = information(1 : N_bit);            

%% ��Ϣӳ��
[SymbolIn, Table] = Mapping(bit_send, Mod);   

%% IQ������������ϲ������ز�
signal_IQ = IQmodulate(PulseShape, N_up, SymbolIn, f0, fs);
signal_IQ = signal_IQ ./ max(abs(signal_IQ));
%% �����ղ����ź� ��LFM��
T_syn = 0.1; B = fh - fl; K = B / T_syn;                              %LFM�źŲ�����B����T����K��Ƶб��
t = 0 : 1/fs : T_syn-1/fs;
signal_measure = cos(2*pi*fl*t + pi*K*t.^2);                  
length_measure = T_syn * fs;
length_GI = 0.1 * fs;                                                %�������
signal_GI = zeros(1, length_GI);

%% �����źŹ���
signal_send = [signal_measure signal_GI signal_IQ signal_GI signal_measure signal_GI];    %�źŽṹ[�����ź� ������� �����ź� ������� �����ź� �������]

%% --------------------�ŵ�����------------------------
%% ������
SNR=15;                                         %����ȣ�dB��
signal_add_noise = BandNoiseAdd(signal_send, SNR, b1 ,length_measure+length_GI, length_measure+length_GI+length_BS);
% signal_add_noise=signal_send;

%% �Ӷ�����
% �����չ��ƾ���
% �������ز��������Ķ����չ��ƾ��ȣ����������ֵ�������������ʱ��resample����ʱ���
dup_precision1 = factor_resample(fs) / fs; 
% ������FFT������ֱ��ʵĶ����վ��ȣ��������������ʱ�����ƻ�������
dup_precision2 = 1 / (length_measure + 2*length_GI + length_BS);  
% ����С��������ͬʱ�����������־�������                                           
dup_precision = dup_precision1;

% ���������
m = 5;                      % ���ƾ��ȵ�������
dup_ori= m * dup_precision;  
fs1 = fs * (1 - dup_ori); 
signal_add_dopper = resample(signal_add_noise, fs1, fs);

%% --------------------���ջ�����------------------------
signal_receive = signal_add_noise;
%% ��ͨ�˲�
signal_bandpass = filter(b1, 1, [signal_add_dopper zeros(1,fix(length(b1)/2))]);
signal_rec_pass = signal_bandpass(fix(length(b1)/2)+1:end);

%% �����ղ���
Res_xcorr = corr_fun(signal_rec_pass, signal_measure);

% ��ط������ź��ײ�β��LFM�źź���λ��
[~, pos1] = max(Res_xcorr(1 : length_measure+signal_GI));             
[~, pos2] = max(Res_xcorr(length_GI+length_BS+1 : end));      
pos2 = pos2 + length_GI + length_BS;                                      

% ��������ź���βLFM������뷢�������Ա�
del_rec = pos2 - pos1;
del_send = length_measure + 2*length_GI + length_BS ;

% ���ü���仯�������ղ���
dup_det = (del_send - del_rec) / del_send;

% �����ز������ж����ղ���
fs2 = fs*(1-dup_det);
fs2 = round(fs2 / factor_resample(fs)) * factor_resample(fs);    %ʹ������resample����Ҫ���ֹ����
signal_rec_dc = resample(signal_rec_pass, fs, fs2);   %dc��Doppler compensation

%% ��ȡ��Ϣ����
signal_rec_nodc_information = signal_rec_pass(length_measure+length_GI+1 : length_measure+length_GI+length_BS);
signal_rec_dc_information = signal_rec_dc(length_measure+length_GI+1 : length_measure+length_GI+length_BS);

%% IQ���+���ز�
[symbol_demodulate_nodc] = IQdemodulate(signal_rec_nodc_information, fs, length_BS, f0, PulseShape, N_up);
[symbol_demodulate_dc] = IQdemodulate(signal_rec_dc_information, fs, length_BS, f0, PulseShape, N_up);

%% �о���ӳ���������
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

%% ��ͼ�����
figure('Name', 'ʱ���ź�');
subplot(2,1,1)
plot(signal_send);
title('�����ź�')
subplot(2,1,2)
plot(signal_receive);
title('�����ź�')

fprintf(['����������ʵ��ֵ��'  num2str(dup_ori) '\n']);
fprintf(['���������Ӳ���ֵ��' num2str(dup_det) '\n']);

fprintf(['�����ղ����������ʣ�'  num2str(BER_nodc) '\n'] );
fprintf(['�����ղ��������ʣ�' num2str(BER_dc) '\n']);

scatterplot(symbol_demodulate_nodc);
title('���ö����ղ���ǰ')
scatterplot(symbol_demodulate_dc);
title('���ö����ղ�����')





