function [SignalAftNoise] = BandNoiseAdd(Signal_in,SNR,BandPass,s_begin,s_end)
% 带限噪声添加函数
% 改进部分：只利用关注部分计算噪声，反正保护间隔参与噪声计算使噪声添加有误差
% Signal_in：输入信号；SNR：信噪比；BandPass：带限滤波器系数；s_begin/s_end：关注部分的开始位置和结束位置

    NumFilter = length(BandPass)-1;
    if isequal(imag(Signal_in), zeros(1, length(Signal_in)))   %实数信号的接收噪声
        noise = normrnd(0, 1, 1, length(Signal_in));      
    else
        noise = sqrt(1/2)*normrnd(0, 1, 1, length(Signal_in)) + 1i*sqrt(1/2)*normrnd(0, 1, 1, length(Signal_in)); %复数信号的接收噪声
    end
    
    % 生成随机噪声
    NoiseAftFilter = filter(BandPass, 1, [noise zeros(1, NumFilter/2)]);           %滤波，1为分母
    NoiseAftFilter = NoiseAftFilter(NumFilter/2+1 : end);
    
    % 功率计算
    EnOfSignal  =  Signal_in(s_begin:s_end) * Signal_in(s_begin:s_end)';
    EnOfNoise = NoiseAftFilter* NoiseAftFilter';
    NorOfNoise = NoiseAftFilter / sqrt(EnOfNoise);          %噪声能量归一化
    
    %SignalAftChannel包含了保护间隔和LFM信号，但实际我们只在乎符号的，所以按symbol的能量扩充到接收信号
    AmpOfNoise = sqrt(10^(-SNR/10)*EnOfSignal*length(Signal_in)/(s_end-s_begin));   
    Noise = NorOfNoise * AmpOfNoise;
    SignalAftNoise = Signal_in + Noise;              %信号加噪声
    snrband = 20 * log10(std(Signal_in(s_begin:s_end))/std(Noise));      %带内信噪比测试    
end

