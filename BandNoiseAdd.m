function [SignalAftNoise] = BandNoiseAdd(Signal_in,SNR,BandPass,s_begin,s_end)
% ����������Ӻ���
% �Ľ����֣�ֻ���ù�ע���ּ��������������������������������ʹ������������
% Signal_in�������źţ�SNR������ȣ�BandPass�������˲���ϵ����s_begin/s_end����ע���ֵĿ�ʼλ�úͽ���λ��

    NumFilter = length(BandPass)-1;
    if isequal(imag(Signal_in), zeros(1, length(Signal_in)))   %ʵ���źŵĽ�������
        noise = normrnd(0, 1, 1, length(Signal_in));      
    else
        noise = sqrt(1/2)*normrnd(0, 1, 1, length(Signal_in)) + 1i*sqrt(1/2)*normrnd(0, 1, 1, length(Signal_in)); %�����źŵĽ�������
    end
    
    % �����������
    NoiseAftFilter = filter(BandPass, 1, [noise zeros(1, NumFilter/2)]);           %�˲���1Ϊ��ĸ
    NoiseAftFilter = NoiseAftFilter(NumFilter/2+1 : end);
    
    % ���ʼ���
    EnOfSignal  =  Signal_in(s_begin:s_end) * Signal_in(s_begin:s_end)';
    EnOfNoise = NoiseAftFilter* NoiseAftFilter';
    NorOfNoise = NoiseAftFilter / sqrt(EnOfNoise);          %����������һ��
    
    %SignalAftChannel�����˱��������LFM�źţ���ʵ������ֻ�ں����ŵģ����԰�symbol���������䵽�����ź�
    AmpOfNoise = sqrt(10^(-SNR/10)*EnOfSignal*length(Signal_in)/(s_end-s_begin));   
    Noise = NorOfNoise * AmpOfNoise;
    SignalAftNoise = Signal_in + Noise;              %�źż�����
    snrband = 20 * log10(std(Signal_in(s_begin:s_end))/std(Noise));      %��������Ȳ���    
end

