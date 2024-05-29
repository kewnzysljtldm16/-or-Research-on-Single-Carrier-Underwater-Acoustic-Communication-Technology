function [symbol_block]= IQdemodulate(signal_in,fs,signal_length,f0,PulseShape,N_up)

    %% IQ解调
    t1 = 0 : 1/fs : (signal_length)/fs-1/fs;
    sigI = signal_in .* cos(2*pi*f0*t1);                                         %I路解调
    sigQ = signal_in .* sin(2*pi*f0*t1);                                         %Q路解调
    sigI_lowpass = filter(PulseShape, 1, [sigI zeros(1,fix(length(PulseShape)/2))]);
    sigQ_lowpass = filter(PulseShape, 1, [sigQ zeros(1,fix(length(PulseShape)/2))]);
    sigI_lowpass = 2*sigI_lowpass(fix(length(PulseShape)/2)+1 : end);
    sigQ_lowpass = 2*sigQ_lowpass(fix(length(PulseShape)/2)+1 : end);
    
    %% 降采样
    N_symbol = signal_length / (N_up);
    for i = 1 : N_symbol
        shuchuI(1,i) =  sigI_lowpass((i - 1) * N_up + 1);                %I路
        shuchuQ(1,i) =  sigQ_lowpass((i - 1) * N_up + 1);                %Q路
    end
    symbol_block=shuchuI + 1i*shuchuQ;
   
end