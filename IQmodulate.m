function [signal] = IQmodulate(p,N_up,SymbolIn,f0,fs)

%%%% p:脉冲成型滤波器系数；N_pn：升采样点数；SymbolIn：映射符号（复数）；f0：载波频率；fs：采样频率；

%%%% IQ两路信号分别升采样，分别过脉冲成型滤波器
%%%% pams_I和pams_Q分别是是两路信号升采样后的信号
%%%% signal调制后信号

    %% 上采样
    pams_I = upsample(real(SymbolIn), N_up);
    pams_Q = upsample(imag(SymbolIn), N_up);
    
    %% 整形滤波器%
    pams_I = [pams_I zeros(1, fix(length(p)/2))];
    pams_Q = [pams_Q zeros(1, fix(length(p)/2))];
    figure
    subplot(2,1,1);
    plot(pams_I);
    title('pams_I');
    subplot(2,1,2);
    plot(pams_Q);
    title('pams_Q');
    
    ynI = filter(p, 1, pams_I);                                      %脉冲成形滤波
    ynQ = filter(p, 1, pams_Q);                                   %脉冲成形滤波
    ynI = ynI(fix(length(p)/2)+1 : end);
    ynQ = ynQ(fix(length(p)/2)+1 : end);
    figure
    subplot(2,1,1);
    plot(ynI);
    title('ynI');
    subplot(2,1,2);
    plot(ynQ);
    title('ynQ');
    
    %% 调制
    t = 0 : 1/fs : length(ynI)/fs-1/fs;
    y_cos = cos(2*pi*f0*t);
    y_sin = sin(2*pi*f0*t);
    signal = ynI.*y_cos + ynQ.*y_sin;                          %%sending signal
    
end



