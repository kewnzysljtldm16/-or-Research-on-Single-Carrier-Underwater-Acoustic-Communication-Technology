function [signal] = IQmodulate(p,N_up,SymbolIn,f0,fs)

%%%% p:��������˲���ϵ����N_pn��������������SymbolIn��ӳ����ţ���������f0���ز�Ƶ�ʣ�fs������Ƶ�ʣ�

%%%% IQ��·�źŷֱ����������ֱ����������˲���
%%%% pams_I��pams_Q�ֱ�������·�ź�����������ź�
%%%% signal���ƺ��ź�

    %% �ϲ���
    pams_I = upsample(real(SymbolIn), N_up);
    pams_Q = upsample(imag(SymbolIn), N_up);
    
    %% �����˲���%
    pams_I = [pams_I zeros(1, fix(length(p)/2))];
    pams_Q = [pams_Q zeros(1, fix(length(p)/2))];
    figure
    subplot(2,1,1);
    plot(pams_I);
    title('pams_I');
    subplot(2,1,2);
    plot(pams_Q);
    title('pams_Q');
    
    ynI = filter(p, 1, pams_I);                                      %��������˲�
    ynQ = filter(p, 1, pams_Q);                                   %��������˲�
    ynI = ynI(fix(length(p)/2)+1 : end);
    ynQ = ynQ(fix(length(p)/2)+1 : end);
    figure
    subplot(2,1,1);
    plot(ynI);
    title('ynI');
    subplot(2,1,2);
    plot(ynQ);
    title('ynQ');
    
    %% ����
    t = 0 : 1/fs : length(ynI)/fs-1/fs;
    y_cos = cos(2*pi*f0*t);
    y_sin = sin(2*pi*f0*t);
    signal = ynI.*y_cos + ynQ.*y_sin;                          %%sending signal
    
end



