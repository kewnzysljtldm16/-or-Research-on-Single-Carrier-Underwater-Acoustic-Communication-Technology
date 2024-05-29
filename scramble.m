% % 加扰函数LFSR
% function scrambledSignal = scramble(inputSignal)
%     msig = [0 0 1 0 1 0 0 1 0 1 0 1 0 0 1];
%     n_m = length(msig);
%     N = length(inputSignal);
%     
%     % Initialization
%     count = 0;
%     a = zeros(1, N);
%     y = zeros(1, N);
%     
%     % Scrambling Operation
%     for i = 1:N
%         a(i) = xor(msig(1), msig(2));
%         
%         for j = 1:n_m - 1
%             msig(j) = msig(j + 1);
%         end
%         
%         msig(n_m) = a(i);
%         
%         y(i) = xor(inputSignal(i), a(i));
%         
%         if a(i) == 1
%             count = count + 1;
%         end
%     end
%     
%     p1 = count / N;
%     
%     % Output the scrambled signal
%     scrambledSignal = y;
% end

function scrambled_signal = scramble(input_signal)
    msig = [0 0 1 0 1 0 0 1 0 1 0 1 0 0 1]; % 高位到底位
    n_m = length(msig);

    % 产生原始序列
    N = length(input_signal);

    % 扰码
    count = 0;
    scrambled_signal = zeros(1, N);

    for i = 1:N
        a = xor(msig(1), msig(2));

        for j = 1:n_m - 1
            msig(j) = msig(j + 1);
        end

        msig(n_m) = a;

        scrambled_signal(i) = xor(input_signal(i), a);

        if a == 1
            count = count + 1;
        end
    end

    p1 = count / N;

%     % 画图
%     t = 1:100;
%     subplot(2, 1, 1);
%     stem(t, input_signal(1:100));
%     ylabel('sig');
%     title('Original Signal');
% 
%     subplot(2, 1, 2);
%     stem(t, scrambled_signal(1:100));
%     ylabel('msig');
%     title('Scrambled Signal');
end