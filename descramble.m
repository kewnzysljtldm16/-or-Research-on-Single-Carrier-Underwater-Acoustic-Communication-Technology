% function descrambledSignal = descramble(scrambledSignal)
%     msig = [0 0 1 0 1 0 0 1 0 1 0 1 0 0 1];
%     n_m = length(msig);
%     N = length(scrambledSignal);
%     
%     % Initialization
%     a = zeros(1, N);
%     descrambledSignal = zeros(1, N);
%     
%     % Descrambling Operation
%     for i = 1:N
%         a(i) = xor(msig(1), msig(2));
%         
%         for j = 1:n_m - 1
%             msig(j) = msig(j + 1);
%         end
%         
%         msig(n_m) = a(i);
%         
%         descrambledSignal(i) = xor(scrambledSignal(i), a(i));
%     end
% end

function descrambled_signal = descramble(scrambled_signal)
    msig = [0 0 1 0 1 0 0 1 0 1 0 1 0 0 1]; % 高位到底位
    n_m = length(msig);

    N = length(scrambled_signal);
    descrambled_signal = zeros(1, N);

    for i = 1:N
        a = xor(msig(1), msig(2));

        for j = 1:n_m - 1
            msig(j) = msig(j + 1);
        end

        msig(n_m) = a;

        descrambled_signal(i) = xor(scrambled_signal(i), a);
    end
end