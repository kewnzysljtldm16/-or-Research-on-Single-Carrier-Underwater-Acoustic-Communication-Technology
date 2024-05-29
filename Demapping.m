function [ Output ] = Demapping( Input , Table , M)
% 功能：解映射
% 输入：
% Input：输入符号
% Table：映射星座图
% M：调制阶数

% 输出：
% Output：输出二进制比特
Q = log2(M);
LenInput = length(Input);
Output = Q * LenInput;

for i = 1 : LenInput
    
    if Q == 1 % BPSK
        Position = find(Input(i) == Table);
        if Position == 1
            Output(i) = 0;
        else
            Output(i) = 1;
        end
        
    elseif Q == 2 % QPSK
        Position = find(Input(i) == Table);
        if Position == 1
            Output(2*i-1) = 0;
            Output(2*i)   = 0;
        elseif Position == 2
            Output(2*i-1) = 1;
            Output(2*i)   = 0;
        elseif Position == 3
            Output(2*i-1) = 0;
            Output(2*i)   = 1;
        else
            Output(2*i-1) = 1;
            Output(2*i)   = 1;
        end
        
    elseif Q == 3 % 8PSK
        Position = find(Input(i) == Table);
        if Position == 1
            Output(3*i-2) = 0;
            Output(3*i-1) = 0;
            Output(3*i)   = 0;
        elseif Position == 2
            Output(3*i-2) = 1;
            Output(3*i-1) = 0;
            Output(3*i)   = 0;
        elseif Position == 3
            Output(3*i-2) = 0;
            Output(3*i-1) = 1;
            Output(3*i)   = 0;
        elseif Position == 4
            Output(3*i-2) = 1;
            Output(3*i-1) = 1;
            Output(3*i)   = 0;
        elseif Position == 5
            Output(3*i-2) = 0;
            Output(3*i-1) = 0;
            Output(3*i)   = 1;
        elseif Position == 6
            Output(3*i-2) = 1;
            Output(3*i-1) = 0;
            Output(3*i)   = 1;
        elseif Position == 7
            Output(3*i-2) = 0;
            Output(3*i-1) = 1;
            Output(3*i)   = 1;
        else
            Output(3*i-2) = 1;
            Output(3*i-1) = 1;
            Output(3*i)   = 1;
        end

    else  % 16QAM
        Position = find(Input(i) == Table);
        if Position == 1
            Output(4*i-3) = 0;
            Output(4*i-2) = 0;
            Output(4*i-1) = 0;
            Output(4*i)   = 0;
        elseif Position == 2
            Output(4*i-3) = 1;
            Output(4*i-2) = 0;
            Output(4*i-1) = 0;
            Output(4*i)   = 0;
        elseif Position == 3
            Output(4*i-3) = 0;
            Output(4*i-2) = 1;
            Output(4*i-1) = 0;
            Output(4*i)   = 0;
        elseif Position == 4
            Output(4*i-3) = 1;
            Output(4*i-2) = 1;
            Output(4*i-1) = 0;
            Output(4*i)   = 0;
        elseif Position == 5
            Output(4*i-3) = 0;
            Output(4*i-2) = 0;
            Output(4*i-1) = 1;
            Output(4*i)   = 0;
        elseif Position == 6
            Output(4*i-3) = 1;
            Output(4*i-2) = 0;
            Output(4*i-1) = 1;
            Output(4*i)   = 0;
        elseif Position == 7
            Output(4*i-3) = 0;
            Output(4*i-2) = 1;
            Output(4*i-1) = 1;
            Output(4*i)   = 0;
        elseif Position == 8
            Output(4*i-3) = 1;
            Output(4*i-2) = 1;
            Output(4*i-1) = 1;
            Output(4*i)   = 0;
        elseif Position == 9
            Output(4*i-3) = 0;
            Output(4*i-2) = 0;
            Output(4*i-1) = 0;
            Output(4*i)   = 1;
        elseif Position == 10
            Output(4*i-3) = 1;
            Output(4*i-2) = 0;
            Output(4*i-1) = 0;
            Output(4*i)   = 1;
        elseif Position == 11
            Output(4*i-3) = 0;
            Output(4*i-2) = 1;
            Output(4*i-1) = 0;
            Output(4*i)   = 1;
        elseif Position == 12
            Output(4*i-3) = 1;
            Output(4*i-2) = 1;
            Output(4*i-1) = 0;
            Output(4*i)   = 1;
        elseif Position == 13
            Output(4*i-3) = 0;
            Output(4*i-2) = 0;
            Output(4*i-1) = 1;
            Output(4*i)   = 1;
        elseif Position == 14
            Output(4*i-3) = 1;
            Output(4*i-2) = 0;
            Output(4*i-1) = 1;
            Output(4*i)   = 1;
        elseif Position == 15
            Output(4*i-3) = 0;
            Output(4*i-2) = 1;
            Output(4*i-1) = 1;
            Output(4*i)   = 1;
        else 
            Output(4*i-3) = 1;
            Output(4*i-2) = 1;
            Output(4*i-1) = 1;
            Output(4*i)   = 1;
        end

    end

end



end

