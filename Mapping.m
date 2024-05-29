function [ Output , Table] = Mapping( InputBit , M )
% 功能：星座图映射
% 输入：
% InputBit：输入二进制比特
% M：调制阶数
% 带格列码映射
% 输出：
% Output：输出符号
% Table：星座图


load TableCon.mat
Q = log2(M);
ReshapeSeq = reshape(InputBit , Q , length(InputBit)/Q).'; 

if Q == 1
    Table = TableBPSK;
    DecadeSeq = ReshapeSeq(: , 1);
elseif Q == 2
    Table = TableQPSK; 
    DecadeSeq = ReshapeSeq(: , 1) + ReshapeSeq(: , 2)*2; 
elseif Q == 3
    Table = Table8PSK;
    DecadeSeq = ReshapeSeq(: , 1) + ReshapeSeq(: , 2)*2 + ReshapeSeq(: , 3)*4; 
else
    Table = Table16QAM;
    DecadeSeq = ReshapeSeq(: , 1) + ReshapeSeq(: , 2)*2 + ReshapeSeq(: , 3)*4 + ReshapeSeq(: , 4)*8; 
end

Output = Table(DecadeSeq+1); 

end




