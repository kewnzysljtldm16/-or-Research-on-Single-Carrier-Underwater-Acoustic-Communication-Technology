function [ Output , Table] = Mapping( InputBit , M )
% ���ܣ�����ͼӳ��
% ���룺
% InputBit����������Ʊ���
% M�����ƽ���
% ��������ӳ��
% �����
% Output���������
% Table������ͼ


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




