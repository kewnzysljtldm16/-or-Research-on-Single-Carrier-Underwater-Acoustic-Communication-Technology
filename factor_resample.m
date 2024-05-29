function [Precision] = factor_resample(fs)
%% 计算fs1的最小resample不报错的可约因数
% fs1和fs2经过约去公因数
% 需要保证fs和fs1被同时约分后的乘积小于resample函数的上线，即q*p<2^31
% 而多普勒下重采样的频率于采样率总体上差距不大 因此可以近似看作q^2<2^31，即q<sqrt(2^31)或p<sqrt(2^31)
% 通常采样率为规整整数，能被2，5整除，因此限制fs2的也要能被2，5整除
% 部分能被3整除，但不具有普遍性，这部分需要考虑采样率的因式
% 能被越小值整除的时候先考虑小数整除，因为这样的精度更高

factor_fs = unique(factor(fs));                 %获取fs的因数
factor_choose = ceil(fs / sqrt(2^31));                     %获取最小的resample不报错约分倍率

%寻找最靠近不报错约分倍数的fs因数，这个因数就是最合适的不报错约分倍率
index = find(factor_fs >= factor_choose, 1, 'first');         % 1表示只对第一个这样的值感兴趣               
Precision = factor_fs(index);
end


