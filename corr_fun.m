function yout=corr_fun(input1,input2)
%%%input1����input2.input1������input2ƽ����input1���
len1 = length(input1);
len2 = length(input2);
yout1 = xcorr(input1, input2);
yout = yout1(abs(len1-len2)+1 : end);
yout = yout ./ sqrt(sum(abs(input1.^2)) * sum(abs(input2.^2)));%%��һ��
end