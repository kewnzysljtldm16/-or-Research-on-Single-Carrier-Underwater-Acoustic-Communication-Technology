function yout=corr_fun(input1,input2)
%%%input1长于input2.input1不动，input2平移与input1相乘
len1 = length(input1);
len2 = length(input2);
yout1 = xcorr(input1, input2);
yout = yout1(abs(len1-len2)+1 : end);
yout = yout ./ sqrt(sum(abs(input1.^2)) * sum(abs(input2.^2)));%%归一化
end