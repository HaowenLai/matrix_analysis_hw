% Matrix Analysis problem No.7
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/03/27

clc;
clear;
load('MA0324_No7_crout.mat');

[L, U] = crout(A);

fprintf('Error ||A-L*U|| = %e\n\n',norm(A-L*U));

