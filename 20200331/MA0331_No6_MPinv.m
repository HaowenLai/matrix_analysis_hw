% Matrix Analysis problem No.6
% Moore-Penrose inverse problem
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/04/03

clear;clc;
load('MA0331_No6_MPinv.mat');

Ap_col   = MPinv_col(A);
Ap_trace = MPinv_trace(A);
Ap_std   = pinv(A); % NOTE: standard MP inverse, used for analysis

fprintf('Error between column iteration and standard methods:\n');
fprintf('   ||Ap_col - Ap_std|| = %e\n',norm(Ap_col-Ap_std));
fprintf('Error between trace and standard methods:\n');
fprintf('   ||Ap_trace - Ap_std|| = %e\n\n',norm(Ap_trace-Ap_std));

% ----------------- END OF FILE --------------------
