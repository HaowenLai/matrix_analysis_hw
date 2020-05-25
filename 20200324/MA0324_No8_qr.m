% Matrix Analysis problem No.8
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/03/27

clc;clear;
load('MA0324_No8_qr.mat');

% two QR methods
[Q_g,R_g]=qr_givens(A);
[Q_h,R_h]=qr_householder(A);

% analysis
disp('-------------Analysis-------------')
disp('Givens method:')
fprintf('  norm ||A-Q_g*R_g||: %e\n',norm(A-Q_g*R_g))
fprintf('  norm ||I-Q_g*Q_g''||: %e\n',norm(eye(length(Q_g))-Q_g*Q_g'))
disp('Householder method:')
fprintf('  norm ||A-Q_h*R_h||: %e\n',norm(A-Q_h*R_h))
fprintf('  norm ||I-Q_h*Q_h''||: %e\n',norm(eye(length(Q_h))-Q_h*Q_h'))

