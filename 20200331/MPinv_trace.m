function [ Ap ] = MPinv_trace( A )
%MPINV_TRACE The Moore-Penrose inverse of A using trace method.
% Param:
%   A: m*n input matrix.
% Return:
%   Ap: Moore-Penrose inverse of A
% 
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/04/03

[~, n] = size(A);
r = rank(A);

% according to each step
% Step 1
B = A'*A;

% Step 2
C = eye(n);

% Step 3
for k = 1:r-1
    C = trace(C*B)/k*eye(n) - C*B;
end

% Step 4
Ap = r/trace(C*B) * C * A';

end

