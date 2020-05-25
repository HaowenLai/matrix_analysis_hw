function [ Ap ] = MPinv_col( A )
%MPINV_COL The Moore-Penrose inverse of A using column iteration method.
% Param:
%   A: m*n input matrix.
% Return:
%   Ap: Moore-Penrose inverse of A
% 
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/04/03

[m, n] = size(A);
Ap = zeros(n, m);

% initial A1p
Ap(1,:) = A(:,1)'/(A(:,1)'*A(:,1));

for k = 2:n
    d = Ap(1:k-1,:)*A(:,k);
    
    a_temp = A(:,k)-A(:,1:k-1)*d; % test condition
    if norm(a_temp) < 1e-10
        b = d' * Ap(1:k-1,:) / (1+d'*d);
    else
        b = a_temp' / (a_temp'*a_temp);
    end
    
    % enlarge Ap
    Ap(1:k,:) = [Ap(1:k-1,:)-d*b;b];
end

end

