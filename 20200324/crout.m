function [ L, U ] = crout( A )
%CROUT Crout decomposition of square matrix A.
%   In Crout decomposition, we find a lower triangle matrix L and a unit
% upper triangle matrix U, s.t. A=L*U.
% param: 
%   A: n*n square matrix.
% return:
%   L: n*n lower triangle matrix;
%   U: n*n unit upper triangle matrix.
% 
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/03/27

% some pre-operation
N = length(A);
L = zeros(N);
U = eye(N);

% statistics
num_add = 0; % numbers of add operation
num_mul = 0; % numbers of multiply operation

% big loop, for the k-th layer:
for k = 1:N
    % first, calculate l_kk
    L(k,k) = A(k,k);
    for i = 1:(k-1)
        L(k,k) = L(k,k) - L(k,i)*U(i,k);
        num_add = num_add + 1;
        num_mul = num_mul + 1;
    end
    
    % second, cauculate u_kj and l_ik
    % we use 'ij' to unite 'i' and 'j' subscriptor
    for ij0 = (k+1):N
        U(k,ij0) = A(k,ij0);
        L(ij0,k) = A(ij0,k);        
        % //
        for ij1 = 1:(k-1)
            U(k,ij0) = U(k,ij0) - L(k,ij1)*U(ij1,ij0);
            L(ij0,k) = L(ij0,k) - L(ij0,ij1)*U(ij1,k);
            num_add = num_add + 2;
            num_mul = num_mul + 2;
        end
        % //
        U(k,ij0) = U(k,ij0)/L(k,k);
        num_mul = num_mul + 1;
    end % ij0
    
end % k

% display statistics
fprintf('Storage used: %d double unit(s)\n',N*N);
fprintf('Operation used: %d add(s) and %d multiply(s)\n',num_add,num_mul);

end
