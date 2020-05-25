function [ Q, R ] = qr_householder( A )
%QR_HOUSEHOLDER QR decomposition using Householder method.
%   This function uses Householder method to find QR decompositon of A, where Q
% is a orthogonal matrix and R is a upper triangle matrix, s.t. A=Q*R.
% 
% param: 
%   A: m*n matrix.
% return:
%   Q: m*n orthogonal matrix;
%   R: n*n upper triangle matrix.
% 
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/03/27

[m,n] = size(A);
R=A;
Q=eye(m);

for i=1:m-1
    % it is a tall matrix
    if i>n,break,end;
    
    % householder transform
    x=R(i:m,i);
    y=[1;zeros(m-i,1)];
    Ht=my_householder(x,y);
    H=blkdiag(eye(i-1),Ht);
    Q=Q*H;
    R=H*R;
end

end


function [H,r] = my_householder(x,y)
% Find Householder matrix H, s.t. Hx=r*y, where r=-sign(x(1))*||x||/||y||
%
% Param:
%   x£ºcolumn vector
%   y£ºcolumn vector
%

x=x(:);
y=y(:);

% compute orientation v and transform H
r=-sign(x(1))*norm(x)/norm(y);
y=r*y;
v=(x-y)/norm(x-y);
I=eye(length(x));
H=I-2*v*v';

end
