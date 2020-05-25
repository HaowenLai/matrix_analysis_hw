function [ Q, R ] = qr_givens( A )
%QR_GIVENS QR decomposition using Givens method.
%   This function uses Given method to find QR decompositon of A, where Q
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

[m,n]=size(A);
R=A;
Q=eye(m);

% for any diagonal element
for i=1:m-1
    % it is a tall matrix
    if i>n,break,end;
    
    % big loop, using givens
    for j=2:m-i+1
        x=R(i:m,i);
        t=my_givens(x,1,j);
        r=blkdiag(eye(i-1),t);
        Q=Q*r';
        R=r*R;
    end
end

end


function [R,y] = my_givens(x,i,j)
% Find the transformation matrix R, s.t. R*x=y, in which 
% y(j)=0£¬y(i)=sqrt(x(i)^2+x(j)^2)
%
% Param:
%   x£ºcolumn vector needed to transform
%   i£ºsqrt(x(i)^2+x(j)^2) subscriptor
%   j£º0 element subscriptor
% Return:
%   R£ºGivens transform matrix
%   y£ºGivens transform result

xi=x(i);
xj=x(j);
r=sqrt(xi^2+xj^2);
cost=xi/r;
sint=xj/r;
R=eye(length(x));

% rotate
R(i,i)=cost;
R(i,j)=sint;
R(j,i)=-sint;
R(j,j)=cost;
y=x(:);
y([i,j])=[r,0];

end

