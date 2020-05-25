function  omegaHat = EspritOnce( xn )
%ESPRITONCE Esprit frequency estimation for one time.
%       omegaHat = EspritOnce( xn )
% Param:
%   xn: input signal, column vector
% Return:
%   omegaHat: estimated normalize angular frequency, row vector.
% 
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/04/30

% self-related matrix
N = length(xn);         % n samples
M = 4;                  % m * m matrix
Xn = zeros(M, N-M);     % Xn = [x(1,2,...,m), x(2,3,...,m+1),...]
for k = 1:(N-M)
    Xn(:,k) = xn(k:k+M-1).';
end
Rxx = Xn(:, 1:end-1)*Xn(:, 1:end-1)' / (N-M-1); % x-x related
Rxy = Xn(:, 1:end-1)*Xn(:, 2:end)'   / (N-M-1); % x-y related

% estimate variances sigma2 = lambda_min(Rxx)
lambda_min = min(svd(Rxx));

% construct Cxx and Cxy
Z = [zeros(1,M);eye(M-1),zeros(M-1,1)];
Cxx = Rxx - lambda_min * eye(M);
Cxy = Rxy - lambda_min * Z;

% eigenvalue decomposition
z = eig(Cxx,Cxy);
omega = angle(z);
[~,index] = sort( abs(abs(z)-1) );   % sort eigenvalue closest to unit circle
omegaHat = sort( omega(index(1:3)) );% output form small to big

end

