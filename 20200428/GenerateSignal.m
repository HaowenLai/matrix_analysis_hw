function xn = GenerateSignal( nSample, sigma2 )
%GENERATESIGNAL Generate signal xn with white Gaussian noise.
%   The signal is nSample long and is defined by constants below, with
% noise that has 0 mean and sigma2 variance:
%     s = [1.31*exp(pi/4*1j), 2.07*exp(pi/3*1j), 1.88*exp(pi/5*1j)]; 
%     omega = [0.12*pi, 0.37*pi, 0.72*pi];
%     x(n) = s*exp(1j * omega' * n) + wgn(nSample,1,10*log10(sigma2))
% 
%   xn = GenerateSignal( nSample, sigma2 ), colunm vector
% 
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/04/30

% some constants
s = [1.31*exp(pi/4*1j), 2.07*exp(pi/3*1j), 1.88*exp(pi/5*1j)];  % signal magnitude
omega = [0.12*pi, 0.37*pi, 0.72*pi];         % normalize angle frequency
xn_pure = @(n)( s*exp(1j * omega' * n) );    % without noise

% generate signal
xn = zeros(nSample, 1);
for n = 1:nSample
    xn(n) = xn_pure(n);
end
xn = xn + wgn(nSample, 1, 10*log10(sigma2), 'complex');  % add noise

end
