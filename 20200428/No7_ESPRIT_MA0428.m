% Matrix Analysis problem No.7
% ESPRIT frequency estimation problem.
% Author: Haowen Lai, lhw19@mails.tsinghua.edu.cn
% Date: 2020/04/30

clear;clc;

omega = [0.12*pi, 0.37*pi, 0.72*pi];  % normalize angle frequency
K = 200;                              % test K times

% ------------------ Question (1) -----------------------
N     = 100;   % sample number of signal
nStep = 50;    % how many points to plot
sigma2s  = linspace(0.01,1.5,nStep);  % list of sigma2
OmegaHat = zeros(K,3);                % each K independent experiment has a row
MSE1     = zeros(1,nStep);            % for each nStep

% Main big loop
disp('------- start processing question (1) ---------');
for n_step = 1:nStep
    for i = 1:K
        OmegaHat(i,:) = EspritOnce(GenerateSignal( N, sigma2s(n_step) ));
    end
    MSE1(n_step) = sum(sum( (OmegaHat - repmat(omega,K,1)).^2 )) / K;
    
    % output process status
    fprintf('  %d point(s) finished. Total %d points\n',n_step,nStep);
end


% ------------------ Question (2) -----------------------
sigma2 = 0.01;                             % variance
N      = round(linspace(40,1000,nStep));   % list of N samples
MSE2   = zeros(1,nStep);                   % for each nStep

% Main big loop
disp('------- start processing question (2) ---------');
for n_step = 1:nStep
    for i = 1:K
        OmegaHat(i,:) = EspritOnce(GenerateSignal( N(n_step), sigma2 ));
    end
    MSE2(n_step) = sum(sum( (OmegaHat - repmat(omega,K,1)).^2 )) / K;
    
    % output process status
    fprintf('  %d point(s) finished. Total %d points\n',n_step,nStep);
end
disp('All done. Ploting...');


% -----------------   Plot   -----------------------
subplot(1,2,1)
plot(sigma2s,MSE1)
title('N=100,  MSE - \sigma^2 plot')
xlabel('\sigma^2')
ylabel('MSE')
% 
subplot(1,2,2)
plot(N,MSE2)
title('\sigma^2=0.01,  MSE - N plot')
xlabel('N')
ylabel('MSE')

% ----------------- END OF FILE --------------------
