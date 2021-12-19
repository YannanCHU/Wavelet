clc;
clear;
close all;

T = 64;
t = (0:2048-1)/T;

% Creat a stream of Diracs
K = 2;
tk = [8, 18];
xk = [1.6, 2.4];
xt = zeros(1, 2048);
% Since K = 2, there should have 2 non-zero diracs
xt(T * tk(1) + 1) = xk(1);
xt(T * tk(2) + 1) = xk(2);

% N > 2K = 4
Ns = 5:1:8;
noise_var = [0.001, 0.01, 0.1, 1];
seed = 0;       rng(seed);      randomNumbers = randn(max(Ns), 1);
sm_noisy_cell = cell(4,1);      % each cell stores sm corresponding to different Ns

% generate the noisy sm
for N = Ns
    waveName = "db"+ N;
    phi = zeros(1, 2048);
    [phi_T, psi_T, xval] = wavefun(waveName, 6);
    phi(1:length(phi_T)) = phi_T;
    L = round(length(phi_T)/T);
    n = T * (0:1:32-L);
    
    Cmn = zeros(4, 32-L+1);
    for m = 0:N-1
        for i = 1:length(n)
            Cmn(m+1, i) = dot((1/T) *  shift(phi, n(i)), t.^m);
        end
    end
    yn = zeros(size(Cmn, 1),1);
    for i = 1:length(n)
        yn(i) = shift(phi, n(i)) * xt';
    end
    
    % Get the sequence s[m]
    sm = Cmn * yn;
    % Add the White Gaussian Noise
    
    sm_hat = zeros(length(sm), length(noise_var));
    for i = 1:1:length(noise_var)
        epsilonm = sqrt(noise_var(i)) * randomNumbers(1:N);
        sm_hat(:, i) = sm + epsilonm;
    end
    
    sm_noisy_cell{N-min(Ns)+1} = sm_hat;
end
%%
for i = 1:1:length(noise_var)
    figure();
    for N = Ns
        j = N-min(Ns)+1;
        sm_hat = sm_noisy_cell{j};
        [tks_TLS, xks_TLS] = annihilating_TLS(K, (sm_hat(:, i))');
        [tks_Cadzow, xks_Cadzow] = annihilating_TLS_Cadzow(K, (sm_hat(:, i))');
        subplot(4,1,j);
        stem(tk, xk, '-o', 'LineWidth', 1);  hold on;
        title("Original and Reconstructed continuous time Diracs (N = " + N + ", \sigma^2 = " + noise_var(i) + ")");
        stem(tks_TLS, xks_TLS, '-*', 'LineWidth', 1);
        stem(tks_Cadzow, xks_Cadzow, '-x', 'LineWidth', 1);
        legend("Original Noiseless signal", "TLS results", "Cadzow results");
        xlabel("Time"); ylabel("Amplitude");
        hold off;
    end
end


function  [tks, xks] = annihilating_TLS_Cadzow(K, sm)
N = length(sm)-1;

% Construct the matrix of taus
smMatrix = ones(N-K, K+1);
for i = 1:1:K+1
    smMatrix(:, i) = sm(K-i+2:N-i+1);
end

% Cadzow
times = 100;
for iter = 1:times
    [U, Lambda,V] = svd(smMatrix);
    LambdaSorted = sort(diag(Lambda), 'descend');
    % retain K singular values while force all others as 0 by using the Kth
    % singular value as a threshold.
    Lambda2 = Lambda .* (Lambda >= LambdaSorted(K));
    smMatrix2 = U * Lambda2 * V';
    
    % take the average along diagonals
    smMatrix3 = zeros(size(smMatrix2));
    for k = 1-size(smMatrix2,1):1:0
        
        diag_ele_num = length(diag(smMatrix2, k));
        diag_k_mean = mean(diag(smMatrix2, k));
        diag_k = diag(diag_k_mean * ones(diag_ele_num,1));
        smMatrix3(-k+1:-k+1+diag_ele_num-1,1:diag_ele_num) = smMatrix3(-k+1:-k+1+diag_ele_num-1,1:diag_ele_num) + diag_k;
    end
    
    for k = 1:1:size(smMatrix2,2)-1
        
        diag_ele_num = length(diag(smMatrix2, k));
        diag_k_mean = mean(diag(smMatrix2, k));
        diag_k = diag(diag_k_mean * ones(diag_ele_num,1));
        smMatrix3(1:diag_ele_num, k+1:k+1+diag_ele_num-1) = smMatrix3(1:diag_ele_num, k+1:k+1+diag_ele_num-1) + diag_k;
    end
    
    smMatrix = smMatrix3;
    if (rank(smMatrix) == K)
%         disp("Iteration Times: " + iter);
        break;
    end
end
% rank(smMatrix)

[~, ~,V] = svd(smMatrix);
H = V(:, end);

% find the distinct tks
tks = (roots(H))';

% Construct the matrix of time
tMatrix = ones(K, K);
for j = 2:K
    tMatrix(j,:) = tks .^ (j-1);
end

% solve the amplitude xks of detected time tks
xks = tMatrix \ (sm(1:K)');
xks = xks';
end

function  [tks, xks] = annihilating_TLS(K, sm)
N = length(sm)-1;

% Construct the matrix of taus
smMatrix = ones(N-K, K+1);
for i = 1:1:K+1
    smMatrix(:, i) = sm(K-i+2:N-i+1);
end

[~, ~,V] = svd(smMatrix);
H = V(:, end);

% find the distinct tks
tks = (roots(H))';

% Construct the matrix of time
tMatrix = ones(K, K);
for j = 2:K
    tMatrix(j,:) = tks .^ (j-1);
end

% solve the amplitude xks of detected time tks
xks = tMatrix \ (sm(1:K)');
xks = xks';
end
