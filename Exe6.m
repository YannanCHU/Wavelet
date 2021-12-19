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

phi = zeros(1, 2048);
[phi_T, psi_T, xval] = wavefun('db5', 6);
phi(1:length(phi_T)) = phi_T;
L = round(length(phi_T)/T);
n = T * (0:1:32-L);

Cmn = zeros(4, 32-L+1);
for m = 0:3
    for i = 1:length(n)
        Cmn(m+1, i) = dot((1/T) *  shift(phi, n(i)), t.^m);
    end
end

% Sampling
yn = zeros(size(Cmn, 1),1);

for i = 1:length(n)
    yn(i) = shift(phi, n(i)) * xt';
end

% Get the sequence s[m]
sm = Cmn * yn;
% Add the White Gaussian Noise
seed = 3;
rng(seed);
noise_var = [0.001, 0.01, 0.1, 1];
sm_hat = zeros(length(sm), length(noise_var));
for i = 1:1:length(noise_var)
    sm_hat(:, i) = sm + sqrt(noise_var(i)) * randn(size(sm));
end

xt_reconstructed = zeros(4, 2048);
for i = 1:1:4
    [tks, xks] = annihilating_func(K, (sm_hat(:, i))');
    xt_reconstructed(i, round(tks*64)+1) = xks;
end

figure(1);      
plot(t, xt, '-o', t, xt_reconstructed, '-*', 'LineWidth', 2);     
title("Original Noiseless and Reconstructed continuous time Diracs (without noise reduction)", 'fontsize', 24);   xlim([0,32]);
legend("Original noiseless signal", ...
       "\sigma^2 = " + noise_var(1), "\sigma^2 = " + noise_var(2), ...
       "\sigma^2 = " + noise_var(3), "\sigma^2 = " + noise_var(4), 'fontsize', 12)
xlabel("Time", 'fontsize', 16);
ylabel("Amplitude", 'fontsize', 16);