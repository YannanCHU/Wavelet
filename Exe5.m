clc;
clear;
close all;

T = 64;
t = (0:2048-1)/T;
K = 2;

load samples.mat

phi = zeros(1, 2048);
[phi_T, psi_T, xval] = wavefun('db4', 6);
phi(1:length(phi_T)) = phi_T;
L = round(length(phi_T)/T);
n = T * (0:1:32-L);

Cmn = zeros(4, 32-L+1);
for m = 0:3
    for i = 1:length(n)
        Cmn(m+1, i) = dot((1/T) *  shift(phi, n(i)), t.^m);
    end
end

% Get the sequence s[m]
sm = Cmn * y_sampled(1:32-L+1)';

[tks, xks] = annihilating_func(K, sm');
display(tks);
display(xks);
xt_reconstructed = zeros(1, 2048);
xt_reconstructed(1, round(tks*64)+1) = xks;

figure(1);
plot(t, xt_reconstructed, '-*');     
title("Reconstructed continuous time Diracs", 'fontsize', 24);   xlim([0,32]);
xlabel("Time", 'fontsize', 12);
ylabel("Amplitude", 'fontsize', 12);
