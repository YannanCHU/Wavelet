clc;
clear;
close all;

T = 64;
t = (0:2048-1)/T;

% Creat a stream of Diracs
K = 2;
xt = zeros(1, 2048);
% Since K = 2, there should have 2 non-zero diracs
xt(T * 8 + 1) = 1.6;
xt(T * 18 + 1) = 2.4;

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

% Sampling
yn = zeros(size(Cmn, 1),1);

for i = 1:length(n)
    yn(i) = shift(phi, n(i)) * xt';
end

% Get the sequence s[m]
sm = Cmn * yn;

% poly_recovered = zeros(4, 2048);
% sm = zeros(size(Cmn, 1),1);
% for m = 0:3
%     for i = 1:length(n)
%         poly_recovered(m+1, :) = poly_recovered(m+1, :) + Cmn(m+1, i) * shift(phi, n(i));
%     end
%     sm(m+1) = dot(xt, poly_recovered(m+1, :));
% end

[tks, xks] = annihilating_func(K, sm');
xt_reconstructed = zeros(1, 2048);
xt_reconstructed(1, round(tks*64)+1) = xks;

figure(1);
subplot(2,1,1);     plot(t, xt, '-*');      
title("Original stream of Diracs", 'fontsize', 24); xlim([0,32]);
xlabel("Time", 'fontsize', 12);
ylabel("Amplitude", 'fontsize', 12);
subplot(2,1,2);    plot(t, xt_reconstructed, '-*');     
title("Reconstructed continuous time Diracs", 'fontsize', 24);   xlim([0,32]);
xlabel("Time", 'fontsize', 12);
ylabel("Amplitude", 'fontsize', 12);
%% functions
function A = shift(A, shift_pos)
    A(1+shift_pos:end) = A(1:end-shift_pos);
    A(1:shift_pos) = 0;
end

function  [tks, xks] = annihilating_func(K, tau)
N = length(tau)-1;

% Construct the matrix of taus
tauMatrix = ones(N-K+1, K);
for i = 1:1:K
    tauMatrix(:, i) = tau(K-i+1:N-i+1);
end

% Construct the filter h by solving tauMatrix * h = tau
% h[0] = h(1) = 1;
h = ones(K+1,1);
h(2: K+1) = tauMatrix \ (- tau(K+1: N+1)');

% find the distinct tks
tks = (roots(h))';

% Construct the matrix of time
tMatrix = ones(K, K);
for j = 2:K
    tMatrix(j,:) = tks .^ (j-1);
end

% solve the amplitude xks of detected time tks
xks = tMatrix \ (tau(1:K)');
xks = xks';
end
