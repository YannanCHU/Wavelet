clc;
clear;
close all;

T = 64;
t = (0:2048-1)/T;

load tau.mat

K = 2;
[tks, xks] = annihilating_func(K, tau)

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

