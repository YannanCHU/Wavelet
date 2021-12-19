clc
clear
close all

syms y z
y = ((1-z)/2) * ((1-z^(-1))/2);
T = 64;
t = (0:2048-1)/T;

% we want to get reproduce polynomials with maximum degree of 3
% Thus, the B-spline with order 4 should be used
p = 4;
% Due to the half-band condition:
% P(z) + P(z^-1) = 2
Pz = 0;
for k=0:1:(p-1)
    Pz = Pz + nchoosek(p+k-1,k) * y^k;
end
Pz = 2 * (1-y)^p * Pz;

% Factorization
Pz_factors = factor(Pz);

betai = bspline(2, 3);
g0 = betai;

G0z = 0;
for i = 1:length(betai)
    G0z = G0z + betai(i) * (z^(-(i-1)));
end

% G0z = ztrans(phi);
H0z = Pz / G0z;
H0z_expand = expand(H0z);
[numH0z, denH0z] = numden(H0z_expand);
numH0zCoef = sym2poly(numH0z);
denH0zCoef = sym2poly(denH0z);
h0 = numH0zCoef / denH0zCoef(1);

% figure(1);
% subplot(2,1,1);      plot((0: length(g0)-1)/2, g0);   title("g_0^{(1)}");
% subplot(2,1,2);      plot((0: length(h0)-1)/2, h0);   title("h_0^{(1)}");

% In order to get the resolution of 1/64, 6 iterations are required (from 0 to 5)
Sa = h0;    % analysis scaling function
g0=[0 0 0 g0 0 0 0];
Ss = g0;    % synthesis sacling function

for iter = 1:5
    h0_up = upsample(h0,2^iter);
    h0_up(end-(2^iter)+2: end) = [];
    Sa = conv(h0_up, Sa);
    
    g0_up = upsample(g0,2^iter);
    g0_up(end-(2^iter)+2: end) = [];
    Ss = conv(g0_up, Ss);
end

figure(1);
subplot(2,1,1);      plot((0: length(Ss)-1)/(T), Ss, 'LineWidth', 2);   title("g_0^{(6)} - scaling function of B-Spline", 'fontsize', 20);
xlabel("Time", 'fontsize', 12);      ylabel("Amplitude", 'fontsize', 12);
subplot(2,1,2);      plot((0: length(Sa)-1)/(T), Sa, 'LineWidth', 2);   title("h_0^{(6)} - dual of scaling function of B-Spline", 'fontsize', 20);
xlabel("Time", 'fontsize', 12);      ylabel("Amplitude", 'fontsize', 12);
%% Signal Reproduction
phi = zeros(1, 2048);
phi(1,1:length(Ss)) = Ss;
phi_dual = zeros(1, 2048);
phi_dual(1,1:length(Sa)) = Sa;

L_dual = floor(length(Sa)/T);
n_dual = T * (0:1:32-L_dual);
L = floor(length(Ss)/T);
n = T * (0:1:32-L);

Cmn = zeros(4, 32);
for m = 0:3
    for i = 1:length(n_dual)
        Cmn(m+1, i) = dot(shift(phi_dual, n_dual(i)), t.^m);
    end
end

poly_recovered = zeros(4, 2048);
kernel_shifted = zeros(4, 2048);
for m = 0:3
    for i = 1:length(n)
        poly_recovered(m+1, :) = poly_recovered(m+1, :) + Cmn(m+1, i) * shift(phi, n(i));
    end
end

for m = 0:3
    figure();
    plot(t, t.^m, 'g--', t, poly_recovered(m+1,:), 'r', 'LineWidth', 2)
    xlim([1 32])
    hold on
    for i = 1:length(n)
        plot(t, Cmn(m+1, i) * shift(phi, n(i)), 'b', 'LineWidth', 1);
    end
    hold off
    legend("Original t^"+m, "Reproduced polynomials", "Scaled and shifted kernels");
    title("Polynomial of degree " + m);
    xlabel("Time")
end

function A = shift(A, shift_pos)
    A(1+shift_pos:end) = A(1:end-shift_pos);
    A(1:shift_pos) = 0;
end

function betai = bspline(T, order)
beta0 = ones(1,T);
betai = beta0;

for i = 1: order
    betai = conv(betai, beta0) / T;
end

% figure();
% plot(betai);
end