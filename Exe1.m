clc;
clear;
close all;

T = 64;
t = (0:2048-1)/T;

phi = zeros(1, 2048);
[phi_T, psi_T, xval] = wavefun('db4', 6);
phi(1:length(phi_T)) = phi_T;
figure();   plot((1:length(phi_T))/T, phi_T, 'LineWidth', 2); title("The scaling function - db4", 'fontsize', 20);
xlabel("Time", 'fontsize', 12); ylabel("Amplitude", 'fontsize', 12)
% Since the sacling function is orthogonal to itself. There is no need to
% calculate its dual.
L = round(length(phi_T)/T);

n = T * (0:1:32-L);

Cmn = zeros(4, 32);
for m = 0:3
    for i = 1:length(n)
        Cmn(m+1, i) = dot((1/T) *  shift(phi, n(i)), t.^m);
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
