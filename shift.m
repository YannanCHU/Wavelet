function A = shift(A, shift_pos)
    A(1+shift_pos:end) = A(1:end-shift_pos);
    A(1:shift_pos) = 0;
end