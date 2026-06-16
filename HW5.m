function [Fx, Fy] = fftreal(x, y)
    x = x(:);   
    y = y(:);
    N = length(x);

    % z = x + j*y
    z = x + 1j * y;

    % ONE N-point FFT
    Z = fft(z);
    idx = [1; (N:-1:2)'];
    Z_conj_mirror = conj( Z(idx) );

    Fx = (Z + Z_conj_mirror) / 2;
    Fy = (Z - Z_conj_mirror) / (2j);

end
% Test parameters
N  = 8;
rng(42);                      
x  = randn(N, 1);
y  = randn(N, 1);

fprintf('fftreal: FFT of two real signals using one FFT\n');

fprintf('Input signals (N = %d):\n', N);
fprintf('  x = '); fprintf('%8.4f  ', x); fprintf('\n');
fprintf('  y = '); fprintf('%8.4f  ', y); fprintf('\n\n');

[Fx, Fy] = fftreal(x, y);
Fx_ref = fft(x);
Fy_ref = fft(y);

fx_ok = max(abs(Fx - Fx_ref)) < 1e-10;
fy_ok = max(abs(Fy - Fy_ref)) < 1e-10;

fprintf('Fx (fftreal):\n');
disp(Fx);
fprintf('Fx (reference fft):\n');
disp(Fx_ref);
fprintf('Fx correct: %d\n\n', fx_ok);

fprintf('Fy (fftreal):\n');
disp(Fy);
fprintf('Fy (reference fft):\n');
disp(Fy_ref);
fprintf('Fy correct: %d\n\n', fy_ok);

if fx_ok && fy_ok
    fprintf('PASS: fftreal matches MATLAB fft reference.\n');
else
    fprintf('FAIL: results do not match.\n');
end
