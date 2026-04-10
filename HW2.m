clc; clear; close all;

k = input('Enter k = ');
N = 2*k + 1;

% Freq Sampling
H = zeros(1, N);
for n = 1:N
    if n == 1
        H(n) = 0;        
    elseif n <= (N+1)/2
        H(n) = -1j;    
    else
        H(n) = 1j;     
    end
end
h = real(ifft(H));
h = fftshift(h);
n = -k:k;

% Plot
figure;
stem(n, h, 'filled');
title(['Impulse Response h[n]']);
xlabel('n'); ylabel('h[n]');
grid on;
w = linspace(-pi, pi, 1024);
H_w = zeros(size(w));
for i = 1:length(w)
    H_w(i) = sum(h .* exp(-1j*w(i)*n));
end
figure;
plot(w/pi, imag(H_w), 'LineWidth', 1.5);
title('Imaginary Part of Frequency Response');
xlabel('\omega / \pi'); ylabel('Imag(H)');
grid on;