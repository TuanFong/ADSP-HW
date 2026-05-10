% Load Image
% Input your link to iamge here
a = ('image_A.jpg');
b = ('image_B.jpg');

try
    ImagA = imread(a);
    ImagB = imread(b);
catch
    error('Image not found');
end

L = 255; k1 = 0.01; k2 = 0.03;
C1 = (k1*L)^2; C2 = (k2*L)^2;

% SSIM calc - Imag A & B have the same size
if size(ImagA,3)==3, A = double(rgb2gray(ImagA)); else, A = double(ImagA); end
if size(ImagB,3)==3, B = double(rgb2gray(ImagB)); else, B = double(ImagB); end

w = fspecial('gaussian', 11, 1.5);
mu1 = filter2(w, A, 'valid');
mu2 = filter2(w, B, 'valid');
sigma1 = filter2(w, A.^2, 'valid') - mu1.^2;
sigma2 = filter2(w, B.^2, 'valid') - mu2.^2;
sigma12 = filter2(w, A.*B, 'valid') - mu1.*mu2;
map = ((2*mu1.*mu2 + C1).*(2*sigma12 + C2)) ./ ...
           ((mu1.^2 + mu2.^2 + C1).*(sigma1 + sigma2 + C2));
result = mean2(map);

% Plot
[h, w, c] = size(ImagA);
padding = 130;
canvas = uint8(255 * ones(h + padding, w * 2, c));
canvas(1:h, :, :) = [ImagA, ImagB];

figure('Name', 'SSIM Analysis', 'NumberTitle', 'off');
imshow(canvas); hold on;

text(20, h + 30, 'Image A (original)', 'FontSize', 12, 'FontWeight', 'bold');
text(w + 20, h + 30, 'Imag B (after reconstructed)', 'FontSize', 12, 'FontWeight', 'bold');

info_str = {sprintf('SSIM value: %.4f', result), ...
            sprintf('Adjust constant c1: %.2f', C1), ...
            sprintf('Adjust constant c2: %.2f', C2)};
text(20, h + 90, info_str, 'FontSize', 11, 'EdgeColor', 'red', ...
     'BackgroundColor', 'white', 'LineWidth', 1.2);

hold off;