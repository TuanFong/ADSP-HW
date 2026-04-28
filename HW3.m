function Main()
    % read image
    input_image = 'peppers.png';
    try
        A = imread(input_image);
    catch
        fprintf('image not found, try with other image!\n');
        return;
    end
   
    [B, size_original, size_comp] = C420(A);
    
    % Plot
    figure('Name', '4:2:0 Image Compression Technique');
   
    subplot(1,2,1); 
    imshow(A); 
    title('Input color image (A)');
    xlabel(['Original Size: ', num2str(size_original), ' Bytes']);
    
    subplot(1,2,2); 
    imshow(B); 
    title('Reconstructed image B = C420(A)');
    xlabel(['Reconstructed Size: ', num2str(size_comp), ' Bytes']);
end

%B = C420(A)
function [B, size_original, size_comp] = C420(A)
    [rows, cols, ~] = size(A);
    
    % Size Calc
    size_original = rows * cols * 3;
    
    A_double = double(A);
    % Change RGB to YCbCr
    R = A_double(:,:,1); G = A_double(:,:,2); B_ch = A_double(:,:,3);
    Y  = 16  + ( 65.481*R + 128.553*G +  24.966*B_ch) / 256;
    Cb = 128 + (-37.797*R -  74.203*G + 112.000*B_ch) / 256;
    Cr = 128 + (112.000*R -  93.786*G -  18.214*B_ch) / 256;
    Cb_sub = Cb(1:2:end, 1:2:end);
    Cr_sub = Cr(1:2:end, 1:2:end);
    
    % Size Calc after compression
    size_comp = numel(Y) + numel(Cb_sub) + numel(Cr_sub);
    
    % Interpolation
    Cb_up = Interpolation(Cb_sub, rows, cols);
    Cr_up = Interpolation(Cr_sub, rows, cols);
    
    % Change YCbCr to RGB
    Y_norm = Y - 16;
    Cb_norm = Cb_up - 128;
    Cr_norm = Cr_up - 128;
    
    R_out = (298.082 * Y_norm + 408.583 * Cr_norm) / 256;
    G_out = (298.082 * Y_norm - 100.291 * Cb_norm - 208.120 * Cr_norm) / 256;
    B_out = (298.082 * Y_norm + 516.411 * Cb_norm) / 256;
    
    R_out = max(0, min(255, R_out));
    G_out = max(0, min(255, G_out));
    B_out = max(0, min(255, B_out));
    
    B = uint8(cat(3, R_out, G_out, B_out));
end

% Interpolation function
function out = Interpolation(img_sub, target_h, target_w)
    [h, w] = size(img_sub);
    out = zeros(target_h, target_w);
    sh = h / target_h;
    sw = w / target_w;
    
    for r = 1:target_h
        for c = 1:target_w
            x = c * sw + 0.5 * (1 - sw);
            y = r * sh + 0.5 * (1 - sh);
            
            x1 = floor(x); x2 = ceil(x);
            y1 = floor(y); y2 = ceil(y);
            
            x1 = max(1, min(x1, w)); x2 = max(1, min(x2, w));
            y1 = max(1, min(y1, h)); y2 = max(1, min(y2, h));
            
            dx = x - x1;
            dy = y - y1;
            
            v1 = img_sub(y1, x1);
            v2 = img_sub(y1, x2);
            v3 = img_sub(y2, x1);
            v4 = img_sub(y2, x2);
            
            out(r, c) = (1-dx)*(1-dy)*v1 + dx*(1-dy)*v2 + (1-dx)*dy*v3 + dx*dy*v4;
        end
    end
end