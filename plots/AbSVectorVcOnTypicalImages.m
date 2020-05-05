function AbSVectorVcOnTypicalImages()

close all; clear all; 
% Lena, or pillars
imageIn = imread('..\images\lena.tiff','tiff');
if size(imageIn,3)>1
    imageIn = rgb2gray(imageIn);
end
imageIn = imresize(imageIn, [512 512]);
w = 512;
%imageIn = 255.*(imageIn>127); % uncomment this line to test thresholded secret image

% Abs
tic
[shares, stacked, imageHtNoVc] = BlkVcAbs(imageIn);
toc
figure;
imshow(stacked,[]); title('AbS');
print -depsc -r200 'fig_AbsVectorVcLena.eps'
imwrite(stacked,'fig_AbsVectorVcLena.tif','tif');

figure(200); 
imshow(imageIn,[]);
figure(201);
imshow(shares(:,:,1),[]);
figure(202);
imshow(shares(:,:,2),[]);
figure(203);
imshow(stacked,[]);

refImage = double(imageIn)/2;
sigma = 2;
hSize = 11;
h = fspecial('gaussian', hSize, sigma);
stackedFilt = imfilter(double((stacked>0)*255), h, 'replicate');
psnrAbs = calPSNR(refImage, stackedFilt); % HPSNR
[mssimAbs, ssim_map] = ssim(refImage, stackedFilt);
psnrAbs
mssimAbs

% AbS with Equalization
imageInEq = histeq(imageIn);
[shares, stacked, imageHtNoVc] = BlkVcAbs(imageInEq);
figure;
imshow(stacked,[]); title('AbS equalized');
print -depsc -r200 'fig_AbsVectorVcEqLena.eps'
imwrite(stacked,'fig_AbsVectorVcEqLena.tif','tif');
% show histogram
sigma = 2;
hSize = 11;
h = fspecial('gaussian', hSize, sigma);
stackedFilt = imfilter(double((stacked>0)*255), h, 'replicate');
figure; imhist(uint8(stackedFilt)); title('AbS equalized');