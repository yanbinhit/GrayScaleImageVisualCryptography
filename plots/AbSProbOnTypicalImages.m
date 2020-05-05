function AbSProbOnTypicalImages()
% Use typical images to compare the visual quality of Yang, Wang, AbS and
% Halftone result.
close all; clear all; 

imageIn = imread('..\images\lena.tiff','tiff');

if size(imageIn,3)>1
    imageIn = rgb2gray(imageIn);
end
imageIn = imresize(imageIn, [512 512]);
w = 512;

[shares, stacked] = ItoAbsFun(uint8(double(imageIn)/2));

%subplot(2,2,3);

figure(200); 
imshow(imageIn,[]);
figure(201);
imshow(shares(:,:,1),[]);
figure(202);
imshow(shares(:,:,2),[]);
figure(203);
imshow(stacked,[]);

figure;
imshow(stacked,[]); 
%title('AbS');
%print('-deps','-r200','fig_probAbsLena.eps')
imwrite(stacked,'fig_probAbsLena.tif','tif');
stackedAbs = stacked;
imageInAbs = imageIn>127;
save dataAbs.mat stackedAbs imageInAbs;

refImage = double(imageIn)/2;
sigma = 2;
hSize = 11;
h = fspecial('gaussian', hSize, sigma);
stackedFilt = imfilter(double((stacked>0)*255), h, 'replicate');
psnrAbs = calPSNR(refImage, stackedFilt); % HPSNR
[mssimAbs, ssim_map] = ssim(refImage, stackedFilt);
psnrAbs
mssimAbs
