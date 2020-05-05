function g = HistEqHalf(f)
% Equalize the histogram of input image to range [0 127]
%<Inputs>
%   f: input gray scale image with 8 bit quantization
%<Outputs>
%   g: output gray scale image, equalized to the range [0, 127]
%
%<Author>
%   Bin Yan, Modified.


if size(f,3)>1
    f = rgb2gray(f);
end
f = double(f); %load image

[mr,mc] = size(f);
m = mr*mc; %total number of pixels
g = zeros(mr,mc); % output image

%start with integer intensities range from 0 to 255
f = floor(f);
L = 256;
gmax = 127; %number of intensity levels

%define histogram
p = zeros(L,1);
for i = 1:mr
   for j = 1:mc
      p(f(i,j)+1) = p(f(i,j)+1) + 1/m;
   end
end

%histogram equalization
for i = 1:mr
   for j = 1:mc
      g(i,j) = gmax*sum(p(1:f(i,j)+1));
   end
end
g = floor(g); %round down to nearest integer
g = uint8(g);
