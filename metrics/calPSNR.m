function psnr = calPSNR(imageOrig, imageAtt)
[Mo, No, So] = size(imageOrig);
[Ma, Na, Sa] = size(imageAtt);
if (Mo ~= Ma) | (No ~= Na) | (So ~= Sa)
    error('The size of the two input images must be the same!');
end

imageOrig = double(imageOrig);
imageAtt = double(imageAtt);

peakPower = max(imageOrig(:))^2;
noisePower = (imageOrig - imageAtt).^2;
noisePowerAvg = sum(noisePower(:))/(Ma*Na*Sa);
psnr = 10*log10((peakPower)/(noisePowerAvg+eps));


    

