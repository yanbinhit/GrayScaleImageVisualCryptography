function [shares, stacked] = ProbAbs33Fun(imageIn)
% AbS based Prob. Scheme for (3,3)-threshold

if size(imageIn, 3)>1
    imageIn = rgb2gray(imageIn);
end

imageIn = round(imageIn/4);  %Gamut mapping

M = 1;  % the size of the block is M-by-N
N = 1;
M0 =[0 0 1 1; 0 1 0 1; 0 1 1 0]; % Basis matrix for (3,3) deterministic scheme
M1 =[1 1 0 0; 1 0 1 0; 1 0 0 1];

[nR, nC] = size(imageIn); 
shares = zeros(nR+2*M, nC+2*N, 3);
imageHT = zeros(nR+2*M, nC+2*N);
imageEq = imageIn;
imageEq = padarray(double(imageEq), [2 2], 'replicate', 'both');


for i = (M+1):M:(nR)
    for j = (N+1):N:(nC)
        
        % Try ED halftoing current pixel (don't distribute error)
        T = 127;
        blkHt(i,j) = double((imageEq(i,j)>=T)*255);      
        
        % See the intensity of trail halftone block
        p = unidrnd(4);
        if blkHt(i,j) == 255
            Mp = M0(:,p);
            shares(i,j,1) = Mp(1,1);
            shares(i,j,2) = Mp(2,1);
            shares(i,j,3) = Mp(3,1);
            Ms = ~(Mp(1,1) | Mp(2,1) | Mp(3,1) );
            imageHT(i,j) = Ms;               % imageHT is the stacked
        else
            Mp = M1(:,p);
            shares(i,j,1) = Mp(1,1);
            shares(i,j,2) = Mp(2,1);
            shares(i,j,3) =  Mp(3,1);
            Ms = ~(Mp(1,1) | Mp(2,1) | Mp(3,1) );
            imageHT(i,j) = Ms;
        end  
        imageHT(i,j) = double(imageHT(i,j) .* 255);
        err = imageEq(i,j) - imageHT(i,j);
        % now diffuse error according to reconstructed pixel
        imageEq(i,j+1) = imageEq(i,j+1) + err * (7/16);
        imageEq(i+1,j-1) = imageEq(i+1,j-1) + err * (3/16);
        imageEq(i+1,j) = imageEq(i+1,j) + err * (5/16);
        imageEq(i+1,j+1) = imageEq(i+1,j+1) + err * (1/16);      
             
    end
end

stacked = uint8(imageHT(M+1:M+nR, N+1:N+nC));
% figure; imshow(uint8(imageEqBackup));
% figure; imshow(imageHT,[]);
% figure; 
% subplot(1,2,1); imshow(shares(:,:,1),[]); 
% subplot(1,2,2); imshow(shares(:,:,2),[]);
