function [sharesOut, imageHTStack, imageHtNoVc] = BlkVcAbs(imageIn)
% Implement the block encoding based VC, with AbS loops. Now implement only
% the (2,2)-threshold scheme. 
% For the (3,3)-threshold scheme, we must use a larger block size, such as
% 2x3 or 3x3. For the 2x2 blocksize, 



if size(imageIn, 3)>1
    imageIn = rgb2gray(imageIn);
end

M = 2;  % the size of the block is M-by-N
N = 2;
M2 = [1 1 0 0; 1 1 0 0]; % 0: black, subscript indicates # whiteness
M1 = [1 1 0 0; 0 1 1 0];
M0 = [1 1 0 0; 0 0 1 1];

%imageEq = HistEqHalf(imageIn); % equalize to range [0, 127]
imageEq = uint8(double(imageIn)/2); % gamut mapping
imageHtNoVc = HalftoningED(imageEq);

[nR, nC] = size(imageIn); 
imageHTStack = zeros(nR+2*M, nC+2*N); % stacked result
imageHtED = zeros(nR+2*M, nC+2*N);    % used in the diffusion process
shares = zeros(nR+2*M, nC+2*N, 2);    % generated shares
sharesOut = zeros(nR, nC, 2);         % clipped shares, having same size as input image
imageEq = padarray(double(imageEq), [2 2], 'replicate', 'both'); 
                                      % sec.image after gamut mapping (no equalization)
imageEqBackup = imageEq;
imageClass = zeros(nR+2*M, nC+2*N);   % record the pixel block class, for testing only

for i = (M+1):M:(M+nR)
    for j = (N+1):N:(M+nC)
        
        % Try ED halftoing current block
        blkIdxR = i:(i+M-1);
        blkIdxC = j:(j+N-1);
        currentBlk = imageEq(blkIdxR, blkIdxC); % get current block
        blkHt = zeros(M,N);
        blkHt(1,1) = double((currentBlk(1,1)>127)*255); % (1,1)
        err11 = currentBlk(1,1)-blkHt(1,1);
        currentBlk(2,1) = currentBlk(2,1)+ err11*5/16;
        currentBlk(1,2) = currentBlk(1,2)+ err11*7/16;
        currentBlk(2,2) = currentBlk(2,2) + err11/16;
        
        blkHt(1,2) = double((currentBlk(1,2)>127)*255); % (1,2)
        err12 = currentBlk(1,2)-blkHt(1,2);
        currentBlk(2,1) = currentBlk(2,1)+ err12*3/16;
        currentBlk(2,2) = currentBlk(2,2) + err12*5/16;
        
        blkHt(2,1) = double((currentBlk(2,1)>127)*255); % (2,1)
        err21 = currentBlk(2,1)-blkHt(2,1);
        currentBlk(2,2) = currentBlk(2,2) + err21*7/16;
        
        blkHt(2,2) = double((currentBlk(2,2)>127))*255; % (2,2)
        
        
        imageClass(blkIdxR, blkIdxC) = ones(2,2).* sum(blkHt(:)/255); 
                                        % record image class for display
        
        % See the intensity of trail halftone block
        
        p = randperm(4);
        switch sum(blkHt(:)/255)
            case 4 % four white pixels
                Mp = M2(:,p);
                Ms = ~(~Mp(1,:) | ~Mp(2,:)); % simulate stacking
                imageHTStack(blkIdxR, blkIdxC) = reshape(Ms,2,2);
                %imageHtED(blkIdxR, blkIdxC) = ones(2,2);
                imageHtED(blkIdxR, blkIdxC) = imageHTStack(blkIdxR, blkIdxC);
            case 3
                Mp = M2(:,p);
                Ms = ~(~Mp(1,:) | ~Mp(2,:)); % simulate stacking
                imageHTStack(blkIdxR, blkIdxC) = reshape(Ms,2,2);
                %imageHtED(blkIdxR, blkIdxC) = blkHt./255;
                imageHtED(blkIdxR, blkIdxC) = imageHTStack(blkIdxR, blkIdxC);
            case 2
                Mp = M1(:,p);
                Ms = ~(~Mp(1,:) | ~Mp(2,:)); % simulate stacking
                imageHTStack(blkIdxR, blkIdxC) = reshape(Ms,2,2);
                imageHtED(blkIdxR, blkIdxC) = imageHTStack(blkIdxR, blkIdxC);
            case 1
                Mp = M0(:,p);
                Ms = ~(~Mp(1,:) | ~Mp(2,:)); % simulate stacking
                imageHTStack(blkIdxR, blkIdxC) = reshape(Ms,2,2);
                imageHtED(blkIdxR, blkIdxC) = imageHTStack(blkIdxR, blkIdxC);
            case 0
                Mp = M0(:,p);
                Ms = ~(~Mp(1,:) | ~Mp(2,:)); % simulate stacking
                imageHTStack(blkIdxR, blkIdxC) = reshape(Ms,2,2);
                imageHtED(blkIdxR, blkIdxC) = imageHTStack(blkIdxR, blkIdxC);
            otherwise
                error('wrong block type in trail halftone block');
        end
        
        shares(blkIdxR, blkIdxC, 1) = reshape(Mp(1,:),2,2);
        shares(blkIdxR, blkIdxC, 2) = reshape(Mp(2,:),2,2);
        
        imageHTStack(blkIdxR, blkIdxC) = double(imageHTStack(blkIdxR, blkIdxC) .* 255);
        imageHtED(blkIdxR, blkIdxC) = double(imageHtED(blkIdxR, blkIdxC) .* 255);
        
        % Now distribute the quantization noise to neighbors (2 steps)
        % Basic idea: If you can't compensate for it, pass it to your neighbor
        errUl = imageEq(i,j)-imageHtED(i,j);
        
        imageEq(i+2,j-2) = imageEq(i+2,j-2) + errUl * (3/16) *(3/16);
        imageEq(i+2,j-1) = imageEq(i+2,j-1) + errUl * (3/16) *(5/16);
        imageEq(i+2,j) = imageEq(i+2,j) + errUl * (3/16) *(1/16);
        imageEq(i+1,j) = imageEq(i+1,j) + errUl * (3/16) * (7/16);
        imageEq(i,j+1) = imageEq(i,j+1) + errUl * (7/16);
        imageEq(i+1,j) = imageEq(i+1,j) + errUl * (5/16);
        imageEq(i+1,j+1) = imageEq(i+1,j+1) + errUl * (1/16);
        
        errUr = imageEq(i,j+1)-imageHtED(i,j+1);
        imageEq(i,j+2) = imageEq(i,j+2) + errUr * (7/16);
        imageEq(i+1,j) = imageEq(i+1,j) + errUr * (3/16);
        imageEq(i+1,j+1) = imageEq(i+1,j+1) + errUr * (5/16);
        imageEq(i+1,j+2) = imageEq(i+1,j+2) + errUr * (1/16);
        
        
        errll = imageEq(i+1,j)- imageHtED(i+1,j);
        imageEq(i+1,j+1) = imageEq(i+1,j+1) + errll * (7/16);
        imageEq(i+2,j-1) = imageEq(i+2,j-1) + errll * (3/16);
        imageEq(i+2,j) = imageEq(i+2,j) + errll * (5/16);
        imageEq(i+2,j+1) = imageEq(i+2,j+1) + errll * (1/16);
        
        errlr = imageEq(i+1,j+1) - imageHtED(i+1,j+1);
        imageEq(i+1,j+2) = imageEq(i+1,j+2) + errlr * (7/16);
        imageEq(i+2,j) = imageEq(i+2,j) + errlr * (3/16);
        imageEq(i+2,j+1) = imageEq(i+2,j+1) + errlr * (5/16);
        imageEq(i+2,j+2) = imageEq(i+2,j+2) + errlr * (1/16);
    end
end

imageHTStack = imageHTStack(M+1:M+nR, N+1:N+nC);
sharesOut(:,:,1) = shares(M+1:M+nR, N+1:N+nC, 1); % shares to return
sharesOut(:,:,2) = shares(M+1:M+nR, N+1:N+nC, 2);



