function AbSProbOnConstantImages()
% Use constant images to test the spectral properties of the (n,n)-AbSProb
% algorithm. For each constant image, Plot the reconstructed image, the
% RAPSD, and the DDF. 
%
%
close all; clear;

% change to current directory


n = 2; % Only illustrate (2,2) in the paper
G = [1/8, 1/4, 1/2, 3/4, 7/8];  % Also need to test <1/2 and >1/2
M = 512; % size of the constant image is M-by-N
N = M; 

w = 132;
stackedAll = ones(5*w, 5*w);

%======================== AbS ====================
for i = 1:length(G)
    img = uint8(ones(M,N).*G(i).*255);
    img = uint8(double(img)/2);
    % 
    [sharesYang, stackedYang, imgHt] = ItoFun(img); % we need imgHt
    [shares, stackedAbs] = ItoAbsFun(img);
    stackedAbs = stackedAbs>0; 
    [P,fr]=rapsd(double(imgHt>0), 128, [], []);
       
    figure(10); % Plot dither patterns
    subplot(5,1,i); imshow(stackedAbs(1:128,1:128),[]);
    stackedAll((i-1)*w+1:(i-1)*w+128, (3*w+1):(3*w+128)) = stackedAbs(1:128,1:128);  
    
    figure(11); % Plot rapsd
    subplot(5,1,i);
    plot(fr,P,'b-.','LineWidth',2); 
    %title('RAPSD','Fontsize',14); 
    axis([min(fr) max(fr) 0 4]);
    xlabel('Radial frequency \rightarrow','Fontsize',14);
    ylabel('RAPSD \rightarrow','Fontsize',14);
    hold on;
    
    [P,fr]=rapsd(double(stackedAbs), 128, [], []);
    plot(fr,P,'k','LineWidth',2); 
    %title('RAPSD','Fontsize',14); 
    xlabel('Radial frequency \rightarrow','Fontsize',14);
    ylabel('RAPSD \rightarrow','Fontsize',14);
    axis([min(fr) max(fr) 0 4]);
    if G(i)/2<=1/2
        hold on; plot(min(max(fr),sqrt(G(i)/2)),0,'dk','LineWidth',2);
    else
        hold on; plot(min(max(fr),sqrt((1-G(i))/2)),0,'dk','LineWidth',2);
    end
    %legend('Ordinary error diffusion','AbS Prob.'); 
    grid on;
    set(gca, 'Fontsize',14);
    

end

figure;
imshow(stackedAll,[]);


 end