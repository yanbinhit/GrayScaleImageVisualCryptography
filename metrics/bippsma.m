%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = bippsma(N, g, M, bp, B_mask)
% BIPPSMA (function)
%       Binary-Pattern-Power-Spectrum-Matching-Algorithm, 
%	used by bnm, to create halftone patterns for a constant
%	gray-scale with blue-noise characteristics.
%
% INPUT ARGUMENT:
%       N      -> Size of desired output pattern.
%	g      -> Gray level.
%	M      -> Initial number of pixel pairs to swap.
%	bp     -> Initial binary pattern.
%	B_mask -> Constraint mask such that only pixels 
%		  corresponding to 1 can be swapped.
% OUTPUT ARGUMENT:
%       y -> The blue noise mask created
%
% EXAMPLE:
%       m = bnm(128);
%
% MICHAEL McCANN
% mccann@ee.udel.edu
%
% August 7, 1996
% Copyright 1996 Michael McCann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_bp=bippsma(N, g, M, bp, B_mask)

if (nargin<2)
        g=0.5;
end;
if (nargin<3)
        M=N;
end;
if (nargin<4)
        bp=init_pattern(N, g);
end;
if (nargin<5)
        B_mask=ones(N,N);
end;
fb=(g<0.5)*sqrt(g)+(g>=0.5)*sqrt(1-g);
c=interp1([0 1/32 1/16 1/8 1/4 3/8 1/2 5/8 3/4 7/8 15/16 31/32 1], ...
          [0.5 0.5  0.9  1.0 1.5 1.25 .5 1.25 1.5 1.0 0.9 0.5 0.5], g);
[F1,F2]=freqspace([N N], 'meshgrid');
%Hd=fftshift(exp(-(F1.^2+F2.^2)/(2*(fb/4.5)^2)));
Hd=fftshift(exp(-(F1.^2+F2.^2)/(c*fb^2)));
MSE_old = inf;                                  %Allows MSE_old to decrease
MSE_min = inf;                                  %Allows MSE_min to improve
new_bp=bp;
while(M >= 1)
  Y = real(ifft2(fft2(bp).*Hd));                %Filters bp
  err_mat = (Y-g).*double(B_mask);                      %Forms error matrix

  ind_ones = find(bp==1);                       %Finds ones
  ind_zero = find(bp==0);                       %Finds zeros
  
  err_ones = err_mat(ind_ones);                 %Finds ones errors
  err_zero = err_mat(ind_zero);                 %Finds zeros errors
 
  [S_EO,I_EO] = sort(err_ones);                 %Sorts ones errors
  [S_EZ,I_EZ] = sort(err_zero);                 %Sorts zeros errors

  for m = 1:min([M length(I_EO) length(I_EZ)])  %Switches ones and zeros
      if (err_mat(ind_zero(I_EZ(m))) < 0 & ...  %Checks for pos and neg errors
          err_mat(ind_ones(I_EO(length(I_EO)-m+1))) > 0)
         bp(ind_zero(I_EZ(m))) = 1;             %Switches zeros to ones
         bp(ind_ones(I_EO(length(I_EO)-m+1))) = 0;
                                                %Switches ones to zeros
      end
  end

  MSE = mean2(err_mat.^2);                      %Determines mean squared error
  if (MSE >= MSE_old)                           %MSE increase
     M = M/2;
  end
  if (MSE < MSE_min)                            %MSE decrease
     MSE_min = MSE;
     new_bp = bp;
  end
  MSE_old = MSE;
% imshow(bp,2); drawnow
end
new_bp=uint8(new_bp);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IP=init_pattern(N, g)

P=rand(N*N,1);
N_pixels=ceil(N*N*g);
[Ps, Pi]=sort(P);
IP=zeros(N,N);
IP(Pi(1:N_pixels))=ones(1,N_pixels);
