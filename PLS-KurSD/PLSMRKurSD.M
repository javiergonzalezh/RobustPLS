function [b,W]=PLSMRKurSD(Xp,yp,K)

%
% [b,W]=PLSMRKurSD(Xp,yp,K)
%
% This functions runs the method PLSMRKurSD for robust partial least squares regression in Gonzalez et. al., 2009. 
%
% Inputs:  Xp, Matrix of dependent variables
%          yp, Vector of independent variables
%          K, Number of latent components
%          mode,   0  use p maxizing and p minimizing directions
%                     iterate until no outliers are detected (default)
%                 >0  use p maxizing and p minimizing directions
%                     iterate "mode" times
%                 <0  use only maximizing directions
%                     iterate "mode" times
% Outputs: b, Vector of estimated parameters
%          W, weigths matrices 
%

% Javier Gonzalez/Daniel Peña/Rosario Romera    16/04/05

y=cen(yp);
X=cen(Xp);

yx = [y,X];
[n,p] = size(yx);

V = kur_main(yx);		
s = V(2:p,1);
S = V(2:p,2:p);
W = s;

if K >= 2,
    for i = 2:K,
        w = s-S*W*inv(W' * S * W)*W'*s;   
	W = [W,w];                           
    end
end

b = W*inv(W'*S*W)*W'*s;                
