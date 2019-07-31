function [L N] = matrixfrominvmap(invmap,p,r)
% MATRIXFROMINVMAP.  Use circle sampling method to get lower-tri coeff fwd map
%
% L = matrixfrominvmap(invmap,p,r) returns L, a (p+1)-by-(p+1) low-triangular
%  matrix mapping the vector of p+1 Taylor coeffs of f(z) about 0 in the z-plane
%  to the p+1 coeffs of f(invmap(w)) about 0 in the w-plane.
%  Collocation on the circle of radius r in the w-plane is used.
%
% [L N] = matrixfrominvmap(invmap,p,r) also returns # colloc pts used.
%
% Called without arguments, performs a self-test.

% Barnett 7/3/19
if nargin==0, test_matrixfrominvmap; return; end

N = 6*p+1;                      % # colloc pts, empirical >3p for emach
wj = r*exp(2i*pi*(0:N-1)/N);    % r-circle
zj = invmap(wj);                % preimages
%figure(2); plot([zj zj(1)],'.-'); axis equal; title('z');
V = ones(N,p+1); for i=1:p, V(:,i+1) = V(:,i).*zj(:); end   % Vandermonde
F = (1/N)*exp(-2i*pi*(0:p)'*(0:N-1)/N);   % rows 0..p of iDFT matrix
L = diag(r.^-(0:p)) * (F*V);    % rescale coeffs
%%%%%%%%%%%%%

% test circle preimage sampling method for getting L matrix...
function test_matrixfrominvmap
b1=-0.7; b2=0.3; b3=-0.1; iw = @(w) b1*w+b2*w.^2+b3*w.^3;  % toy inv map p=3
L4exact = [1 0 0 0;0 b1 0 0; 0 b2 b1^2 0; 0 b3 2*b1*b2 b1^3];
p = 5;    % max degree of poly to fit in L (>=3)
rs = [4 2 1 1e-1 1e-2 1e-3];    % w-radius
for i=1:numel(rs), r = rs(i);
  L = matrixfrominvmap(iw,p,r);
  err = norm(L(1:4,1:4) - L4exact);
  fprintf('r=%.3g:     \terr=%.3g\n',r,err)
end
