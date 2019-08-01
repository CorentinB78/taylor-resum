% Test Taylor series resummation via general Schwarz-Christoffel conformal map.
% Bertrand et al arxiv:1903.11646v3, App. A, toy func going to back side of pole
% Barnett 7/31/19
clear; setupsc; verb=1;    % verbosity: 0 just text, 1 figs, 2 more tests...

f = @(z) 1./log(1i*(1-z)+1);     % Bertrand toy func
p=15;                            % # available coeffs - fails if too large
r = 0.3; N=100; zj = r*exp(2i*pi*(0:N-1)'/N);  % r-circle pts to extract Taylor
F = (1/N)*exp(-2i*pi*(0:p)'*(0:N-1)/N);   % rows 0..p of iDFT matrix
an = diag(r.^-(0:p)) * (F*f(zj));         % rescale to get f's Taylor coeffs
if verb>1, x = 0.3; fe = f(x); fa = sum(an.'.*x.^(0:p));   % check trunc Taylor
  fprintf('rel err in Taylor at sensible x=%.3g:  %.3g\n',x,(fe-fa)./fe), end

ztarg = 2.0+0i;                  % a desired z target (see App. A)
ftrue = f(ztarg)                 % true ans, for validation

a = 0.1; po = polygon([1+(-1-1i)*a, 1, 1+(1-1i)*a, inf]+1i*a, [1.25,1.5,1.25,-2]);  % exterior of hairpin
%po=polygon([inf 1 1-1i],[-2 2 1]);   % exterior domain of vertical slit, fails
%po=polygon([inf 1 2-1i],[-2 1.75 1.25]);   % 2-point hairpin, worse
g = center(diskmap(po),0);    % set conformal center (preimage of origin)
if verb>1, figure; plot(g); title('SC poly and map, in z plane'); end
iw = @(w) eval(g,w);     % z(w), our inv map, inline func
w = @(z) evalinv(g,z);   % w(z), our fwd map. Note SC evalinv is slow
  
if verb, z = 0.4+0.6i; fprintf('check inv map good: %.3g\n',abs(iw(w(z))-z)),end
wtarg = w(ztarg);                    % target in w-plane

% get L...
r = 0.3; [L Nc] = matrixfrominvmap(iw,p,r);   % r<1 but not too small

if verb, figure(2); imagesc(0:p,0:p,log10(abs(L)));   % check L
  axis equal tight; colorbar; title('coeff matrix: log_{10} L_{mn}');
  v=caxis; v(1)=-2; caxis(v); xlabel('n'); ylabel('m');
  r2 = 1.2*r; L2 = matrixfrominvmap(iw,p,r2);   % check stability wrt r
  fprintf('est L rel acc: %.3g\n',norm(L-L2)/norm(L)) % p=10: L acc to 1e-13
end

if verb    % plot both planes, first z then w...
  figure(1); clf; subplot(1,2,1);
  warning('off','MATLAB:hg:EraseModeIgnored');   % since toolbox is too old
  plot(1,0,'r*'); hold on; plot([1 1],[-1 -10],'r*-');  % pole & branch cut
  plot(0,0,'ko'); axis equal; axis(3*[-1 1 -1 1]);
  plot(g);    % overloaded SC map and polygon plot method
  plot(real(ztarg),imag(ztarg),'k.'); title('z plane');
  text(real(ztarg),imag(ztarg),'z_{targ}');
  plot(iw(r*exp(2i*pi*(1:Nc)/Nc)),'m.');   % colloc pts used to get L
  %
  subplot(1,2,2); plot(exp(1i*2*pi*(1:100)/100),'r-');  % check is unit circ
  hold on; plot(real(wtarg),imag(wtarg),'k.');
  text(real(wtarg),imag(wtarg),'w(z_{targ})');
  plot(r*exp(2i*pi*(1:Nc)/Nc),'m.'); axis tight equal; title('w plane');
end

cn = L*an;                        % get w-plane Taylor coeffs
%cn = cn .* erfc([-p/2:p/2]'/p*5.0)/2;  % roll off the series? - fails
disp('w-plane Taylor coeff mags c_n (should decay!):'); disp(abs(cn))
ftarg = sum(cn.'.*wtarg.^(0:p))   % do the eval
e = abs((ftarg-ftrue)/ftrue);     % validate error
fprintf('ztarg=(%.3g,%.3g): f(ztarg) rel err = %.3g\n',real(ztarg),imag(ztarg),e)
