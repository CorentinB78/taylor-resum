% Simple test Taylor series resummation via conformal map, for Parcollet/CCQ
% Barnett 7/2/19. Use SCtoolbox 7/29/19.
clear; setupsc; verb=1;    % verbosity: 0 just text, 1 figs, 2 more tests...

f = @(z) log(1+z); fsing = -1;   % func to eval in z-plane, & its sing loc
p=15;                            % # available coeffs
an = [0 -(-1).^(1:p)./(1:p)]';   % data: col vec of Taylor coeffs of f about 0
if verb>1, x = 0.3; fe = f(x); fa = sum(an.'.*x.^(0:p));   % check trunc Taylor
  fprintf('rel err in eval at sensible x=%.3g:  %.3g\n',x,(fe-fa)./fe), end

ztarg = 3+0i;                    % desired z target
ftrue = f(ztarg)                 % true ans, for validation

map = 'p';                       % map type: 'a' RHP parabola, 'p' SC polygon
                                 % (both work for log(1+z) )

if map=='a'  % map right-opening parabola interior to unit disc...
  b = 0.9;   % minus x-intersection of parabola (also scales it, fixing w(0)=0)
  w = @(z) -tan(pi*sqrt(-z/b)/4).^2;        % fwd map z -> w
  iw = @(w) -b*((4/pi)*atan(sqrt(-w))).^2;  % inv map w -> z
elseif map=='p'  % arb polygon via SC toolbox
  %xmax = 10; po = polygon([-1-1i, xmax-1i, xmax+1i, -1+1i]); % rect
  %x0 = 10; po = polygon([-1-1i, x0-x0*1i, x0+x0*1i, -1+1i]); % trapezoid
  a=1.25; po = polygon([-1+1i, -1-1i, inf], [a,a, -2*a+1]);  % a=ang, like RHP
  g = center(diskmap(po),0);    % set conformal center (preimage of origin)
  if verb>1, figure; plot(g); title('SC poly and map, in z plane'); end
  iw = @(w) eval(g,w);     % z(w), our inv map, inline func
  w = @(z) evalinv(g,z);   % w(z), our fwd map. Note SC evalinv is slow
end
  
if verb, z = 0.4+0.6i; fprintf('check inv map good: %.3g\n',abs(iw(w(z))-z)),end
wtarg = w(ztarg);                    % target in w-plane

% get L...
r = 0.3; [L Nc] = matrixfrominvmap(iw,p,r);  % r<1 but not too small - how pick?

if verb, figure(2); imagesc(0:p,0:p,log10(abs(L)));   % check L
  axis equal tight; colorbar; title('coeff matrix: log_{10} L_{mn}');
  v=caxis; v(1)=-2; caxis(v); xlabel('n'); ylabel('m');
  r2 = 1.3*r; L2 = matrixfrominvmap(iw,p,r2);   % check stability wrt r
  fprintf('est L rel acc: %.3g\n',norm(L-L2)/norm(L)) % p=10: L acc to 1e-13
end

if verb    % plot both planes, first z then w...
  figure(1); clf; subplot(1,2,1);
  plot(real(fsing)+[0 -10],imag(fsing)+[0 0],'r*-'); hold on;  % branch cut
  plot(0,0,'ko'); axis equal; axis(3*[-1 1 -1 1]);
  if map=='a'
    y=-10:0.1:10; x = -1+y.^2/4; z = b*(x+1i*y);   %  unit circ preimage in z
    plot(z,'.-'); hold on; plot([-b 0],[0,2*b],'ko');  % check parabola params
  elseif map=='p'
    warning('off','MATLAB:hg:EraseModeIgnored');   % since toolbox is too old
    plot(g); hold on    % overloaded SC map and polygon plot method
  end
  plot(real(ztarg),imag(ztarg),'k.'); title('z plane');
  text(real(ztarg),imag(ztarg),'z_{targ}');
  plot(iw(r*exp(2i*pi*(1:Nc)/Nc)),'m.');   % colloc pts used to get L
  %
  subplot(1,2,2); plot(exp(1i*2*pi*(1:100)/100),'r-');  % check is unit circ
  hold on; plot(real(w(fsing)),imag(w(fsing)),'r*');
  if map=='a', plot(w(z),'.-'); plot([w(0) w(-b) w(2i*b)],'ko'); end   % some pts and parabola
  plot(real(wtarg),imag(wtarg),'k.');
  text(real(wtarg),imag(wtarg),'w(z_{targ})');
  plot(r*exp(2i*pi*(1:Nc)/Nc),'m.'); axis tight equal; title('w plane');
end

cn = L*an;                        % get w-plane Taylor coeffs
ftarg = sum(cn.'.*wtarg.^(0:p))   % do the eval
e = abs((ftarg-ftrue)/ftrue);     % validate error
fprintf('ztarg=(%.3g,%.3g): f(ztarg) rel err = %.3g\n',real(ztarg),imag(ztarg),e)
