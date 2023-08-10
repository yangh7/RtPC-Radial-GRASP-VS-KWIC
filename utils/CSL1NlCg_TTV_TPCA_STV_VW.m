function x = CSL1NlCg_TTV_TPCA_STV_VW(x0,param)
% 
% res = CSL1NlCg(param)
%
% Compressed sensing reconstruction of undersampled k-space MRI data
%
% L1-norm minimization using non linear conjugate gradient iterations
% 
% Given the acquisition model y = E*x, and the sparsifying transform W, 
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda1 * ||W*x||_1 + lambda2 * TV(x) 
%
% Based on the paper: Sparse MRI: The application of compressed sensing for rapid MR imaging. 
% Lustig M, Donoho D, Pauly JM. Magn Reson Med. 2007 Dec;58(6):1182-95.
%
% Ricardo Otazo 2008
%
 
% starting point
x=x0;

% line search parameters
maxlsiter = 150 ;
gradToll = 1e-3 ;
param.l1Smooth = 1e-15;	
alpha = 0.01;  
beta = 0.6;
t0 = 1 ; 
k = 0;

% lambda calculation
tmp_max = max( abs(x(:)) );
param.TPCAWeight = param.TPCAw * tmp_max;
param.TTVWeight = param.TTVw * tmp_max;
param.STVWeight = param.STVw * tmp_max;

% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% iterations
while(1)
    tic
    % backtracking line-search
	f0 = objective(x,dx,0,param);
	t = t0;
    f1 = objective(x,dx,t,param);
	lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
        if lsiter >= 50; sprintf( 'obj = %f,  L-S iteration = %d\n', f1, lsiter ); end
	end

	if lsiter == maxlsiter
		disp('Error - line search ...');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end

    % update x
	x = (x + t*dx);
	% print some numbers for debug purposes	
    if param.display,
        disp(sprintf('%3d:  obj= %f,  Max. signal= %f,  L-S iter= %d', k, f1, tmp_max, lsiter));
    end

    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	% stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end

    tmp_max = max( abs(x(:)) );
    if strcmp( param.VW, 'on' )
        % lambda update
        param.TPCAWeight = param.TPCAw * tmp_max;
        param.TTVWeight = param.TTVw * tmp_max;
        param.STVWeight = param.STVw * tmp_max;
    end
    toc
end
return;

function res = objective(x,dx,t,param) %**********************************

% L2-norm part
w=param.E*(x+t*dx)-param.y;
L2Obj=w(:)'*w(:);

% L1: temp PCA part
if param.TPCAWeight
   w = param.TPCA*(x+t*dx); 
   w=w.data;
   TPCAObj = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
else
   TPCAObj = 0;
end


% L1: temp TV part
if param.TTVWeight
   w = param.TTV*(x+t*dx); 
   TTVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TTVObj=0;
end


% L1: spatial TV part
if param.STVWeight
   w = param.STV*(x+t*dx); 
   STVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    STVObj=0;
end


% objective function
res = L2Obj + param.TPCAWeight*TPCAObj + param.TTVWeight*TTVObj + param.STVWeight*STVObj;

function g = grad(x,param)%***********************************************

% L2-norm part
L2Grad = 2.*(param.E'*(param.E*x-param.y));

% L1: temp PCA part
if param.TPCAWeight
    ww = param.TPCA*x;
    w=ww.data;
    www=(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    ww.data=www;
    TPCAGrad = param.TPCA'*ww;
else
    TPCAGrad = 0;
end

% L1: temp TV part
if param.TTVWeight
    w = param.TTV*x;
    TTVGrad = param.TTV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TTVGrad=0;
end

% L1: spatial TV part
if param.STVWeight
    w = param.STV*x;
    STVGrad = param.STV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    STVGrad=0;
end


% complete gradient
g = L2Grad + param.TPCAWeight*TPCAGrad + param.TTVWeight*TTVGrad + param.STVWeight*STVGrad;
