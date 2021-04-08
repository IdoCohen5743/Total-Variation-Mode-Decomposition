function [u,px,py]=proj_tvl2(f,lam,iter,dt,p0x,p0y)
%% Chambolle's projection method for solving ROF-TV
%% input: image, lam, #iterations, dt, p0
%% p0x,y - optional, if there is a good approximation of p0 (vector)
%% output: processed image
%% example: J=proj_tvl2(I,lam,iter,dt,p0)

[ny,nx]=size(f); 

if ~exist('iter')
	iter=30;
end
if ~exist('dt')
	dt=0.125;
end
if ~exist('p0x')
	p0x=zeros(ny,nx); p0y=p0x;
end

px=p0x; py=p0y;
lami=1/lam;  % here the algorithm works with inverse of lambda
% for i=1:iter,  %% do iterations
while(true)
   %% compute projection % eq. 9 in Chambolle (dual problem for p)
   [Gx,Gy] = grad(div(px,py) - f/lami);
   aG = sqrt(Gx.^2+Gy.^2); %abs(G)
   pxnew = (px + dt*Gx)./(1+dt*aG);
   pynew = (py + dt*Gy)./(1+dt*aG);
   diffx = max(abs(pxnew(:)-px(:)));
   diffy = max(abs(pynew(:)-py(:)));
   if((diffx<1e-8)&(diffy<1e-8))
       break;
   else
       px = pxnew;
       py = pynew;
   end
end % for i
%% return image
v = lami*div(px,py);
u=f-v;
u = u+mean(f(:))-mean(u(:));  % keep mean unchanged

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions              %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient (forward difference)
function [fx,fy] = grad(P)
	fx = P(:,[2:end end])-P;
	fy = P([2:end end],:)-P;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Divergence (backward difference)

% corrected version 24.5.17
function fd = div(Px,Py)
	fx = Px-Px(:,[1 1:end-1]);   
    fx(:,1) = Px(:,1); fx(:,end) = -Px(:,end-1);
	fy = Py-Py([1 1:end-1],:);
    fy(1,:) = Py(1,:); fy(end,:) = -Py(end-1,:);
    fd = fx+fy;
      
% % version with error
% function fd = div(Px,Py)
% 	fx = Px-Px(:,[1 1:end-1]);
% 	fy = Py-Py([1 1:end-1],:);
%     fd = fx+fy;



