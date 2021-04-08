function [S,T,Phi,f_r,U,P] = ss_freq_tv_evolve(f, Max_time, dt, Method)
% private function by Guy Gilboa (June 2013)
% Computes TV spectrum of an image, returns S, time interval T, Phi(T)
% and residual image f_r (to be added in the reconstruction).
% U, P - the flow and subgradient along time, if desired.
% Method is an optional struct specifying numerical method and params
% Example: [S,T,Phi] = ss_freq_tv_evolve(f, Max_time, dt)
% Based on: [1] G. Gilboa, “A spectral approach to total variation”, A. Kuijper et al. (Eds.): SSVM 2013, LNCS 7893, pp. 36--47. Springer, Heidelberg, 2013.
%           [2] G. Gilboa, “A total variation spectral framework for scale and texture analysis”, CCIT Report 833, Dept. of Electrical Engineering, Technion , June 2013.


if exist('Method','var')
    Num_method = Method.Num_method;
    if (strcmp(Num_method,'split_breg')) 
        tol  = Method.tol;
    else
        dt_proj = Method.dt_proj;        
        iter_proj = Method.iter_proj; 
    end    
else
    Num_method = 'proj'; % Projection algorithm, default        
    dt_proj=0.2; 
    iter_proj=500; 
end    
%Num_method = 'proj'; % projection algorithm
%Num_method = 'split_breg'; % split bregman

% Split Bregman params
if (strcmp(Num_method,'split_breg')) 
    % Split Bregman params:
    addpath SplitBregman_Rice    
    %tol = 0.000001;
end


% Check if user wants also u and p
if (nargout >4)
    up_flag = true;
else
    up_flag = false;
end

%iter_proj = 80; % debug
%iter_proj = 100; % debug

mu = 1/(2*dt);
NumIter = round(Max_time/dt);

S = zeros(1,NumIter); 
Phi = zeros(size(f,1),size(f,2),NumIter);
if up_flag
    U=Phi; P=Phi;
end
T = (1:NumIter)*dt;

u0 = f;
if (strcmp(Num_method,'proj')) 
   [u1,px,py]=proj_tvl2(u0,mu,iter_proj,dt);
   [u2,px,py]=proj_tvl2(u1,mu,iter_proj,dt,px,py);
   %u1 = tvl2_primal_dual(u0,2*dt,iter_proj);  %# debug
   %u2 = tvl2_primal_dual(u1,2*dt,iter_proj);
   %u1 = tvl2_primal_dual_sym(u0,2*dt,iter_proj);  %# debug
   %u2 = tvl2_primal_dual_sym(u1,2*dt,iter_proj);

else %split breg
    u1 = splitBregmanROF(u0,mu,tol);
    u2 = splitBregmanROF(u1,mu,tol);
end

for i=1:NumIter,
    ddu = (u0+u2-2*u1)/(dt*dt);  % one/two more iter
    t = i*dt;
    phi = ddu*t;
    Phi(:,:,i) = phi;
    if up_flag
       U(:,:,i) = u1; P(:,:,i) = (u0-u1)/dt;  
    end
    S(i) = sum(abs(phi(:)));
    if (i<NumIter) % not last iteration
        i
        
        NumIter
        u0=u1;
        u1=u2;
        if (strcmp(Num_method,'proj')) 
            [u2,px,py]=proj_tvl2(u2,mu,iter_proj,dt,px,py);
            %u2 = tvl2_primal_dual(u2,2*dt,iter_proj);  %# debug
            %u2 = tvl2_primal_dual_sym(u2,2*dt,iter_proj);  %# debug
        else %split breg
            u2 = splitBregmanROF(u2,mu,tol);
        end
    end

end % for i

f_r = (NumIter+1)*u1-NumIter*u2;  % residual image

end

