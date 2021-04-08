clear; close all; clc;

loadData = 1;   % 1 = Load saved TV data ; 0 = compute from scratch

%% Compute TV flow on circles image
tic;
addpath('../')
if ~loadData
    zebra = im2double(rgb2gray(imread('zebra_media_gmu.jpg')));
    Max_time = 3; dt = 1e-3; % 3 peaks
    
    Method.Num_method = 'proj'; % currently single method
    Method.dt_proj = 0.2;
    Method.iter_proj = 10000;  % very accurate
    
    
    
    N=100;
    f1d = zeros(1,N);
    
    
    % 3 peaks
    f1d(20:21)=1;
    f1d(38:41)=1;
    f1d(70:80)=1;
    f1d = f1d + 0.01*randn(size(f1d));
    h202 = figure(202);
    plot(f1d, 'k', 'LineWidth', 6); grid on;
    set(gca, 'fontsize', 360, 'TickLabelInterpreter', 'latex');
    ystart = min(f1d);
    ystop = max(f1d);
    ygap = ystop - ystart;
    set(gca, 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);
    write_pdf_New_Image('peaks_src',h202,6000,4000)
    
    
    
    f1d = f1d - mean(f1d);
    N = length(f1d);
    
    
    f1d = f1d - mean(f1d);  % subtract mean
    figure(1); movegui('northwest'); plot(f1d,'k','LineWidth',2); axis([1 N -0.5 1.5]);
    title('f')
    
    f = [f1d; f1d; f1d];
    
    scale = 1;
    f = f*scale;
    
    [S,T,Phi,f_r,utmp,ptmp] = ss_freq_tv_evolve(f, Max_time, dt, Method);  % evolve image
    
    
else
    load('TV_data1D.mat', 'S', 'T', 'Phi', 'f_r', 'utmp','ptmp','dt','f1d');
    
    
    
end
t_chambolle = toc;

% Draw initial condition

h203 = figure(203);
plot(f1d, 'k','LineWidth', 6); grid on;
set(gca, 'fontsize', 360, 'TickLabelInterpreter', 'latex');
ystart = min(f1d);
ystop = max(f1d);
ygap = ystop - ystart;
set(gca, 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);
write_pdf_New_Image('peaks_src',h203,6000,4000)



zebra_flag = 0;
figure(); plot(T,S)
figure();plot(squeeze(ptmp(2,:,200))); hold on;
plot(squeeze(ptmp(2,:,500)));
plot(squeeze(ptmp(2,:,1000)));
plot(squeeze(ptmp(2,:,1200)));
plot(squeeze(ptmp(2,:,1600)));
plot(squeeze(ptmp(2,:,2500)));hold off;
utmpMain = utmp;
utx = reshape(utmp(2,:,:),size(utmp,2),size(utmp,3));
h70 = figure(70); mesh(utx); % title('u(x,t)')
grid on;
h70.Children.CameraPosition = 1.0e+04 *[ 1.4607   -0.0627    0.0005 ];
set(gca,'FontSize',240, 'TickLabelInterpreter', 'latex');

timePoints = [0.65,0.95,1.2,1.65,2.5];
band1 = find(T>=0 & T<timePoints(1));
band2 = find(T>=timePoints(1) & T<timePoints(2));
band3 = find(T>=timePoints(2) & T<timePoints(3));
band4 = find(T>=timePoints(3)& T<timePoints(4));
band5 = find(T>=timePoints(4)& T<timePoints(5));
band6 = find(T>=timePoints(5));

imgBand1 = dt*sum(Phi(:,:,band1),3);
imgBand2 = dt*sum(Phi(:,:,band2),3);
imgBand3 = dt*sum(Phi(:,:,band3),3);
imgBand4 = dt*sum(Phi(:,:,band4),3);
imgBand5 = dt*sum(Phi(:,:,band5),3);
imgBand6 = dt*sum(Phi(:,:,band6),3) + f_r(2,:);
figure(700);plot(imgBand1(2,:)+imgBand2(2,:)+imgBand3(2,:)+imgBand4(2,:)+imgBand5(2,:)+imgBand6(2,:))
figure; plot(imgBand1(2,:))
figure; plot(imgBand2(2,:))
figure; plot(imgBand3(2,:))
figure; plot(imgBand4(2,:))
figure; plot(imgBand5(2,:))
figure; plot(imgBand6(2,:))


%%
tic;
pf=-squeeze(ptmp(2,:,1)); pf=pf(:);
f = f1d(:);

ui = f;pui = pf;
delta = 1e-3;
epsi = 1e-16;
[uix,uiy] = gradIdo(ui); absGrad = abs(uix+uiy*1i);
maskUi = absGrad<=epsi;
CC = bwconncomp(maskUi);
for jjj=1:1:CC.NumObjects
    ind = CC.PixelIdxList{jjj};
    if(jjj<CC.NumObjects)
        ind(end+1) = ind(end)+1;
    end
    pui(ind) = sum(pui(ind))/length(ind);
end
puiT = pui;
T = [];
t = 0;

while norm(pui)>delta
    dt = findT(ui,pui);
    
    ui = ui+ dt*pui; % update ui
    t = t+dt;
    [uix,uiy] = gradIdo(ui); absGrad = abs(uix+uiy*1i);
    maskUi = absGrad<=epsi;
    CC = bwconncomp(maskUi);
    for jjj=1:1:CC.NumObjects
        ind = CC.PixelIdxList{jjj};
        if(jjj<CC.NumObjects)
            ind(end+1) = ind(end)+1;
        end
        pui(ind) = sum(pui(ind))/length(ind);
    end
    
    puiT = [puiT,pui];
    T = [T;t];
    
    
    
    figure(54);plot(-pui);
    figure(53);plot(ui);
    
    
end
fast_TV_time = toc;

phIdo = puiT(:,2:1:end)-puiT(:,1:1:end-1);
phIdo = phIdo * diag(T);
timePoints = [0.65,0.95,1.2,1.65,2.5];



%% Compute time scale for exponantial decay, resample uniformly
utmp = utmpMain(:,:,1:1:2800);
f1d = squeeze(utmp(2,:,1));
[m, n, z] = size(utmp);
Psi_vecs = squeeze(utmp(2,:,:));
Psi_diff = Psi_vecs(:,2:end) - Psi_vecs(:,1:end-1);
Psi_diff_sq_norms = diag(Psi_diff' * Psi_diff);
Psi_diffs_dot_Psi = diag(Psi_diff' * Psi_vecs(:,1:end-1));

dt_tildas = - Psi_diff_sq_norms ./ Psi_diffs_dot_Psi;
sum(dt_tildas < 0)
dt_tildas = [0;dt_tildas];
sum(isnan(dt_tildas))
t_tilda = cumsum(dt_tildas);

[Psi_uniform_new, t_new] = resample(Psi_vecs(:,1:end)', t_tilda);
Psi_uniform_new = Psi_uniform_new';
utx=Psi_uniform_new;

utx=Psi_uniform_new;
h71 = figure(71); mesh(utx);
grid on;
h71.Children.CameraPosition = 1.0e+04 *[ 1.4607   -0.0627    0.0005 ];
set(gca,'FontSize',240, 'TickLabelInterpreter', 'latex');

%% Compute DMD + Diagonalize into Eigenvalues and Mods

DMD_mods1T = [];
DMD_mods2T = [];
alphas1 = [];
alphas2 = [];
starts = [1,500,900,1100,1400,2200]; ends = [400,800,1000,1250,2000,2800];
t_tildaStarts = t_tilda(starts);t_tildaEnds=t_tilda(ends);
for iii=1:1:length(starts)
    band = find(t_new>t_tildaStarts(iii)&t_new<t_tildaEnds(iii));
    V1 = Psi_uniform_new(:,band(1):1:band(end)-1);
    V2 = Psi_uniform_new(:,band(1)+1:1:band(end));
    
    num_mods = 2;                           % Set number of DMD mods to compute
    [U, Sigma, W_T] = svd(V1, 'econ');
    diag(Sigma)'
    U_r = U(:,1:num_mods);Sigma_r = Sigma(1:1:num_mods,1:1:num_mods);W_T_r=W_T(:,1:1:num_mods);
    figure();imagesc(V1-U_r*Sigma_r*W_T_r')
    X = U_r' * V1;
    Y = U_r' * V2;
    F = sylvester(X*X', X*X', X*Y'+Y*X');   %  F = U' A U
    [y, mu,~] = eigs(F);
    DMD_mods = U_r * real(y);
    DMD_mods1T = [DMD_mods1T,DMD_mods(:,1)];
    DMD_mods2T = [DMD_mods2T,DMD_mods(:,2)];
    alphas1 = [alphas1, Psi_uniform_new(:, band(1))' * DMD_mods(:,1)];
    alphas2 = [alphas2, Psi_uniform_new(:, band(1))' * DMD_mods(:,2)];
end




scl = 1 + zebra_flag*3;
tick_jmp = 25*(1 + zebra_flag);
tick_num = 3 + 2*zebra_flag;
fsz = 125;
lnw = 8;


%---- rDMD based reconstruction ------------------

reconstructed_flow_DMD = zeros(size(DMD_mods2T));
for iii=1:1:length(starts)
    if(iii<length(starts))
        comp = DMD_mods2T(:,iii)-(DMD_mods2T(:,iii+1)'*DMD_mods2T(:,iii)/norm(DMD_mods2T(:,iii+1))^2)*DMD_mods2T(:,iii+1);
    else
        comp = DMD_mods2T(:,end);
    end
    comp = comp/norm(comp);
    alpha = f1d*comp;
    reconstructed_flow_DMD(:, iii) = alpha*comp;
end


filter = zeros(size(T));
band1 = find(T>=0 & T<timePoints(1));
filter1 = filter;
filter1(band1) = 1;
phIdoBand1 = phIdo*filter1(:);
figure(800); subplot(2,6,1);plot(phIdoBand1, 'b' ,'linewidth', lnw);hold on; plot(reconstructed_flow_DMD(:,1), ':k' ,'linewidth', 2*lnw); plot(imgBand1(2,:), '--r' ,'linewidth', lnw/2); hold off; grid on;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min(imgBand1(2,:));
ystop = max(imgBand1(2,:));
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


band2 = find(T>=timePoints(1) & T<timePoints(2));
filter2 = filter;
filter2(band2) = 1;
phIdoBand2 = phIdo*filter2(:);
figure(800); subplot(2,6,2); plot(phIdoBand2, 'b' ,'linewidth', lnw);hold on; plot(reconstructed_flow_DMD(:,2), ':k' ,'linewidth', 2*lnw); plot(imgBand2(2,:), '--r' ,'linewidth', lnw/2); hold off; grid on;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min(imgBand2(2,:));
ystop = max(imgBand2(2,:));
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


band3 = find(T>=timePoints(2) & T<timePoints(3));
filter3 = filter;
filter3(band3) = 1;
phIdoBand3 = phIdo*filter3(:);
figure(800); subplot(2,6,3); plot(phIdoBand3, 'b' ,'linewidth', lnw);hold on; plot(reconstructed_flow_DMD(:,3), ':k' ,'linewidth', 2*lnw);plot(imgBand3(2,:), '--r', 'linewidth', lnw/2); hold off; grid on;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min(imgBand3(2,:));
ystop = max(imgBand3(2,:));
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


Lgnd = legend('Fast TV decomposition','DMD based TV decomposition', 'Reference Decomposition','FontSize', fsz, 'Location', 'South', 'Orientation', 'horizontal', 'interpreter', 'latex');
Lgnd.Units = 'points';
Lgnd.Position(1) = -720;
Lgnd.Position(2) = 380;

band4 = find(T>=timePoints(3)& T<timePoints(4));
filter4 = filter;
filter4(band4) = 1;
phIdoBand4 = phIdo*filter4(:);
figure(800); subplot(2,6,4); plot(phIdoBand4, 'b' ,'linewidth', lnw);hold on; plot(reconstructed_flow_DMD(:,4), ':k' ,'linewidth', 2*lnw); plot(imgBand4(2,:), '--r' ,'linewidth', lnw/2); hold off; grid on;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min(imgBand4(2,:));
ystop = max(imgBand4(2,:));
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);

band5 = find(T>=timePoints(4)& T<timePoints(5));
filter5 = filter;
filter5(band5) = 1;
phIdoBand5 = phIdo*filter5(:);
figure(800); subplot(2,6,5); plot(phIdoBand5, 'b' ,'linewidth', lnw);hold on; plot(reconstructed_flow_DMD(:,5), ':k' ,'linewidth', 2*lnw); plot(imgBand5(2,:), '--r' ,'linewidth', lnw/2);  hold off; grid on;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min(imgBand5(2,:));
ystop = max(imgBand5(2,:));
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);

band6 = find(T>=timePoints(5));
filter6 = filter;
filter6(band6) = 1;
phIdoBand6 = phIdo*filter6(:);
figure(800); subplot(2,6,6); plot(phIdoBand6, 'b' ,'linewidth', lnw);hold on;  plot(reconstructed_flow_DMD(:,6), ':k' ,'linewidth', 2*lnw); plot(imgBand6(2,:), '--r' ,'linewidth', lnw/2); hold off; grid on;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min(imgBand6(2,:));
ystop = max(imgBand6(2,:));
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);





err1 = phIdoBand1 - imgBand1(2,:)';
err1_dmd = reconstructed_flow_DMD(:,1) - imgBand1(2,:)';
figure(900); subplot(2,6,1);plot(err1, 'b' ,'linewidth', lnw); hold on; plot(err1_dmd, ':k' ,'linewidth', 2*lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(err1_dmd) min(err1)]);
ystop = max([max(err1_dmd) max(err1)]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


err2 = phIdoBand2 - imgBand2(2,:)';
err2_dmd = reconstructed_flow_DMD(:,2) - imgBand2(2,:)';
figure(900); subplot(2,6,2);plot(err2, 'b' ,'linewidth', lnw); hold on; plot(err2_dmd, ':k' ,'linewidth', 2*lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(err2_dmd) min(err2)]);
ystop = max([max(err2_dmd) max(err2)]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


err3 = phIdoBand3 - imgBand3(2,:)';
err3_dmd = reconstructed_flow_DMD(:,3) - imgBand3(2,:)';
figure(900); subplot(2,6,3);plot(err3, 'b' ,'linewidth', lnw);hold on; plot(err3_dmd, ':k' ,'linewidth', 2*lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(err3_dmd) min(err3)]);
ystop = max([max(err3_dmd) max(err3)]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);

Lgnd = legend('Fast TV error','DMD based TV error', 'FontSize', fsz, 'Location', 'South', 'Orientation', 'horizontal', 'interpreter', 'latex');
Lgnd.Units = 'points';
Lgnd.Position(1) = 470;
Lgnd.Position(2) = 360;

err4 = phIdoBand4 - imgBand4(2,:)';
err4_dmd = reconstructed_flow_DMD(:,4) - imgBand4(2,:)';
figure(900); subplot(2,6,4);plot(err4, 'b' ,'linewidth', lnw); hold on; plot(err4_dmd, ':k' ,'linewidth', 2*lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(err4_dmd) min(err4)]);
ystop = max([max(err4_dmd) max(err4)]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


err5 = phIdoBand5 - imgBand5(2,:)';
err5_dmd = reconstructed_flow_DMD(:,5) - imgBand5(2,:)';
figure(900); subplot(2,6,5);plot(err5, 'b' ,'linewidth', lnw); hold on; plot(err5_dmd, ':k' ,'linewidth', 2*lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(err5_dmd) min(err5)]);
ystop = max([max(err5_dmd) max(err5)]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);

err6 = phIdoBand6 - imgBand6(2,:)';
err6_dmd = reconstructed_flow_DMD(:,6) - imgBand6(2,:)';
figure(900); subplot(2,6,6);plot(err6, 'b' ,'linewidth', lnw);grid on; hold on; plot(err6_dmd, ':k' ,'linewidth', 2*lnw); hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(err6_dmd) min(err6)]);
ystop = max([max(err6_dmd) max(err6)]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);






alpha1_mtx = repmat(alphas1, n, 1);
alpha2_mtx = repmat(alphas2, n, 1);

DMD_mods1T = DMD_mods1T.*alpha1_mtx;
DMD_mods2T = DMD_mods2T.*alpha2_mtx;


% ----- Plot Both mods 1 and 2 together ------

figure(999); subplot(2,6,1);plot(DMD_mods1T(:,1), 'Color',[0 0.5451 0.5451] ,'linewidth', lnw); hold on; plot(DMD_mods2T(:,1), '-.', 'Color',[0.8500 0.3250 0.0980] ,'linewidth', lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(DMD_mods1T(:,1)) min(DMD_mods2T(:,1))]);
ystop = max([max(DMD_mods1T(:,1)) max(DMD_mods2T(:,1))]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


figure(999); subplot(2,6,2);plot(DMD_mods1T(:,2), 'Color',[0 0.5451 0.5451] ,'linewidth', lnw); hold on; plot(DMD_mods2T(:,2), '-.','Color',[0.8500 0.3250 0.0980] ,'linewidth', lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(DMD_mods1T(:,2)) min(DMD_mods2T(:,2))]);
ystop = max([max(DMD_mods1T(:,2)) max(DMD_mods2T(:,2))]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


figure(999); subplot(2,6,3);plot(DMD_mods1T(:,3), 'Color',[0 0.5451 0.5451] ,'linewidth', lnw); hold on; plot(DMD_mods2T(:,3), '-.', 'Color',[0.8500 0.3250 0.0980] ,'linewidth', lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(DMD_mods1T(:,3)) min(DMD_mods2T(:,3))]);
ystop = max([max(DMD_mods1T(:,3)) max(DMD_mods2T(:,3))]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);

Lgnd = legend('$\xi_1^k$','$\xi_2^k$', 'FontSize', 0.8*fsz, 'Location', 'South', 'Orientation', 'horizontal', 'interpreter', 'latex');
Lgnd.Units = 'points';
Lgnd.Position(1) = 170;
Lgnd.Position(2) = 600;


figure(999); subplot(2,6,4);plot(DMD_mods1T(:,4), 'Color',[0 0.5451 0.5451] ,'linewidth', lnw); hold on; plot(DMD_mods2T(:,4), '-.', 'Color',[0.8500 0.3250 0.0980] ,'linewidth', lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(DMD_mods1T(:,4)) min(DMD_mods2T(:,4))]);
ystop = max([max(DMD_mods1T(:,4)) max(DMD_mods2T(:,4))]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


figure(999); subplot(2,6,5);plot(DMD_mods1T(:,5), 'Color',[0 0.5451 0.5451] ,'linewidth', lnw); hold on; plot(DMD_mods2T(:,5), '-.', 'Color',[0.8500 0.3250 0.0980] ,'linewidth', lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(DMD_mods1T(:,5)) min(DMD_mods2T(:,5))]);
ystop = max([max(DMD_mods1T(:,5)) max(DMD_mods2T(:,5))]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);


figure(999); subplot(2,6,6);plot(DMD_mods1T(:,6), 'Color',[0 0.5451 0.5451] ,'linewidth', lnw); hold on; plot(DMD_mods2T(:,6), '-.', 'Color',[0.8500 0.3250 0.0980] ,'linewidth', lnw); grid on; hold off;
set(gca,'FontSize',fsz, 'TickLabelInterpreter', 'latex');
ystart = min([min(DMD_mods1T(:,6)) min(DMD_mods2T(:,6))]);
ystop = max([max(DMD_mods1T(:,6)) max(DMD_mods2T(:,6))]);
ygap = ystop - ystart;
y_scale = yticks;
set(gca,'xtick',scl*[0:tick_jmp:100], 'ylim', [ystart - 0.1*ygap, ystop + 0.1*ygap]);%, 'ytick', [y_scale(1) (y_scale(1)+y_scale(end))/2 y_scale(end)]);



h1=figure(800);
write_pdf_New_Image('fast_TV_pulse_flow',h1,8000,2000)

h2=figure(900);
write_pdf_New_Image('fast_TV_pulse_error',h2,8000,2000)


h9=figure(999);
write_pdf_New_Image('DMD_mods_1_and_2',h9,8000,2000)

write_pdf_New_Image('pulse_flow_mesh',h70,6000,4000)
write_pdf_New_Image('pulse_flow_dmd_mesh',h71,6000,4000)

fprintf('done!\n');
close all
%%%%%%%%%%%%%%%%%%%%%%%%%% find t %%%%%%%%%%%%%%%%%%%%%%%%%
function [dt] = findT(f,pf)
[~,dpfy] = gradIdo(pf);
[~,dfy] = gradIdo(f);
dt = -dfy(:)./dpfy(:);
dt = dt(find(dt>0));
dt = min(dt);
dt = dt(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient (forward difference)
function [fx,fy] = gradIdo(P)
fx = P(:,[2:end end])-P;
fy = P([2:end end],:)-P;
end