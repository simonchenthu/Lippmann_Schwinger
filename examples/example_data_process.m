% This code post-process the data
clear;
%% starting up the path
LS_startup()
%% Frequency
omega = 2.0^11;

%% Domain parameters
a = 1;
b = 1;

h = 2.0^(-12);

x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

X = repmat(x',1,m);
Y = repmat(y,n,1);

%% Boundary
% define the boundary parameters
xb = 0; yb = 0; rb = 0.3;
BD = boundary(xb,yb,rb);

%% the heterogeneity
xn = 0; yn = 0; An = -2^(-7); rn = 0.25; MDName = 2;
nu = @(x,y) An*exp_bump(sqrt((x-xn).^2+(y-yn).^2)/rn);
gradx_nu = @(x,y) -An*2*(x-xn)/rn^2/((1-((x-xn).^2+(y-yn).^2)/rn^2)^2).*...
                    exp_bump(sqrt((x-xn).^2+(y-yn).^2)/rn);
grady_nu = @(x,y) -An*2*(y-xn)/rn^2/((1-((x-xn).^2+(y-yn).^2)/rn^2)^2).*...
                    exp_bump(sqrt((x-xn).^2+(y-yn).^2)/rn);
MD = medium(nu,rn,xn,yn,An,gradx_nu,grady_nu);

%%  Source term
% width in direction
sigma0 = 2^(-5);

%% define husimi transform class
Hk = husimi(omega,h,BD);

%% Liouville trajectory class
TR = trajectory(BD,MD);
ode_options = odeset('RelTol',1e-10); T = 3;

%% Position parameters
dtheta = pi/30; 
% center of source
theta_s = pi/4; 
% relative Incident direction
theta_i = -pi/2+dtheta:dtheta:pi/2-dtheta; Ni = length(theta_i);
% relative receiver position
theta_r = 0:dtheta:2*pi-dtheta; Nr = length(theta_r);
% relative outgoing direction
theta_o = -pi/2+dtheta:dtheta:pi/2-dtheta; No = length(theta_o);

%% H-transform
Hu_bd = zeros(No,Nr,Ni);
% load(fullfile('data',['Hu_medium',int2str(MDName),'_An',num2str(An),...
%         '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_thetas',num2str(rad2deg(theta_s)),'.mat']),...
%         'Hu_bd');
for i = 1:Ni
    %% Load data
    load(fullfile('data',['u_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'_thetai',num2str(rad2deg(theta_i(i))),'.mat']),...
        'u');
    
    %% Husimi transform
    Hu_bd(:,:,i) = Hk.transform(X,Y,u,theta_s,theta_r,theta_o);
    
    %% Liouville limit
    [theta_rL,theta_oL] = TR.outgoing(theta_i(i),theta_s,T,ode_options);
    
    %% Plot
    figure(333)
    Hk.plot(Hu_bd(:,:,i),theta_r,theta_o)
    title(['Hu, \theta_s = ',num2str(rad2deg(theta_s)),...
        ', \theta_i = ',num2str(rad2deg(theta_i(i)))]);
    hold on;
    scatter(rad2deg(theta_rL),rad2deg(theta_oL),...
        100,'x','LineWidth',2,'MarkerEdgeColor','r'); pause(0.1);
    hold off;
    
    print(fullfile('plot',['Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'_thetai',num2str(rad2deg(theta_i(i))),'.pdf']),...
        '-dpdf');
    
    savefig(fullfile('plot',['Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'_thetai',num2str(rad2deg(theta_i(i))),'.fig']));
    
end

save(fullfile('data',['Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'.mat']),...
        'Hu_bd');
    
%% Integral in theta_o
% compute int Hu dtheta_o
int_thetao_Hu = permute(sum(Hu_bd,1)*dtheta,[3,2,1]);
% plot
figure(334)
imagesc(rad2deg(theta_r),rad2deg(theta_i),int_thetao_Hu); 
ylabel('\theta_i'); xlabel('\theta_r'); colorbar('southoutside');
axis('equal');
ylim([-90,90]); yticks(-90:30:90); set(gca,'YDir','normal');
% Liouville limit
dtheta_plot = pi/100;
theta_i_plot = -pi/2+dtheta_plot:dtheta_plot:pi/2-dtheta_plot;
[theta_rL_plot,~] = TR.outgoing(theta_i_plot,theta_s,T,ode_options);
% plot
hold on;
plot(rad2deg(theta_rL_plot),rad2deg(theta_i_plot),'LineWidth',2,'Color','r');
hold off;
% Save
print(fullfile('plot',['int_o_Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'.pdf']),'-dpdf');
savefig(fullfile('plot',['int_o_Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'.fig']));

%% Integral in theta_r

% Extend theta_o to thetaor
theta_or = 0:dtheta:2*pi-dtheta; Nor = length(theta_or);
Hu_bd_e = cat(1,Hu_bd,zeros(Nor-No,Nr,Ni));
% convert to theta_o relative to receiver
for q = 1:Nr
    Hu_bd_e(:,q,:) = circshift(Hu_bd_e(:,q,:),-(No-1)/2+q-1,1);
end
int_thetar_Hu = permute(sum(Hu_bd_e,2)*dtheta,[3,1,2]);

% plot
figure(335)
imagesc(rad2deg(theta_or),rad2deg(theta_i),int_thetar_Hu); 
ylabel('\theta_i'); xlabel('\theta_{or}'); colorbar('southoutside');
axis('equal');
ylim([-90,90]); yticks(-90:30:90); set(gca,'YDir','normal');

% Liouville limit
dtheta_plot = pi/100;
theta_i_plot = -pi/2+dtheta_plot:dtheta_plot:pi/2-dtheta_plot;
[theta_rL_plot,theta_oL_plot] = TR.outgoing(theta_i_plot,theta_s,T,ode_options);
theta_orL_plot = wrapTo2Pi(theta_rL_plot+theta_oL_plot);

% plot
hold on;
plot(rad2deg(theta_orL_plot),rad2deg(theta_i_plot),'LineWidth',2,'Color','r');
hold off;

% Save
print(fullfile('plot',['int_r_Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'.pdf']),'-dpdf');
savefig(fullfile('plot',['int_r_Hu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'.fig']));
