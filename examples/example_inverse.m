% This code runs an example to show the setup of our inverse problem
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
xn = 0; yn = 0; An = -2^(-7); rn = 0.25; 
nu = @(x,y) An*exp_bump(sqrt((x-xn).^2+(y-yn).^2)/rn);
% nu = @(x,y) -0.3*exp(-320*(x.^2 + y.^2)).*(abs(x)<0.48).*(abs(y)<0.48);
gradx_nu = @(x,y) -An*2*(x-xn)/rn^2/((1-((x-xn).^2+(y-yn).^2)/rn^2)^2).*...
                    exp_bump(sqrt((x-xn).^2+(y-yn).^2)/rn);
grady_nu = @(x,y) -An*2*(y-xn)/rn^2/((1-((x-xn).^2+(y-yn).^2)/rn^2)^2).*...
                    exp_bump(sqrt((x-xn).^2+(y-yn).^2)/rn);
% medium class
MD = medium(nu,rn,xn,yn,An,gradx_nu,grady_nu); MDName = 2;
% plot
% figure(666); MD.plot(x,y,BD); title('Medium');

%% the trajectory
TR = trajectory(BD,MD);
    
%% we define the Lippmann-Schwinger operator
LS = LippmannSchwinger_precompute(x,y,omega,nu);

%%  Source term
% width in direction
sigma0 = 2^(-5);
% Source class
SRC = source(BD,sigma0);

% center of source/Relative Incident direction
theta_s = pi/4; 
dtheta = pi/30;
for theta_i = -pi/2+dtheta:dtheta:pi/2-dtheta
    % figure(667); SRC.plot(omega,theta_s,theta_i,x,y,'r',BD); title("Source");
    
    tic
    %% Building the incident wave
    u_inc = LS.apply_Green(SRC.source_func(omega,theta_s,theta_i,X,Y));
    
    %% building the right hand-side
    rhsDual = -omega^2*nu(X,Y).*u_inc;
    
    %% solving the Lippmann-Schwinger equation
    sigma = LS\rhsDual(:);
    
    %% computing the wavefield
    u_sca = LS.apply_Green(sigma);
    u = u_inc + u_sca;
    
    toc
    
    %% Plot
    u = reshape(u,n,m);
    
    figure(777)
    LS.plot(abs(u).^2,BD); title("|u|^2");
    T = 0:0.02:3;
    hold on; TR.plot(theta_s,theta_i,T); hold off;
    
    figure(778)
    LS.plot(real(u),BD); title("Re(u)")
    T = 0:0.02:3;
    hold on; TR.plot(theta_s,theta_i,T); hold off;
    
    %% Save
    
    save(fullfile('data',['u_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...%'_a',num2str(a),'_b',num2str(b),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'_thetai',num2str(rad2deg(theta_i)),'.mat']),...
        'u','-v7.3');
    
    figure(777)
    print(fullfile('plot',['u_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...%'_a',num2str(a),'_b',num2str(b),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'_thetai',num2str(rad2deg(theta_i)),'.pdf']),...
        '-dpdf');
    
    figure(778)
    print(fullfile('plot',['Reu_medium',int2str(MDName),'_An',num2str(An),...
        '_k',int2str(log2(omega)),'_h',int2str(-log2(h)),...%'_a',num2str(a),'_b',num2str(b),...
        '_sigma',int2str(-log2(sigma0)),...
        '_thetas',num2str(rad2deg(theta_s)),'_thetai',num2str(rad2deg(theta_i)),'.pdf']),...
        '-dpdf');

end
