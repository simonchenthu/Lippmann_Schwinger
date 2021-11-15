classdef husimi
    % This class define the Husimi transform
    properties (SetAccess = public)
        % Coherent state
        phi
        % frequency
        omega
        % step size
        h
        % boundary information
        xb
        yb
        rb
    end
    
    methods
        % constructor
        function Hk = husimi(omega,h,BD)
            % Prototype function
            Hk.phi = @(omega,x,y,thetax,thetay) 0.5*(omega/pi)^(3/2)*...
                exp(-omega/2*(x.^2+y.^2) + 1i*omega*(thetax*x+thetay*y));
            % frequency
            Hk.omega = omega;
            % step size
            Hk.h = h;
            % boundary information
            Hk.xb = BD.xb;
            Hk.yb = BD.yb;
            Hk.rb = BD.rb;
        end
        
        % Husimi Transform
        function Hu = transform(Hk,X,Y,u,theta_s,theta_r,theta_o)
            Nr = length(theta_r); No = length(theta_o);
            
            Hu = zeros(No,Nr);
            for p = 1:Nr*No
                [i,j] = ind2sub([No,Nr],p);    % use linear indices to avoid multiple loop
                
                theta_ra = theta_s + theta_r(j);
                x = Hk.xb +  Hk.rb*cos(theta_ra);
                y = Hk.yb +  Hk.rb*sin(theta_ra);
                
                theta_oa = theta_ra + theta_o(i);
                thetax = cos(theta_oa); thetay = sin(theta_oa);
                
                Hu(i,j) = abs(sum(u.*Hk.phi(Hk.omega,x-X,y-Y,thetax,thetay),'all')*Hk.h^2).^2;
            end
        end
        
        % Plot
        function plot(~,Hu,theta_r,theta_o)
            imagesc(rad2deg(theta_r),rad2deg(theta_o),Hu); 
            ylabel('\theta_o'); xlabel('\theta_r'); colorbar('southoutside');
            axis('equal');
            ylim([-90,90]); yticks(-90:30:90); set(gca,'YDir','normal');
        end
        
    end
end