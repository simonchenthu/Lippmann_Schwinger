classdef source
    % This class contain the information of the medium
    
    properties (SetAccess = public)
        % concentration
        sigma
        % radius of boundary
        rb
        % center of boundary
        xb
        yb
        % Prototype function
        S
        % Normalization constant
        Cd
    end
     
    methods
        % constructor
        function SRC = source(BD,sigma)
            SRC.rb = BD.rb;
            SRC.xb = BD.xb;
            SRC.yb = BD.yb;
            % Width
            SRC.sigma = sigma;
            % Prototype function
            SRC.S = @(x,y,thetax,thetay) exp(-sigma^2/2*(x.^2+y.^2) + 1i*(thetax*x+thetay*y));
            % Normalization constant
            SRC.Cd = sqrt(2)*(sigma/sqrt(pi))^(3/2);
        end
        
        % output the source matrix
        function f = source_func(SRC,omega,theta_s,theta_i,X,Y)
            % position of source
            xs = SRC.xb + SRC.rb*cos(theta_s); ys = SRC.yb + SRC.rb*sin(theta_s);
            % incident direction (absolute direction)
            theta_ia = theta_i + pi + theta_s;
            theta = [cos(theta_ia),sin(theta_ia)];
            % Source (Delta u + omega^2 = -f)
            f = omega^(5/2)*SRC.Cd*SRC.S(omega*(X-xs),omega*(Y-ys),theta(1),theta(2));
        end
        
        % plot function
        function plot(SRC,omega,theta_s,theta_i,x,y,mode,BD)
            f = SRC.source_func(omega,theta_s,theta_i,x,y');
            if mode == 'a'
                f = abs(f);
            elseif mode == 'r'
                f = real(f);
            elseif mode == 'i'
                f = imag(f);
            else
                error(message('unknownMode'));
            end
            
            % plot source
            imagesc(x,y,f); pbaspect([1 1 1]);
            ylabel('y'); xlabel('x'); set(gca,'YDir','normal'); colorbar;
            
            if nargin >=8
                % plot boundary
                rectangle('Position',[BD.xb-BD.rb BD.yb-BD.rb 2*BD.rb 2*BD.rb],...
                    'Curvature',[1,1],'LineWidth',2);
                hline_b = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color','k');
                legend(hline_b,'$\partial\Omega$','Interpreter','latex');
            end
        end
        
    end
end