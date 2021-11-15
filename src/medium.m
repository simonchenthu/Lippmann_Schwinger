classdef medium
    % This class contain the information of the medium
    
    properties (SetAccess = public)
        % heterogeneity (function)
        nu
        % radius of support
        rn
        % center of heterogeneity
        xn
        yn
        % amplitude
        An
        % gradient of heterogeneity (function)
        gradx_nu
        grady_nu
    end
     
    methods
        % constructor
        function MD = medium(nu,rn,xn,yn,An,gradx_nu,grady_nu)
            MD.nu = nu;
            MD.rn = rn;
            MD.xn = xn;
            MD.yn = yn;
            MD.An = An;
            MD.gradx_nu = gradx_nu;
            MD.grady_nu = grady_nu;
        end
        
        % plot function
        function plot(MD,x,y,BD)
            % compute refractive index
            m = 1 + MD.nu(x,y');
            % plot medium
            imagesc(x,y,m); pbaspect([1 1 1]);
            ylabel('y'); xlabel('x'); set(gca,'YDir','normal'); colorbar;
            % plot support
            rectangle('Position',[MD.xn-MD.rn MD.yn-MD.rn 2*MD.rn 2*MD.rn],...
                'Curvature',[1,1],'LineWidth',2,'LineStyle','--'); 
            
            if nargin >=4
                hline_n = line(NaN,NaN,'LineWidth',2,'LineStyle','--','Color','k');
                % plot boundary
                hline_b = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color','k');
                rectangle('Position',[BD.xb-BD.rb BD.yb-BD.rb 2*BD.rb 2*BD.rb],...
                    'Curvature',[1,1],'LineWidth',2);
                % legend
                legend([hline_b,hline_n],'$\partial\Omega$','supp$(n-1)$','Interpreter','latex');
            end
        end
        
    end
end