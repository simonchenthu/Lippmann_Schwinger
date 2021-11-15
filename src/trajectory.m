classdef trajectory
    % This class contains the information of Liouville trajectory
    properties (SetAccess = public)
        % boundary information
        xb
        yb
        rb
        % ode equation
        dxdp
    end
    
    methods
        % constructor
        function TR = trajectory(BD,MD)
            TR.dxdp = @(t,xp) [xp(3:4);...
                                MD.gradx_nu(xp(1),xp(2))/2;...
                                MD.grady_nu(xp(1),xp(2))/2];
            TR.xb = BD.xb;
            TR.yb = BD.yb;
            TR.rb = BD.rb;
        end
        
        % Wrap up the solver of trajectory
        function [x,p] = trajectory_solver(TR,t,x_i,p_i,ode_options)
            % t is the time vector (at least with initial and terminal t)
            % x_i, p_i is initial position and momentum
            % The 1st/2nd columns in x, p are the time sequence in the
            % 1st/2nd dimension
            xp_i = [x_i;p_i];
            [~,xp] = ode23(@(t,xp)TR.dxdp(t,xp),t,xp_i,ode_options);
            x = xp(:,1:2); p = xp(:,3:4);
        end
        
        % Plot trajectory
        function plot(TR,theta_s,theta_i,T)
            % T is a sequence of time for plot
            
            % initial condition
            theta_ia = theta_i + pi + theta_s;
            x_i = [TR.xb+TR.rb*cos(theta_s);TR.yb+TR.rb*sin(theta_s)];
            p_i = [cos(theta_ia);sin(theta_ia)];
            
            % Solve trajectory
            [x,~] = TR.trajectory_solver(T,x_i,p_i,[]);
            
            % Plot trajectory
            plot(x(:,1),x(:,2),'LineWidth',1.5);
        end
        
        % Compute outgoing position and direction
        function [theta_r,theta_o] = outgoing(TR,theta_i_list,theta_s_list,T,ode_options)
            % T is a guess of intersection time
            % theta_i_list and theta_s_list are vectors
            % theta_r and theta_o are matrices of 
            % the size (length(theta_i_list),length(theta_s_list))
            
            Ns = length(theta_s_list); Ni = length(theta_i_list);
            theta_r = zeros(Ni,Ns); theta_o = zeros(Ni,Ns);
            for q = 1:Ns*Ni
                [m,n] = ind2sub([Ni,Ns],q);
                theta_i = theta_i_list(m); theta_s = theta_s_list(n);
                
                % initial position and velocity
                theta_ia = theta_i + pi + theta_s;
                x_i = [TR.xb+TR.rb*cos(theta_s);TR.yb+TR.rb*sin(theta_s)];
                p_i = [cos(theta_ia);sin(theta_ia)];
                
                % solve for t when the two curves intersect
                t_o = fsolve(@(t)TR.BD_intersect(t,x_i,p_i,ode_options),T);
                if t_o < 1e-3
                    fprintf(['OutgoingTimeTooCloseToZero at theta_s = ',rad2deg(theta_s),...
                        ' theta_i = ',rad2deg(theta_i)])
                end
                
                % compute the outgoing position and velocity
                T_o = [0,t_o];
                [x,p] = TR.trajectory_solver(T_o,x_i,p_i,ode_options);
                
                % convert x/p to relative angle
                theta_ra = atan2(x(end,2),x(end,1));
                theta_oa = atan2(p(end,2),p(end,1));
                
                theta_r(m,n) = wrapTo2Pi(theta_ra-theta_s);
                theta_o(m,n) = wrapToPi(theta_oa-theta_ra);
            end
        end
    end
    
    methods (Access = private)
        % function to be solved for the outgoing point
        function f = BD_intersect(TR,t,x_i,p_i,ode_options)
            
            % solve position x
            [x,~] = TR.trajectory_solver([0,t],x_i,p_i,ode_options);
            
            % plug it into the boundary equation
            f = (x(end,1)-TR.xb)^2+(x(end,2)-TR.yb)^2-TR.rb^2;
        end
    end
    
end