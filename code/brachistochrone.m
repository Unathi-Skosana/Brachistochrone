function brachistochrone
    % Adaptive step control constants
    sft  = 0.9;
    exp_bot = 0.25;
    exp_top = 0.2;
    eps = 1.0e-10;
    
    % % % % Set up % % % %
    dy(1) = -10000; % Initial value of derivative;
    y(1) = -0.001; % Starting point in y;
    x(1) = (2/3)*(y(1)/dy(1)); % Starting point in x;
    x_e  = 100; % End point 
    dx = 0.001;  % Initial step size
    g = -9.8;  % Gravitional constant
    
    k = 1; % Index
    while (x(k) < x_e)
        afdeq = deq(x(k), y(k), dy(k));
        dela0 = eps * (abs(y(k))  + 2.0 * dx * abs(dy(k)) + 1.0e-30);
        delb0 = eps * (abs(dy(k)) + 2.0 * dx * abs(afdeq) + 1.0e-30);
        
        % % %  Adaptive step size control % % %
      
        % first small step
        [x_1, y_1, dy_1] = rk42(x(k), dx, y(k), dy(k));
        
        % second small step
        [x_s, y_s, dy_s] = rk42(x_1, dx, y_1, dy_1);
        
        % single large step
        [x_f, y_f, dy_f] = rk42(x(k), 2.0 * dx, y(k), dy(k));
        
        % step size control %
        dela = y_s - y_f;
        delb = dy_s - dy_f;
        delta = max(abs(dela/dela0),abs(delb/delb0));
        
         if (delta > 1.0)
             dx = sft * dx * (1.0/delta)^exp_bot;
         else
              x(k+1)  = x_f;
              y(k+1)  = y_s + dela/15.0;
              dy(k+1) = dy_s + delb/15.0;
              dx1 = sft * dx * (1.0/delta)^exp_top;
              dx = max(100.0 * dx1, dx1);
              k = k + 1;
         end
         
         % Prevent overshooting %
         if (x(k) + 2 * dx > x_e)
             dx = 0.5 * (x_e - x(k));
         end
    end
     
    %%% Straight line connecting A and B. %%%
    m  = (y(k) - y(1))/(x(k) - x(1));
    c = y(k) - m*x(k); 
    for l=1:k
        % Straight line %
        str_l(l) = m*x(l) + c;
        
        % Variation function 1%  % x^2 - c*x %
        vy(l)  = 0.1*(x(l)^2 - x(k)*x(l));
        
        % Derivative of variation function % % 2x - c% 
        dvy(l) = 0.1*(2*x(l)- x(k));
        
        % Variation function 2 % % sin(pi*x/x_e) %
        vvy(l) = -50*sin(pi*x(l)/x(k));
        
        % Derivative of variation function 2 % (pi/xe)cos(pi*x(l)/x_e) %
        dvvy(l) = -(50*pi/x(k))*cos(pi*x(l)/x(k));
       
    end
       
        
    intrgl = 0.0;
    str_intrgl = 0.0;
    var_intrgl_1 = 0.0;
    var_intrgl_2 = 0.0;
    
    % General trapezium rule %
    for i = 1:k - 1
         h   = x(i+1) - x(i);
         
         ym      = 0.5 * (y(i+1) +  y(i));
         varym   = ym + 0.5 * (vy(i+1) + vy(i));
         vvarym  = ym + 0.5 * (vvy(i+1) + vvy(i));
         y_strl  = 0.5 * (str_l(i+1) + str_l(i));
         
         dym_strl = 0.5 * (m + m); % Constant derivative.
         dym = 0.5 * (dy(i+1) + dy(i));
         vardym = dym + 0.5 * (dvy(i+1) + dvy(i));
         vvardym = dym + 0.5 * (dvvy(i+1) + dvvy(i));
         
         intrgl = intrgl + h * sqrt((1 + dym^2) /(2*ym*g));
         str_intrgl = str_intrgl + h * sqrt((1 + dym_strl^2) /(2*y_strl*g));
         var_intrgl_1 = var_intrgl_1 + h * sqrt((1 + vardym^2) /(2*varym*g));
         var_intrgl_2 = var_intrgl_2 + h * sqrt((1 + vvardym^2) /(2*vvarym*g));
         
    end
    
    % Printing to console
    fprintf('time taken for cycloid: %.5f seconds\n', intrgl);
    fprintf('time taken for a straight line: %.5f seconds\n', str_intrgl);
    fprintf('time taken for a small variation1: %.5f seconds\n', var_intrgl_1);
    fprintf('time taken for a small variation2: %.5f seconds\n', var_intrgl_2);
    
    % Plotting 
    plot(x, str_l, '-b', x, y, '-r', x, vy + y, '-g', x , y  +vvy, '-k');
    title('Shortest descent candidate trajectories.');
    xlabel('x (m)');
    ylabel('y (m)');
    legend('straight line', 'cycloid', 'small variation1', 'small variation2');
    
    % % % End of main % % % 
   
    % Fourth order runge-kutta algorithm %
    function [x2,y2,dy2] = rk42(x1,dx,y1,dy1)
        x2 = x1 + dx;
        q1 = dx * dy1;
        p1 = dx * deq(x1,y1,dy1);
        q2 = dx * (dy1 + 0.5 * p1);
        p2 = dx * deq(x1 + 0.5 * dx, y1 + 0.5 * q1, dy1 + 0.5 * p1);
        q3 = dx * (dy1 + 0.5 * p2);
        p3 = dx * deq(x1 + 0.5 * dx, y1 + 0.5 * q2, dy1 + 0.5 * p2);
        q4 = dx * (dy1 + p3);
        p4 = dx * deq(x1 + dx, y1 + q3, dy1 + p3);
        y2 = y1 + (q1 + 2.0 * q2 + 2.0 * q3 + q4)/6.0;
        dy2 = dy1 + (p1 + 2.0 * p2+2.0 * p3 + p4)/6.0;
    end

    function [y] = deq(~, y ,dy)
        y = -(1 + dy^2)/(2*y);
    end
end