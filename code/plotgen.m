function plotgen
    % Function for generating plots with different values of r,
    % on the same graph.
    rho0 = 23;
    [t,rx, ry, rz] = lorenz(rho0, 10001);
 %   [t1, rx1, ry1, rz1] = lorenz(rho1, 5001);
%    [t2, rx2, ry2, rz2] = lorenz(rho2, 5001); 

    plot(t, rx, '-r', t, ry, '-b', t , rz, '-g');
    %plot(t, rx, '-r', t1, rx1, '-b', t2, rx2, '-g');
    xlabel('time (s)');
    ylabel('system state (x,y,z)');
    legend("x", "y", "z");
    title(strcat('rho =',num2str(rho0)," sigma = 10", " beta = 8/3 "));
   % legend(strcat("rho=",num2str(rho0)), strcat("rho=",num2str(rho1)), strcat("rho=",num2str(rho2)));
    
    function [t,rx,ry,rz] = lorenz(rho, num)
        % Constants.
        sigma = 10;
        beta = 8/3;

        % Initial conditions.
        t(1)  = 0;
        r(1)  = 0;
        r(2)  = 1;
        r(3)  = 0;

        % Storing the functions for later access.
        rx(1) = r(1);
        ry(1) = r(2);
        rz(1) = r(3);

        % Maximum time for running.
        tmax = 500;

        % Step size.
        dt = tmax/(num-1);

        % Differentiating using the Runge-Kutta method.
        for i=1:num-1
        [t(i+1), r] = rk42(t(i), dt, r);
        rx(i+1) = r(1);
        ry(i+1) = r(2);
        rz(i+1) = r(3);
        end

        % Plotting the function.
        %plot3(ry,rz,rx, '-b');
        %xlabel('x');
        %ylabel('y');
        %ylabel('z');

        % Runge-Kutta Method.
        function [t2, r2] = rk42(t1, dt, r1)
        t2  = t1 + dt;
        yk1 = dt * fdeq(t1, r1);
        yk2 = dt * fdeq(t1 + 0.5 * dt, r1 + 0.5 * yk1);
        yk3 = dt * fdeq(t1 + 0.5 * dt, r1 + 0.5 * yk2);
        yk4 = dt * fdeq(t1 + dt, r1 + yk3);
        r2  = r1 + (yk1 + 2.0 * yk2 + 2.0 * yk3 + yk4) / 6.0;
        end

        % Differential equation.
        function [ff] = fdeq(~, r)
            ff(1) = sigma*(r(2) - r(1));
            ff(2) = r(1)*(rho - r(3)) - r(2);
            ff(3) = r(1)*r(2) - beta*r(3);
        end
    end
end
