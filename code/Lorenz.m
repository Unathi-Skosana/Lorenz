function Lorenz
    % Constants.
    sigma = 10;
    rho = 23;
    beta = 8/3;
    
    % Equilibrium points.
    
    cx1 = sqrt(beta*(rho - 1));
    cy1 = cx1;
    cz1 = rho - 1;
    
    cx2 = -cx1;
    cy2 = cx2;
    cz2 = cz1;
    

    % Initial conditions.
    t(1)  = 0;
    r(1)  = 2;
    r(2)  = 1;
    r(3)  = 10;

    % Storing the functions for later access.
    rx(1) = r(1);
    ry(1) = r(2);
    rz(1) = r(3);

    % Number of points.
    num = 35001;

    % Maximum time for running.
    tmax = 900;

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
    plot3(rx,ry,rz, 'r');
    hold on;
     
    %Plot equilibrium points.
    plot3(cx1,cy1,cz1, '-o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5,0.5,0.5]);
    plot3(cx2,cy2,cz2, '-o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5,0.5,0.5]);
    plot3(0,0,0, '-o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5,0.5,0.5]);
    plot3(rx(1),ry(1),rz(1), '-s', 'MarkerSize', 5, 'MarkerFaceColor', [0.5,0.5,0.5]);
    hold off;
    
    % Label axes.
    xlabel('y');
    ylabel('z');
    ylabel('x');
    
    % Label the graph.
    title(strcat(strcat('rho =',num2str(rho)), strcat(' sigma=', num2str(sigma)), strcat(" beta=", num2str(beta))));

    % Test
    test();
    
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
    
    % Reliability test
    function [ss] = test()
        if (rho > 1 && rho <= 24.74)
           fprintf("Reliability test...\n\n");
          fprintf("The final point is (%d,%d,%d) corresponds to either equilibrium point (%d, %d, %d) or (%d, %d, %d)\n", rx(num-1),ry(num-1),rz(num-1),cx1,cy1,cy1,cx2,cy2,cz2);
        end
        if (rho < 1.0)
          fprintf("Reliability test...\n\n");
          fprintf("The final point is (%d, %d, %d) corresponds to the null equilibrium (%d, %d, %d)\n", rx(num-1), ry(num-1), rz(num-1), 0, 0, 0);
        end
    end
end
