function result = nonlinFit(y, fs, xStart, xMax, disp)
    if ~exist('disp','var')
        disp = 0;
    end
    
    % Fitting functions
    yP = @(x,k,a,b) (1 - k) * gampdf(x,a,b) / max(gampdf(x,a,b));
    yR = @(x,k,a,b) k * gamcdf(x,a,b);
    bolus = @(x,p) p(1) * (yP(x - p(2),p(3),p(4),p(5)) + yR(x - p(2),p(3),p(4),p(5)));
    
    % Recalculation of peak time
    xMax = xMax - xStart + 1;
    
    % Baseline removal
    yMin = min(y);
    y = y - yMin;
    yMax = y(xMax);
    x = [0:length(y)-1]' / fs;
    
    % Initial guesses and bounds
    p0 = [yMax;         0;          0.9;    7.5;    1.0];
    lb = [0.1 * yMax;   0;          0;      0;      0];
    ub = [10 * yMax;    xMax/fs;    1;      30;     25];
    
    cost = @(p) bolus(x,p) - y;
    
    opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective', 'Display', 'off');
    tic
    [p, ~, res] = lsqnonlin(cost, p0, lb, ub, opts);
    TC = toc;
    NE = rms((res .^ 2) / length(y));
    
    result.yP =  p(1) * yP(x - p(2),p(3),p(4),p(5));
    result.yR =  p(1) * yR(x - p(2),p(3),p(4),p(5));
    result.t = x;
    result.A = p(1);
    result.t0 = p(2) + (xStart - 1) / fs;
    result.ttp = x(xMax) - p(2);
    result.tpeak = result.t0 + result.ttp;
    result.k = p(3);
    result.alpha = p(4);
    result.beta = p(5);
    result.b = yMin;
    result.NE = NE;
    result.TC = TC;
    
    if disp
        figure; 
        hold on;
        plot(x, y + yMin, 'LineWidth', 2);
        plot(x, bolus(x,p) + yMin, 'LineWidth', 2);
        
        plot(x, p(1) * yP(x - p(2),p(3),p(4),p(5)) + yMin, 'LineWidth', 2);
        plot(x, p(1) * yR(x - p(2),p(3),p(4),p(5)) + yMin, 'LineWidth', 2); 
        
        plot(x, y - bolus(x,p) + yMin, 'LineWidth', 2);
        
        scatter(x(xMax), yMax + yMin, 'filled');
        
        hold off;
        legend('Data', 'Total model', 'First-pass', 'Recirculation', 'Residual', 'Peak');
        title(strjoin(["RMSE = " NE],""));
        xlabel("Time [s]"); ylabel("Conductivity change [-]");
        grid on;
    end
end

