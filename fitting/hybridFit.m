function result = hybridFit(y, fs, xStart, xMax, disp)
    if ~exist('disp','var')
        disp = 0;
    end
 
    % Recalculation of peak time
    xMax = xMax - xStart + 1;
    
    % Baseline removal
    yMin = min(y);
    y = y - yMin;
    x = [0:length(y)-1]' / fs;
    
    % tPrime
    xPrime = (x+eps)/((xMax - 1)/fs);
 
    % initial guess
    alpha0 = 5;

    % gamma formulation (first pass)
    yP = @(tPrime,a) tPrime .^ a .* exp(a .* (1 - tPrime));
    % Integral of gamma - recirculation
    yR = @(tPrime,a) cumsum(yP(tPrime,a)) / sum(yP(tPrime,a));

    % setting up X and C matrices
    X = @(tPrime,a) [yR(tPrime,a) - yP(tPrime,a), yP(tPrime,a)];
    C = @(tPrime,a,y) X(tPrime,a)\y;
    
    % signal Y = X * C (after initial alpha guess)
    yT = @(tPrime,a,y) X(tPrime,a) * C(tPrime,a,y);
    
    % cost function
    nonlincost = @(a) yT(xPrime,a,y) - y;
    opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','Display','off');

    % nonlinear optimization to get alpha
    tic
    [alpha, ~, res] = lsqnonlin(nonlincost, alpha0, [], [], opts);
    TC = toc;
    NE = rms((res .^ 2) / length(y));

    Cmat = C(xPrime,alpha,y); 
    yMax = Cmat(2);
    k = Cmat(1) / Cmat(2);
    
    result.yP = yMax * (1 - k) * yP(xPrime,alpha);
    result.yR = yMax * k * yR(xPrime,alpha);
    result.A = yMax;
    result.k = k;
    result.beta = x(xMax) / alpha;
    result.alpha = alpha - 1;
    result.t = x;
    result.t0 = (xStart - 1) / fs;
    result.ttp = x(xMax);
    result.tpeak = result.t0 + result.ttp;
    result.b = yMin;
    result.NE = NE;
    result.TC = TC;
    
    if disp
        figure; 
        hold on;
        P = @(tPrime, a) yMax * (1 - k) * yP(tPrime,a) + yMin;
        R = @(tPrime, a) yMax * k * yR(tPrime,a) + yMin;
        plot(x, y  + yMin, 'LineWidth', 2);
        plot(x, yT(xPrime,alpha,y)  + yMin, 'LineWidth', 2);
        plot(x, P(xPrime, alpha), 'LineWidth', 2);
        plot(x, R(xPrime, alpha), 'LineWidth', 2);
        plot(x, y - yT(xPrime,alpha,y), 'LineWidth', 2);
        scatter(x(xMax), y(xMax)  + yMin, 'filled');
        hold off;
        legend('Data', 'Total model', 'First-pass', 'Recirculation', 'Residual', 'Peak');
        title(strjoin(["Hybrid fitting, RMSE = " NE],""));
        xlabel("Time [s]"); ylabel("Conductivity change [-]");
        grid on;
        xlim([min(x) max(x)]);
    end
end