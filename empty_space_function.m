%% Integrate F(r)

lambda = 0.4492;
alpha = 1.558;
% Gauss DPP
K_0 = @(x) lambda*exp(-(x/alpha)^2);

% SobolPoints
M = 2^11; 
P = sobolset(2); % set dimension
scramble(P,'MatousekAffineOwen');

%% Integral
% Generate points
r = 0.1:0.1:2;
F = zeros(length(r),1);

num_integrals = 20;

parfor i =1:length(r)
    X_full = net(P,N);
    for n=1:num_integrals
        % Calculate front constant
        c = (-1)^(n-1)*(2*r(i))^(2*n)/factorial(n);
        if isnan(c)
            c = 0; % factorial is larger than 4^(2n)
            fprint('c=0\n');
            continue;
        end
        % Integrals summation
        integral_result = InnerIntegral(r(i), K_0, M, X_full, n);
        F(i) = F(i) + c*integral_result;
        fprintf('n = %d\n',n);
    end
end

plot(r,F);

