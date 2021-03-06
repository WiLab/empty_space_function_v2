%% F(r) integration

% SobolPoints
N = 2^11; % Number of points to sum over to estimate inner integral
P = sobolset(2); % set dimension
X_full = net(P,N); % Generate points

%% Integrate
r = 0:0.1:2;
%r = (0.1:0.1:2).*1000;

F = zeros(length(r),1);
Fc = zeros(length(r),1);
Fin = zeros(length(r),1);

% The smaller num_integrals, larger integral_result
% num_integrals>12  -> integral_result = 0
% increasing num_integrals provides more stables results for F(r)
num_integrals = 256; % > 256 unstable aka c=NaN

% calculate F(r)
parfor i =1:length(r)
    
    % Outer most summation
    for n=1:num_integrals
        % Calculate front constant
        c = (-1)^(n-1)*(2*r(i))^(2*n)/factorial(n);
        if isnan(c)
            disp('NAN');
        end
        % Integrals summation
        integral_result = InnerIntegral(r(i), N, X_full, n);
        F(i) = F(i) + c*integral_result;
        Fc(i) = Fc(i) + c;
        Fin(i) = Fin(i) + integral_result;
        fprintf('n = %d\n',n);
    end
end

figure(1);plot(r,F);
axis([0 2.5 0 2]);
figure(2);plot(r,Fc);
figure(3);plot(r,Fin);

