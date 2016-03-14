%% Integrate F(r)

lambda = 0.4492;
alpha = 1.558;
% Gauss DPP
K_0 = @(x) lambda.*exp(-(x./alpha).^2);
N = 2^9; % Num SobolPoints
P = sobolset(2);
X_full = net(P,N);
%% Integral
% Generate points
r = 0.1:0.1:2;
F = zeros(length(r),1);
indexes = 1:N;

for i =1:length(r)
    
    for n=1:N
        c = (-1)^(n-1)*(2*r(i))^(2*n)/factorial(n);
        if isnan(c)
            c = 0;
            continue;
        end
        % Generate [0,1]^n pairs for n dimensional unit squares
        %X = net(P,n_int);
        indexes_shuffled = indexes(randperm(N)); % Shuffle indexes
        X = X_full(indexes_shuffled(1:5),:);
        
        % Product
        a = bsxfun(@minus,X,[0.5,0.5]);
        norma = sqrt(a(:,1).^2 + a(:,2).^2);
        Product = prod(double(norma<=1/2));
        
        % Calculate K_0
        K = PairWiseDifferences(X,2*r(i));
        K_0_val = K_0(K);
        
        % Sum integral
        integral_result = det(K_0_val)*Product;
        F(i) = F(i) + 1/n*c*integral_result;
        fprintf('n = %d\n',n);
    end
end

plot(r,F);

