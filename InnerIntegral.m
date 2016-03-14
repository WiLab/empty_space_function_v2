function integral_result = InnerIntegral(r, K_0, N, X_full, num_integrals  )

indexes = 1:N;
lambda = 0.4492;
alpha = 1.558;

% SobolPoints
P = sobolset(2); % set dimension
scramble(P,'MatousekAffineOwen');

% Integrals summation
integral_result = 0;
for n_int=1:N
    % Select ([0,1]x[0,1])^n pairs for n dimensional unit squares
    %X = X_full(1:num_integrals,:);
    % Shuffle sampling values
    indexes_shuffled = indexes(randperm(N));
    X = X_full(indexes_shuffled(1:num_integrals),:);
    
    % Product
    a = bsxfun(@minus,X,[0.5,0.5]);
    norma = sqrt(a(:,1).^2 + a(:,2).^2);
    Product = prod(double(norma<=1/2));
    
    % Calculate determinenant if product nonzero
    if (Product>0)
        %disp('Nonzero Product');
        K = PairWiseDifferences(X,2*r);
        %K_0_val = K_0(K);
        K_0_val = lambda.*exp( -(K./alpha).^2);
        integral_result = integral_result + det(K_0_val);
    end
end

integral_result = 1/N*integral_result;

end