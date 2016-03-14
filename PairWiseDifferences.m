function differencesWithSelf = PairWiseDifferences(P,scaler)

% Input is Nx2
N = size(P,1);
%differences = zeros(N,N-1);
differencesWithSelf = zeros(N,N);

for i = 1:N
    % Get pairwise difference of vectors
    tmp = scaler.*bsxfun(@minus, P , P(i,:));
    % Get norm
    tmpNorm = sqrt(tmp(:,1).^2 + tmp(:,2).^2);
    % Can compare with self
    differencesWithSelf(i,:) = tmpNorm.';
end


end