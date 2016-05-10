function alpha = IIA( y, K, theta_center, theta0, R, eta, epsilon )
%IIA	Cortes Interpolated Iterative Algorithm
%
% Input:
%		K               3-way tensor of kernel matrix
%		theta0          initial theta
%		R               a factor control theta changes
%		eta             parameter controls alpha changes
%		epsilon         stop criterion
%       theta_center    the center of ball where theta concentrate
% Output:
%		alpha

MAX_ITER = 1000;
theta = NaN(MAX_ITER,1);
alpha = NaN(MAX_ITER,1);

for iter = 1 : MAX_ITER
    % find v
    for p = 1:P
        v(p) = alpha(p)'*K(:,:,p)*alpha(p);
    end
    
    % find theta
    theta(iter) = theta_center + R * v / norm(v);
    
    % update alpha
    alpha(iter) = eta*alpha(iter-1) + (1-eta)* (K + mu*eye(N)) \ y;
    
    % check stop criterion
    if norm( alpha(iter) - alpha(iter-1) ) < epsilon
        break;
    end
end


end

