%% rbvel2twist(xi_hat) is the inverse mapping of twist2rbvel. 

function xi = rbvel2twist(xi_hat)
    % TODO: construct the 6-vector twist xi from xi_hat
    xi = zeros(6, 1);
    xi = [xi_hat(1,4); xi_hat(2,4); xi_hat(3,4); xi_hat(3,2); xi_hat(1,3); xi_hat(2,1)]
end
