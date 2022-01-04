%% angvel2skew(w) maps the 3-vector w to the 3x3 skew-symmetric matrix w_hat
function w_hat=angvel2skew(w)
    % TODO: construct the 3x3 skew-symmetric matrix w_hat from w
    w_hat = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0]
end


