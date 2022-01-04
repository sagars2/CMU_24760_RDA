%% The function should return w, t, rotm and rotm_matlab, and rotm and rotm_matlab should be the same.

function [w, t, rotm, rotm_matlab] = compare_rotm()
    % TODO: generate a random unit vector w and rotation amount t
    w = zeros(3, 1);
    t = 0;
    w = rand(3,1)
    w = w/norm(w) % Unit vector
    
    
    % TODO: Compute the 3x3 rotation matrix rotm_matlab generated with 
    % Matlab built-in function axang2rotm 
    rotm_matlab = zeros(3, 3);
    axang = [transpose(w) t]
    rotm_matlab = axang2rotm(axang)
    
    % TODO: Compute the 3x3 rotation matrix rotm generated with generated 
    % with your function angvel2skew and the Matlab function expm
    rotm = zeros(3, 3);
    rotm = expm(angvel2skew(w)*t)
end