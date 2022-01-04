%% The function should return V, g, V_s_Ad_g, V_s_tform, V_b_Ad_g and V_b_tform.

function [V, g, V_s_Ad_g, V_s_tform, V_b_Ad_g, V_b_tform] = compare_twist()
    % TODO: generate a random unit vector w, rotation amount t, and the
    % twist pitch p
    w = zeros(3, 1);
    w = rand(3,1);
    w = w/norm(w);
    t = zeros();
    t = rand()
    p = zeros();
    p = rand();
    % TODO: compute the velocity of 3-vector v = wp
    v = zeros(3, 1);
    v = w*p
    % TODO: using w, t, and v, construct a random twist V and compute the 
    % 4x4 rigid body transformation matrix g
    V = zeros(6, 1);
    g = zeros(4, 4);
    V = [v;w]
    g = expm(twist2rbvel(V)*t)
    % TODO: first, treating V as a body velocity, compute the conversion to 
    % spatial velocity V_s_Ad_g using tform2adjoint and the conversion in 
    % homogeneous coordinates V_s_tform using twist2rbvel and rbvel2twist
    Ad_g = zeros(6, 6);
    V_s_Ad_g = zeros(6, 1);
    V_s_tform = zeros(6, 1);
    Ad_g = tform2adjoint(g);
    V_s_Ad_g = Ad_g*V;
    V_s_tform = rbvel2twist(twist2rbvel(V_s_Ad_g))

    % TODO: repeat the test the other way, treating V as a spatial velocity 
    % and converting to body velocity and compute V_b_Ad_g and V_b_tform
    V_b_Ad_g = zeros(6, 1);
    V_b_tform = zeros(6, 1);
    Ad_g = tform2adjoint(g);
    V_b_Ad_g = inv(Ad_g)*V;
    V_b_tform = rbvel2twist(twist2rbvel(V_b_Ad_g))

end