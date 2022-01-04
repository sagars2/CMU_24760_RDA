%% The function should return w, t, p, g_w, g_w_matlab, g_v, g_v_matlab, g and g_matlab.

function [w, t, p, g_w, g_w_matlab, g_v, g_v_matlab, g, g_matlab] = compare_tform()
    % TODO: generate a random unit vector w, rotation amount t, and the
    % twist pitch p
    t = zeros();
    p = zeros();
    v = zeros(3,1)
    w = rand(3,1)
    w = w/norm(w);
    t = rand()
    p = rand()
    % TODO: using w and t, compute the 4x4 rigid body transformation matrix 
    % g_w_matlab generated with axang2tform, and g_w generated with your 
    % function twist2rbvel and the Matlab function expm
    xi = [v; w]
    g_w_matlab = zeros(4, 4);
    g_w = zeros(4, 4);
    axang = [transpose(w) t];
    g_w_matlab = axang2tform(axang)
    g_w = expm(twist2rbvel(xi)*t)
    % TODO: compute the velocity of 3-vector v = wp
    v = zeros(3, 1);
    v = w*p

    % TODO: using the pure translation with a velocity of v and the amount 
    % t, compute the 4x4 rigid body transformation matrix g_v_matlab 
    % generated with trvec2tform, and g_v generated with your function 
    % twist2rbvel and the Matlab function expm
    g_v_matlab = zeros(4, 4);
    g_v = zeros(4, 4);
    v1 = transpose(v);
    g_v_matlab = trvec2tform(v1*t);
    xi = [v;0;0;0];
    g_v = expm(twist2rbvel(xi)*t);
    
    % TODO: using w, t, and v, compute the 4x4 rigid body transformation 
    % matrix g generated with your function twist2rbvel and the Matlab 
    % function expm and g_matlab generated with the composition of 
    % axang2tform and trvec2tform
    g_matlab = zeros(4, 4);
    g = zeros(4, 4);
    
    l = [g_w_matlab(1,1) g_w_matlab(1,2) g_w_matlab(1,3); g_w_matlab(2,1) g_w_matlab(2,2) g_w_matlab(2,3); g_w_matlab(3,1) g_w_matlab(3,2) g_w_matlab(3,3); 0 0 0]
    m = [g_v_matlab(:,4)]
    
     g_matlab = [l m];

     s = [g_w(1,1) g_w(1,2) g_w(1,3); g_w(2,1) g_w(2,2) g_w(2,3); g_w(3,1) g_w(3,2) g_w(3,3); 0 0 0]
     e = g_v(:,4);
     g = [s e]; 
end