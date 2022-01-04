clear;clc;close all;


%block, dynamics,reset map, constraints 


% Define start and stop times, set a dt to keep outputs at constant time
% interval
tstart = 0;
tfinal = 5;
dt = 0.001;

% Initialize state and contact mode
q0 = [1;5];
dq0 = [0;0];
x0 = [q0;dq0];
contactMode = [];
dq_p = [0,0];
%% Main Loop

% Tell ode45 what event function to use and set max step size to make sure
% we don't miss a zero crossing
options = odeset('Events', @(t,x)guardFunctions(t,x,contactMode),'MaxStep',dt);

% Initialize output arrays
x = x0;
tout = [];
xout = x0';
teout = [];
xeout = [];
ieout = [];
ie = [];
t = [];
te = [];
xe = [];

% Main simulation loop
while tstart < tfinal
    options = odeset('Events', @(t,x)guardFunctions(t,x,contactMode),'MaxStep',dt);
    % Initialize simulation time vector
    tspan = [tstart:dt:tfinal];
    
    % Simulate with ode45
    [t,x,te,xe,ie] = ode45(@(t,x)Dynamics(t,x,contactMode),tspan,x0,options);
    te
    contactMode

    % Sometimes the events function will record a nonterminal event if the
    % initial condition is a zero. We want to ignore this, so we will only
    % use the last row in the terminal state, time, and index.
    if ~isempty(ie)
        te = te(end,:);
        xe = xe(end,:);
        ie = ie(end,:);
    end
    
    % Log output
    nt = length(t);
    tout = [tout; t(2:nt)];
    xout = [xout; x(2:nt,:)];
    teout = [teout; te];
    xeout = [xeout; xe];
    ieout = [ieout; ie];
    
    % Quit if simulation completes
    if isempty(ie) 
        disp('Final time reached');
        break; % abort if simulation has completed
    end
    
    % If flag was caused by a_i < 0 (i not in contact mode), compute the
    % proper contact mode via IV complemetarity
    
    
    % Check to see if there should be liftoff (positive lambda), if so
    % compute the proper contact mode via FA complementarity
%     disp(contactMode)  
    contactMode = compIV(xe,contactMode);
     [dq_p,~]=resetMap(xe,contactMode);
     xe(3) = dq_p(1);
     xe(4) = dq_p(2);

     if isequal(contactMode,[3])
           contactMode = compFA(xe,contactMode);
     end

        % Update initial conditions for next iteration
        x0 = xe';
        tstart = t(end);
        
        % Stop if the particle comes to a rest
        if all(abs(dq_p)<1e-6)
            break;
        end

end

% This function shows animation for the simulation, don't modify it
animateHW10(xout, dt);