clear;clc;close all;
%% Initialization

% Define start and stop times, set an h to keep outputs at constant time
% interval. You should only change h for different requirements.
t = 0;
tfinal = 3;
% h = 0.04;%  seconds
% h = 0.02; % seconds
h = 0.01;% seconds

tspan = t:h:tfinal;

% Initialize state and contact mode
q = [0.2;1];
dq = [0;0];

% Initialize arrays for logging data
xout = []; % collection of state vectors at each time
lambdaout = []; % collection of contact forces at each time


% Initialize contact mode
oldContactMode = zeros(0,1);
ContactMode = zeros(0,1);
disp(['Initialize in mode {', num2str(oldContactMode'), '}.']);
%% Main Loop

for i = 1:length(tspan)
    % Log state and contact forces
    % Your code here
    [q_next,dq_next,lambda] = eom_solver(q,dq,h);

    % Check new contact mode and determine if there has been a transition
    % Your code here, and display the following message when appropriate
    
    % Reset data
    % Your code here
    var = -lambda>0;
    if var == [0;0]
        ContactMode = 0;
    elseif var == [1;0]
        ContactMode = 1;
    elseif var == [0;1]
        ContactMode = 2;
    elseif var == [1;1]
        ContactMode = 12;
    end
    
    if oldContactMode ~= ContactMode
            disp(['Transition from mode {', num2str(oldContactMode), '} to mode {', num2str(ContactMode), '} at t = ', num2str(tspan(i)), '.']);
    end
    oldContactMode = ContactMode;
    xout = [xout;q'];
    q = q_next;
    dq = dq_next;
end
% disp(['Terminate in mode {', num2str(), '} at t = ', num2str(), '.']);

% This function shows animation for the simulation, don't modify it
animateHW11(xout, h);