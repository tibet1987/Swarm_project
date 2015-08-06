% Testing rule:
%  1. Move towards the midpoint of those agents you see in your surrounding
%   and
%  2. Move away if someone get too close to you
%
% Scenario: 
%  - Two agents move towards each other
%  - then the agents see each other (but nothing else happens)
%  - then the agents get too close and start to repell each other via a
%    virtual spring-damper
%  - when the agents get close to the space walls they also get repelled
%    again via a virutal spring-damper

% clear all;
close all;
clc;

%% initialize parameters
global space_lims       % size of space where agents move in
space_lims = 20*[1,1,1]'; %[x,y,z]
global step_size
step_size = 0.001;


% agent parameters
view_dist = 8;  % [m] how far can the agent see
too_close_dist = 2; % if agents get closer than this radius, the agent moves away
mass = 0.3;  % [kg], mass of agent
k_dist = 20; % stiffness of virtual spring between too close neighbors and agent
d_dist = 0.05; % virtual damping between too close neighbors and agent
k_midpnt = 10;
d_midpnt = 0.5;


%% Initialize test
fprintf('Generating agents...')
% creating agents
myagent(1) = swarm_agent(-3*[1,1,0]',5*[1,0.95,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,...
                'k_midpnt',k_midpnt,'d_midpnt',d_midpnt,...
                'too_close_dist',too_close_dist);
%%%%%%%%%%%%%%%
myagent(2) = swarm_agent(3*[1,1,0]',-5*[1,0.8,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,...
                'k_midpnt',k_midpnt,'d_midpnt',d_midpnt,...
                'too_close_dist',too_close_dist);
%%%%%%%%%%%%%%%         
myagent(3) = swarm_agent(3*[1,-1,0]',-5*[1,-0.8,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,...
                'k_midpnt',k_midpnt,'d_midpnt',d_midpnt,...
                'too_close_dist',too_close_dist);
            
num_agents = numel(myagent);
fprintf(' done!\n')


%% preparing simulation 
T_end = 30; % simulation end time in seconds
global agent_list
for i=1:num_agents
    recorder(i).agentPos = zeros(3,T_end/step_size);
end

%% simulating
fprintf('Simulating swarm...')
for i=1:T_end/step_size
    updateSwarm();
    for k=1:num_agents
        recorder(k).agentPos(:,i) = agent_list(k).handle.getPos;
    end

	if mod(i,T_end/step_size/10) == 0
        fprintf('.')
	end
end
fprintf(' done!\n')

%% now plotting the results
close all
fprintf('Plotting simulated data swarm... ')
plotSimulatedData(recorder);
fprintf('done!\n')

