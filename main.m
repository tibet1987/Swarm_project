% Modelling of a swarm based on agents

% Three basic laws for every agent are:
%  1. Move towards the midpoint of those agents you see in your surrounding
%  2. Move away if someone get too close to you
%  3. Move approximately in the same direction as your neighbors

% NOTE:
%  - every agent is a infinitesimally small point
%  - every agent has a position [x,y,z] and a direction with a certain
%       velocity [vx,vy,vz]

clear all, close all;
clc;

% general parameters
num_agents = 4;
gen_veloctiy = 1;
global space_lims 
space_lims = 100*[1,1,1]; %[x,y,z]


% agent parameters
agents.view_dist = 5;  % [cm] how far can the agent see
agents.mass = 0.1;  % [kg], mass of agent
k_neighbor_dist = 10; % stiffness of virtual spring between too close neighbors and agent
d_neighbir_dist = 1; % virtual damping between too close neighbors and agent

% testing
myagent(1) = swarm_agent([0,0,0]',gen_veloctiy*[1,1,0]',agents.view_horizon);
myagent(2) = swarm_agent(10*[1,0,0]',gen_veloctiy*[1,1,0]',agents.view_horizon);

myagent(1).calc_distance(myagent(2))
myagent(1).rate_distance_to_neighbor(myagent(2))
clc

for i=1:50
    myagent.plot_swarm(myagent);
    myagent.update();
    drawnow;
    pause(0.05)
end


