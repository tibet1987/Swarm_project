% Modelling of a swarm based on agents

% Three basic laws for every agent are:
%  1. Move towards the midpoint of those agents you see in your surrounding
%  2. Move away if someone get too close to you
%  3. Move approximately in the same direction as your neighbors

% NOTE:
%  - every agent has a position [x,y,z] and a direction and a certain
%       velocity [vx,vy,vz]

clear classes;
% clear all;
close all;
clc;

%% initialize parameters
global space_lims       % size of space where agents move in
space_lims = 100*[1,1,1]; %[x,y,z]
global step_size
step_size = 0.02;
num_agents = 3;     % number of agents


% agent parameters
view_dist = 0.5;  % [m] how far can the agent see
mass = 0.3;  % [kg], mass of agent
k_neighbor_dist = 10; % stiffness of virtual spring between too close neighbors and agent
d_neighbor_dist = 1; % virtual damping between too close neighbors and agent


%% Initialize test

% creating agents
myagent(1) = swarm_agent(   [0,0,0]',0  *[0,0,0]','view_dist',0.5,...
                'mass',mass,'k_neighbor_dist',k_neighbor_dist,...
                'd_neighbor_dist',d_neighbor_dist);
%%%%%%%%%%%%%%%
myagent(2) = swarm_agent(0.8*[1,0,0]',-0.1*[1,0,0]','view_dist',0.5,...
                'mass',mass,'k_neighbor_dist',k_neighbor_dist,...
                'd_neighbor_dist',d_neighbor_dist);
%%%%%%%%%%%%%%%
myagent(3) = swarm_agent(0.2*[0,1,1]',0.1*[1,1,0]','view_dist',0.5,...
                'mass',mass,'k_neighbor_dist',k_neighbor_dist,...
                'd_neighbor_dist',d_neighbor_dist);

            
% testing some agent operations
% agent_distance = myagent(1).calcDistance(myagent(2))
% myagent.getNeighborIDs()


% myagent.updateNeighborList()

% for i=1:50
%     myagent.plotSwarm(myagent);
%     myagent.update();
%     drawnow;
%     pause(0.05)
% end


