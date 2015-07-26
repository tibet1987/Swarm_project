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
space_lims = 5*[1,1,1]; %[x,y,z]
global step_size
step_size = 0.01;


% agent parameters
view_dist = 1;  % [m] how far can the agent see
mass = 0.3;  % [kg], mass of agent
k_dist = 1; % stiffness of virtual spring between too close neighbors and agent
d_dist = 1; % virtual damping between too close neighbors and agent


%% Initialize test

% creating agents
myagent(1) = swarm_agent(-1*[1,0,0]',2  *[1,0,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist);
%%%%%%%%%%%%%%%
myagent(2) = swarm_agent(1*[1,0,0]',-2*[1,0,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist);
            
% testing some agent operations
% agent_distance = myagent(1).calcDistance(myagent(2))
% myagent.getNeighborIDs()


% myagent.updateNeighborList()

for i=1:200
    plotSwarm();
    updateSwarm();
    drawnow;
%     pause(0.005)
end


