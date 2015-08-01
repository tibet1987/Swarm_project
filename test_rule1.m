% Modelling of a swarm based on agents

% Three basic laws for every agent are:
%  1. Move towards the midpoint of those agents you see in your surrounding
%  2. Move away if someone get too close to you
%  3. Move approximately in the same direction as your neighbors

% NOTE:
%  - every agent has a position [x,y,z] and a direction and a certain
%       velocity [vx,vy,vz]

% clear all;
close all;
clc;


%% initialize parameters
global space_lims       % size of space where agents move in
space_lims = 5*[1,1,1]; %[x,y,z]
global step_size
step_size = 0.02;


% agent parameters
view_dist = 3;  
too_close_dist = 1; 
mass = 0.3;  % [kg], mass of agent
k_dist = 5; % stiffness of virtual spring between too close neighbors and agent
d_dist = 0.5; % virtual damping between too close neighbors and agent


%% Initialize test

% creating agents
myagent(1) = swarm_agent(-2*[1,0,0]',2  *[1,0,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,'too_close_dist',too_close_dist);
%%%%%%%%%%%%%%%
myagent(2) = swarm_agent(2*[1,0,0]',-2*[1,0,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,'too_close_dist',too_close_dist);
%%%%%%%%%%%%%%%
% myagent(3) = swarm_agent(1*[0,1,0]',1*[0,1,0]','view_dist',view_dist,...
%                 'mass',mass,'k_dist',k_dist,'d_dist',d_dist,'too_close_dist',too_close_dist);
            
% testing some agent operations
% agent_distance = myagent(1).calcDistance(myagent(2))
% myagent.getNeighborIDs()


% myagent.updateNeighborList()
global agent_list 

for i=1:300
    plotSwarm();
    updateSwarm();
    
    
    [agents_see_each_other,agents_are_too_close] = checkIfAgentsSeeEachOther(agent_list(1).handle,agent_list(2).handle);
    
    if agents_see_each_other
        if agents_are_too_close
            disp('Agent1 & Agent2 see each other and are too close')
        else
            disp('Agent1 & Agent2 see each other')
        end
    else
        disp('Agent1 & Agent2 do not see each other')
    end
    
    drawnow;
    pause(0.01)
end


