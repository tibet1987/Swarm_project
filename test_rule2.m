% Testing rule:
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
space_lims = 5*[1,1,1]'; %[x,y,z]
global step_size
step_size = 0.001;


% agent parameters
view_dist = 3;  % [m] how far can the agent see
too_close_dist = 1; % if agents get closer than this radius, the agent moves away
mass = 0.3;  % [kg], mass of agent
k_dist = 20; % stiffness of virtual spring between too close neighbors and agent
d_dist = 0.5; % virtual damping between too close neighbors and agent


%% Initialize test
fprintf('Generating agents...')
% creating agents
myagent(1) = swarm_agent(-3*[1,0,0]',5  *[1,0,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,'too_close_dist',too_close_dist);
%%%%%%%%%%%%%%%
myagent(2) = swarm_agent(3*[1,0,0]',-5*[1,0,0]','view_dist',view_dist,...
                'mass',mass,'k_dist',k_dist,'d_dist',d_dist,'too_close_dist',too_close_dist);
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
% 	[agents_see_each_other,agents_are_too_close] = checkIfAgentsSeeEachOther(agent_list(1).handle,agent_list(2).handle);
%     if agents_see_each_other
%         if agents_are_too_close
%             disp('Agent1 & Agent2 see each other and are too close')
%         else
%             disp('Agent1 & Agent2 see each other')
%         end
%     else
%         disp('Agent1 & Agent2 do not see each other')
%     end
	if mod(i,T_end/step_size/10) == 0
        fprintf('.')
	end
end
fprintf(' done!\n')

%% now plotting the results
close all
fprintf('Plotting simulated data swarm\n')
plotSimulatedData(recorder);

