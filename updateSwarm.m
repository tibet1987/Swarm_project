function updateSwarm()
    global agent_list
    % looking for new neighbors nearby
    for i=1:numel(agent_list)
        agent_list(i).handle.updateNeighborList();
    end

    % updating states of movement dynamics
    global step_size
    for i=1:numel(agent_list)
        
        forces = agent_list(i).handle.calcAllForces(); % returns a force vector
        % update dynamic (mechanical) states
        agent_list(i).handle.updateDynamics(forces);

    end
end