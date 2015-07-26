function updateSwarm()
    global agent_list
    % looking for new neighbors nearby
    for i=1:numel(agent_list)
        agent_list(i).handle.updateNeighborList();
    end

    % updating states of movement dynamics
    global step_size
    for i=1:numel(agent_list)
        % IMPORTANT: 
        %   FIRST:  Calculate respecive forces based on current states
        %   SECOND: Update dynamic states based on those forces
        forces(:,i) = agent_list(i).handle.calcAllForces(); % returns a force vector
    end
    
    for i=1:numel(agent_list)
    % update dynamic (mechanical) states
        agent_list(i).handle.updateDynamics(forces(:,i));
    end
end