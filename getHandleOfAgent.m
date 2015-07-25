function obj_handle = getHandleOfAgent( neighborID )
% obj_handle = getHandleOfAgent( neighborID )
% This function returns the handle of the object whose ID is specified with
% 'neighborID'
% Required: An active list of swarm agents (created when the first swarm
%           agent lives
%

    global agent_list
    obj_handle = [];
    for i=1:numel(agent_list)
        if agent_list(i).handle.getID() == neighborID
            obj_handle = agent_list(i).handle;
        end
    end
end