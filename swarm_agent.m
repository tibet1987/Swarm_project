classdef swarm_agent < handle % is a subclass of the 'Handle' class --> objects can be created
    % This class creates agents that live in a swarm
    % The agents follow three basic rules:
    %  1. Move towards the midpoint of those agents you see in your surrounding
    %  2. Move away if someone gets too close to you
    %  3. Move approximately in the same direction as your neighbors
    % (C) 2015 Timo Benzel
    
    properties (Access = private)
        % dynamic values (change during simulation)
        %   these values are accessible via 'get...' functions
        pos    % current position vector of the agent
        vel    % current velocity vector of the agent
        viewable_neighbor_ID % viewable_neighbor_ID --> neighbors in viewable distance
        too_close_neighbor_ID % too_close_neighbor_ID --> neighbors that are too close (Law2)
        
        % parameters (static during simulation)
        ID           % unique ID number
        view_dist    % [m] agents view horizon
        too_close_dist % [m] <= view_dist, distance to viewable neighbor that is too close
        m            % [kg] mass
        k_dist % [N/m] stiffness of virtual spring to neighbor
        d_dist % [Ns/m]damping of virtual damper to neighbor
        k_midpnt % [N/m] stiffness of virtual spring to neighbor midpoint
        d_midpnt % [Ns/m]damping of virtual damper to neighbor midpoint
        
        % non-states
        %neighbor_midpoint_force; % vector of force that pulls agent towards midpoint of neighbors (rule 1)
        %neighbor_dist_force;     % vectors of the forces created by too close neighbors (rule 2)
        %neighbor_dir_force;      % vector of force that pulls agent into heading direction of neighbors (rule 3)
    end
    
    %%%%%%%%%%%%%%
    methods
        
        function obj = swarm_agent( pos, vel, varargin )  % CONSTRUCTOR
        % obj = swarm_agent( pos, vel )
        % obj = swarm_agent( pos, vel, 'valueName1', value1, ...)
        %
        % Constructor for class 'swarm_agent'
        % input:
        %  pos          [3 x 1] double - position column vector
        %  vel          [3 x 1] double - veloctiy column vector
        %
        %  optional input "varargin":
        %    valueName       value           description
        %    -------------------------------------------------------------
        %    'ID'            integer         unique number for every agent
        %    'view_dist'     double [m]      agent can see this far in
        %                                    every direction
        %    'm'             double [kg]     mass of agent
        %    'k_dist'        double [N/m]    stiffness of virtual spring 
        %                                    between agent and neighbors
        %    'd_dist'        double [Ns/m]   virtual damping between 
        %                                    neighbors and agent
        % output:
        %  object of class 'swarm_agent'
        %
        % (C) 2015 Timo Benzel
        global agent_list active_swarm_laws
        
            % default settings for agent
            obj.ID = numel(agent_list)+1;
            obj.view_dist = 2;
            obj.too_close_dist = 1;
            obj.m = 0.3;
            obj.k_dist = 1;
            obj.d_dist = 0.08;
            obj.k_midpnt = 0.1;
            obj.d_midpnt = 0.01;
            obj.viewable_neighbor_ID = [];
            
            if isempty(active_swarm_laws)
                % if not defined otherwise: enable all swarm laws
                active_swarm_laws = [1,1,1];  
            end
            
            % checking varargin
            if mod(nargin,2) ~= 0
                error('class ''swarm_agent'' constructor: Initial parameters need to be in pairs of the form "...,''variableName'',variableValue"');
            end
            % reading in additional agent parameters
            for i=1:numel(varargin)/2   
                if strcmpi(varargin{2*i-1},'view_dist')
                    obj.view_dist = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'too_close_dist')
                    obj.too_close_dist = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'mass')
                    obj.m = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'k_dist')
                    obj.k_dist = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'d_dist')
                    obj.d_dist = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'k_midpnt')
                    obj.k_dist = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'d_midpnt')
                    obj.d_dist = varargin{2*i};
                end
            end
            
            % initializing position vector
            if ~isequal(size(pos),[3,1]) && ~isequal(size(pos),[1,3])
                error('swarm_agent constructor: wrong dimensions for ''pos''');
            end
            obj.pos = pos;
            
            % initializing velocity vector
            if ~isequal(size(vel),[3,1]) && ~isequal(size(vel),[1,3])
                error('swarm_agent constructor: wrong dimensions for ''vel''');
            end
            obj.vel = vel; 
            
            % adding new object to agent_list
            if isempty(agent_list)
                agent_list(1).handle = handle(obj); 
            else
                agent_list(end+1).handle = handle(obj);
                % looking for neighbors that are already around
                obj.updateNeighborList();
            end
        end
        
        function pos = getPos( obj )
        % return position vector of this swarm agent
            pos = obj.pos;
        end
        
        function vel = getVel( obj ) 
        % return velocity vector of this swarm agent
            vel = obj.vel;
        end
        
        function ID = getID( obj ) 
        % return unique ID of this swarm agent
            ID = obj.ID;
        end
        
        function m = getMass( obj ) 
        % return mass of this swarm agent
            m = obj.m;
        end
        
        function isNeighbor = checkIfIsNeighbor( obj , neighborID ) 
        % check if a an agent (with ID) is a neighbor
            isNeighbor = any(obj.viewable_neighbor_ID == neighborID);
        end
        
        function appendNeighbor( obj, neighborID )
        % add a new neighbor ID to agents list of viewable neighbors
            obj.viewable_neighbor_ID(1,end+1) = neighborID;
        % sorting neighbor IDs in ascending order
%             obj.viewable_neighbor_ID = sort(obj.viewable_neighbor_ID);
        end
        
        function found = removeNeighbor( obj, neighborID )
        % remove a neighbor ID from agents list of viewable neighbors
            idx = find(obj.viewable_neighbor_ID == neighborID);
            if isempty(idx)
                found = 0;
                disp(['Method ''removeNeighbor'' of object in class',...
                    '''swarm_agent'': Neighbor not found, could not ',...
                    'remove neighbor']);
            else
                found = 1;
                obj.viewable_neighbor_ID(idx) = [];
                % sorting entries in ascending order
                obj.viewable_neighbor_ID = sort(obj.viewable_neighbor_ID);
            end
        end
        
        function addNeighborFromCloseList( obj, neighborID )
        % add a new neighbor ID to agents list of too close neighbors
            obj.too_close_neighbor_ID(1,end+1) = neighborID;
        end
        
        function isTooClose = isNeighborTooClose( obj, neighborID )
        % add a new neighbor ID to agents list of too close neighbors
            isTooClose = any(obj.too_close_neighbor_ID == neighborID);
        end
        
        function found = removeNeighborFromCloseList( obj, neighborID )
        % remove a neighbor ID from agents list of too close neighbors
            idx = find(obj.too_close_neighbor_ID == neighborID);
            if isempty(idx)
                found = 0;
                disp(['Method ''removeNeighbor'' of object in class',...
                    '''swarm_agent'': Neighbor not found, could not ',...
                    'remove neighbor']);
            else
                found = 1;
                obj.too_close_neighbor_ID(idx) = [];
                % sorting entries in ascending order
                obj.too_close_neighbor_ID = sort(obj.too_close_neighbor_ID);
            end
        end
        
        function updateNeighborList( obj )
            global agent_list

            % looking through agent list if someone is in viewable range
            for i=1:numel(agent_list)

                % how far away is this agent
                dist = obj.calcDistance(agent_list(i).handle);
                
                % is this agent already in neighbor list?
                alreadyInNeighborList = obj.checkIfIsNeighbor( agent_list(i).handle.getID() );
                isInTooCloseList = obj.isNeighborTooClose( agent_list(i).handle.getID() ); 
                object_is_not_myself = obj.ID ~= agent_list(i).handle.getID();
                
                if dist <= obj.view_dist  && object_is_not_myself && ~alreadyInNeighborList
                    % other agent is within viewable range
                    % he is not myself and he is not yet in neighbor list 
                    % --> put him on neighbor list
                    obj.appendNeighbor( agent_list(i).handle.getID() );
                    
                elseif dist <= obj.view_dist && object_is_not_myself && alreadyInNeighborList
                    % other agent is within viewable range
                    % he is not myself and 
                    % he is in the neighbor list 
                    %   --> check if agent is too close
                    if dist <= obj.too_close_dist && ~isInTooCloseList
                        % if other agent closer than 'too_close_dist', if he is
                        % not myself and if he is a neighbor --> he is too
                        % close
                        obj.addNeighborFromCloseList(  agent_list(i).handle.getID() );
                    elseif dist > obj.too_close_dist && isInTooCloseList
                        % if other agent not closer than 'too_close_dist' 
                        % anymore, if he is not myself and if he is a neighbor 
                        % --> set him to "not close"
                        obj.removeNeighborFromCloseList(  agent_list(i).handle.getID() );
                    end
                    
                elseif dist > obj.view_dist  &&  alreadyInNeighborList
                    % if other agent is not in viewable range but he is in
                    % the neighbor list --> remove him from the list
                    obj.removeNeighbor( agent_list(i).handle.getID() );
                    
                end

            end

        end
        
        function dist = calcDistance( obj, neighbor )
        % distance = scalar product of vector pointing from obj to neighbor
            pointer_to_neighbor = obj.getPos() - neighbor.getPos();
            dist = sqrt(  dot(pointer_to_neighbor,pointer_to_neighbor)  ); 
        end
        
        function force_vec = calcAllForces( obj )
            % calculating forces of the 3 laws:
            global active_swarm_laws
            global space_lims
            
            % -------------------------------------------
            % Law 1: midpoint force
            neighbor_midpoint_force =  zeros(3,1); % init
            
            if active_swarm_laws(1) == 1 % only evaluate if Law1 is active
                if ~isempty(obj.viewable_neighbor_ID)
                    midpoint_sum = zeros(3,1);
                    for i=1:numel(obj.viewable_neighbor_ID)
                        neighbor = getHandleOfAgent(obj.viewable_neighbor_ID(i));
                        midpoint_sum = midpoint_sum+neighbor.getPos();
                    end
                    midpoint = midpoint_sum/numel(obj.viewable_neighbor_ID);
                    neighbor_midpoint_force = obj.k_midpnt * (midpoint - obj.pos);
                end
            end
            
            % -------------------------------------------
            % Law 2: neighbor distance force
            neighbor_dist_force = zeros(3,1); % init
            
            if active_swarm_laws(2) == 1
                if ~isempty(obj.too_close_neighbor_ID)
                    for i=obj.too_close_neighbor_ID % only those neighbors that are too close
                        % if he is in viewable range and if he is too close -->
                        % take evasive action
                        neighbor = getHandleOfAgent(i);

                        pos_diff_vector = neighbor.getPos() - obj.pos;
                        abs_spring_deflection = obj.too_close_dist - norm(pos_diff_vector);
                        deflection_vector = abs_spring_deflection * pos_diff_vector/norm(pos_diff_vector);

                        velocity_difference = obj.vel - neighbor.getVel();
                        
                        % Hooke's Law with damping: F = -k*x + D*x_dot
                        neighbor_dist_force = neighbor_dist_force - obj.k_dist * deflection_vector - obj.d_dist * velocity_difference;
                    end
%                 else
%                     if obj.ID == 1
%                         obj.vel
%                     end
                end
                
                % ... also move away from the borders of the defined space
                xlim_violated = abs( abs(obj.pos(1)) - space_lims(1) ) < obj.too_close_dist;
                ylim_violated = abs( abs(obj.pos(2)) - space_lims(2) ) < obj.too_close_dist;
                zlim_violated = abs( abs(obj.pos(3)) - space_lims(3) ) < obj.too_close_dist;
                if xlim_violated || ylim_violated || zlim_violated
                    % agent too close to space border
                    viol_vector = [xlim_violated; ylim_violated; zlim_violated];
                    dist_to_space_border = - sign(obj.pos).* (space_lims-obj.too_close_dist - abs(obj.pos)) .* viol_vector;
                    neighbor_dist_force = neighbor_dist_force - obj.k_dist * dist_to_space_border - obj.d_dist * obj.vel;
%                 else
%                     if obj.ID == 1
%                         obj.vel
%                     end
                end
                    
                % and for objects that are in the way
                
            end
            
            

            
            
            % -------------------------------------------
            % Law 3: neighbor direction/heading force
            neighbor_dir_force =  zeros(3,1);
            urge_force = zeros(3,1); % maybe another force for urges --> food, shelter, sleep
            
            if active_swarm_laws(3) == 0
                neighbor_dir_force = zeros(3,1);
            end
            
            % vector sum of all forces pulling/pushing on the agent
            force_vec = neighbor_midpoint_force + neighbor_dist_force...
                        + neighbor_dir_force + urge_force;
                    
%         % Gaussian approach (3D-Gaussian) to rate 
%             c = 0; % center is always zero
%             sigma_norm = sqrt(1/(2*log(100))); % normalizing sigma to gaussmf(1) = 0.01;
%             sig = sigma_norm * obj.view_dist;
%             rating = gaussmf(dist,[sig c]);
        end
        
        function updateDynamics( obj, force )
            global step_size
            % evaluating the dynamic equations
            % 'x' and 'u' are 2-dimensional column vectors
            for coord=1:3 % iterate over x,y,z direction
                x = [obj.pos(coord) ; obj.vel(coord)];
    %                [x,y,z , dx,dy,dz]
                u = force(coord);
                A = [0, 1;
                     0, 0];
                B = [0; 1/obj.getMass()];
                x_dot = A*x + B*u;

                % solving dynamics with euler forward (explicit)
                % x_(k+1) = x(k) + h * f( x(k) , t(k) ) 
                obj.pos(coord) = obj.pos(coord) + step_size*x_dot(1);
                obj.vel(coord) = obj.vel(coord) + step_size*x_dot(2);
            end
        end
        
    end
    %%%%%%%%%%%%%%%%
    methods(Static)

    end
end

