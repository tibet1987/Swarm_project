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
        viewable_neighbor_ID % list of neighbors in viewable distance
        
        % parameters (static during simulation)
        ID           % unique ID number
        view_dist    % [m] agents view horizon
        m            % [kg] mass
        k_dist % [N/m] stiffness of virtual spring-damper to neighbor
        d_dist % [Ns/m]damping of virtual spring-damper to neighbor
        
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
            obj.view_dist = 0.5;
            obj.m = 0.3;
            obj.k_dist = 1;
            obj.d_dist = 0.08;
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
                elseif strcmpi(varargin{2*i-1},'mass')
                    obj.m = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'k_dist')
                    obj.k_dist = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'d_dist')
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
        % return velocity vector of this swarm agent
            ID = obj.ID;
        end
        
        function m = getMass( obj ) 
        % return velocity vector of this swarm agent
            m = obj.m;
        end
        
        function ID = getNeighborIDs( obj ) 
        % return velocity vector of this swarm agent
            ID = obj.viewable_neighbor_ID;
        end
        
        function appendNeighbor( obj, neighborID )
        % add a new neighbor ID to agents list of viewable neighbors
            obj.viewable_neighbor_ID(end+1) = neighborID;
        % sorting neighbor IDs in ascending order
            obj.viewable_neighbor_ID = sort(obj.viewable_neighbor_ID);
        end
        
        function found = removeNeighbor( obj, neighbor_ID )
        % remove a neighbor ID from agents list of viewable neighbors
            idx = find(obj.viewable_neighbor_ID == neighbor_ID);
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
        
        function updateNeighborList( obj )
            global agent_list

            % looking through agent list if someone is in viewable range
            for i=1:numel(agent_list)

                % how far away is this agent
                dist = obj.calcDistance(agent_list(i).handle);
                
                % is this agent already in neighbor list?
                alreadyInList = ~isempty( find(obj.getNeighborIDs() == agent_list(i).handle.getID(),1,'first') );
                object_is_not_myself = obj.ID ~= agent_list(i).handle.getID();
                
                if dist <= obj.view_dist  && object_is_not_myself && ~alreadyInList
                    obj.appendNeighbor( agent_list(i).handle.getID() );
                elseif dist > obj.view_dist  &&  alreadyInList
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
            
            % Law 1: midpoint force
            neighbor_midpoint_force =  zeros(3,1); % init
            if active_swarm_laws(1) == 0
                neighbor_midpoint_force = zeros(3,1);
            end
            
            % Law 2: neighbor distance force
            neighbor_dist_force = zeros(3,1); % init
            neighbor_IDs = obj.viewable_neighbor_ID;

            for i=1:numel(neighbor_IDs)
                neighbor = getHandleOfAgent(neighbor_IDs(i));
                
                pos_diff_vector = neighbor.getPos() - obj.pos;
                abs_spring_deflection = obj.view_dist-norm(pos_diff_vector);
                deflection_vector = abs_spring_deflection * pos_diff_vector/norm(pos_diff_vector);
                
                velocity_difference = neighbor.getVel() - obj.vel;
                
                neighbor_dist_force = neighbor_dist_force - obj.k_dist * deflection_vector + obj.d_dist * velocity_difference; % Hooke's Law with damping: F = -k*x
            end
            
            % ... this also counts for objects that are in the way
            
            % and for the borders of the defined space

            if active_swarm_laws(2) == 0
                neighbor_dist_force = zeros(3,1);
            end
            
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

                % solving dynamics with euler forward
                obj.pos(coord) = obj.pos(coord) + step_size*x_dot(1);
                obj.vel(coord) = obj.vel(coord) + step_size*x_dot(2);
            end
        end
        
    end
    %%%%%%%%%%%%%%%%
    methods(Static)

    end
end

