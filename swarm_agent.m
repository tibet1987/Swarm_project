classdef swarm_agent < handle % is a subclass of the 'Handle' class --> objects can be created
    % This class creates agents that live in a swarm
    % The agents follow three basic rules:
    %  1. Move towards the midpoint of those agents you see in your surrounding
    %  2. Move away if someone get too close to you
    %  3. Move approximately in the same direction as your neighbors
    % (C) 2015 Timo Benzel
    
    properties (Access = private)
        % dynamic values (change during simulation)
        %   these values are accessible via 'get...' functions
        pos    % current position vector of the agent
        vel    % current velocity vector of the agent
        viewable_neighbor % list of neighbors in viewable distance
        
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
        
            % default settings for agent
            ID = 1;
            view_dist = 0.5;
            m = 0.3;
            k_dist = 10;
            d_dist = 1;
            viewable_neighbor = [];
            
            % checking varargin
            if mod(nargin,2) ~= 0
                error('class ''swarm_agent'' constructor: Initial parameters need to be in pairs of the form "...,''variableName'',variableValue"');
            end
            for i=1:numel(varargin)/2
                if     strcmpi(varargin{2*i-1},'ID')
                    obj.ID = varargin{2*i};
                elseif strcmpi(varargin{2*i-1},'view_dist')
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
            if ~isequal(size(pos),[3,1])
                if ~isequal(size(pos),[1,3])
                    error('swarm_agent constructor: wrong dimensions for ''pos''');
                end
            end
            obj.pos = pos; 
            
            % initializing velocity vector
            if ~isequal(size(vel),[3,1])
                if ~isequal(size(vel),[1,3])
                    error('swarm_agent constructor: wrong dimensions for ''vel''');
                end
            end
            obj.vel = vel; 
            
            % looking for neighbors that are already around
            obj.updateNeighborList();
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
        
        function ID = getNeighborIDs( obj ) 
        % return velocity vector of this swarm agent
            ID = obj.viewable_neighbor;
        end
        
        function appendNeighbor( obj, neighbor_ID )
        % add a new neighbor ID to agents list of viewable neighbors
            obj.viewable_neighbor(end+1) = neighbor_ID;
        % sorting entries in ascending order
            obj.viewable_neighbor = sort(obj.viewable_neighbor);
        end
        
        function found = removeNeighbor( obj, neighbor_ID )
        % remove a neighbor ID from agents list of viewable neighbors
            idx = find(obj.viewable_neighbor == neighbor_ID);
            if isempty(idx)
                found = 0;
                disp(['Method ''removeNeighbor'' of object in class',...
                    '''swarm_agent'': Neighbor not found, could not ',...
                    'remove neighbor']);
            else
                found = 1;
                obj.viewable_neighbor(idx) = [];
                % sorting entries in ascending order
                obj.viewable_neighbor = sort(obj.viewable_neighbor);
            end
        end
        
        function dist = calcDistance(obj, neighbor)
        % distance = scalar product of vector pointing from obj to neighbor
            pointer_to_neighbor = obj.getPos() - neighbor.getPos();
            dist = sqrt(  dot(pointer_to_neighbor,pointer_to_neighbor)  ); 
        end
        
        function force = calcNeighborDistForce(obj, neighbor)

            dist = obj.calc_distance(neighbor);
            
            force = obj.k_dist * dist; % Hooke's Law: F = k*x
            
%         % Gaussian approach (3D-Gaussian) to rate 
%             c = 0; % center is always zero
%             sigma_norm = sqrt(1/(2*log(100))); % normalizing sigma to gaussmf(1) = 0.01;
%             sig = sigma_norm * obj.view_dist;
%             rating = gaussmf(dist,[sig c]);
        end
        
        function updateNeighborList(obj)
            global agent_list
            % looking through agent list if someone is in viewable range
            for i=1:numel(agent_list)
                dist = obj.calcDistance(agent_list(i).handle);
                if dist < obj.view_dist
                    obj.appendNeighbor(agent_list(i).handle);
                end
            end
        end
            
        function obj = update( obj )
            % looking for new neighbors nearby
            
            
            global step_size
            for i=1:numel(obj)
                % calculating forces of the 3 laws
                neighbor_midpoint_force =  zeros(3,1);
                neighbor_dist_force = obj.calcNeighborDistForce();
                neighbor_dir_force =  zeros(3,1);
                urge_force = zeros(3,1);
                
                % vector sum of all forces pulling/pushing on the agent
                forces = neighbor_midpoint_force + neighbor_dist_force...
                            + neighbor_dir_force + urge_force;
                        
                % solving dynamics with euler forward
                x_dot = dynamics(obj(i).pos,obj(i).vel,forces,obj(i).m);
                obj(i).pos = obj(i).pos + step_size*x_dot(1);
                obj(i).vel = obj(i).vel + step_size*x_dot(2);
                
            end
        end

    end
    %%%%%%%%%%%%%%%%
    methods(Static)
        
        function x_dot = dynamics(pos,vel,forces,m)
            % this function contains the dynamic equations of the system
            % 'x' and 'u' are 2-dimensional column vectors
            x = [pos; vel];
            u = forces;
            A = [0, 1;
                 0, 0];
            B = [0; 1/m];
            x_dot = A*x + B*u;
        end
        
        function plot_swarm( obj )
            global space_lims 
            for i=1:numel(obj)
%                 scatter3(obj(i).pos(1),obj(i).pos(2),obj(i).pos(3),'*','MarkerSize',20)
                scatter3(obj(i).pos(1),obj(i).pos(2),obj(i).pos(3),1000,'.')
                hold all
            end
            hold off
            xlim([-space_lims(1),space_lims(1)])
            ylim([-space_lims(2),space_lims(2)])
            zlim([-space_lims(3),space_lims(3)])
            view(0,90); % view from above
            
            set(gcf,'position',[10,40,1600,900]);
        end
        
    end
end

