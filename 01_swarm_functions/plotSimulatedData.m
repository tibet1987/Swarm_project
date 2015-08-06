function plotSimulatedData(recorder)
    global space_lims step_size
    if isempty(recorder)
        error('Function ''plotSimulatedData'': There is no recorded data available');
    end
    
    num_sim_steps = size(recorder(1).agentPos,2);
    n = 1;
    %% Test run to determine plot duration
	for i=1:5 % five test samples 
        tic;
%                 scatter3(obj(i).pos(1),obj(i).pos(2),obj(i).pos(3),'*','MarkerSize',20)
        for k = 1:numel(recorder) % number of active agents in swarm
            position = recorder(k).agentPos(:,i);
            scatter3(position(1),position(2),position(3),1000,'.')
            hold all
        end
        
        xlim([-space_lims(1),space_lims(1)])
        ylim([-space_lims(2),space_lims(2)])
        zlim([-space_lims(3),space_lims(3)])
        view(0,90); % view from above
%             view(45,45); % view from the side
%             set(gcf,'position',[10,40,1600,900]);
        set(gcf,'position',[800,300,800,600]);
        drawnow;
%         hold on
%         velocity = agent_list(i).handle.getVel();
%         quiver(position(1),position(2),velocity(1)*step_size,velocity(2)*step_size,50,'linewidth',4,'color','c')
        hold off
        plotDuration(n) = toc;
        n = n+1;
    end
    % not considering the first time sample because the first sample takes a little longer
    meanPlotDuration = mean(plotDuration(2:end)); 
    
    
    %% Main run
    
    % the simulation time stamp difference has to be larger than 
    % meanPlotDuration in order for the plot to look "real-time"
    downsample = ceil(meanPlotDuration/step_size);
    n = 1;
    
	for i=1:downsample:num_sim_steps
%                 scatter3(obj(i).pos(1),obj(i).pos(2),obj(i).pos(3),'*','MarkerSize',20)
        for k = 1:numel(recorder) % number of active agents in swarm
            position = recorder(k).agentPos(:,i);
            scatter3(position(1),position(2),position(3),1000,'.')
            hold all
        end
        
        xlim([-space_lims(1),space_lims(1)])
        ylim([-space_lims(2),space_lims(2)])
        zlim([-space_lims(3),space_lims(3)])
        view(0,90); % view from above
%             view(45,45); % view from the side
%             set(gcf,'position',[10,40,1600,900]);
        set(gcf,'position',[800,300,800,600]);
        drawnow;
%         hold on
%         velocity = agent_list(i).handle.getVel();
%         quiver(position(1),position(2),velocity(1)*step_size,velocity(2)*step_size,50,'linewidth',4,'color','c')
        hold off
        plotDuration(n) = toc;
	end
    meanPlotDuration = mean(plotDuration);
end
