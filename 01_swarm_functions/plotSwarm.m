function plotSwarm()
    global space_lims agent_list step_size
    
    for i=1:numel(agent_list)
%                 scatter3(obj(i).pos(1),obj(i).pos(2),obj(i).pos(3),'*','MarkerSize',20)
        position = agent_list(i).handle.getPos();
        scatter3(position(1),position(2),position(3),1000,'.')
        if i==1
            xlim([-space_lims(1),space_lims(1)])
            ylim([-space_lims(2),space_lims(2)])
            zlim([-space_lims(3),space_lims(3)])
            view(0,90); % view from above
%             view(45,45); % view from the side

%             set(gcf,'position',[10,40,1600,900]);
            set(gcf,'position',[800,300,800,600]);
        end
%         hold on
%         velocity = agent_list(i).handle.getVel();
%         quiver(position(1),position(2),velocity(1)*step_size,velocity(2)*step_size,50,'linewidth',4,'color','c')
        hold all
    end
    hold off

    
end
