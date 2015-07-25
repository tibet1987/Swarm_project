function plotSwarm()
    global space_lims agent_list
    for i=1:numel(agent_list)
%                 scatter3(obj(i).pos(1),obj(i).pos(2),obj(i).pos(3),'*','MarkerSize',20)
        position = agent_list(i).handle.getPos();
        scatter3(position(1),position(2),position(3),1000,'.')
        hold all
    end
    hold off
    xlim([-space_lims(1),space_lims(1)])
    ylim([-space_lims(2),space_lims(2)])
    zlim([-space_lims(3),space_lims(3)])
    view(0,90); % view from above

    set(gcf,'position',[10,40,1600,900]);
end
