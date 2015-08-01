function [ agentsSeeEachOther,agents_are_too_close ] = checkIfAgentsSeeEachOther( agent1, agent2 )
%CHECKIFAGENTSSEEEACHOTHER Summary of this function goes here
%   Detailed explanation goes here

agent1_sees_agent2 = agent1.checkIfIsNeighbor(agent2.getID);
agent2_sees_agent1 = agent2.checkIfIsNeighbor(agent1.getID);

agentsSeeEachOther = agent1_sees_agent2 && agent2_sees_agent1;

agent1_is_close_to_agent2 = agent1.isNeighborTooClose(agent2.getID);
agent2_is_close_to_agent1 = agent2.isNeighborTooClose(agent1.getID);

agents_are_too_close = agent1_is_close_to_agent2 && agent2_is_close_to_agent1;


end