% Modelling of a swarm based on agents

% Three basic laws for every agent are:
%  1. Move towards the midpoint of those agents you see in your surrounding
%  2. Move away if someone get too close to you
%  3. Move approximately in the same direction as your neighbors

% NOTE:
%  - every agent is a infinitesimally small point
%  - every agent has a position [x,y,z] and a direction with a certain
%       velocity [vx,vy,vz]

clear all, close all;
clc;

%% setting paths
addpath 01_swarm_functions

%% First test --> Law 2: "Move away if someone get too close to you"
% test_rule2;

%% Second test --> Adding Law 1: "Move towards the midpoint of those agents
%%                                you see in your surrounding"
test_rule1;