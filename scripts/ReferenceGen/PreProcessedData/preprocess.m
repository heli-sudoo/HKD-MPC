clear all;

% Gait ="RunJump/";
% Gait = "MixedHopping/";
% Gait = "RunJump_ICRA23/";
Gait = "Trot/";

body_states = readmatrix(Gait + "body_state.csv");
contacts = readmatrix(Gait + "contact.csv");
foot_placements = readmatrix(Gait + "ee_pos.csv");
qJs = readmatrix(Gait + "jnt.csv");
t = readmatrix(Gait + "time.csv", "Delimiter",",");