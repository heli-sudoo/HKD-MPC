clear all;

% Gait ="RunJump/";
% Gait = "MixedHopping/";
% Gait = "RunJump_ICRA23/";
% Gait = "Trot/";
name = 'SingleHop_dqJ';
Gait = "MIP_Hopping/"+name+"/";

body_states = readmatrix(Gait + "body_state.csv");
contacts = readmatrix(Gait + "contact.csv");
foot_placements = readmatrix(Gait + "ee_pos.csv");
qJs = readmatrix(Gait + "jnt.csv");
t = readmatrix(Gait + "time.csv", "Delimiter",",");

center_point = readmatrix(Gait+"center_point.csv");
plane_coefficients = readmatrix(Gait+"plane_coefficients.csv");

grfs = readmatrix(Gait+'grfs.csv');

qJds = readmatrix(Gait+'djnt.csv');


save(Gait+"Data.mat",'body_states','contacts','foot_placements','qJs','t','center_point','plane_coefficients','grfs','qJds');