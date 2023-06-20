clear all;

% Gait ="RunJump/";
% Gait = "MixedHopping/";
% Gait = "RunJump_ICRA23/";
% Gait = "Trot/";
name = 'SingleHopHold';
Gait = "MIP_Hopping/"+name+"/";

body_states = readmatrix(Gait + "body_state.csv");
body_states(:,1:3) = flip(body_states(:,1:3),2);
body_states(:,7:9) = flip(body_states(:,7:9),2);
contacts = readmatrix(Gait + "contact.csv");
foot_placements = readmatrix(Gait + "ee_pos.csv");
qJs = readmatrix(Gait + "jnt.csv");
t = readmatrix(Gait + "time.csv", "Delimiter",",");

center_point = readmatrix(Gait+"center_point.csv");
plane_coefficients = readmatrix(Gait+"plane_coefficients.csv");

grfs = readmatrix(Gait+'grfs.csv');

qJds = readmatrix(Gait+'djnt.csv');



name = 'TestOrientation';
Gait = "MIP_Hopping/"+name+"/";
length_traj = size(body_states,1);
for i = 1:length_traj
    body_states(i,4:6) = body_states(1,4:6);
    foot_placements(i,:) = foot_placements(1,:);
    qJs(i,:) = qJs(1,:);
    center_point(i,:) = center_point(1,:);
    plane_coefficients(i,:) = plane_coefficients(1,:);
    grfs(i,:) = [0,0,22.05,0,0,22.05,0,0,22.05,0,0,22.05];
end
contacts = ones(size(contacts));
body_states(:,1:3) = 0*body_states(:,1:3);
body_states(:,7:12) = 0*body_states(:,7:12);

turn = pi/8;
for i = 1:length_traj
    for j = 1:7
        if (i < j/7 * length_traj)
            if j == 2
                body_states(i,1) = turn;
            elseif j == 3
                body_states(i,1) = -turn;
            elseif j == 4
                body_states(i,2) = turn;
            elseif j == 5
                body_states(i,2) = -turn;
            elseif j == 6
                body_states(i,3) = turn;
            elseif j == 7
                body_states(i,3) = -turn;
            end
            break
        end
    end
end

save(Gait+"Data.mat",'body_states','contacts','foot_placements','qJs','t','center_point','plane_coefficients','grfs','qJds');