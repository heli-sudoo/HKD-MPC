clear, clc, close all

%% read in reference data
name = 'DoubleJumpHold';
Gait = "../PreProcessedData/MIP_Hopping/"+name+"/";

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

%% read in lcm logs
load('lcmlog_2023_06_24_00.mat')

% find start of trajectory
for i = 1:size(leg_control_command.f_ff,1)
    if (abs(leg_control_command.f_ff(i,1)) > 0)
        start_traj = i;
        break
    end
end

length_traj = 200;%length(body_states(:,4));
time = 0:0.002:length_traj*0.01;
scale_ind = round(0.01/0.002);
%% plot COM trajectories

subplot(3,1,1)
plot(t(1:length_traj),body_states(1:length_traj,4),'red')
hold on
plot(time,state_estimator.p(start_traj:start_traj+scale_ind*length_traj,1),'blue')

subplot(3,1,2)
plot(t(1:length_traj),body_states(1:length_traj,5),'red')
hold on
plot(time,state_estimator.p(start_traj:start_traj+scale_ind*length_traj,2),'blue')

subplot(3,1,3)
plot(t(1:length_traj),body_states(1:length_traj,6),'red')
hold on
plot(time,state_estimator.p(start_traj:start_traj+scale_ind*length_traj,3),'blue')

%% read in lcm logs with 10x cost on tracking
load('lcmlog_2023_06_24_01.mat')

% find start of trajectory
for i = 1:size(leg_control_command.f_ff,1)
    if (abs(leg_control_command.f_ff(i,1)) > 0)
        start_traj = i;
        break
    end
end

length_traj = 200;%length(body_states(:,4));
time = 0:0.002:length_traj*0.01;
scale_ind = round(0.01/0.002);
%% plot COM trajectories
figure(1)

subplot(3,1,1)
plot(time,state_estimator.p(start_traj:start_traj+scale_ind*length_traj,1),'green')

subplot(3,1,2)
plot(time,state_estimator.p(start_traj:start_traj+scale_ind*length_traj,2),'green')

subplot(3,1,3)
plot(time,state_estimator.p(start_traj:start_traj+scale_ind*length_traj,3),'green')

legend('reference','1x','10x')

%% plot GRFs
figure(2)

f_ff_w = zeros(size(leg_control_command.f_ff(start_traj:start_traj+scale_ind*length_traj,1:3)));
for i = 1:scale_ind*length_traj
    R = eul2rotm(flip(state_estimator.rpy(start_traj + i,:)));
    f_ff_w(i,:) =  -R * leg_control_command.f_ff(start_traj + i,1:3)';
end

plot(time,f_ff_w(:,3),'red');
hold on
plot(time,simulator_state.f_foot(start_traj:start_traj+scale_ind*length_traj,1,3),'blue');

legend('cmd','sim')


