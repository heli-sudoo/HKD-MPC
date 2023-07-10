clear, clc, close all;
%% mpc logs
mpc_time = 0.2;
iter = mpc_time / 0.01 + 1;
mpc_soln = readmatrix("../../../HKDMPC/log/state_log.txt");
time = mpc_time:0.01:mpc_time+0.01*(size(mpc_soln(:,1))-1);

%% reference
name = 'AngledBoxRight';
Gait = "../PreProcessedData/MIP_Hopping/"+name+"/";

body_states = readmatrix(Gait + "body_state.csv");
body_states(:,1:3) = flip(body_states(:,1:3),2);

plot(time,mpc_soln(:,3))
hold on
plot(time,body_states(iter:iter+length(time)-1,3));

legend('MPC Solution', 'Reference')