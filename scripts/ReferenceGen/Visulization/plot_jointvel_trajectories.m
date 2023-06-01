function plot_jointvel_trajectories(qd, qdr, dT, dt)
N = size(qd, 2);
Nr = size(qdr, 2);
t = dT*(0:N-1);
tr = dt*(0:Nr-1);
figure
subplot(2,3,1)
plot(t, qd(1,:));
hold on
plot(tr, qdr(1,:));
plot(t, qd(4,:));
plot(tr, qdr(4,:));
legend('abad_r','abad_r_r','abad_l','abad_l_r');

subplot(2,3,2)
plot(t, qd(2,:));
hold on
plot(tr, qdr(2,:));
plot(t, qd(5,:));
plot(tr, qdr(5,:));
legend('hip_r','hip_r_r','hip_l','hip_l_r');


subplot(2,3,3)
plot(t, qd(3,:));
hold on
plot(tr, qdr(3,:));
plot(t, qd(6,:));
plot(tr, qdr(6,:));
legend('knee_r','knee_r_r','knee_l','knee_l_r');

subplot(2,3,4)
plot(t, qd(7,:));
hold on
plot(tr, qdr(7,:));
plot(t, qd(10,:));
plot(tr, qdr(10,:));
% legend('abad_r','abad_r_r','abad_l','abad_l_r');


subplot(2,3,5)
plot(t, qd(8,:));
hold on
plot(tr, qdr(8,:));
plot(t, qd(11,:));
plot(tr, qdr(11,:));
% legend('hip_r','hip_r_r','hip_l','hip_l_r');


subplot(2,3,6)
plot(t, qd(9,:));
hold on
plot(tr, qdr(9,:));
plot(t, qd(12,:));
plot(tr, qdr(12,:));
% legend('knee_r','knee_r_r','knee_l','knee_l_r');


end