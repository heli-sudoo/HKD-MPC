function plot_joint_trajectories(q, qr, dT, dt)
N = size(q, 2);
Nr = size(qr, 2);
t = dT*(0:N-1);
tr = dt*(0:Nr-1);

figure
subplot(2,3,1)
plot(t, q(1,:));
hold on
plot(tr, qr(1,:));
plot(t, q(4,:));
plot(tr, qr(4,:));
legend('abad_r','abad_r_r','abad_l','abad_l_r');

subplot(2,3,2)
plot(t, q(2,:));
hold on
plot(tr, qr(2,:));
plot(t, q(5,:));
plot(tr, qr(5,:));
legend('hip_r','hip_r_r','hip_l','hip_l_r');


subplot(2,3,3)
plot(t, q(3,:));
hold on
plot(tr, qr(3,:));
plot(t, q(6,:));
plot(tr, qr(6,:));
legend('knee_r','knee_r_r','knee_l','knee_l_r');

subplot(2,3,4)
plot(t, q(7,:));
hold on
plot(tr, qr(7,:));
plot(t, q(10,:));
plot(tr, qr(10,:));
% legend('abad_r','abad_r_r','abad_l','abad_l_r');


subplot(2,3,5)
plot(t, q(8,:));
hold on
plot(tr, qr(8,:));
plot(t, q(11,:));
plot(tr, qr(11,:));
% legend('hip_r','hip_r_r','hip_l','hip_l_r');


subplot(2,3,6)
plot(t, q(9,:));
hold on
plot(tr, qr(9,:));
plot(t, q(12,:));
plot(tr, qr(12,:));
% legend('knee_r','knee_r_r','knee_l','knee_l_r');


end