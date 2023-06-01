function plot_torque_trajectories(u, ur, dT, dt)
N = size(u, 2);
Nr = size(ur, 2);
t = dT*(0:N-1);
tr = dt*(0:Nr-1);
figure
subplot(2,3,1)
plot(t, u(1,:));
hold on
plot(tr, ur(1,:));
plot(t, u(4,:));
plot(tr, ur(4,:));
legend('abad_r','abad_r_r','abad_l','abad_l_r');


subplot(2,3,2)
plot(t, u(2,:));
hold on
plot(tr, ur(2,:));
plot(t, u(5,:));
plot(tr, ur(5,:));
legend('hip_r','hip_r_r','hip_l','hip_l_r');

subplot(2,3,3)
plot(t, u(3,:));
hold on
plot(tr, ur(3,:));
plot(t, u(6,:));
plot(tr, ur(6,:));
legend('knee_r','knee_r_r','knee_l','knee_l_r');

subplot(2,3,4)
plot(t, u(7,:));
hold on
plot(tr, ur(7,:));
plot(t, u(10,:));
plot(tr, ur(10,:));

subplot(2,3,5)
plot(t, u(8,:));
hold on
plot(tr, ur(8,:));
plot(t, u(11,:));
plot(tr, ur(11,:));

subplot(2,3,6)
plot(t, u(9,:));
hold on
plot(tr, ur(9,:));
plot(t, u(12,:));
plot(tr, ur(12,:));

end