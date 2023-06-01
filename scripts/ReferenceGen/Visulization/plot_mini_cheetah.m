function plot_mini_cheetah(config)
robot = importrobot("mini_cheetah_mesh.urdf");
robot = float_base_quadruped(robot);
robot.DataFormat = 'column';

figure(198)
show(robot, config, "Frames","off");
hold on
end