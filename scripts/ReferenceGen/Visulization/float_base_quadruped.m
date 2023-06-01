function float_base_robot = float_base_quadruped(robot)
% Create a float-base quadruped from a fix-base quadruped
% robot -> a rigidBodyTree converted from urdf file

% Modify the fixed-based quadruped to floating base
float_base_robot = rigidBodyTree;

body_x = rigidBody('vbody_x');
jnt_x = rigidBodyJoint('jnt_x', 'prismatic');
jnt_x.HomePosition = 0;
jnt_x.JointAxis = [1, 0, 0];
setFixedTransform(jnt_x, trvec2tform([0, 0, 0]));
body_x.Joint = jnt_x;
addBody(float_base_robot, body_x, 'base');

body_y = rigidBody('vbody_y');
jnt_y = rigidBodyJoint('jnt_y', 'prismatic');
jnt_y.HomePosition = 0;
jnt_y.JointAxis = [0, 1, 0];
setFixedTransform(jnt_y, trvec2tform([0, 0, 0]));
body_y.Joint = jnt_y;
addBody(float_base_robot, body_y, 'vbody_x');

body_z = rigidBody('vbody_z');
jnt_z = rigidBodyJoint('jnt_z', 'prismatic');
jnt_z.HomePosition = 0;
jnt_z.JointAxis = [0, 0, 1];
setFixedTransform(jnt_z, trvec2tform([0, 0, 0]));
body_z.Joint = jnt_z;
addBody(float_base_robot, body_z, 'vbody_y');

body_yaw = rigidBody('vbody_yaw');
jnt_yaw = rigidBodyJoint('vjnt_yaw', 'revolute');
jnt_yaw.HomePosition = 0;
jnt_yaw.JointAxis = [0, 0, 1];
setFixedTransform(jnt_yaw, trvec2tform([0, 0, 0]));
body_yaw.Joint = jnt_yaw;
addBody(float_base_robot, body_yaw, 'vbody_z');

body_pitch = rigidBody('vbody_pitch');
jnt_pitch = rigidBodyJoint('jnt_pitch', 'revolute');
jnt_pitch.HomePosition = 0;
jnt_pitch.JointAxis = [0, 1, 0];
setFixedTransform(jnt_pitch, trvec2tform([0, 0, 0]));
body_pitch.Joint = jnt_pitch;
addBody(float_base_robot, body_pitch, 'vbody_yaw');

body_roll = copy(robot.Base);
jnt_roll = rigidBodyJoint('jnt_roll', 'revolute');
jnt_roll.HomePosition = 0;
jnt_roll.JointAxis = [1, 0, 0];
setFixedTransform(jnt_roll, trvec2tform([0, 0, 0]));
body_roll.Joint = jnt_roll;
addBody(float_base_robot, body_roll, 'vbody_pitch');

% Modify the order of legs
% Mini Cheetah (in Cheetah Software) orders its legs as follows
% https://github.com/heli-sudoo/Cheetah-Software/blob/master/documentation/getting_started.md
%   Front
%   1  0    Right
%   3  2
%   Back
% URDF file orders its legs as follows
%   Front
%   0  1    Right
%   2  3
%   Back

leg_fl = removeBody(robot, 'abduct_fl');
leg_fr = removeBody(robot, 'abduct_fr');
leg_hl = removeBody(robot, 'abduct_hl');
leg_hr = removeBody(robot, 'abduct_hr');

addSubtree(float_base_robot, body_roll.Name, leg_fr);
addSubtree(float_base_robot, body_roll.Name, leg_fl);
addSubtree(float_base_robot, body_roll.Name, leg_hr);
addSubtree(float_base_robot, body_roll.Name, leg_hl);

end