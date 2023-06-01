classdef MiniCheetahGraphics < handle
    properties
        state_traj
        time
        state_ref
        time_ref
    end

    methods
        function G = MiniCheetahGraphics()
        end

        function G = addTrajectory(G,time_in,state_traj_in)
            G.time = time_in;
            G.state_traj = state_traj_in;
        end

        function G = addReference(G, time_ref_in, state_ref_in)
            G.time_ref = time_ref_in;
            G.state_ref = state_ref_in;     
            assert(length(G.time)==length(time_ref_in),"dimensions of reference and trajectory do not match");
        end
        
        function G = plot_base_pos(G)
            base_traj = G.state_traj(4:6, :);
            base_ref = [];
            if ~isempty(G.state_ref)
                base_ref = G.state_ref(4:6, :);
            end
            figure
            subplot(3,1,1)
            plot(G.time, base_traj(1,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(1,:));
            end
            ylabel('x (m)');
            subplot(3,1,2)
            plot(G.time, base_traj(2,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(2,:));
            end
            ylabel('y (m)');
            subplot(3,1,3)
            plot(G.time, base_traj(3,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(3,:));
            end
            ylabel('z (m)');
            xlabel('time (s)');
        end

        function G = plot_base_vel(G)
            base_traj = G.state_traj(10:12, :);
            base_ref = [];
            if ~isempty(G.state_ref)
                base_ref = G.state_ref(10:12, :);
            end
            figure
            subplot(3,1,1)
            plot(G.time, base_traj(1,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(1,:));
            end
            ylabel('vx (m/s)');
            legend('vel', 'ref');
            subplot(3,1,2)
            plot(G.time, base_traj(2,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(2,:));
            end
            ylabel('vy (m/s)');
            subplot(3,1,3)
            plot(G.time, base_traj(3,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(3,:));
            end
            ylabel('vz (m/s)');
            xlabel('time (s)');
        end

        function G = plot_base_euler(G)
            base_traj = G.state_traj(1:3, :);
            base_ref = [];
            if ~isempty(G.state_ref)
                base_ref = G.state_ref(1:3, :);
            end
            figure
            subplot(3,1,1)
            plot(G.time, base_traj(1,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(1,:));
            end
            ylabel('yaw (rad)');
            subplot(3,1,2)
            plot(G.time, base_traj(2,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(2,:));
            end
            ylabel('pitch (rad)');
            subplot(3,1,3)
            plot(G.time, base_traj(3,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(3,:));
            end
            ylabel('roll (rad)');
            xlabel('time (s)');
        end

        function G = plot_base_angvel(G)
            base_traj = G.state_traj(7:9, :);
            base_ref = [];
            if ~isempty(G.state_ref)
                base_ref = G.state_ref(7:9, :);
            end
            figure
            subplot(3,1,1)
            plot(G.time, base_traj(1,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(1,:));
            end
            ylabel('wx (rad/s)');
            subplot(3,1,2)
            plot(G.time, base_traj(2,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(2,:));
            end
            ylabel('wy (rad/s)');
            subplot(3,1,3)
            plot(G.time, base_traj(3,:));
            hold on
            if ~isempty(base_ref)
                plot(G.time_ref, base_ref(3,:));
            end
            ylabel('wz (rad/s)');
            xlabel('time (s)');
        end
        function G = plot_abad_angles(G)
            abad_traj = G.state_traj([13,16,19,22], :);
            abad_ref = [];
            if ~isempty(G.state_ref)
                abad_ref = G.state_ref([13,16,19,22], :);
            end
            figure
            subplot(2,2,1)
            plot(G.time, abad_traj(1,:));
            hold on
            if ~isempty(abad_ref)
                plot(G.time_ref, abad_ref(1,:));
            end
            ylabel('angle (rad)');
            subplot(2,2,2)
            plot(G.time, abad_traj(2,:));
            hold on
            if ~isempty(abad_ref)
                plot(G.time_ref, abad_ref(2,:));
            end
            subplot(2,2,3)
            plot(G.time, abad_traj(3,:));
            hold on
            if ~isempty(abad_ref)
                plot(G.time_ref, abad_ref(3,:));
            end
            ylabel('angle (rad)');
            xlabel('time (s)');
            subplot(2,2,4)
            plot(G.time, abad_traj(4,:));
            hold on
            if ~isempty(abad_ref)
                plot(G.time_ref, abad_ref(4,:));
            end
            xlabel('time (s)');
        end
        function G = plot_hip_angles(G)            
            hip_traj = G.state_traj([14,17,20,23], :);
            hip_ref = [];
            if ~isempty(G.state_ref)
                hip_ref = G.state_ref([14,17,20,23], :);
            end
            figure
            subplot(2,2,1)
            plot(G.time, hip_traj(1,:));
            hold on
            if ~isempty(hip_ref)
                plot(G.time_ref, hip_ref(1,:));
            end
            ylabel('angle (rad)');
            subplot(2,2,2)
            plot(G.time, hip_traj(2,:));
            hold on
            if ~isempty(hip_ref)
                plot(G.time_ref, hip_ref(2,:));
            end
            subplot(2,2,3)
            plot(G.time, hip_traj(3,:));
            hold on
            if ~isempty(hip_ref)
                plot(G.time_ref, hip_ref(3,:));
            end
            ylabel('angle (rad)');
            xlabel('time (s)');
            subplot(2,2,4)
            plot(G.time, hip_traj(4,:));
            hold on
            if ~isempty(hip_ref)
                plot(G.time_ref, hip_ref(4,:));
            end
            xlabel('time (s)');
        end
        function G = plot_knee_angles(G)
            knee_traj = G.state_traj([15,18,21,24], :);
            knee_ref = [];
            if ~isempty(G.state_ref)
                knee_ref = G.state_ref([15,18,21,24], :);
            end
            figure
            subplot(2,2,1)
            plot(G.time, knee_traj(1,:));
            hold on
            if ~isempty(knee_ref)
                plot(G.time_ref, knee_ref(1,:));
            end
            ylabel('angle (rad)');
            subplot(2,2,2)
            plot(G.time, knee_traj(2,:));
            hold on
            if ~isempty(knee_ref)
                plot(G.time_ref, knee_ref(2,:));
            end
            subplot(2,2,3)
            plot(G.time, knee_traj(3,:));
            hold on
            if ~isempty(knee_ref)
                plot(G.time_ref, knee_ref(3,:));
            end
            ylabel('angle (rad)');
            xlabel('time (s)');
            subplot(2,2,4)
            plot(G.time, knee_traj(4,:));
            hold on
            if ~isempty(knee_ref)
                plot(G.time_ref, knee_ref(4,:));
            end
            xlabel('time (s)');
        end
    end
end