
import numpy as np    
from pybullet_utils.bullet_client import BulletClient
import pybullet
import pybullet_data as pd



class MiniCheetah:
    DEFAULT_JOINT_POSE = [0, 0.8, -1.6, 
                        0, 0.8, -1.6, 
                        0, 0.8, -1.6, 
                        0, 0.8, -1.6]        
    INIT_POS = [0,0,0.25]
    INIT_QUAT = [0,0,0,1]
    EE_ID = [3, 7, 11, 15] # FL, FR, HL, HR
    JOINT_DAMPING = [0.1, 0.05, 0.01,
                     0.1, 0.05, 0.01,
                     0.1, 0.05, 0.01,
                     0.1, 0.05, 0.01]  
    def __init__(self, urdf_file=None):        
        # damping used by numerical IK solver
        self.pb = BulletClient(pybullet.DIRECT)
        self.pb.setAdditionalSearchPath(pd.getDataPath())        
        if urdf_file==None:
            self.urdf_file = 'mini_cheetah/mini_cheetah.urdf'
        else:
            self.urdf_file = urdf_file
        self.robot = self.pb.loadURDF(self.urdf_file)
        self.set_pose(self.INIT_POS, self.INIT_QUAT, self.DEFAULT_JOINT_POSE)

    def ik(self, 
           base_pos = None, 
           base_quat = None, 
           foot_pos = None):
        """
            Compute inverse kinmatics
            args: 
                base_pos: base position nadrray of dim 3
                base_quat: base orientation in quaterion
                foot_pos: foot positions in world frame ndarray of dim 12            
        """
        # reset base position and orientation to solve IK in world frame
        self.pb.resetBasePositionAndOrientation(self.robot, base_pos, base_quat)
        joint_lim_low, joint_lim_high = self.get_joint_limits()
        joint_angles = self.pb.calculateInverseKinematics2(self.robot, self.EE_ID,
                                                    foot_pos.reshape([4,3]),
                                                    jointDamping = self.JOINT_DAMPING,
                                                    lowerLimits=joint_lim_low,
                                                    upperLimits=joint_lim_high,
                                                    solver=0)
        joint_angles = np.array(joint_angles)
        return joint_angles.reshape(12)    
        
    def print_link_jnt_info(self):
        num_joints = self.pb.getNumJoints(self.robot)
        for j in range(num_joints):
            jnt_info = self.pb.getJointInfo(self.robot, j)
            jnt_name = jnt_info[1]
            link_name = jnt_info[12]
            print("joint {} name".format(j), jnt_name)
            print("link of joint {} name".format(j), link_name)
    
    def getURDFFile(self):
        return self.urdf_file
    
    def set_pose(self, root_pos, root_quat, jnt_pose):
        num_joints = self.pb.getNumJoints(self.robot) 
        self.pb.resetBasePositionAndOrientation(self.robot, root_pos, root_quat)
        for j in range(num_joints):
            j_info = self.pb.getJointInfo(self.robot, j)
            j_state = self.pb.getJointStateMultiDof(self.robot, j)

            j_pose_idx = j_info[3]
            j_pose_size = len(j_state[0])
            j_vel_size = len(j_state[1])

            if (j_pose_size > 0):
                j_pose_idx = j_pose_idx - 7 # shift joint index by 7 to exclude base position and orientation
                j_pose = jnt_pose[j_pose_idx:(j_pose_idx + j_pose_size)]
                j_vel = np.zeros(j_vel_size)
                self.pb.resetJointStateMultiDof(self.robot, j, j_pose, j_vel)

    
