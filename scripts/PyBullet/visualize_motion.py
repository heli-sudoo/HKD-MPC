
import time
from turtle import pd
import numpy as np


from scripts.PyBullet.mini_cheetah import MiniCheetah
from scripts.PyBullet.animator import Animator

import lcm
from lcmtypes.python.visualize_quadState_lcmt import visualize_quadState_lcmt
from lcmtypes.python.visualize_quadTraj_lcmt import visualize_quadTraj_lcmt
from lcmtypes.python.MHPC_Command_lcmt import MHPC_Command_lcmt
import threading

urdf_filename =  "urdf/mini_cheetah_simple_correctedInertia.urdf"
robot = MiniCheetah(urdf_file=urdf_filename)    
animator = Animator(robot)
animator.initialization()

def MHPC_COMMAND_lcm_handler(channel, data):
    print("received MHPC Command")
    msg = MHPC_Command_lcmt.decode(data)
    eul = np.array(msg.eul[0])
    pos = np.array(msg.pos[0])
    rpy = eul[[2,1,0]]    
    quat = np.array(robot.pb.getQuaternionFromEuler(rpy))
    qJ = np.array(msg.qJ[0])

    animator.set_pose(pos, quat, qJ)
 
def visualize_config_lcm_handler(channel, data):
    print("received visualization lcm message")
    msg = visualize_quadState_lcmt.decode(data)
        
    eul = np.array(msg.eul)
    pos = np.array(msg.pos)
    rpy = eul[[2,1,0]]
    ee_pos = np.array(msg.pFoot)
    quat = np.array(robot.pb.getQuaternionFromEuler(rpy))
    qJ = np.array(msg.qJ)
    qJd = np.array(msg.qJd)
    pFoot = np.array(msg.pFoot)    

    animator.set_pose(pos, quat, qJ)

def visualize_motion_lcm_handler(channel, data):
    print("received visualization lcm message")    
    msg = visualize_quadTraj_lcmt.decode(data)

    print("Number of states received: ", msg.len)    
    t, t_WB, t_SRB = [],[],[]
    WB_len = int(np.rint(msg.WB_plan_dur/msg.WB_dt+1)) 
    SRB_len = int(np.rint(msg.SRB_plan_dur/msg.SRB_dt))
        
    if msg.WB_plan_dur > 0:
        t_WB = np.linspace(0.0, msg.WB_plan_dur, WB_len+1)
    if msg.SRB_plan_dur > 0:
        t_SRB = np.linspace(0.0, msg.SRB_plan_dur, SRB_len+1)
        if len(t_WB) > 0:
            t_SRB = t_SRB + t_WB[-1]
    t = np.hstack((t_WB, t_SRB))

    eul = np.asarray([np.array(eul_k) for eul_k in msg.eul])
    pos = np.asarray([np.array(pos_k) for pos_k in msg.pos])    
    torque = np.asarray([np.array(torque_k) for torque_k in msg.torque])        
    GRF = np.asarray([np.array(GRF_k) for GRF_k in msg.grf])    

    # animator.plot_eul(t, eul)
    # animator.plot_pos(t, pos)
    # animator.plot_joint_torques(t, torque)
    # animator.plot_GRF(t, GRF)
    # animator.plot_feasibility(t, msg.feas)

    tau_sz = msg.len
    for k in range(tau_sz-1):
        eul_k = np.array(msg.eul[k])
        rpy_k = eul_k[[2,1,0]]
        quat_k = np.array(robot.pb.getQuaternionFromEuler(rpy_k))
        pos_k = np.array(msg.pos[k])        
        qJ_k = np.array(msg.qJ[k]) 

        animator.set_pose(pos_k, quat_k, qJ_k)        
        if k * msg.WB_dt <= msg.WB_plan_dur + 0.05:
            animator.recover_legs()
            time.sleep(2*msg.WB_dt)
        else:
            animator.hide_legs()
            time.sleep(msg.SRB_dt)
    

lc = lcm.LCM()
subscription1 = lc.subscribe("visualize_mc_state", visualize_config_lcm_handler)
subscription2 = lc.subscribe("visualize_mc_motion", visualize_motion_lcm_handler)
subscription3 = lc.subscribe("MHPC_COMMAND", MHPC_COMMAND_lcm_handler)
render_thread = threading.Thread(target=animator.render)
render_thread.start()
# robot.print_link_jnt_info()
try:
    while True:        
        lc.handle()
except KeyboardInterrupt:
    pass