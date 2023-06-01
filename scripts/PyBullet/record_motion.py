import numpy as np

import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
os.sys.path.insert(0, parentdir)


from mini_cheetah import MiniCheetah

import lcm
from lcmtypes.mc_towr_data_t import mc_towr_data_t

robot = MiniCheetah()    

def towr_data_lcm_handler(channel, data):
    print("received towr lcm message")
    towr_msg = mc_towr_data_t.decode(data)
    base_traj = []
    jnt_traj = []
    contact_traj = []
    ee_pos_traj = []

    traj_len = towr_msg.len
    time_traj = [t/1e6 for t in towr_msg.microtime]
    for k in range(traj_len):
        eul_k = np.array(towr_msg.eul[k])
        pos_k = np.array(towr_msg.base_pos[k])
        ee_pos_k = np.array(towr_msg.ee_pos[k])
        quat_k = np.array(robot.pb.getQuaternionFromEuler(eul_k))
        jnt_k = robot.ik(pos_k, quat_k, ee_pos_k)
        contact_k = np.array(towr_msg.contact[k])
        eulrate_k = np.array(towr_msg.eulrate[k])
        angrate_k = eulrate2angrate(eul_k, eulrate_k)
        vel_k = np.array(towr_msg.base_vel[k])

        base_traj.append(np.hstack((eul_k, pos_k, angrate_k, vel_k)))
        jnt_traj.append(jnt_k)
        contact_traj.append(contact_k)
        ee_pos_traj.append(ee_pos_k)
        
    write_traj_to_file(time_traj, base_traj, jnt_traj, ee_pos_traj, contact_traj)

def write_traj_to_file(time_traj, base_traj, jnt_traj, ee_pos_traj, ctact_traj):
    np.savetxt("data/time.csv", np.asarray(time_traj), delimiter=",", fmt='%8.4f')
    np.savetxt("data/body_state.csv", base_traj, delimiter=",", fmt='%8.4f')
    np.savetxt("data/ee_pos.csv", np.asarray(ee_pos_traj), delimiter=",", fmt='%8.4f')
    np.savetxt("data/jnt.csv", np.asarray(jnt_traj), delimiter=",", fmt='%8.4f')
    np.savetxt("data/contact.csv", np.asarray(ctact_traj), delimiter=",", fmt='%u')

def eulrate2angrate(eul, eulrate):
    b = eul[1]
    r = eul[2]
    T = np.matrix([[-np.sin(b), 0,  1], 
        [np.cos(b)*np.sin(r), np.cos(r), 0],
        [np.cos(b)*np.cos(r), -np.sin(r),  0]])
    angrate = np.matmul(T, eulrate)
    angrate = np.asarray(angrate)
    return angrate.reshape(3)

lc = lcm.LCM()
subscription = lc.subscribe("TOWR", towr_data_lcm_handler)
robot.print_link_jnt_info()

try:
    while True:
        lc.handle()
except KeyboardInterrupt:
    pass