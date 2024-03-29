/* LCM type definition class file
 * This file was automatically generated by lcm-gen
 * DO NOT MODIFY BY HAND!!!!
 */

package lcmtypes;
 
import java.io.*;
import java.util.*;
import lcm.lcm.*;
 
public final class MHPC_Command_lcmt implements lcm.lcm.LCMEncodable
{
    public int N_mpcsteps;
    public float mpc_times[];
    public float torque[][];
    public float eul[][];
    public float pos[][];
    public float qJ[][];
    public float vWorld[][];
    public float eulrate[][];
    public float qJd[][];
    public float feedback[][];
    public int contacts[][];
    public float statusTimes[][];
    public float solve_time;
 
    public MHPC_Command_lcmt()
    {
        mpc_times = new float[4];
        torque = new float[4][12];
        eul = new float[4][3];
        pos = new float[4][3];
        qJ = new float[4][12];
        vWorld = new float[4][3];
        eulrate = new float[4][3];
        qJd = new float[4][12];
        feedback = new float[4][432];
        contacts = new int[4][4];
        statusTimes = new float[4][4];
    }
 
    public static final long LCM_FINGERPRINT;
    public static final long LCM_FINGERPRINT_BASE = 0xadd72c9ceb3b5e16L;
 
    static {
        LCM_FINGERPRINT = _hashRecursive(new ArrayList<Class<?>>());
    }
 
    public static long _hashRecursive(ArrayList<Class<?>> classes)
    {
        if (classes.contains(lcmtypes.MHPC_Command_lcmt.class))
            return 0L;
 
        classes.add(lcmtypes.MHPC_Command_lcmt.class);
        long hash = LCM_FINGERPRINT_BASE
            ;
        classes.remove(classes.size() - 1);
        return (hash<<1) + ((hash>>63)&1);
    }
 
    public void encode(DataOutput outs) throws IOException
    {
        outs.writeLong(LCM_FINGERPRINT);
        _encodeRecursive(outs);
    }
 
    public void _encodeRecursive(DataOutput outs) throws IOException
    {
        outs.writeInt(this.N_mpcsteps); 
 
        for (int a = 0; a < 4; a++) {
            outs.writeFloat(this.mpc_times[a]); 
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 12; b++) {
                outs.writeFloat(this.torque[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                outs.writeFloat(this.eul[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                outs.writeFloat(this.pos[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 12; b++) {
                outs.writeFloat(this.qJ[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                outs.writeFloat(this.vWorld[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                outs.writeFloat(this.eulrate[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 12; b++) {
                outs.writeFloat(this.qJd[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 432; b++) {
                outs.writeFloat(this.feedback[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                outs.writeInt(this.contacts[a][b]); 
            }
        }
 
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                outs.writeFloat(this.statusTimes[a][b]); 
            }
        }
 
        outs.writeFloat(this.solve_time); 
 
    }
 
    public MHPC_Command_lcmt(byte[] data) throws IOException
    {
        this(new LCMDataInputStream(data));
    }
 
    public MHPC_Command_lcmt(DataInput ins) throws IOException
    {
        if (ins.readLong() != LCM_FINGERPRINT)
            throw new IOException("LCM Decode error: bad fingerprint");
 
        _decodeRecursive(ins);
    }
 
    public static lcmtypes.MHPC_Command_lcmt _decodeRecursiveFactory(DataInput ins) throws IOException
    {
        lcmtypes.MHPC_Command_lcmt o = new lcmtypes.MHPC_Command_lcmt();
        o._decodeRecursive(ins);
        return o;
    }
 
    public void _decodeRecursive(DataInput ins) throws IOException
    {
        this.N_mpcsteps = ins.readInt();
 
        this.mpc_times = new float[(int) 4];
        for (int a = 0; a < 4; a++) {
            this.mpc_times[a] = ins.readFloat();
        }
 
        this.torque = new float[(int) 4][(int) 12];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 12; b++) {
                this.torque[a][b] = ins.readFloat();
            }
        }
 
        this.eul = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                this.eul[a][b] = ins.readFloat();
            }
        }
 
        this.pos = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                this.pos[a][b] = ins.readFloat();
            }
        }
 
        this.qJ = new float[(int) 4][(int) 12];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 12; b++) {
                this.qJ[a][b] = ins.readFloat();
            }
        }
 
        this.vWorld = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                this.vWorld[a][b] = ins.readFloat();
            }
        }
 
        this.eulrate = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 3; b++) {
                this.eulrate[a][b] = ins.readFloat();
            }
        }
 
        this.qJd = new float[(int) 4][(int) 12];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 12; b++) {
                this.qJd[a][b] = ins.readFloat();
            }
        }
 
        this.feedback = new float[(int) 4][(int) 432];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 432; b++) {
                this.feedback[a][b] = ins.readFloat();
            }
        }
 
        this.contacts = new int[(int) 4][(int) 4];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                this.contacts[a][b] = ins.readInt();
            }
        }
 
        this.statusTimes = new float[(int) 4][(int) 4];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                this.statusTimes[a][b] = ins.readFloat();
            }
        }
 
        this.solve_time = ins.readFloat();
 
    }
 
    public lcmtypes.MHPC_Command_lcmt copy()
    {
        lcmtypes.MHPC_Command_lcmt outobj = new lcmtypes.MHPC_Command_lcmt();
        outobj.N_mpcsteps = this.N_mpcsteps;
 
        outobj.mpc_times = new float[(int) 4];
        System.arraycopy(this.mpc_times, 0, outobj.mpc_times, 0, 4); 
        outobj.torque = new float[(int) 4][(int) 12];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.torque[a], 0, outobj.torque[a], 0, 12);        }
 
        outobj.eul = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.eul[a], 0, outobj.eul[a], 0, 3);        }
 
        outobj.pos = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.pos[a], 0, outobj.pos[a], 0, 3);        }
 
        outobj.qJ = new float[(int) 4][(int) 12];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.qJ[a], 0, outobj.qJ[a], 0, 12);        }
 
        outobj.vWorld = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.vWorld[a], 0, outobj.vWorld[a], 0, 3);        }
 
        outobj.eulrate = new float[(int) 4][(int) 3];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.eulrate[a], 0, outobj.eulrate[a], 0, 3);        }
 
        outobj.qJd = new float[(int) 4][(int) 12];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.qJd[a], 0, outobj.qJd[a], 0, 12);        }
 
        outobj.feedback = new float[(int) 4][(int) 432];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.feedback[a], 0, outobj.feedback[a], 0, 432);        }
 
        outobj.contacts = new int[(int) 4][(int) 4];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.contacts[a], 0, outobj.contacts[a], 0, 4);        }
 
        outobj.statusTimes = new float[(int) 4][(int) 4];
        for (int a = 0; a < 4; a++) {
            System.arraycopy(this.statusTimes[a], 0, outobj.statusTimes[a], 0, 4);        }
 
        outobj.solve_time = this.solve_time;
 
        return outobj;
    }
 
}

