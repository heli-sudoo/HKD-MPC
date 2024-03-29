/* LCM type definition class file
 * This file was automatically generated by lcm-gen
 * DO NOT MODIFY BY HAND!!!!
 */

package lcmtypes;
 
import java.io.*;
import java.util.*;
import lcm.lcm.*;
 
public final class MHPC_Data_lcmt implements lcm.lcm.LCMEncodable
{
    public boolean reset_mpc;
    public boolean MS;
    public double mpctime;
    public float pos[];
    public float eul[];
    public float qJ[];
    public float vWorld[];
    public float eulrate[];
    public float qJd[];
 
    public MHPC_Data_lcmt()
    {
        pos = new float[3];
        eul = new float[3];
        qJ = new float[12];
        vWorld = new float[3];
        eulrate = new float[3];
        qJd = new float[12];
    }
 
    public static final long LCM_FINGERPRINT;
    public static final long LCM_FINGERPRINT_BASE = 0xfcc6a78db3f0ef2bL;
 
    static {
        LCM_FINGERPRINT = _hashRecursive(new ArrayList<Class<?>>());
    }
 
    public static long _hashRecursive(ArrayList<Class<?>> classes)
    {
        if (classes.contains(lcmtypes.MHPC_Data_lcmt.class))
            return 0L;
 
        classes.add(lcmtypes.MHPC_Data_lcmt.class);
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
        outs.writeByte( this.reset_mpc ? 1 : 0); 
 
        outs.writeByte( this.MS ? 1 : 0); 
 
        outs.writeDouble(this.mpctime); 
 
        for (int a = 0; a < 3; a++) {
            outs.writeFloat(this.pos[a]); 
        }
 
        for (int a = 0; a < 3; a++) {
            outs.writeFloat(this.eul[a]); 
        }
 
        for (int a = 0; a < 12; a++) {
            outs.writeFloat(this.qJ[a]); 
        }
 
        for (int a = 0; a < 3; a++) {
            outs.writeFloat(this.vWorld[a]); 
        }
 
        for (int a = 0; a < 3; a++) {
            outs.writeFloat(this.eulrate[a]); 
        }
 
        for (int a = 0; a < 12; a++) {
            outs.writeFloat(this.qJd[a]); 
        }
 
    }
 
    public MHPC_Data_lcmt(byte[] data) throws IOException
    {
        this(new LCMDataInputStream(data));
    }
 
    public MHPC_Data_lcmt(DataInput ins) throws IOException
    {
        if (ins.readLong() != LCM_FINGERPRINT)
            throw new IOException("LCM Decode error: bad fingerprint");
 
        _decodeRecursive(ins);
    }
 
    public static lcmtypes.MHPC_Data_lcmt _decodeRecursiveFactory(DataInput ins) throws IOException
    {
        lcmtypes.MHPC_Data_lcmt o = new lcmtypes.MHPC_Data_lcmt();
        o._decodeRecursive(ins);
        return o;
    }
 
    public void _decodeRecursive(DataInput ins) throws IOException
    {
        this.reset_mpc = ins.readByte()!=0;
 
        this.MS = ins.readByte()!=0;
 
        this.mpctime = ins.readDouble();
 
        this.pos = new float[(int) 3];
        for (int a = 0; a < 3; a++) {
            this.pos[a] = ins.readFloat();
        }
 
        this.eul = new float[(int) 3];
        for (int a = 0; a < 3; a++) {
            this.eul[a] = ins.readFloat();
        }
 
        this.qJ = new float[(int) 12];
        for (int a = 0; a < 12; a++) {
            this.qJ[a] = ins.readFloat();
        }
 
        this.vWorld = new float[(int) 3];
        for (int a = 0; a < 3; a++) {
            this.vWorld[a] = ins.readFloat();
        }
 
        this.eulrate = new float[(int) 3];
        for (int a = 0; a < 3; a++) {
            this.eulrate[a] = ins.readFloat();
        }
 
        this.qJd = new float[(int) 12];
        for (int a = 0; a < 12; a++) {
            this.qJd[a] = ins.readFloat();
        }
 
    }
 
    public lcmtypes.MHPC_Data_lcmt copy()
    {
        lcmtypes.MHPC_Data_lcmt outobj = new lcmtypes.MHPC_Data_lcmt();
        outobj.reset_mpc = this.reset_mpc;
 
        outobj.MS = this.MS;
 
        outobj.mpctime = this.mpctime;
 
        outobj.pos = new float[(int) 3];
        System.arraycopy(this.pos, 0, outobj.pos, 0, 3); 
        outobj.eul = new float[(int) 3];
        System.arraycopy(this.eul, 0, outobj.eul, 0, 3); 
        outobj.qJ = new float[(int) 12];
        System.arraycopy(this.qJ, 0, outobj.qJ, 0, 12); 
        outobj.vWorld = new float[(int) 3];
        System.arraycopy(this.vWorld, 0, outobj.vWorld, 0, 3); 
        outobj.eulrate = new float[(int) 3];
        System.arraycopy(this.eulrate, 0, outobj.eulrate, 0, 3); 
        outobj.qJd = new float[(int) 12];
        System.arraycopy(this.qJd, 0, outobj.qJd, 0, 12); 
        return outobj;
    }
 
}

