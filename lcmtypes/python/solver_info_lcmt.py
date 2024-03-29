"""LCM type definitions
This file automatically generated by lcm.
DO NOT MODIFY BY HAND!!!!
"""

try:
    import cStringIO.StringIO as BytesIO
except ImportError:
    from io import BytesIO
import struct

class solver_info_lcmt(object):
    __slots__ = ["n_iter", "cost", "dyn_feas", "eqn_feas", "ineq_feas"]

    __typenames__ = ["int32_t", "float", "float", "float", "float"]

    __dimensions__ = [None, ["n_iter"], ["n_iter"], ["n_iter"], ["n_iter"]]

    def __init__(self):
        self.n_iter = 0
        self.cost = []
        self.dyn_feas = []
        self.eqn_feas = []
        self.ineq_feas = []

    def encode(self):
        buf = BytesIO()
        buf.write(solver_info_lcmt._get_packed_fingerprint())
        self._encode_one(buf)
        return buf.getvalue()

    def _encode_one(self, buf):
        buf.write(struct.pack(">i", self.n_iter))
        buf.write(struct.pack('>%df' % self.n_iter, *self.cost[:self.n_iter]))
        buf.write(struct.pack('>%df' % self.n_iter, *self.dyn_feas[:self.n_iter]))
        buf.write(struct.pack('>%df' % self.n_iter, *self.eqn_feas[:self.n_iter]))
        buf.write(struct.pack('>%df' % self.n_iter, *self.ineq_feas[:self.n_iter]))

    def decode(data):
        if hasattr(data, 'read'):
            buf = data
        else:
            buf = BytesIO(data)
        if buf.read(8) != solver_info_lcmt._get_packed_fingerprint():
            raise ValueError("Decode error")
        return solver_info_lcmt._decode_one(buf)
    decode = staticmethod(decode)

    def _decode_one(buf):
        self = solver_info_lcmt()
        self.n_iter = struct.unpack(">i", buf.read(4))[0]
        self.cost = struct.unpack('>%df' % self.n_iter, buf.read(self.n_iter * 4))
        self.dyn_feas = struct.unpack('>%df' % self.n_iter, buf.read(self.n_iter * 4))
        self.eqn_feas = struct.unpack('>%df' % self.n_iter, buf.read(self.n_iter * 4))
        self.ineq_feas = struct.unpack('>%df' % self.n_iter, buf.read(self.n_iter * 4))
        return self
    _decode_one = staticmethod(_decode_one)

    _hash = None
    def _get_hash_recursive(parents):
        if solver_info_lcmt in parents: return 0
        tmphash = (0xc0fe49bc4d9f56eb) & 0xffffffffffffffff
        tmphash  = (((tmphash<<1)&0xffffffffffffffff) + (tmphash>>63)) & 0xffffffffffffffff
        return tmphash
    _get_hash_recursive = staticmethod(_get_hash_recursive)
    _packed_fingerprint = None

    def _get_packed_fingerprint():
        if solver_info_lcmt._packed_fingerprint is None:
            solver_info_lcmt._packed_fingerprint = struct.pack(">Q", solver_info_lcmt._get_hash_recursive([]))
        return solver_info_lcmt._packed_fingerprint
    _get_packed_fingerprint = staticmethod(_get_packed_fingerprint)

