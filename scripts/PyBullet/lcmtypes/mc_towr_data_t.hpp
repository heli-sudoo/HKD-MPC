/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 **/

#ifndef __mc_towr_data_t_hpp__
#define __mc_towr_data_t_hpp__

#include <lcm/lcm_coretypes.h>

#include <vector>


class mc_towr_data_t
{
    public:
        int16_t    len;

        std::vector< int32_t > microtime;

        std::vector< std::vector< float > > base_pos;

        std::vector< std::vector< float > > eul;

        std::vector< std::vector< float > > ee_pos;

        std::vector< std::vector< int16_t > > contact;

    public:
        /**
         * Encode a message into binary form.
         *
         * @param buf The output buffer.
         * @param offset Encoding starts at thie byte offset into @p buf.
         * @param maxlen Maximum number of bytes to write.  This should generally be
         *  equal to getEncodedSize().
         * @return The number of bytes encoded, or <0 on error.
         */
        inline int encode(void *buf, int offset, int maxlen) const;

        /**
         * Check how many bytes are required to encode this message.
         */
        inline int getEncodedSize() const;

        /**
         * Decode a message from binary form into this instance.
         *
         * @param buf The buffer containing the encoded message.
         * @param offset The byte offset into @p buf where the encoded message starts.
         * @param maxlen The maximum number of bytes to read while decoding.
         * @return The number of bytes decoded, or <0 if an error occured.
         */
        inline int decode(const void *buf, int offset, int maxlen);

        /**
         * Retrieve the 64-bit fingerprint identifying the structure of the message.
         * Note that the fingerprint is the same for all instances of the same
         * message type, and is a fingerprint on the message type definition, not on
         * the message contents.
         */
        inline static int64_t getHash();

        /**
         * Returns "mc_towr_data_t"
         */
        inline static const char* getTypeName();

        // LCM support functions. Users should not call these
        inline int _encodeNoHash(void *buf, int offset, int maxlen) const;
        inline int _getEncodedSizeNoHash() const;
        inline int _decodeNoHash(const void *buf, int offset, int maxlen);
        inline static uint64_t _computeHash(const __lcm_hash_ptr *p);
};

int mc_towr_data_t::encode(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;
    int64_t hash = getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

int mc_towr_data_t::decode(const void *buf, int offset, int maxlen)
{
    int pos = 0, thislen;

    int64_t msg_hash;
    thislen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &msg_hash, 1);
    if (thislen < 0) return thislen; else pos += thislen;
    if (msg_hash != getHash()) return -1;

    thislen = this->_decodeNoHash(buf, offset + pos, maxlen - pos);
    if (thislen < 0) return thislen; else pos += thislen;

    return pos;
}

int mc_towr_data_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t mc_towr_data_t::getHash()
{
    static int64_t hash = static_cast<int64_t>(_computeHash(NULL));
    return hash;
}

const char* mc_towr_data_t::getTypeName()
{
    return "mc_towr_data_t";
}

int mc_towr_data_t::_encodeNoHash(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;

    tlen = __int16_t_encode_array(buf, offset + pos, maxlen - pos, &this->len, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    if(this->len > 0) {
        tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->microtime[0], this->len);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < this->len; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->base_pos[a0][0], 3);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < this->len; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->eul[a0][0], 3);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < this->len; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->ee_pos[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < this->len; a0++) {
        tlen = __int16_t_encode_array(buf, offset + pos, maxlen - pos, &this->contact[a0][0], 4);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

int mc_towr_data_t::_decodeNoHash(const void *buf, int offset, int maxlen)
{
    int pos = 0, tlen;

    tlen = __int16_t_decode_array(buf, offset + pos, maxlen - pos, &this->len, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    if(this->len) {
        this->microtime.resize(this->len);
        tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->microtime[0], this->len);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    try {
        this->base_pos.resize(this->len);
    } catch (...) {
        return -1;
    }
    for (int a0 = 0; a0 < this->len; a0++) {
        if(3) {
            this->base_pos[a0].resize(3);
            tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->base_pos[a0][0], 3);
            if(tlen < 0) return tlen; else pos += tlen;
        }
    }

    try {
        this->eul.resize(this->len);
    } catch (...) {
        return -1;
    }
    for (int a0 = 0; a0 < this->len; a0++) {
        if(3) {
            this->eul[a0].resize(3);
            tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->eul[a0][0], 3);
            if(tlen < 0) return tlen; else pos += tlen;
        }
    }

    try {
        this->ee_pos.resize(this->len);
    } catch (...) {
        return -1;
    }
    for (int a0 = 0; a0 < this->len; a0++) {
        if(12) {
            this->ee_pos[a0].resize(12);
            tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->ee_pos[a0][0], 12);
            if(tlen < 0) return tlen; else pos += tlen;
        }
    }

    try {
        this->contact.resize(this->len);
    } catch (...) {
        return -1;
    }
    for (int a0 = 0; a0 < this->len; a0++) {
        if(4) {
            this->contact[a0].resize(4);
            tlen = __int16_t_decode_array(buf, offset + pos, maxlen - pos, &this->contact[a0][0], 4);
            if(tlen < 0) return tlen; else pos += tlen;
        }
    }

    return pos;
}

int mc_towr_data_t::_getEncodedSizeNoHash() const
{
    int enc_size = 0;
    enc_size += __int16_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, this->len);
    enc_size += this->len * __float_encoded_array_size(NULL, 3);
    enc_size += this->len * __float_encoded_array_size(NULL, 3);
    enc_size += this->len * __float_encoded_array_size(NULL, 12);
    enc_size += this->len * __int16_t_encoded_array_size(NULL, 4);
    return enc_size;
}

uint64_t mc_towr_data_t::_computeHash(const __lcm_hash_ptr *)
{
    uint64_t hash = 0x6a164ca064407912LL;
    return (hash<<1) + ((hash>>63)&1);
}

#endif
