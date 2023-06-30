/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 **/

#ifndef __hkd_command_lcmt_hpp__
#define __hkd_command_lcmt_hpp__

#include <lcm/lcm_coretypes.h>



class hkd_command_lcmt
{
    public:
        int32_t    N_mpcsteps;

        double     mpc_times[25];

        float      hkd_controls[25][24];

        float      des_body_state[25][12];

        int32_t    contacts[25][4];

        double     statusTimes[25][4];

        float      foot_placement[12];

        float      feedback[25][12][12];

        float      solve_time;

        float      qJ_ref[25][12];

        float      qJd_ref[25][12];

        float      terrain_info[25][6];

        float      foot_placement_rel_com[12];

        float      vcom_td[3];

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
         * Returns "hkd_command_lcmt"
         */
        inline static const char* getTypeName();

        // LCM support functions. Users should not call these
        inline int _encodeNoHash(void *buf, int offset, int maxlen) const;
        inline int _getEncodedSizeNoHash() const;
        inline int _decodeNoHash(const void *buf, int offset, int maxlen);
        inline static uint64_t _computeHash(const __lcm_hash_ptr *p);
};

int hkd_command_lcmt::encode(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;
    int64_t hash = getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

int hkd_command_lcmt::decode(const void *buf, int offset, int maxlen)
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

int hkd_command_lcmt::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t hkd_command_lcmt::getHash()
{
    static int64_t hash = static_cast<int64_t>(_computeHash(NULL));
    return hash;
}

const char* hkd_command_lcmt::getTypeName()
{
    return "hkd_command_lcmt";
}

int hkd_command_lcmt::_encodeNoHash(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->N_mpcsteps, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __double_encode_array(buf, offset + pos, maxlen - pos, &this->mpc_times[0], 25);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->hkd_controls[a0][0], 24);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->des_body_state[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->contacts[a0][0], 4);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __double_encode_array(buf, offset + pos, maxlen - pos, &this->statusTimes[a0][0], 4);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->foot_placement[0], 12);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < 25; a0++) {
        for (int a1 = 0; a1 < 12; a1++) {
            tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->feedback[a0][a1][0], 12);
            if(tlen < 0) return tlen; else pos += tlen;
        }
    }

    tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->solve_time, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->qJ_ref[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->qJd_ref[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->terrain_info[a0][0], 6);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->foot_placement_rel_com[0], 12);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __float_encode_array(buf, offset + pos, maxlen - pos, &this->vcom_td[0], 3);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

int hkd_command_lcmt::_decodeNoHash(const void *buf, int offset, int maxlen)
{
    int pos = 0, tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->N_mpcsteps, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __double_decode_array(buf, offset + pos, maxlen - pos, &this->mpc_times[0], 25);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->hkd_controls[a0][0], 24);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->des_body_state[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->contacts[a0][0], 4);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __double_decode_array(buf, offset + pos, maxlen - pos, &this->statusTimes[a0][0], 4);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->foot_placement[0], 12);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < 25; a0++) {
        for (int a1 = 0; a1 < 12; a1++) {
            tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->feedback[a0][a1][0], 12);
            if(tlen < 0) return tlen; else pos += tlen;
        }
    }

    tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->solve_time, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->qJ_ref[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->qJd_ref[a0][0], 12);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    for (int a0 = 0; a0 < 25; a0++) {
        tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->terrain_info[a0][0], 6);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->foot_placement_rel_com[0], 12);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __float_decode_array(buf, offset + pos, maxlen - pos, &this->vcom_td[0], 3);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

int hkd_command_lcmt::_getEncodedSizeNoHash() const
{
    int enc_size = 0;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __double_encoded_array_size(NULL, 25);
    enc_size += 25 * __float_encoded_array_size(NULL, 24);
    enc_size += 25 * __float_encoded_array_size(NULL, 12);
    enc_size += 25 * __int32_t_encoded_array_size(NULL, 4);
    enc_size += 25 * __double_encoded_array_size(NULL, 4);
    enc_size += __float_encoded_array_size(NULL, 12);
    enc_size += 25 * 12 * __float_encoded_array_size(NULL, 12);
    enc_size += __float_encoded_array_size(NULL, 1);
    enc_size += 25 * __float_encoded_array_size(NULL, 12);
    enc_size += 25 * __float_encoded_array_size(NULL, 12);
    enc_size += 25 * __float_encoded_array_size(NULL, 6);
    enc_size += __float_encoded_array_size(NULL, 12);
    enc_size += __float_encoded_array_size(NULL, 3);
    return enc_size;
}

uint64_t hkd_command_lcmt::_computeHash(const __lcm_hash_ptr *)
{
    uint64_t hash = 0x8943aaf62325041bLL;
    return (hash<<1) + ((hash>>63)&1);
}

#endif
