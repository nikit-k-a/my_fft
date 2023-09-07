#include "fft.h"

#include <algorithm> // std::transform
#include <cassert>

namespace fft {

constexpr uint32_t BIT_LENGTH = 32;

MyFFT::MyFFT (const std::vector<complexD> &vec):
    data_(std::move(vec)),
    vec_sz_(zero_padding())
    {}

void MyFFT::bit_reverse_swaps () {
    uint32_t new_index = 0;

    const unsigned int power = std::log2(vec_sz_);

    for (uint32_t i = 0; i < vec_sz_; i++) {
        new_index = bit_reverse_index(i, power);

        if (new_index > i) {
            assert(new_index < vec_sz_);
            std::swap(data_[i], data_[new_index]);
        }
    }
}

uint32_t MyFFT::bit_reverse_index(const uint32_t index, const unsigned int power) const {
    uint32_t reversed = index;
    //swap odd and even bits
    reversed = ((reversed >> 1) & 0x55555555) | ((reversed & 0x55555555) << 1);
    // swap consecutive pairs
    reversed = ((reversed >> 2) & 0x33333333) | ((reversed & 0x33333333) << 2);
    // swap nibbles
    reversed = ((reversed >> 4) & 0x0F0F0F0F) | ((reversed & 0x0F0F0F0F) << 4);
    //swap bytes
    reversed = ((reversed >> 8) & 0x00FF00FF) | ((reversed & 0x00FF00FF) << 8);
    //swap 2-byte long pairs
    reversed = (( reversed >> 16) | ( reversed << 16)) >> (BIT_LENGTH - power);

    return reversed;
}

void MyFFT::fft_func() {
    MyFFT::bit_reverse_swaps();
    for (int i = 2; i <= vec_sz_; i <<= 1) {

        double angle = -2 * M_PI / i;
        complexD root(cos(angle), sin(angle));

        for (int j = 0; j < vec_sz_; j += i) {
            complexD tmp(1, 0);
            for (int k = j; k < j + i/2; k++) {
                assert(k + i/2 < vec_sz_);
                complexD twid = tmp * data_[k + i/2];
                complexD left = data_[k];

                data_[k] = left + twid;
                data_[k + i/2] = left - twid;

                tmp = tmp * root;
            }
        }
    }
}

void MyFFT::ifft_func() {
    // complex conjugation
    idft();

    // FFT
    fft_func();

    // normalization
    std::transform(data_.begin(), data_.end(), data_.begin(),
              [vec_sz_by_val = vec_sz_](complexD &c){
                  return complexD(c.real()/vec_sz_by_val, c.imag()/vec_sz_by_val);
              });

    // complex conjugation
    idft();
}

uint32_t MyFFT::roundup_power2 (const uint32_t value) const {
    uint32_t tmp = value;

    tmp--;
    tmp |= tmp >> 1;
    tmp |= tmp >> 2;
    tmp |= tmp >> 4;
    tmp |= tmp >> 8;
    tmp |= tmp >> 16;
    tmp++;

    return tmp;
}

uint32_t MyFFT::zero_padding () {
    const uint32_t data_sz = data_.size();
    assert(data_sz > 0);
    const uint32_t next_power_of2 = roundup_power2(data_sz);

    if (next_power_of2 - data_sz > 0)
        data_.resize(next_power_of2, complexD (0, 0));

    return next_power_of2;
}

void MyFFT::idft() {
    for (auto& x : data_)
    x = std::conj(x);
}

} // namespace fft
