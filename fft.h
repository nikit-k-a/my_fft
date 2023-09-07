#ifndef MYFFT_H

#define MYFFT_H

#include <complex>
#include <vector>

namespace fft {

typedef std::complex<double> complexD;

/**
 * Fast Fourier Transform (FFT) class.
 *
 * Perfofms forward and inverse FFT.
 * Uses iterative radix-2 Cooley–Tukey algorithm with bit-reversal permutation.
 * Input data with size which is not equal to the power of two,
 * is padded with zeros. The data size should be 32-bit. Input data is invalid
 * afterwards, it is moved.
 */
class MyFFT {
public:
    explicit MyFFT(const std::vector<complexD> &vec);

    MyFFT(const MyFFT&) = delete;
    MyFFT& operator=(const MyFFT&) = delete;

    MyFFT(MyFFT&&) = default;
    MyFFT& operator=(MyFFT&&) = delete;

    // Computes FFT.
    void fft_func();

    // Computes inverse FFT.
    void ifft_func();

    // returns data_
    std::vector<complexD> get_data() const {return data_;}

private:
    // data vector size
    std::vector<complexD> data_;
    const uint32_t vec_sz_;

    /**
     * Bit reversal permutation.
     *
     * Mapping each item to the item whose
     * representation has the same bits in the reversed order.
     */
    void bit_reverse_swaps ();

    /**
     * Reversing bits in a 32-bit word.
     *
     * Adapted from: https://graphics.stanford.edu/~seander/bithacks.html
     *
     * @param [in] index index to modify.
     * @param [in] power log2 of the padded vector size.
     * @return new index, has same bits in reverse order.
     */
    uint32_t bit_reverse_index(const uint32_t index, const unsigned int power) const;

    /**
     * Calculates the next power of two.
     *
     * @param value
     * @return next (closest greater) power of two to the value.
     */
    uint32_t roundup_power2 (const uint32_t value) const;

    /**
     * Padding the data with zeros.
     *
     * If the data size is not equal to the power of two,
     * it is padded up to the next power of two with zeros to perform
     * radix-2 Cooley–Tukey algorithm. Radix-2 version of this algorithm is
     * higly efficient.
     *
     * @return size of the padded data vector.
     */
    uint32_t zero_padding ();

    /**
     * Complex conjugation of the all elements of the vector.
     */
    void idft();
};

}
#endif /* end of include guard: MYFFT_H */
