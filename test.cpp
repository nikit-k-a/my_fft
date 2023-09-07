#include "fft.h"

#include <iostream>
#include <complex>
#include <vector>
#include <ctime>
#include <algorithm> // max_element
#include <numeric> // accumulate

// Fills vector with random complex double values.
void generate_data (std::vector<fft::complexD> &data);

// Simple visual test of the FFT. Used mainly for debug purposes.
void test0_visual_debug (const std::vector<fft::complexD> &vec);

/**
 * Tests FFT error.
 *
 * Create vector of random complex double values. Performs FFT and IFFT, then
 * averages maximum and average error.
 *
 * @param [in] input_data_sz vector length.
 * @param [in] averaging.
 */
void test1 (const uint32_t input_data_sz, const size_t averaging = 1);

// Calculates the average of the vector vec
double average(const std::vector<double> &vec);

int main() {

    test0_visual_debug({0, 1, 2, 3});
    test0_visual_debug({1, 2, 3, 4, 0});
    // test0_visual({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
    // test0_visual({fft::complexD( 6, 0), fft::complexD(-2,  2),
    //               fft::complexD(-2, 0), fft::complexD(-2, -2)});

    test1(5, 100);
    test1(16, 100);
    test1(333, 100);
    test1(1000, 100);
    test1(4096, 100);
    test1(65536, 100);

    return 0;
}

void generate_data (std::vector<fft::complexD> &data) {
    std::srand(unsigned(std::time(nullptr)));
    std::generate(data.begin(), data.end(),
                  []()->fft::complexD
                  {
                      return fft::complexD(std::rand(), std::rand());
                  }
                 );
}

void test1 (const uint32_t input_data_sz, const size_t averaging) {
    std::cout << "-----Test----\n";
    std::cout << "Data size = " << input_data_sz << "\n";
    std::cout << "Averaging = " << averaging << "\n";
    std::vector<double> avg_errors;//(averaging);
    std::vector<double> max_errors;//(averaging);

    for (size_t i = 0; i < averaging; i++) {
        std::vector<fft::complexD> input_data (input_data_sz);
        generate_data(input_data);

        // copy input data to check error later
        std::vector<fft::complexD> check (input_data);


        fft::MyFFT test_fft (input_data);

        // FFT
        test_fft.fft_func();

        //Inderse FFT
        test_fft.ifft_func();
        std::vector<fft::complexD> after_ifft_vec = test_fft.get_data();

        // Calculating error
        std::vector<double> diff (input_data_sz);
        for (int i = 0; i < input_data_sz; i++) {
            diff[i] = static_cast<double>(abs(after_ifft_vec[i] - check[i]));
        }

        // Maximum error
        auto max_error = max_element(std::begin(diff), std::end(diff));
        max_errors.push_back(*max_error);
        // Average error
        avg_errors.push_back(average(diff));
    }

    std::cout << "Average max error = " << average(max_errors) << "\n";
    std::cout << "Average average error = " << average(avg_errors) << "\n";
    std::cout << "-----End Test-----\n";
}

void test0_visual_debug(const std::vector<fft::complexD> &vec) {
    std::cout << "-----Test_visual----\n";
    std::vector<fft::complexD> input_data(vec);

    fft::MyFFT test_fft (input_data);

    std::cout << "Input padded data:\n";
    for (auto val: input_data){
        std::cout << val << "\n";
    }
    std::cout << "--------------------\n";

    // FFT
    test_fft.fft_func();
    std::vector<fft::complexD> after_fft_vec = test_fft.get_data();

    std::cout << "After FFT:\n";
    for (auto val: after_fft_vec) {
        std::cout << val << "\n";
    }
    std::cout << "--------------------\n";

    // Inverse FFT
    test_fft.ifft_func();
    std::vector<fft::complexD> after_ifft_vec = test_fft.get_data();
    std::cout << "After IFFT:\n";
    for (auto val: after_ifft_vec) {
        std::cout << val << "\n";
    }
    std::cout << "----------End Test_visual----------\n";
}

double average(const std::vector<double> &vec){
    if(vec.empty()){
        return 0;
    }

    auto const count = static_cast<double>(vec.size());
    return accumulate(vec.begin(), vec.end(), 0.0) / count;
}
