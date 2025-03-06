#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <brotli/encode.h>
#include <brotli/decode.h>
#include "Literal.h" 


// Binary search approximation function
float binary_approximation_with_error(float num, float e, const std::string& error_bound_mode) {

    float lower_bound = std::floor(num);
    float upper_bound = lower_bound + 1;

    // If num is within the error bounds of the lower or upper bound, return directly
    if ((error_bound_mode == "A" && std::abs(num - lower_bound) <= e) ||
        (error_bound_mode == "R" && std::abs(num - lower_bound) <= e * std::abs(num))) {
        return lower_bound;
    }
    if ((error_bound_mode == "A" && std::abs(num - upper_bound) <= e) ||
        (error_bound_mode == "R" && std::abs(num - upper_bound) <= e * std::abs(num))) {
        return upper_bound;
    }

    // Use binary search to find the best alternative value
    while (true) {
        float mid = (lower_bound + upper_bound) / 2.0f; // Calculate the midpoint of the lower and upper bounds
        
        if ((error_bound_mode == "A" && std::abs(num - mid) <= e) ||
            (error_bound_mode == "R" && std::abs(num - mid) <= e * std::abs(num))) {
            return mid;
        }

        // Adjust the search range based on the error
        if (num < mid) {
            upper_bound = mid;
        }
        else {
            lower_bound = mid;
        }
    }
}


std::vector<float> approximate_representation(const std::vector<float>& data, const std::string& error_bound_mode, float max_error) {
    std::vector<float> approximated_data;
    approximated_data.reserve(data.size());

    // The first value is directly approximated using binary search
    approximated_data.push_back(binary_approximation_with_error(data[0], max_error, error_bound_mode));

    for (size_t i = 1; i < data.size(); ++i) {
        float previous_value = approximated_data.back();
        float error = std::abs(data[i] - previous_value);

        // Determine if the previous value can be used for approximation
        if ((error_bound_mode == "A" && error <= max_error) ||
            (error_bound_mode == "R" && error <= max_error * std::abs(data[i]))) {
            approximated_data.push_back(previous_value);
        }
        else {
            // Otherwise, use binary search to approximate
            approximated_data.push_back(binary_approximation_with_error(data[i], max_error, error_bound_mode));
        }
    }

    return approximated_data;
}


union FloatBits {
    float f;
    uint32_t bits;
};


std::vector<float> xor_with_previous(const std::vector<float>& data) {
    if (data.empty()) {
        return std::vector<float>();
    }

    std::vector<float> result;
    result.reserve(data.size());

    // The first value remains unchanged
    result.push_back(data[0]);

    FloatBits current, previous;

    // Start from the second value, XOR each value with the previous value
    for (size_t i = 1; i < data.size(); ++i) {
        current.f = data[i];
        previous.f = data[i - 1];

        uint32_t xor_result = current.bits ^ previous.bits;

        FloatBits result_union;
        result_union.bits = xor_result;

        result.push_back(result_union.f);
    }

    return result;
}



bool is_little_endian() {
    uint16_t test_value = 0x1;
    return *reinterpret_cast<uint8_t*>(&test_value) == 0x1;
}


std::vector<uint8_t> vertical_concat_little_endian(const std::vector<float>& xor_result) {
    if (xor_result.empty()) {
        return {};
    }

    // Get the current machine's byte order
    bool little_endian = is_little_endian();

    const size_t FLOAT_BYTES = sizeof(float);

    std::vector<uint8_t> concatenated_bytes;
    concatenated_bytes.reserve(xor_result.size() * FLOAT_BYTES);

    std::vector<uint8_t> temp_buffer(xor_result.size() * FLOAT_BYTES);

    for (size_t i = 0; i < xor_result.size(); ++i) {
        std::memcpy(&temp_buffer[i * FLOAT_BYTES], &xor_result[i], FLOAT_BYTES);
    }

    // If the machine is big-endian, adjust the byte order to little-endian
    if (!little_endian) {
        for (size_t i = 0; i < xor_result.size(); ++i) {
            std::reverse(temp_buffer.begin() + i * FLOAT_BYTES, temp_buffer.begin() + (i + 1) * FLOAT_BYTES);
        }
    }

    for (size_t byte_index = 0; byte_index < FLOAT_BYTES; ++byte_index) {
        for (size_t float_index = 0; float_index < xor_result.size(); ++float_index) {
            concatenated_bytes.push_back(temp_buffer[float_index * FLOAT_BYTES + byte_index]);
        }
    }

    return concatenated_bytes;
}


std::vector<uint8_t> compress_with_brotli(const std::vector<uint8_t>& input) {
    size_t compressed_size = BrotliEncoderMaxCompressedSize(input.size());
    std::vector<uint8_t> compressed_data(compressed_size);

    BrotliEncoderCompress(BROTLI_DEFAULT_QUALITY, BROTLI_DEFAULT_WINDOW, BROTLI_MODE_GENERIC,
        input.size(), input.data(),
        &compressed_size, compressed_data.data());

    compressed_data.resize(compressed_size);
    return compressed_data;
}


std::vector<uint8_t> Approximation_XOR_Brotli_Encoding(const std::vector<float>& Literal, const std::string& error_bound_mode, float max_error) {

    if (Literal.empty()) {
        return std::vector<uint8_t>{};
    }

    // Approximate representation
    std::vector<float> approximate_Literal = approximate_representation(Literal, error_bound_mode, max_error);

    // XOR operation
    std::vector<float> xor_result = xor_with_previous(approximate_Literal);

    // Vertical Concatenation in Little-Endian Byte Order
    std::vector<uint8_t> concatenated_bytes = vertical_concat_little_endian(xor_result);

    // Brotli compression
    return compress_with_brotli(concatenated_bytes);
}


std::vector<uint8_t> decompress_with_brotli(const std::vector<uint8_t>& compressed_data) {
    BrotliDecoderState* decoder_state = BrotliDecoderCreateInstance(nullptr, nullptr, nullptr);
    if (!decoder_state) {
        throw std::runtime_error("Failed to create Brotli decoder instance");
    }

    std::vector<uint8_t> decompressed_data;
    size_t buffer_size = 4096;
    std::vector<uint8_t> buffer(buffer_size);

    const uint8_t* next_in = compressed_data.data();
    size_t available_in = compressed_data.size();

    BrotliDecoderResult result;
    do {
        uint8_t* next_out = buffer.data();
        size_t available_out = buffer_size;

        result = BrotliDecoderDecompressStream(
            decoder_state,
            &available_in, &next_in,
            &available_out, &next_out,
            nullptr);

        size_t bytes_decompressed = buffer_size - available_out;

        if (bytes_decompressed > 0) {
            decompressed_data.insert(
                decompressed_data.end(),
                buffer.begin(),
                buffer.begin() + bytes_decompressed
            );
        }

        if (result == BROTLI_DECODER_RESULT_ERROR) {
            BrotliDecoderDestroyInstance(decoder_state);
            throw std::runtime_error("Brotli decompression error");
        }
    } while (result == BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT);

    BrotliDecoderDestroyInstance(decoder_state);

    return decompressed_data;
}


std::vector<float> restore_floats_from_bytes(const std::vector<uint8_t>& interleaved_bytes) {
    const size_t FLOAT_BYTES = sizeof(float);

    if (interleaved_bytes.size() % FLOAT_BYTES != 0) {
        throw std::runtime_error("Invalid interleaved bytes size");
    }

    size_t float_count = interleaved_bytes.size() / FLOAT_BYTES;
    std::vector<float> restored_floats(float_count);

    // Get the current machine's byte order
    bool little_endian = []() {
        uint16_t test_value = 0x1;
        return *reinterpret_cast<uint8_t*>(&test_value) == 0x1;
    }();

    for (size_t i = 0; i < float_count; ++i) {
        uint8_t temp[FLOAT_BYTES] = { 0 };
        for (size_t byte_index = 0; byte_index < FLOAT_BYTES; ++byte_index) {
            temp[byte_index] = interleaved_bytes[i + byte_index * float_count];
        }

        if (!little_endian) {
            std::reverse(std::begin(temp), std::end(temp));
        }

        std::memcpy(&restored_floats[i], temp, FLOAT_BYTES);
    }

    return restored_floats;
}


std::vector<float> decompress_literal(const std::vector<uint8_t>& compressed_data) {
    std::vector<uint8_t> interleaved_bytes = decompress_with_brotli(compressed_data);

    std::vector<float> restored_floats = restore_floats_from_bytes(interleaved_bytes);

    std::vector<float> restored_literals;
    restored_literals.reserve(restored_floats.size());

    if (!restored_floats.empty()) {
        restored_literals.push_back(restored_floats[0]);
        for (size_t i = 1; i < restored_floats.size(); ++i) {
            FloatBits current, previous;
            current.f = restored_floats[i];
            previous.f = restored_literals.back();

            FloatBits original;
            original.bits = current.bits ^ previous.bits;

            restored_literals.push_back(original.f);
        }
    }

    return restored_literals;
}
