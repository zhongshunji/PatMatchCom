#ifndef LITERAL_HPP
#define LITERAL_HPP

#include <vector>
#include <cstdint>

std::vector<uint8_t> Approximation_XOR_Brotli_Encoding(const std::vector<float>& Literal, const std::string& error_bound_mode, float max_error);
std::vector<float> decompress_literal(const std::vector<uint8_t>& compressed_data);

#endif // LITERAL_HPP
