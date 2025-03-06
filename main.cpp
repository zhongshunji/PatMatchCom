#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <array>
#include <chrono>  
#include "huffman.h"
#include <iomanip>

#include "Literal.h"


using namespace std;



// Data type for time series
#define DATA_TYPE float

//Maximum search times
#define MAX_SEARCH_TIMES 800

// Minimum match length: only when the match length ≥ MIN_MATCH, it is represented as a (length, index) pair; otherwise, it is output as a Literal
#define MIN_LENGTH 3

// Maximum match length
#define MAX_LENGTH  1938

// Maximum value for Index
#define MAX_INDEX 1048703

// Number of intervals for Literal and Length
#define L_INTE_NUM 29

// Number of intervals for Index
#define INDEX_INTE_NUM 19

// Length of the extra bits for Literal
#define LITERAL_EXTRA_BITS_LEN 32


// The starting position of each interval for Literal/Length
unsigned l_inte_starts[L_INTE_NUM] = { 0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,27,35,51,83,147,275,403,915 };
// The length of the extra bits for each interval for Literal/Length
unsigned l_inte_extra_bits_lens[L_INTE_NUM] = { LITERAL_EXTRA_BITS_LEN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10 };
// The starting position of each interval for Index
unsigned index_inte_starts[INDEX_INTE_NUM] = { 0,64,128,384,640,896,1152,2176,3200,4224,5248,6272,8320,16512, 32896, 65664,131200, 262272, 524416 };
// The length of the extra bits for each interval for Index
unsigned index_inte_extra_bits_lens[INDEX_INTE_NUM] = { 6,6,8,8,8,8,10,10,10,10,10,11,13,14,15,16,17,18,19 };


// The following declared functions are used during decompression
void decompress(const std::string& dictionary_file, const std::string& input_file, const std::string& decompressed_file);
std::vector<bool> byteVec_to_bitVec(const std::vector<unsigned char>& byteVec);
std::vector<int> decode_CL3(code_t& encodedBits, unsigned char CL3_MaxValue, unsigned char CL3Size);
std::vector<int> decode_CL1_or_CL2(const std::map<char, code_t>& code_table, code_t& encodedResult, unsigned char CLSize);
vector<std::pair<float, float>> decode_get_matchResults(ct_data s1, ct_data s2, code_t& encodedBits, const std::vector<float>& decompressed_literal);
int decode_interval(code_t::iterator& biter, const code_t::iterator& biter_end, ct_data s1, code_t& encodedVec);
int calculate_Offset(const std::vector<bool>& encodedVec, std::vector<bool>::iterator& biter, int extraBit_len);

std::vector<DATA_TYPE> decode_matchResults(const std::vector<std::pair<float, float>>& matchResults, const std::vector<DATA_TYPE>& D);
void write_data_to_csv(const std::string& csvFileName, const std::vector<DATA_TYPE>& rssiValues);



// In scientific notation, a value can be represented as: value = a*10^b, where the value of |a| is in the range [1,10), where a is called the mantissa and b is the exponent.
// The following function is to obtain the mantissa (a value) corresponding to value
float getMantissa(float value) {
    if (value == 0) return 0;

    while (fabs(value) >= 10) {
        value /= 10;
    }
    while (fabs(value) < 1) {
        value *= 10;
    }
    return value;
}


unordered_map<int, vector<int>> build_hash_bucket1(const vector<float>& D) {
    unordered_map<int, vector<int>> hash_bucket;
    int DSize = D.size();

    for (int i = 0; i < DSize - 2; i++) {

        int key = round((getMantissa(D[i]) + getMantissa(D[i + 1]) + getMantissa(D[i + 2])) * 10);

        hash_bucket[key].push_back(i);
    }

    return hash_bucket;
}


unordered_map<int, vector<int>> build_hash_bucket(const vector<float>& D) {
    unordered_map<int, vector<int>> hash_bucket;
    int DSize = D.size();

    if (DSize < 3) return hash_bucket;

    float d1 = getMantissa(D[0]);
    float d2 = getMantissa(D[1]);
    float d3 = getMantissa(D[2]);

    for (int i = 0; i < DSize - 2; i++) {

        int key = round((d1 + d2 + d3) * 10);
        hash_bucket[key].push_back(i);

        if (i + 3 < DSize) {
            d1 = d2;
            d2 = d3;
            d3 = getMantissa(D[i + 3]);
        }
    }

    return hash_bucket;
}


std::vector<DATA_TYPE> read_file(const std::string& csvFileName) {
    std::vector<DATA_TYPE> Values;
    FILE* fp = fopen(csvFileName.c_str(), "r");

    if (fp == NULL) {
        std::cerr << "Error: Unable to open the file '" << csvFileName << "'" << std::endl;
        return Values; 
    }

    DATA_TYPE value;
    while (fscanf(fp, "%f", &value) == 1) {  
        Values.push_back(value);
    }

    fclose(fp);
    return Values;
}


//Finding the longest matching pattern under the absolute error bound mode.
pair<int, int> find_longest_match_ABS(int i, const vector<DATA_TYPE>& D, const vector<DATA_TYPE>& T, unordered_map<int, vector<int>>& hash_bucket, float max_absolute_error) {

    int search_times = 0;
    int length, index;
    int best_length = 0, best_index = 0;
    vector<int> positions;
    int Tsize = T.size();
    int Dsize = D.size();
    
    // Handling out-of-bounds issues
    if (i >= Tsize - 2) {  
        return std::make_pair(0, 0);
    }

    int key = round((getMantissa(T[i]) + getMantissa(T[i + 1]) + getMantissa(T[i + 2])) * 10);  // Compute the hash key value at position j in T

    int key_threshold = getMantissa(T[i])*10 / T[i] * max_absolute_error * 3;

    for (int offset = 0; offset <= key_threshold; offset++) {
        if (hash_bucket.count(key + offset) > 0) {
            vector<int>& indices = hash_bucket[key + offset];
            for (int k = 0; k < indices.size() && search_times < MAX_SEARCH_TIMES; k++) {
                search_times++;
                index = indices[k];
                length = 0;
                while (i + length < Tsize && index + length < Dsize && abs(T[i + length] - D[index + length]) <= max_absolute_error && length < MAX_LENGTH) {
                    length++;
                }
                if (length > best_length) {
                    best_length = length;
                    best_index = index;
                }
            }
        }

        if (offset != 0 && hash_bucket.count(key - offset) > 0) {
            vector<int>& indices = hash_bucket[key - offset];
            for (int k = 0; k < indices.size() && search_times < MAX_SEARCH_TIMES; k++) {
                search_times++;
                index = indices[k];
                length = 0;

                while (i + length < Tsize && index + length < Dsize && abs(T[i + length] - D[index + length]) <= max_absolute_error && length < MAX_LENGTH) {
                    length++;
                }
                if (length > best_length) {
                    best_length = length;
                    best_index = index;
                }
            }
        }
    }

    return make_pair(best_length, best_index);
}


//Finding the longest matching pattern under the relative error bound mode.
pair<int, int> find_longest_match_REL(int i, const vector<DATA_TYPE>& D, const vector<DATA_TYPE>& T, const vector<DATA_TYPE>& T_mul_mre, unordered_map<int, vector<int>>& hash_bucket, float max_relative_error) {

    int search_times = 0;
    int length, index;
    int best_length = 0, best_index = 0;
    vector<int> positions;
    int Tsize = T.size();
    int Dsize = D.size();

    if (i >= Tsize - 2) {  
        return std::make_pair(0, 0);
    }

    int key = round((getMantissa(T[i]) + getMantissa(T[i + 1]) + getMantissa(T[i + 2])) * 10);  
    int key_threshold = abs(key) * max_relative_error;  

    for (int offset = 0; offset <= key_threshold; offset++) {
        if (hash_bucket.count(key + offset) > 0) {
            vector<int>& indices = hash_bucket[key + offset];
            for (int k = 0; k < indices.size() && search_times < MAX_SEARCH_TIMES; k++) {
                search_times++;
                index = indices[k];
                length = 0;
                
                while (i + length < Tsize && index + length < Dsize && abs(T[i + length] - D[index + length]) <= abs(T_mul_mre[i + length]) && length < MAX_LENGTH) {
                    length++;
                }
                if (length > best_length) {
                    best_length = length;
                    best_index = index;
                }
            }
        }

        if (offset != 0 && hash_bucket.count(key - offset) > 0) {
            vector<int>& indices = hash_bucket[key - offset];
            for (int k = 0; k < indices.size() && search_times < MAX_SEARCH_TIMES; k++) {
                search_times++;
                index = indices[k];
                length = 0;

                while (i + length < Tsize && index + length < Dsize && abs(T[i + length] - D[index + length]) <= abs(T_mul_mre[i + length]) && length < MAX_LENGTH) {
                    length++;
                }
                if (length > best_length) {
                    best_length = length;
                    best_index = index;
                }
            }
        }
    }

    return make_pair(best_length, best_index);
}



size_t calculate_num_bytes(size_t value) {
    size_t num_bytes = 0;
    do {
        num_bytes++;
        value >>= 8; 
    } while (value > 0);
    return num_bytes;
}


void compress(const std::string& dictionary_file, const std::string& input_file, const std::string& compressed_file, const std::string& error_bound_mode,float max_error) {
    
    
    vector<DATA_TYPE> D = read_file(dictionary_file);
    if (D.size() > MAX_INDEX) {
        D.resize(MAX_INDEX);
    }

    unordered_map<int, vector<int>> hash_bucket = build_hash_bucket(D);
   
    vector<DATA_TYPE> T = read_file(input_file);

    vector<int> l_inte_count(L_INTE_NUM, 0); 
    vector<int> index_inte_count(INDEX_INTE_NUM, 0); 
    
    vector<std::pair<float, float>> matchResults;//Save the results of Pattern Matching Encoding
    vector<float> literal_data;
    int l_interval;//literal or length interval number
    int index_interval;//index interval number

    // Pattern Matching Encoding
    if (error_bound_mode == "A") {  // absolute error bound mode
        float max_absolute_error = max_error;
        for (int j = 0; j < T.size(); j++) {
            
            pair<int, int> match = find_longest_match_ABS(j, D, T, hash_bucket, max_absolute_error);
            if (match.first >= MIN_LENGTH) {
                matchResults.push_back(make_pair(match.first, match.second));//The result of pattern matching is (length,index), match.first is length, match.second is index
                l_interval = get_length_interval(match.first);//Get the interval number of length
                index_interval = get_index_interval(match.second);//Get the interval number of index
                l_inte_count[l_interval]++;//The frequency of the interval number of length is increased by 1
                index_inte_count[index_interval]++;//The frequency of the interval number of index is increased by 1
                j += match.first - 1;
            }
            else {
                l_interval = 0;//The interval number of literal is 0
                matchResults.push_back(make_pair(l_interval, T[j]));//The result of pattern matching is literal, save the result as (0,literal)
                literal_data.push_back(T[j]);
                l_inte_count[0]++;//The number of times the literal in the 0 interval is increased by 1
            }
        }
    }
    else if (error_bound_mode == "R") {  // relative error bound mode
        float max_relative_error = max_error;

        std::vector<float> T_mul_mre(T.size());
        
        std::transform(T.begin(), T.end(), T_mul_mre.begin(), [max_relative_error](float x) {
            return x * max_relative_error;
            });

        for (int j = 0; j < T.size(); j++) {
            
            pair<int, int> match = find_longest_match_REL(j, D, T, T_mul_mre, hash_bucket, max_relative_error);
            if (match.first >= MIN_LENGTH) {
                matchResults.push_back(make_pair(match.first, match.second));
                l_interval = get_length_interval(match.first);
                index_interval = get_index_interval(match.second);
                l_inte_count[l_interval]++; 
                index_inte_count[index_interval]++;
                j += match.first - 1;
                
            }
            else {
                l_interval = 0;
                matchResults.push_back(make_pair(l_interval, T[j]));
                literal_data.push_back(T[j]);
                l_inte_count[0]++;  
            }
        }
    }

    
    
    //========================== Pattern-Matched-Results Encoding======================================
    ct_data s1;
    ct_data s2;
    init_ct_data(s1, L_INTE_NUM, l_inte_starts, l_inte_extra_bits_lens);
    init_ct_data(s2, INDEX_INTE_NUM, index_inte_starts, index_inte_extra_bits_lens);

    s1.inte_freq_table = make_freq_table(l_inte_count);// Count the occurrence frequency of each interval of literal/length
    s2.inte_freq_table = make_freq_table(index_inte_count);// Count the occurrence frequency of each interval of index

    s1.huff_tree = build_huff_tree(s1.inte_freq_table); // Build the Huffman tree for literal/length
    s2.huff_tree = build_huff_tree(s2.inte_freq_table); // Build the Huffman tree for index

    // Get the Huffman codebook
    s1.codebook = get_huff_codebook(s1.huff_tree);
    s2.codebook = get_huff_codebook(s2.huff_tree);

    // Get the codeword length sequence
    s1.CL = get_code_length(s1.codebook);
    s2.CL = get_code_length(s2.codebook);

    // Construct a canonical Huffman tree based on the codeword length sequence
    s1.codebook = get_can_huff_codebook(s1.CL);
    s2.codebook = get_can_huff_codebook(s2.CL);

    std::vector<bool> encodedBits;
    std::vector<std::pair<char, unsigned>> CL3_freq_table = make_CL_freq_table(s1.CL, s2.CL);
    HuffmanTree* CL3_htree = build_huff_tree(CL3_freq_table);
    codetable CL3_code_table = get_huff_codebook(CL3_htree);
    std::vector<int> CL3 = get_code_length(CL3_code_table);
    CL3_code_table = get_can_huff_codebook(CL3);

    unsigned char CL3_MaxValue = find_max_value(CL3);// Find the maximum value in CL3

    encode_CL3(CL3, CL3_MaxValue, encodedBits);// Encode the codeword length sequence CL3

    encode_CL1_or_CL2(s1.CL, CL3_code_table, encodedBits);// Encode the codeword length sequence CL1
    encode_CL1_or_CL2(s2.CL, CL3_code_table, encodedBits);// Encode the codeword length sequence CL2

    encode_matchResults(matchResults, s1, s2, encodedBits);  // Encode the pattern matching result (Length,Index) pairs or Literals

    vector<unsigned char> encodedBytes = bitVec_to_byteVec(encodedBits);    // Convert the encoded bits to bytes

    //======================Write the compressed data to the compressed file==========================================
    std::ofstream binaryFile(compressed_file, std::ios::out | std::ios::binary);

    if (!binaryFile) {
        std::cerr << "Failed to open the binary file for writing" << std::endl;
    }

    std::vector<uint8_t> literal_compressed_data = Approximation_XOR_Brotli_Encoding(literal_data, error_bound_mode, max_error);
    size_t literal_compressed_data_size = literal_compressed_data.size();  //literal_compressed_data_size表示literal_compressed_data的大小（字节）
    unsigned char num_bytes = calculate_num_bytes(literal_compressed_data_size);

    binaryFile.write(reinterpret_cast<const char*>(&num_bytes), 1);
    binaryFile.write(reinterpret_cast<const char*>(&literal_compressed_data_size), num_bytes);
    binaryFile.write(reinterpret_cast<const char*>(literal_compressed_data.data()), literal_compressed_data_size);

    unsigned char CL3Size = CL3.size();
    unsigned char CL2Size = s2.CL.size();
    unsigned char CL1Size = s1.CL.size();
    binaryFile.write(reinterpret_cast<const char*>(&CL3_MaxValue), 1);
    binaryFile.write(reinterpret_cast<const char*>(&CL3Size), 1);
    binaryFile.write(reinterpret_cast<const char*>(&CL2Size), 1);
    binaryFile.write(reinterpret_cast<const char*>(&CL1Size), 1);

    binaryFile.write(reinterpret_cast<const char*>(encodedBytes.data()), encodedBytes.size());// Write the data in encodedBytes to the file

    binaryFile.close();

}


// ============================================The following functions are used for decompression======================================================================================

// Decompression function
void decompress(const std::string& dictionary_file, const std::string& input_file, const std::string& decompressed_file) {

    std::ifstream binaryFile(input_file, std::ios::in | std::ios::binary);
    if (!binaryFile) {
        std::cerr << "Failed to open the binary file for reading" << std::endl;
    }

    unsigned char num_bytes;
    size_t literal_compressed_data_size=0;
    std::vector<uint8_t> literal_compressed_data;

    unsigned char CL3_MaxValue;
    unsigned char CL3Size, CL2Size, CL1Size;

    binaryFile.read(reinterpret_cast<char*>(&num_bytes), 1);
    binaryFile.read(reinterpret_cast<char*>(&literal_compressed_data_size), num_bytes);


    literal_compressed_data.resize(literal_compressed_data_size);
    binaryFile.read(reinterpret_cast<char*>(literal_compressed_data.data()), literal_compressed_data_size);
    std::vector<float> decompressed_literal = decompress_literal(literal_compressed_data);


    binaryFile.read(reinterpret_cast<char*>(&CL3_MaxValue), 1);
    binaryFile.read(reinterpret_cast<char*>(&CL3Size), 1);
    binaryFile.read(reinterpret_cast<char*>(&CL2Size), 1);
    binaryFile.read(reinterpret_cast<char*>(&CL1Size), 1);

    std::vector<unsigned char> encodedBytes;
    unsigned char temp_byte;
    while (binaryFile.read(reinterpret_cast<char*>(&temp_byte), 1)) {
        encodedBytes.push_back(temp_byte);
    }
    binaryFile.close();

    std::vector<bool> encodedBits = byteVec_to_bitVec(encodedBytes);

    std::vector<int> CL3 = decode_CL3(encodedBits, CL3_MaxValue, CL3Size);
    codetable CL3_code_table = get_can_huff_codebook(CL3);

    ct_data s1;
    ct_data s2;
    init_ct_data(s1, L_INTE_NUM, l_inte_starts, l_inte_extra_bits_lens);
    init_ct_data(s2, INDEX_INTE_NUM, index_inte_starts, index_inte_extra_bits_lens);

    s1.CL = decode_CL1_or_CL2(CL3_code_table, encodedBits, CL1Size);
    s2.CL = decode_CL1_or_CL2(CL3_code_table, encodedBits, CL2Size);
    s1.codebook = get_can_huff_codebook(s1.CL);
    s2.codebook = get_can_huff_codebook(s2.CL);
    s1.huff_tree = build_can_huff_tree_from_codebook(s1.codebook);
    s2.huff_tree = build_can_huff_tree_from_codebook(s2.codebook);

    vector<std::pair<float, float>> matchResults = decode_get_matchResults(s1, s2, encodedBits, decompressed_literal);

    vector<DATA_TYPE> D = read_file(dictionary_file);

    vector<DATA_TYPE> decoded_T = decode_matchResults(matchResults, D);

    write_data_to_csv(decompressed_file, decoded_T);

}


// This function converts encoded bytes to encoded bits
std::vector<bool> byteVec_to_bitVec(const std::vector<unsigned char>& encodedBytes) {
    std::vector<bool> encodedBits;

    unsigned char hangBits = encodedBytes[0];//The first byte represents the number of hanging bits
    encodedBits.reserve((encodedBytes.size() - 1) * 8); 

    for (size_t i = 1; i < encodedBytes.size(); i++) {
        unsigned char byte = encodedBytes[i];
        for (int boff = 0; boff < 8; boff++) {
            bool bit = (byte >> boff) & 1;
            encodedBits.push_back(bit);
        }
    }

    if (hangBits > 0) {
        encodedBits.resize(encodedBits.size() - (8 - hangBits));
    }

    return encodedBits;
}


// This function decodes the codeword  length sequence of CL3
std::vector<int> decode_CL3(code_t& encodedBits, unsigned char CL3_MaxValue, unsigned char CL3Size) {
    std::vector<int> CL3;

    int numBits = 0; // Calculate the number of bits for each value through CL3_MaxValue
    while (CL3_MaxValue > 0) {
        CL3_MaxValue >>= 1;
        numBits++;
    }

    int currentPosition = 0; 
    int value = 0; 

    for (int i = 0; i < CL3Size; ++i) {
        for (int j = numBits - 1; j >= 0; --j) {
            bool bit = encodedBits[currentPosition];
            value |= (bit << j);
            currentPosition++;
        }
        CL3.push_back(value);
        value = 0;
    }

    encodedBits.erase(encodedBits.begin(), encodedBits.begin() + currentPosition);
    return CL3;
}


// This function decodes the codeword length sequence CL1 or CL2
std::vector<int> decode_CL1_or_CL2(const std::map<char, code_t>& code_table, code_t& encodedBits, unsigned char CLSize) {
    std::vector<int> decodedCL;

    code_t currentCode;
    size_t currentPosition = 0;

    for (bool bit : encodedBits) {
        currentCode.push_back(bit);

        for (const auto& entry : code_table) {
            const code_t& code = entry.second;
            if (currentCode == code) {
                decodedCL.push_back(static_cast<int>(entry.first));
                currentCode.clear();
                break;
            }
        }

        currentPosition++;

        if (decodedCL.size() == CLSize) {
            encodedBits.erase(encodedBits.begin(), encodedBits.begin() + currentPosition);
            break;
        }
    }

    return decodedCL;
}


vector<std::pair<float, float>> decode_get_matchResults(ct_data s1, ct_data s2, code_t& encodedBits, const std::vector<float>& decompressed_literal) {
    vector<std::pair<float, float>> matchResults;//Save the decoded pattern matching results

    int l, index;//l represents the value of literal or length
    float literal;
    int l_interval;//The interval number of literal/length
    int index_interval;//The interval number of index
    int extraBit_len;//The length of the extra bits
    int inte_start;//The starting position of the interval
    int offset; // offset represents the offset of literal / length / index relative to the interval

    code_t::iterator biter = encodedBits.begin();
    size_t literal_index = 0;
    while (true) {
        try {
            l_interval = decode_interval(biter, encodedBits.end(), s1, encodedBits);//获得Literal或Length的区间
            if (l_interval == 0) {// Pattern matching result is a Literal
                
                literal = decompressed_literal[literal_index++];
                matchResults.push_back(make_pair(0, literal));
            }
            else {// Pattern matching result is (Length,Index) pair
                // Decoding Length
                inte_start = s1.inte_starts[l_interval];//Interval starting position
                extraBit_len = s1.inte_extra_bits_lens[l_interval];// The length of the extra bits for the interval
                offset = calculate_Offset(encodedBits, biter, extraBit_len);
                l = inte_start + offset;//literal or length value = interval starting position + offset relative to the interval

                // Decoding Index
                index_interval = decode_interval(biter, encodedBits.end(), s2, encodedBits);
                inte_start = s2.inte_starts[index_interval];
                extraBit_len = s2.inte_extra_bits_lens[index_interval];
                offset = calculate_Offset(encodedBits, biter, extraBit_len);
                index = inte_start + offset;

                matchResults.push_back(make_pair(l, index));

            }

        }
        catch (const std::out_of_range& oor) {
            break;
        }
    }
    return matchResults;
}


// This function decodes the interval number of literal, length or interval
int decode_interval(code_t::iterator& biter, const code_t::iterator& biter_end, ct_data s, code_t& encodedVec) {
    const HuffmanTree* node = s.huff_tree;

    while (true) {
        if (!node->left) {
            break;
        }
        if (biter == biter_end) {
            throw std::out_of_range("No more bits");
        }
        if (*biter) {
            node = node->right;
        }
        else {
            node = node->left;
        }
        ++biter;
    }

    return node->c;
}


// This function calculates the offset of literal/length/index relative to the starting position of the interval
int calculate_Offset(const std::vector<bool>& encodedVec, std::vector<bool>::iterator& biter, int extraBit_len) {
    int offset = 0; // offset represents the offset of literal/length/index relative to the interval

    for (int i = 0; i < extraBit_len; i++) {
        bool bit = *biter;
        offset = offset * 2 + bit;
        ++biter;
    }

    return offset;
}


// This function decodes the pattern matching result (length,index) or literal to decompressed data
std::vector<DATA_TYPE> decode_matchResults(const std::vector<std::pair<float, float>>& matchResults, const std::vector<DATA_TYPE>& D) {
    std::vector<DATA_TYPE> decoded_T;
    int length, index;

    for (const auto& match : matchResults) {

        if (match.first == 0) {
            float literal = match.second;
            decoded_T.push_back(literal);

        }
        else {
            length = match.first;
            index = match.second;
            for (int i = 0; i < length; ++i) {
                decoded_T.push_back(D[index + i]);
            }

        }
    }

    return decoded_T;
}



//Write data into a CSV file, writing to 1 column
void write_data_to_csv(const std::string& csvFileName, const std::vector<DATA_TYPE>& data) {
    std::ofstream outputFile(csvFileName);

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file '" << csvFileName << "' for writing" << std::endl;
        return;
    }

    for (const DATA_TYPE value : data) {
        outputFile << value << std::endl;
    }

    outputFile.close();
}


std::streamsize getFileSize(const std::string& filename) {
    std::ifstream file(filename, std::ifstream::ate | std::ifstream::binary);

    if (!file.is_open()) {
        return -1;
    }
    std::streamsize size = file.tellg();
    file.close();

    return size;
}


//Calculate the size of the CSV file, each point occupies 4 bytes
std::streamsize getCSVDataSize(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return -1;
    }

    std::string line;
    size_t numElements = 0;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            if (!value.empty()) {
                numElements++;
            }
        }
    }

    file.close();

    return static_cast<std::streamsize>(numElements * 4);
}



// Calculate the compression ratio (CR)
double calculate_compression_ratio(const std::string& originalFile, const std::string& compressedFile) {
    std::streamsize originalSize = getCSVDataSize(originalFile);

    std::streamsize compressedSize = getFileSize(compressedFile);

    if (originalSize == -1 || compressedSize == -1 || compressedSize == 0) {
        std::cerr << "Unable to calculate compression ratio because the file cannot be opened or the compressed file size is 0." << std::endl;
        return -1.0;
    }
    std::cout << "Original file size: " << originalSize << " bytes" << std::endl;
    std::cout << "Compressed file size: " << compressedSize << " bytes" << std::endl;

    return static_cast<double>(originalSize) / static_cast<double>(compressedSize);
}



int countDataPoints(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return -1; 
    }

    int count = 0;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            count++;
        }
    }

    file.close();
    return count;
}

void countDataPointsInFiles(const std::string& history_file, const std::string& dict_file) {
    int history_count = countDataPoints(history_file);
    int dict_count = countDataPoints(dict_file);

    if (history_count >= 0 && dict_count >= 0) {
        std::cout << "Number of data points in history file: " << history_count << std::endl;
        std::cout << "Number of data points in dictionary file: " << dict_count << std::endl;
    }
    else {
        std::cerr << "Failed to count data points in one or both files." << std::endl;
    }
}

int main(int argc, char* argv[]) {

    auto start = std::chrono::high_resolution_clock::now();

    if (argc < 3) {
        std::cerr << "Parameter input error! Please check the input format." << std::endl;
        return 1;
    }

    std::vector<std::string> args(argv, argv + argc);
      
    //-c indicates compression
    if (args[1] == "-c" && argc == 6) { 
        std::string dict_file = argv[2];
        std::string input_file = argv[3]; 

        std::string error_bound_mode = argv[4]; //Error bound mode: when error_bound_mode='A', it represents the absolute error bound, when error_bound_mode='R', it represents the relative error bound
        float error_bound = std::stof(argv[5]);// The maximum error, including the maximum absolute error and the maximum relative error
        
        std::string compressed_file = input_file.substr(0, input_file.rfind('.')) + ".bin";

        compress(dict_file, input_file, compressed_file, error_bound_mode, error_bound);

        std::cout << "Compression successful, the compressed file: " << compressed_file << std::endl;
        double compressionRatio = calculate_compression_ratio(input_file, compressed_file);
        if (compressionRatio != -1.0) {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "Compression ratio: " << compressionRatio << std::endl;
        }
    }
    //-d indicates decompression
    else if (args[1] == "-d" && argc == 4) { 
        std::string dict_file = argv[2]; 
        std::string input_file = argv[3]; 
        
        std::string decompressed_file = input_file.substr(0, input_file.rfind('.')) + "_decompressed" + ".csv";
        decompress(dict_file, input_file, decompressed_file);
        std::cout << "Decompression successful, the decompressed file: " << decompressed_file << std::endl;
    }

    else {
        std::cerr << "Parameter input error! Please check the input format." << std::endl;
        return 1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Total program execution time: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}

