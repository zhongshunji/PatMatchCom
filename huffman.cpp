#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <string>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include"huffman.h"
#include <cstdint>


//Get the interval of length
int get_length_interval(int _Length) {
    if (_Length == 3) {
        return 1;
    }
    else if (_Length == 4) {
        return 2;
    }
    else if (_Length == 5) {
        return 3;
    }
    else if (_Length == 6) {
        return 4;
    }
    else if (_Length == 7) {
        return 5;
    }
    else if (_Length == 8) {
        return 6;
    }
    else if (_Length == 9) {
        return 7;
    }
    else if (_Length == 10) {
        return 8;
    }
    else if (_Length == 11) {
        return 9;
    }
    else if (_Length == 12) {
        return 10;
    }
    else if (_Length == 13) {
        return 11;
    }
    else if (_Length == 14) {
        return 12;
    }
    else if (_Length == 15) {
        return 13;
    }
    else if (_Length == 16) {
        return 14;
    }
    else if (_Length == 17) {
        return 15;
    }
    else if (_Length == 18) {
        return 16;
    }
    else if (_Length == 19) {
        return 17;
    }
    else if (_Length == 20) {
        return 18;
    }
    else if (21 <= _Length && _Length <= 22) {
        return 19;
    }
    else if (23 <= _Length && _Length <= 26) {
        return 20;
    }
    else if (27 <= _Length && _Length <= 34) {
        return 21;
    }
    else if (35 <= _Length && _Length <= 50) {
        return 22;
    }
    else if (51 <= _Length && _Length <= 82) {
        return 23;
    }
    else if (83 <= _Length && _Length <= 146) {
        return 24;
    }
    else if (147 <= _Length && _Length <= 274) {
        return 25;
    }
    else if (275 <= _Length && _Length <= 402) {
        return 26;
    }
    else if (403 <= _Length && _Length <= 914) {
        return 27;
    }
    else if (915 <= _Length && _Length <= 1938) {
        return 28;
    }
    else{
        return -1;
    }
}



//Get the interval of Index
int get_index_interval(int _Index) {
    if (0 <= _Index && _Index <= 63) {
        return 0;
    }
    else if (64 <= _Index && _Index <= 127) {
        return 1;
    }
    else if (128 <= _Index && _Index <= 383) {
        return 2;
    }
    else if (384 <= _Index && _Index <= 639) {
        return 3;
    }
    else if (640 <= _Index && _Index <= 895) {
        return 4;
    }
    else if (896 <= _Index && _Index <= 1151) {
        return 5;
    }
    else if (1152 <= _Index && _Index <= 2175) {
        return 6;
    }
    else if (2176 <= _Index && _Index <= 3199) {
        return 7;
    }
    else if (3200 <= _Index && _Index <= 4223) {
        return 8;
    }
    else if (4224 <= _Index && _Index <= 5247) {
        return 9;
    }
    else if (5248 <= _Index && _Index <= 6271) {
        return 10;
    }
    else if (6272 <= _Index && _Index <= 8319) {
        return 11;
    }
    else if (8320 <= _Index && _Index <= 16511) {
        return 12;
    }
    else if (16512 <= _Index && _Index <= 32895) {
        return 13;
    }
    else if (32896 <= _Index && _Index <= 65663) {
        return 14;
    }
    else if (65664 <= _Index && _Index <= 131199) {
        return 15;
    }
    else if (131200 <= _Index && _Index <= 262271) {
        return 16;
    }
    else if (262272 <= _Index && _Index <= 524415) {
        return 17;
    }
    else if (524416 <= _Index && _Index <= 1048703) {
        return 18;
    }
    else{
        return -1;
    }
}


//Initialize the code tree structure
void init_ct_data(ct_data& s, unsigned num, unsigned* inte_starts, unsigned* inte_extra_bits_lens) {
    s.inte_num = num; //Initialize the number of intervals
    s.inte_starts = inte_starts;//Initialize the starting position of each interval
    s.inte_extra_bits_lens = inte_extra_bits_lens;//The length of the extra bits of each interval
}



//This function counts the occurrences of each interval
vector< pair<char, unsigned> > make_freq_table(const vector<int>& counts) {
    vector< pair<char, unsigned> > cfvec;

    for (int i = 0; i < counts.size(); ++i) {
        if (counts[i] != 0) {
            cfvec.push_back(make_pair(i, static_cast<unsigned>(counts[i])));
        }
    }

    return cfvec;
}


HuffmanTree* build_huff_tree(vector< pair<char, unsigned> >& freq_table) {
    priority_queue<HuffmanTree*, vector<HuffmanTree*>, HuffmanTree::Compare > alph_heap;

    for (vector< pair<char, unsigned> >::iterator it = freq_table.begin(); it != freq_table.end(); ++it) {
        HuffmanTree* leaf = new HuffmanTree(it->first, it->second);
        alph_heap.push(leaf);
    }

    //If there is only one node in the frequency table, special processing
    if (alph_heap.size() == 1) {
        HuffmanTree* onlyNode = alph_heap.top();
        alph_heap.pop();

        //Create a virtual parent node (weight is the frequency of the original node), the left node is onlyNode, and the right node is nullptr (null pointer)
        HuffmanTree* root = new HuffmanTree(0, onlyNode->cfreq, onlyNode, nullptr);
        return root;
    }

    //HuffmanTree algorithm: Merge the two leaf nodes with the smallest weight until only one node (root node) is left
    HuffmanTree* root = NULL;
    while (alph_heap.size() > 1) {
        HuffmanTree* l = alph_heap.top();
        alph_heap.pop();
        HuffmanTree* r = alph_heap.top();
        alph_heap.pop();

        //Create a new branch node, the weight is the sum of the weights of the two child nodes, and the character is set to the empty character
        root = new HuffmanTree(0, l->cfreq + r->cfreq, l, r);
        alph_heap.push(root);
    }

    return root;
}


//This function obtains the corresponding Huffman codebook according to the Huffman tree
map<char, code_t> get_huff_codebook(HuffmanTree* htree) {
    codetable codebook; 

    // If the Huffman tree is empty, return an empty codebook directly.
    if (!htree) {
        return codebook;
    }

    // If there is only one node (the left child of the virtual parent node is the unique character node), set the encoding of the unique character node to "0"
    if (htree->left && !htree->right) {
        codebook.insert(make_pair(htree->left->c, code_t(1, 0))); // Assign the encoding "0" to the unique character node
        return codebook;
    }

    deque< pair<HuffmanTree*, code_t> > q;

    q.push_back(make_pair(htree, code_t()));
    while (!q.empty()) {
        HuffmanTree* node, * lc, * rc;
        code_t code;
        node = q.front().first; 
        code = q.front().second;
        q.pop_front();  
        lc = node->left;
        rc = node->right;

        if (lc) {
            code_t code_cp(code);
            q.push_back(make_pair(lc, (code.push_back(0), code)));
            q.push_back(make_pair(rc, (code_cp.push_back(1), code_cp)));
        }
        else {
            codebook.insert(make_pair(node->c, code));
        }
    }

    return codebook;
}


//This function generates the corresponding code length (Code Length,CL) sequence according to the Huffman codebook
std::vector<int> get_code_length(const codetable& code_table) {//typedef std::map<char, code_t> codetable;
    char maxChar = 0; 

    for (const auto& entry : code_table) {
        char currentChar = entry.first;
        if (currentChar > maxChar) {
            maxChar = currentChar;
        }
    }
    int max_value = static_cast<int>(maxChar);

    std::vector<int> code_length(max_value+1, 0);

    unsigned inte;
    for (const auto& entry : code_table) {
        inte = entry.first;
        int code_len = entry.second.size();
        code_length[inte] = code_len;
    }
    return code_length;
}




//This function counts the frequency of each value in the code length sequence CL1 and CL2
std::vector<std::pair<char, unsigned>> make_CL_freq_table(const std::vector<int>& CL1, const std::vector<int>& CL2) {
    std::unordered_map<char, unsigned> CL_count;

    for (int num : CL1) {
        char key = static_cast<char>(num); 
        CL_count[key]++;
    }

    for (int num : CL2) {
        char key = static_cast<char>(num); 
        CL_count[key]++;
    }

    std::vector<std::pair<char, unsigned>> result;
    for (const auto& entry : CL_count) {
        result.push_back(std::make_pair(entry.first, entry.second));
    }

    return result;
}



//Find the maximum value in CL3
int find_max_value(const std::vector<int>& CL3) {
    if (CL3.empty()) {
        return 0;
    }

    int maxVal = CL3[0];
    for (int value : CL3) {
        if (value > maxVal) {
            maxVal = value;
        }
    }

    return maxVal;
}



//This function encodes the code length sequence CL3, the encoding rule is: using fixed-length bit encoding, determine how many bits each value needs according to the maximum value in CL3
void encode_CL3(const std::vector<int>& CL3, int CL3_MaxValue, code_t& encodedBits) {
    int numBits = 0;
    while (CL3_MaxValue > 0) {
        CL3_MaxValue >>= 1;
        numBits++;
    }

    for (int value : CL3) {
        for (int i = numBits - 1; i >= 0; --i) {
            bool bit = (value >> i) & 1;
            encodedBits.push_back(bit);
        }
    }
}



//This function encodes the code length sequence CL1 and CL2 using canonical Huffman encoding
void encode_CL1_or_CL2(const std::vector<int>& CL, const std::map<char, code_t>& code_table, code_t& encodedBits) {
    for (int value : CL) {

        auto it = code_table.find(static_cast<char>(value));

        if (it != code_table.end()) {
            const code_t& code = it->second;
            encodedBits.insert(encodedBits.end(), code.begin(), code.end());
        }
    }
}


std::vector<bool> get_extra_bits(int num, int start, int bits_len) {

    int offset = num - start;
    std::vector<bool> bits_vec(bits_len, false); 

    for (int i = bits_len - 1; i >= 0; i--) {
        bits_vec[i] = (offset & 1) == 1;
        offset >>= 1;
    }

    return bits_vec;
}



#include <iostream>
#include <bitset>
#include <vector>

union FloatUnion {
    float f;
    uint32_t i;
};

//Convert a float number to its 32-bit binary representation
std::vector<bool> float_to_binary_vec32(float value) {
    FloatUnion fu;
    fu.f = value;

    std::bitset<32> bits(fu.i);
    std::vector<bool> binaryVec(32);

    for (int i = 0; i < 32; ++i) {
        binaryVec[i] = bits[i];
    }
    return binaryVec;
}


//Convert a 32-bit binary vector (std::vector<bool>) to the corresponding float type floating-point number.
float binary_vec32_to_float(const std::vector<bool>& binaryVec) {
    std::bitset<32> bits;
    for (int i = 0; i < 32; ++i) {
        bits[i] = binaryVec[i];
    }
   
    FloatUnion fu;
    fu.i = bits.to_ulong();

    return fu.f;
}


//This function encodes the result of pattern matching (length,index) or literal using the canonical Huffman tree. 
//The length and literal are encoded using the same canonical Huffman tree, and the index is encoded using another canonical Huffman tree.
//Encoded bits = Interval Code + Extra Bits
void encode_matchResults(vector<std::pair<float, float>> matchResults, ct_data s1, ct_data s2, code_t& encodedBits) {
    int length, index;
    float literal;
    int l_inte;//  Interval number of literal/length
    int index_inte;  //Interval number of index

    int l_extra_bits_len, index_extra_bits_len;//Extra bits length
    int l_inte_start, index_inte_start;
    std::vector<bool> l_extra_bits, index_extra_bits;//Extra bits
    std::vector<bool> l_inte_code, index_inte_code;//Interval code

    for (const auto& match : matchResults) {
        if (match.first == 0) { //If match.first == 0, it means this is a tuple (0,literal)
            literal = match.second;
            l_inte = 0;//literal is in the 0th interval

            l_inte_code = s1.codebook[l_inte];//literal's interval code

            encodedBits.insert(encodedBits.end(), l_inte_code.begin(), l_inte_code.end());

        }
        else {//If the match result is (length,index)
            length = match.first;
            index = match.second;

            l_inte = get_length_interval(length);//Get the interval where length is located
            index_inte = get_index_interval(index);//Get the interval where index is located

            l_inte_start = s1.inte_starts[l_inte];//The starting position of the interval of length
            index_inte_start = s2.inte_starts[index_inte];//The starting position of the interval of index

            l_extra_bits_len = s1.inte_extra_bits_lens[l_inte];//The length of the extra bits of length
            index_extra_bits_len = s2.inte_extra_bits_lens[index_inte];//The length of the extra bits of index

            l_inte_code = s1.codebook[l_inte];//The interval code of length
            index_inte_code = s2.codebook[index_inte];//The interval code of index
            l_extra_bits = get_extra_bits(length, l_inte_start, l_extra_bits_len);//The extra bits of length
            index_extra_bits = get_extra_bits(index, index_inte_start, index_extra_bits_len);//The extra bits of index

            encodedBits.insert(encodedBits.end(), l_inte_code.begin(), l_inte_code.end());
            encodedBits.insert(encodedBits.end(), l_extra_bits.begin(), l_extra_bits.end());
            encodedBits.insert(encodedBits.end(), index_inte_code.begin(), index_inte_code.end());
            encodedBits.insert(encodedBits.end(), index_extra_bits.begin(), index_extra_bits.end());

        }      
    }

}


vector<unsigned char> bitVec_to_byteVec(std::vector<bool>& encodedBits) {
    vector<unsigned char> encodedBytes;
    unsigned char hangBits;
    hangBits = encodedBits.size() & 7;

    encodedBytes.push_back(static_cast<unsigned char>(hangBits));

    unsigned char byte = 0;
    for (unsigned i = 0; i < encodedBits.size(); i++) {
        unsigned boff = i % 8;
        byte |= encodedBits[i] << boff;
        if (boff == 7) {
            encodedBytes.push_back(byte);
            byte = 0;
        }
    }
    if (hangBits) {
        encodedBytes.push_back(byte);
    }

    return encodedBytes;
}





//======The following is the content of canonical Huffman encoding, only the code length (Code Length,CL) sequence is needed to build the canonical Huffman code table====================
static std::unordered_map<int, int> count_code_length(const std::vector<int>& code_length) {
    std::unordered_map<int, int> cl_count;

    int max_value = *std::max_element(code_length.begin(), code_length.end());

    for (int i = 1; i <= max_value; ++i) {
        cl_count[i] = 0;
    }

    for (int value : code_length) {
        if (value != 0) {
            cl_count[value]++;
        }
    }

    return cl_count;
}


static std::unordered_map<int, int> find_smallest_code(std::unordered_map<int, int>& bl_count) {
    int code = 0;
    bl_count[0] = 0;
    std::unordered_map<int, int> next_code;
    int MAX_BITS = 0;

    for (const auto& pair : bl_count) {
        if (pair.first > MAX_BITS) {
            MAX_BITS = pair.first;
        }
    }

    for (int bits = 1; bits <= MAX_BITS; bits++) {
        code = (code + bl_count[bits - 1]) << 1;
        next_code[bits] = code;
    }

    return next_code;
}


//This function builds the canonical Huffman codebook according to the code length sequence
codetable get_can_huff_codebook(const std::vector<int>& code_length) {
    std::unordered_map<int, int> cl_count = count_code_length(code_length);
    std::unordered_map<int, int> next_code = find_smallest_code(cl_count);

    codetable code_book;

    for (int i = 0; i < code_length.size(); ++i) {
        int len = code_length[i];

        if (len != 0) {
            int code = next_code[len];
            code_t code_vector(len, false);

            for (int i = 0; i < len; ++i) {
                code_vector[len - i - 1] = (code & (1 << i)) != 0;
            }

            code_book[i] = code_vector;
            next_code[len]++;
        }
    }

    return code_book;
}


//This function builds the canonical Huffman tree according to the given code table,
HuffmanTree* build_can_huff_tree_from_codebook(const codetable& codebook) {
    HuffmanTree* root = new HuffmanTree(0, 0, nullptr, nullptr);

    for (const auto& entry : codebook) {
        char character = entry.first;
        const code_t& code = entry.second;

        HuffmanTree* currentNode = root;

        for (bool bit : code) {
            if (bit == 0) {
                if (!currentNode->left) {
                    currentNode->left = new HuffmanTree(0, 0, nullptr, nullptr);
                }
                currentNode = currentNode->left;
            }
            else {
                if (!currentNode->right) {
                    currentNode->right = new HuffmanTree(0, 0, nullptr, nullptr);
                }
                currentNode = currentNode->right;
            }
        }

        currentNode->c = character;
    }

    return root;
}


