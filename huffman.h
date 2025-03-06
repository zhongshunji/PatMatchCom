#ifndef HUFFMAN_H
#define HUFFMAN_H

#include <vector>
#include <map>
#include <string>
#include <array>
#include <unordered_map>
using namespace std;

typedef std::vector<bool> code_t;
typedef std::map<char, code_t> codetable;


// A Huffman Tree Node
struct HuffmanTree {
    char c; // character in an alphabet
    int cfreq; // frequency of c.
    struct HuffmanTree* left;
    struct HuffmanTree* right;
    HuffmanTree(char c, int cfreq, struct HuffmanTree* left = NULL, struct HuffmanTree* right = NULL) :
        c(c), cfreq(cfreq), left(left), right(right) {}

    ~HuffmanTree() {
        delete left, delete right;
    }
    class Compare {
    public:
        bool operator()(HuffmanTree* a, HuffmanTree* b) {
            return a->cfreq > b->cfreq;
        }
    };
};


typedef struct ct_data_s {//这里的"ct"是"code tree"的缩写
    unsigned inte_num;//Number of intervals
    vector< pair<char, unsigned> > inte_freq_table;  //Frequency of each interval
    HuffmanTree* huff_tree;// Huffman tree obtained by normal Huffman coding for the interval
    codetable codebook;//Codebook
    std::vector<int> CL;//codeword length sequence (CL) 
    unsigned* inte_starts;//the starting position of each interval
    unsigned* inte_extra_bits_lens;//the length of the extra bits for each interval
}ct_data;



int get_length_interval(int _Length);
int get_index_interval(int _Index);

void init_ct_data(ct_data& s, unsigned num, unsigned* start, unsigned* extraBitLen);
vector< pair<char, unsigned> > make_freq_table(const vector<int>& counts);
HuffmanTree* build_huff_tree(std::vector<std::pair<char, unsigned>>& alph);
codetable get_huff_codebook(HuffmanTree* htree);
std::vector<int> get_code_length(const codetable& code_table);
std::vector<std::pair<char, unsigned>> make_CL_freq_table(const std::vector<int>& CL1, const std::vector<int>& CL2);
int find_max_value(const std::vector<int>& CL3);
void encode_CL3(const std::vector<int>& CL3, int maxValue, code_t& result);
void encode_CL1_or_CL2(const std::vector<int>& CL, const std::map<char, code_t>& code_table, code_t& result);
float binary_vec32_to_float(const std::vector<bool>& binaryVec);
void encode_matchResults(vector<std::pair<float, float>> matchResults, ct_data s1, ct_data s2, code_t& result);

std::vector<unsigned char> bitVec_to_byteVec(std::vector<bool>& bitvec);


/*==============================================================以下是范式哈夫曼编码的内容====================================================================================*/

codetable get_can_huff_codebook(const std::vector<int>& code_length);
HuffmanTree* build_can_huff_tree_from_codebook(const codetable& codeTable);


#endif  // HUFFMAN_H
