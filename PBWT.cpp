#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>

/**
 * PBWT Implementation in C++
 the Richard Durbin's algorithm for positional Burrows-Wheeler transform.
 */
class PBWT {
public:
    struct Match {
        int h1, h2, start, end;
        Match(int _h1, int _h2, int _start, int _end) 
            : h1(_h1), h2(_h2), start(_start), end(_end) {}
    };

    int N, M;
    std::vector<std::vector<char> > X; 
    std::vector<std::vector<int> > a;
    std::vector<std::vector<int> > d;

    PBWT(const std::vector<std::vector<char> >& matrix) : X(matrix) {
        N = X.size();
        M = (N > 0) ? X[0].size() : 0;
    }

    void build() {
        if (N == 0 || M == 0) return;

        std::vector<int> a_curr(N);
        std::vector<int> d_curr(N, 0);
        for (int i = 0; i < N; ++i) a_curr[i] = i;

        a.push_back(a_curr);
        d.push_back(d_curr);

        for (int k = 0; k < M; ++k) {
            std::vector<int> zeros_idx, ones_idx;
            std::vector<int> zeros_div, ones_div;
            int p = k + 1, q = k + 1;

            for (int i = 0; i < N; ++i) {
                int idx = a_curr[i];
                int div = d_curr[i];
                if (div > p) p = div;
                if (div > q) q = div;

                if (X[idx][k] == 0 || X[idx][k] == '0') {
                    zeros_idx.push_back(idx);
                    zeros_div.push_back(p);
                    p = 0;
                } else {
                    ones_idx.push_back(idx);
                    ones_div.push_back(q);
                    q = 0;
                }
            }

            std::vector<int> a_next = zeros_idx;
            a_next.insert(a_next.end(), ones_idx.begin(), ones_idx.end());
            std::vector<int> d_next = zeros_div;
            d_next.insert(d_next.end(), ones_div.begin(), ones_div.end());

            a.push_back(a_next);
            d.push_back(d_next);
            a_curr = a_next;
            d_curr = d_next;
        }
    }

    std::vector<Match> reportMatches(int minLen) {
        std::vector<Match> matches;
        for (int k = 1; k <= M; ++k) {
            for (int i = 1; i < N; ++i) {
                if (d[k][i] <= k - minLen) {
                    matches.push_back(Match(a[k][i], a[k][i - 1], d[k][i], k));
                }
            }
        }
        return matches;
    }
};

std::vector<std::vector<char> > parseCSV(const std::string& path) {
    std::vector<std::vector<char> > matrix;
    std::ifstream file(path.c_str());
    if (!file) throw std::runtime_error("file is missing there" + path);

    std::string line;
    bool header = true;
    while (std::getline(file, line)) {
        if (header) { header = false; continue; }
        std::vector<char> row;
        std::stringstream ss(line);
        std::string val;
        bool first = true;
        while (std::getline(ss, val, ',')) {
            if (first) { first = false; continue; }
            row.push_back((char)std::atoi(val.c_str()));
        }
        if (!row.empty()) matrix.push_back(row);
    }
    return matrix;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "PBWT Haplotype Matcher" << std::endl;
        std::cout << "Usage: " << argv[0] << " input.csv min_length" << std::endl;
        return 1;
    }

    std::string inputPath = argv[1];
    int minL = std::atoi(argv[2]);

    try {
        std::cout << "Step - 1 " << inputPath << "..." << std::endl;
        std::vector<std::vector<char> > matrix = parseCSV(inputPath);
        
        PBWT pbwt(matrix);
        std::cout << "Step - 2 Building PBWT (" << pbwt.N << " haps, " << pbwt.M << " sites)..." << std::endl;
        pbwt.build();

        std::cout << "Step - 3 Analyzing matches (min length: " << minL << ")..." << std::endl;
        std::vector<PBWT::Match> matches = pbwt.reportMatches(minL);

        if (!matches.empty()) {
            std::string outName = "matched_haplotypes.csv";
            std::ofstream out(outName.c_str());
            out << "hap1,hap2,start,end,length\n";
            for (size_t i = 0; i < matches.size(); ++i) {
                out << matches[i].h1 << "," << matches[i].h2 << "," << matches[i].start 
                    << "," << matches[i].end << "," << (matches[i].end - matches[i].start) << "\n";
            }
            std::cout << "The Algorithm Found " << matches.size() << " matches. Results: " << outName << std::endl;
        } else {
            std::cout << "No matches found." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}