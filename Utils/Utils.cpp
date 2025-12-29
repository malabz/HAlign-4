#include "Utils.hpp"
#include "Pseudo.hpp"
#include "Fasta.hpp"
#include "Arguments.hpp"

#include <stdio.h>
#include <regex>
#include <iostream>
#include <limits>
#include <iomanip>
#include <list>
#include <fstream>
#include <cstring>
#if defined(_WIN32)
#include <io.h>
#include <direct.h>
#elif defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/stat.h>
#include <unistd.h>
#endif
char chars[8] = { 'N','A','C','G','T','N','N','-' };
int my_mk_dir(std::string output_dir)
{
#if defined(_WIN32)
    if (0 != access(output_dir.c_str(), 0))
    {
        return mkdir(output_dir.c_str());
    }
    else
        return 0;
#elif defined(__unix__) || defined(__unix) || defined(unix)
    if (access(output_dir.c_str(), F_OK) == -1)
    {
        return mkdir(output_dir.c_str(), S_IRWXU | S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    }
    else
        return 0;
#else
    return -1;
#endif
    return -1;
}
void cout_cur_time()
{
    auto now = std::chrono::system_clock::now();
    uint64_t dis_millseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count()
                               - std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count() * 1000;
    time_t tt = std::chrono::system_clock::to_time_t(now);
    auto time_tm = localtime(&tt);
    char strTime[25] = { 0 };
    sprintf(strTime, "%d-%02d-%02d %02d:%02d:%02d | ", time_tm->tm_year + 1900,
            time_tm->tm_mon + 1, time_tm->tm_mday, time_tm->tm_hour,
            time_tm->tm_min, time_tm->tm_sec);
    std::cout << strTime;
}

size_t getPeakRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L;        /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t)0L;        /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif
#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;            /* Unsupported. */
#endif
}

#if defined(_WIN32)
void getFiles_win(std::string path, std::vector<std::string>& files)
{

    intptr_t hFile = 0;

    struct _finddata_t fileinfo;
    std::string p;
    if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
    {
        do
        {

            if ((fileinfo.attrib & _A_SUBDIR))
            {
                if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
                    getFiles_win(p.assign(path).append("\\").append(fileinfo.name), files);
            }
            else
            {
                files.emplace_back(p.assign(path).append(fileinfo.name));
            }
        } while (_findnext(hFile, &fileinfo) == 0);
        _findclose(hFile);
    }
}
#elif defined(__unix__) || defined(__unix) || defined(unix)
void getFiles_linux(std::string path, std::vector<std::string>& filenames)
{
    DIR* pDir = NULL;
    struct dirent* ptr=NULL;
    if (!(pDir = opendir(path.c_str()))) {
        std::cout << "Folder doesn't Exist!" << std::endl;
        return;
    }
    while ((ptr = readdir(pDir)) != 0) {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0) {
            filenames.emplace_back(path + ptr->d_name);
        }
    }
    closedir(pDir);
}
#endif

int NSCORE = 0;
int HOXD70[6][6] = { {},{0,91,-114,-31,-123,NSCORE},{0,-114,100,-125,-31,NSCORE},{0,-31,-125,100,-114,NSCORE},
                     {0,-123,-31,-114,91,NSCORE},{0,NSCORE,NSCORE,NSCORE,NSCORE,NSCORE} };
int cs[8] = { 0,91,100,100,91,0,0,0 };
int d = 400, e = 30;
int stop_g = 5;
std::vector<unsigned char> ACGT = { nucleic_acid_pseudo::A ,nucleic_acid_pseudo::C ,nucleic_acid_pseudo::G ,nucleic_acid_pseudo::T };

std::string utils::remove_white_spaces(const std::string& str)
{
    static const std::regex white_spaces("\\s+");
    return std::regex_replace(str, white_spaces, "");
}

std::vector<unsigned char> utils::to_pseudo(const std::string& str)
{
    std::vector<unsigned char> pseu;
    pseu.reserve((size_t)(str.size() * 1.3));
    unsigned char c;
    for (size_t i=0;i<str.size();i++)
    {
        c = to_pseudo(str[i]);
        if (c != nucleic_acid_pseudo::N)
            pseu.emplace_back(c);
        else
            pseu.emplace_back(ACGT[i % 4]);
    }
    return pseu;
}

std::string utils::from_pseudo(const std::vector<unsigned char>& pseu)
{
    static constexpr char map[nucleic_acid_pseudo::NUMBER]{ '-', 'c', 'g', 'a', 't', 'n' };

    std::string str;
    str.reserve(pseu.size());

    for (auto i : pseu) str.push_back(map[i]);
    return str;
}

unsigned char* _get_map()
{
    using namespace nucleic_acid_pseudo;

    static unsigned char map[std::numeric_limits<unsigned char>::max()];
    memset(map, N, sizeof(map));

    // map['-'] = GAP; // we could not process sequences with '-'
    map['c'] = map['C'] = C;
    map['g'] = map['G'] = G;
    map['a'] = map['A'] = A;
    map['t'] = map['T'] = map['u'] = map['U'] = T;
    map['N'] = map['R'] = map['Y'] = map['M'] = map['K'] = map['S'] = map['W'] = map['H'] = map['B'] = map['V'] = map['D'] = N;
    map['n'] = map['r'] = map['y'] = map['m'] = map['k'] = map['s'] = map['w'] = map['h'] = map['b'] = map['v'] = map['d'] = N;
    return map;
}

static const unsigned char* _map = _get_map();

unsigned char utils::to_pseudo(char ch)
{
    return _map[ch];
}

std::vector<std::vector<unsigned char>> utils::read_to_pseudo(std::istream& is, std::string& center_name, int& II, int& center_)
{
    std::vector<std::vector<unsigned char>> sequences;

    std::string each_line;
    std::string each_sequence;
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13))
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_line.erase(std::remove_if(each_line.begin(), each_line.end(), [](char c) {
                return (c == ' ' || c == '\t');
            }), each_line.end());
            if (center_ == -1 && each_line == center_name)
                center_ = II;
            II++;
            if (flag)
            {
                sequences.emplace_back(to_pseudo(each_sequence));
                each_sequence.clear();
            }
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }

    sequences.emplace_back(to_pseudo(each_sequence));
    return sequences;
}

// add for viral msa --- star---
void utils::read_single_fasta_sequence(std::istream& in,
                                       std::string& id,
                                       std::vector<unsigned char>& sequence)
{
    std::string line, seq;
    bool found_id = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (found_id) break;  // 只读取一条序列
            id = line.substr(1);  // 去掉 '>'
            found_id = true;
        } else {
            seq += line;
        }
    }
    sequence.assign(seq.begin(), seq.end());
}

std::vector<std::vector<unsigned char>> utils::read_to_pseudo_for_viral(std::istream& is)
{
    std::vector<std::vector<unsigned char>> sequences;

    std::string each_line;
    std::string each_sequence;
    bool flag = false;

    while (std::getline(is, each_line))
    {
        // 跳过空行或只有回车
        if (each_line.empty() || (each_line.size() == 1 && (int)each_line[0] == 13))
            continue;

        if (each_line[0] == '>')
        {
            // 如果前面已经读入过一条序列，则保存
            if (flag && !each_sequence.empty())
            {
                std::vector<unsigned char> seq_uc;
                seq_uc.reserve(each_sequence.size());
                for (char ch : each_sequence)
                    seq_uc.emplace_back(std::toupper(static_cast<unsigned char>(ch)));

                sequences.emplace_back(std::move(seq_uc));
                each_sequence.clear();
            }
            flag = true;  // 表示开始读取一条新的序列
        }
        else if (flag)
        {
            // 累积序列内容（去除行末回车）
            if (!each_line.empty() && (int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_sequence += each_line;
        }
    }

    // 处理最后一条序列
    if (!each_sequence.empty())
    {
        std::vector<unsigned char> seq_uc;
        seq_uc.reserve(each_sequence.size());
        for (char ch : each_sequence)
            seq_uc.emplace_back(std::toupper(static_cast<unsigned char>(ch)));

        sequences.emplace_back(std::move(seq_uc));
    }

    return sequences;
}

void utils::write_to_fasta_for_vrial_msa(std::ostream& os, std::istream& is, std::vector<std::vector<unsigned char>>& aligned_sequences)
{
    std::string line;
    std::vector<std::string> ids;

    // 提取 ID（以 '>' 开头的行）
    for (; std::getline(is, line); ) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            ids.push_back(line);
        }
    }

    // 检查数量是否一致
    if (ids.size() != aligned_sequences.size()) {
        std::cerr << "ERROR: Number of IDs (" << ids.size()
                  << ") does not match number of aligned sequences ("
                  << aligned_sequences.size() << ").\n";
        return;
    }

    // 输出每条序列，支持 60 字符换行
    for (size_t i = 0; i < ids.size(); ++i) {
        os << ids[i] << "\n";

        const auto& seq = aligned_sequences[i];
        size_t len = seq.size();
        for (size_t pos = 0; pos < len; pos += 60) {
            size_t chunk_size = std::min((size_t)60, len - pos);
            for (size_t j = 0; j < chunk_size; ++j)
                os << seq[pos + j];
            os << "\n";
        }
    }
}

// bool utils::read_batch_of_n_sequences_with_ids(std::istream& in,
//     std::vector<std::pair<std::string, std::vector<unsigned char>>>& batch,
//     size_t N) {
//
//     batch.clear();
//     std::string line, seq, current_id;
//     while (std::getline(in, line)) {
//         if (line.empty()) continue;
//         if (line[0] == '>') {
//             if (!seq.empty()) {
//                 batch.emplace_back(current_id, std::vector<unsigned char>(seq.begin(), seq.end()));
//                 seq.clear();
//                 if (batch.size() >= N) return true;
//             }
//             current_id = line.substr(1);  // 去掉 >
//         } else {
//             seq += line;
//         }
//     }
//     if (!seq.empty()) {
//         batch.emplace_back(current_id, std::vector<unsigned char>(seq.begin(), seq.end()));
//     }
//     return !batch.empty();
// }

// bool utils::read_batch_of_n_sequences_with_ids(
//     std::istream& in,
//     std::vector<std::pair<std::string, std::vector<unsigned char>>>& batch,
//     size_t N)
// {
//     batch.clear();
//     std::string line, seq, current_id;

//     while (std::getline(in, line)) {
//         // 去除末尾可能的 '\r'
//         if (!line.empty() && line.back() == '\r') {
//             line.pop_back();
//         }

//         if (line.empty()) continue;

//         if (line[0] == '>') {
//             if (!seq.empty()) {
//                 // 转为大写并存储
//                 std::vector<unsigned char> seq_vec;
//                 seq_vec.reserve(seq.size());
//                 for (char c : seq) {
//                     seq_vec.emplace_back(std::toupper(static_cast<unsigned char>(c)));
//                 }

//                 batch.emplace_back(current_id, std::move(seq_vec));
//                 if (current_id.empty()) {
//                     std::cerr << "[DEBUG] Empty current_id 1. \n";
//                 }
//                 seq.clear();

//                 if (batch.size() >= N) return true;
//             }
//             current_id = line.substr(1);  // 去掉 '>'
//         } else {
//             seq += line;
//         }
//     }

//     // 文件最后一条序列也别漏掉
//     if (!seq.empty()) {
//         std::vector<unsigned char> seq_vec;
//         seq_vec.reserve(seq.size());
//         for (char c : seq) {
//             seq_vec.emplace_back(std::toupper(static_cast<unsigned char>(c)));
//         }

//         batch.emplace_back(current_id, std::move(seq_vec));
//         if (current_id.empty()) {
//                     std::cerr << "[DEBUG] Empty current_id 2. \n";
//         }
//     }

//     return !batch.empty();
// }

bool utils::read_batch_of_n_sequences_with_ids(
    std::istream& in,
    std::vector<std::pair<std::string, std::vector<unsigned char>>>& batch,
    size_t N)
{
    batch.clear();
    std::string line, seq, current_id;
    std::streampos last_pos;  // <== 记录上一个位置

    while (true) {
        last_pos = in.tellg();
        if (!std::getline(in, line)) break;

        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                std::vector<unsigned char> seq_vec;
                seq_vec.reserve(seq.size());
                for (char c : seq)
                    seq_vec.emplace_back(std::toupper(static_cast<unsigned char>(c)));

                if (current_id.empty())
                    std::cerr << "[WARN] Missing ID before storing a sequence.\n";

                batch.emplace_back(current_id, std::move(seq_vec));
                seq.clear();

                // ⚠️ 达到上限：回退文件指针到当前 header
                if (batch.size() >= N) {
                    in.seekg(last_pos);
                    return true;
                }
            }

            current_id = (line.size() > 1) ? line.substr(1) : "";
        } else {
            seq += line;
        }
    }

    // 最后一条序列别漏
    if (!seq.empty()) {
        std::vector<unsigned char> seq_vec;
        seq_vec.reserve(seq.size());
        for (char c : seq)
            seq_vec.emplace_back(std::toupper(static_cast<unsigned char>(c)));
        batch.emplace_back(current_id, std::move(seq_vec));
    }

    return !batch.empty();
}


// add for viral msa --- end---


void insert_others(utils::Insertion2 That, utils::Insertion2 This, std::vector<std::vector<utils::Insertion2>>& more_insertions,
                   int k, int sequence_number, int* ii)
{
    for (int m = 0; m < sequence_number; m++)
    {
        if (m == k)
        {
            while ((ii[m] < more_insertions[m].size()) && (more_insertions[m][ii[m]].index < That.index)) ii[m]++;
            if (ii[m] == more_insertions[m].size())
                more_insertions[m].emplace_back(That);
            else if (more_insertions[m][ii[m]].index == That.index)
            {
                more_insertions[m][ii[m]].n_num += That.n_num;
                more_insertions[m][ii[m]].gap_num += That.gap_num;
            }
            else
                more_insertions[m].insert(more_insertions[m].begin() + (ii[m]++), That);
        }
        else
        {
            while ((ii[m] < more_insertions[m].size()) && (more_insertions[m][ii[m]].index < This.index)) ii[m]++;
            if (ii[m] == more_insertions[m].size())
                more_insertions[m].emplace_back(This);
            else if (more_insertions[m][ii[m]].index == This.index)
            {
                more_insertions[m][ii[m]].n_num += This.n_num;
                more_insertions[m][ii[m]].gap_num += This.gap_num;
            }
            else
                more_insertions[m].insert(more_insertions[m].begin() + (ii[m]++), This);
        }
    }
}


void utils::insert_and_write_file(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences,
                                  std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions,
                                  std::vector<std::string>& name, std::vector<bool>& sign)
{

    std::string each_sequence, each_line;
    const size_t sequence_number = insertions.size();
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    int i = 0, j = 0, pre = 0, diff, mi, k;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;

    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
    //变insertions
    for (int k = 0; k < sequence_number; k++)
    {
        i = 0; j = 0; g_num = 0;
        for (int kk = 0; kk < sequence_number; kk++) ii[kk] = 0;
        std::vector<Insertion2>().swap(i_all_insertions);
        while (i < insertions[k].size() && j < N_insertions[k].size())
        {
            if (insertions[k][i].index == N_insertions[k][j].index)
            {
                g_num += insertions[k][i].number;
                if (N_insertions[k][j].number > insertions[k][i].number)
                {
                    i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index ,insertions[k][i].number, 0 }));
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
                }
                else
                {
                    i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index ,N_insertions[k][j].number ,insertions[k][i].number - N_insertions[k][j].number }));
                }
                i++;
                j++;
            }
            else if (insertions[k][i].index > N_insertions[k][j].index)
            {
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
                j++;
            }
            else
            {
                g_num += insertions[k][i].number;
                i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
                i++;
            }
        }
        while (j < N_insertions[k].size())
        {
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size())
        {
            g_num += insertions[k][i].number;
            i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.emplace_back(i_all_insertions);
        /*for (int m = 0; m < sequence_number; m++)
        {
            for (int d = 0; d < insertions[m].size(); d++)
            {
                std::cout << insertions[m][d].index << " " << insertions[m][d].number << "\n";
            }
            std::cout << "\n";
        }
            std::cout << k<<" \n";*/
    }

    /*for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < insertions[k].size(); d++)
        {
            std::cout << insertions[k][d].index << " " << insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < N_insertions[k].size(); d++)
        {
            std::cout << N_insertions[k][d].index << " " << N_insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < all_insertions[k].size(); d++)
        {
            std::cout << all_insertions[k][d].index << " " << all_insertions[k][d].n_num<<" "<< all_insertions[k][d].gap_num << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < more_insertions[k].size(); d++)
        {
            std::cout << more_insertions[k][d].index << " " << more_insertions[k][d].n_num << " " << more_insertions[k][d].gap_num << "\n";
        }
        std::cout << k<<"\n";
    }*/

    int min_size, min_i;
    std::vector<Insertion> multi;
    std::vector<Insertion> tmp;
    for (int k = 0; k < sequence_number; k++)
        if (k == 0 || more_insertions[k].size() < min_size)
        {
            min_size = more_insertions[k].size();
            min_i = k;
        }
    i = 0; j = 0;
    if (min_i == 0) diff = 1;
    else diff = 0;
    while (i < more_insertions[diff].size() && j < more_insertions[min_i].size())
    {
        if (more_insertions[diff][i].index == more_insertions[min_i][j].index)
        {
            if (more_insertions[diff][i].gap_num == 0 || more_insertions[min_i][j].gap_num == 0);
            else if (more_insertions[diff][i].gap_num > more_insertions[min_i][j].gap_num)
                multi.emplace_back(Insertion({ more_insertions[diff][i].index ,more_insertions[min_i][j].gap_num }));
            else multi.emplace_back(Insertion({ more_insertions[diff][i].index ,more_insertions[diff][i].gap_num }));
            i++; j++;
        }
        else if (more_insertions[diff][i].index > more_insertions[min_i][j].index) j++;
        else i++;
    }

    if (min_i == 0)
        for (int k = 2; k < sequence_number; k++)
        {
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.emplace_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.emplace_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }
    else
        for (int k = 1; k < sequence_number; k++)
        {
            if (k == min_i) continue;
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.emplace_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.emplace_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }


    for (i = 0; i < sequence_number; i++)
    {
        if (sign[i]) os << "> " << name[i] << " + " << "\n";
        else os << "> " << name[i] << " - " << "\n";
        std::ofstream tmpo(arguments::tmp_file_name, std::ios::binary | std::ios::out);
        k = 0;
        for (j = 0; j < all_insertions[i].size(); j++)
        {
            while (k < all_insertions[i][j].index)
                tmpo << chars[sequences[i][k++]];
            for (int p = 0; p < all_insertions[i][j].n_num; p++)
                tmpo << 'N';
            for (int p = 0; p < all_insertions[i][j].gap_num; p++)
                tmpo << '-';
        }while (k < sequences[i].size())tmpo << chars[sequences[i][k++]];
        tmpo.close();
        //read
        std::ifstream tmpi(arguments::tmp_file_name, std::ios::binary | std::ios::in);
        while (std::getline(tmpi, each_line))
        {
            if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13))
                continue;
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
        tmpi.close();
        sequences[i].assign(each_sequence.begin(), each_sequence.end());
        each_sequence.clear();
        each_line.clear();
        //second insert
        mi = 0;
        k = 0;
        for (j = 0; j < more_insertions[i].size(); j++)
        {
            while (k < more_insertions[i][j].index)
                os << sequences[i][k++];
            for (int p = 0; p < more_insertions[i][j].n_num; p++)
                os << 'N';
            if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            {
                for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                    os << '-';
                mi++;
            }
            else
                for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                    os << '-';
        }while (k < sequences[i].size()) os << sequences[i][k++];
        std::vector<unsigned char>().swap(sequences[i]);
        os << "\n";
    }os << "\n";
    remove(arguments::tmp_file_name.data());

    /*std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;*/
}


int* utils::vector_insertion_gap_N(std::vector<std::vector<unsigned char>>& sequences,
                                   std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions) //写回比对结果
{
    const size_t sequence_number = insertions.size();
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    std::vector<unsigned char> tmp_vector;
    int i = 0, j = 0, pre = 0, diff, mi;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;
    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
    for (i = 0; i < sequence_number; i++)
        for (j = 0; j < N_insertions[i].size(); j++)
            len_sequences[i] += N_insertions[i][j].number;

    for (int k = 0; k < sequence_number; k++)
    {
        i = 0; j = 0; g_num = 0;
        for (int kk = 0; kk < sequence_number; kk++) ii[kk] = 0;
        std::vector<Insertion2>().swap(i_all_insertions);
        while (i < insertions[k].size() && j < N_insertions[k].size())
        {
            if (insertions[k][i].index == N_insertions[k][j].index)
            {
                g_num += insertions[k][i].number;
                if (N_insertions[k][j].number > insertions[k][i].number)
                {
                    i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index ,insertions[k][i].number, 0 }));
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
                }
                else
                {
                    i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index ,N_insertions[k][j].number ,insertions[k][i].number - N_insertions[k][j].number }));
                }
                i++;
                j++;
            }
            else if (insertions[k][i].index > N_insertions[k][j].index)
            {
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
                j++;
            }
            else
            {
                g_num += insertions[k][i].number;
                i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
                i++;
            }
        }
        while (j < N_insertions[k].size())
        {
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size())
        {
            g_num += insertions[k][i].number;
            i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.emplace_back(i_all_insertions);
        /*for (int m = 0; m < sequence_number; m++)
        {
            for (int d = 0; d < insertions[m].size(); d++)
            {
                std::cout << insertions[m][d].index << " " << insertions[m][d].number << "\n";
            }
            std::cout << "\n";
        }
            std::cout << k<<" \n";*/
    }

    /*for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < insertions[k].size(); d++)
        {
            std::cout << insertions[k][d].index << " " << insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < N_insertions[k].size(); d++)
        {
            std::cout << N_insertions[k][d].index << " " << N_insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < all_insertions[k].size(); d++)
        {
            std::cout << all_insertions[k][d].index << " " << all_insertions[k][d].n_num<<" "<< all_insertions[k][d].gap_num << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < more_insertions[k].size(); d++)
        {
            std::cout << more_insertions[k][d].index << " " << more_insertions[k][d].n_num << " " << more_insertions[k][d].gap_num << "\n";
        }
        std::cout << k<<"\n";
    }*/

    //消除  一列多个N导致de多个插空
    int min_size, min_i;
    std::vector<Insertion> multi;
    std::vector<Insertion> tmp;
    for (int k = 0; k < sequence_number; k++)
        if (k == 0 || more_insertions[k].size() < min_size)
        {
            min_size = more_insertions[k].size();
            min_i = k;
        }
    i = 0; j = 0;
    if (min_i == 0) diff = 1;
    else diff = 0;
    while (i < more_insertions[diff].size() && j < more_insertions[min_i].size())
    {
        if (more_insertions[diff][i].index == more_insertions[min_i][j].index)
        {
            if (more_insertions[diff][i].gap_num == 0 || more_insertions[min_i][j].gap_num == 0);
            else if (more_insertions[diff][i].gap_num > more_insertions[min_i][j].gap_num)
                multi.emplace_back(Insertion({ more_insertions[diff][i].index ,more_insertions[min_i][j].gap_num }));
            else multi.emplace_back(Insertion({ more_insertions[diff][i].index ,more_insertions[diff][i].gap_num }));
            i++; j++;
        }
        else if (more_insertions[diff][i].index > more_insertions[min_i][j].index) j++;
        else i++;
    }

    if (min_i == 0)
        for (int k = 2; k < sequence_number; k++)
        {
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.emplace_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.emplace_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }
    else
        for (int k = 1; k < sequence_number; k++)
        {
            if (k == min_i) continue;
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.emplace_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.emplace_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }
    //插入到原串
    int more = 0, all_size = 0, ti = 0, k = 0;
    for (i = 0; i < sequence_number; i++)
    {
        if (i == 0)
        {
            more = sequences[i].size();
            for (j = 0; j < all_insertions[i].size(); j++)
                more += (all_insertions[i][j].n_num + all_insertions[i][j].gap_num);
            all_size = more;
            mi = 0;
            for (j = 0; j < more_insertions[i].size(); j++)
            {
                all_size += more_insertions[i][j].n_num;
                if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
                    all_size -= multi[mi++].number;
                all_size += more_insertions[i][j].gap_num;
            }
            tmp_vector.resize(more);
        }

        ti = 0;
        k = 0;
        for (j = 0; j < all_insertions[i].size(); j++)
        {
            while (k < all_insertions[i][j].index)
                tmp_vector[ti++] = sequences[i][k++];
            for (int p = 0; p < all_insertions[i][j].n_num; p++)
                tmp_vector[ti++] = '\5';
            for (int p = 0; p < all_insertions[i][j].gap_num; p++)
                tmp_vector[ti++] = '\7';
        }while (k < sequences[i].size())tmp_vector[ti++] = sequences[i][k++];

        sequences[i].resize(all_size);
        mi = 0;
        k = 0;
        ti = 0;
        for (j = 0; j < more_insertions[i].size(); j++)
        {
            while (ti < more_insertions[i][j].index)
                sequences[i][k++] = tmp_vector[ti++];
            for (int p = 0; p < more_insertions[i][j].n_num; p++)
                sequences[i][k++] = '\5';
            if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            {
                for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                    sequences[i][k++] = '\7';
                mi++;
            }
            else
                for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                    sequences[i][k++] = '\7';
        }while (ti < more) sequences[i][k++] = tmp_vector[ti++];

    }


    std::vector<unsigned char>().swap(tmp_vector);
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;
    return len_sequences;
}



void utils::insert_and_write_fasta(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences,
                                   std::vector<std::vector<Insertion>>& insertions, std::vector<std::vector<Insertion>>& N_insertions,
                                   std::vector<std::string>& name,bool TU)

{
    if (!TU) chars[4] = 'U';
    //std::cout << "0memory usage: " << getPeakRSS() << " B" << std::endl;//输出内存耗费

    const size_t sequence_number = insertions.size();
    //std::cout << "insertions " << sequence_number << " " << insertions[0].size() << " " << insertions[1].size() << "\n";
    int score = 0, length = 0, name_len = 0;

    for (int i = 0; i < name.size(); i++)
        if (name_len < name[i].size())
            name_len = name[i].size();
    const auto align_start1 = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    int* score_two = NULL;
    std::vector<unsigned char> tmp_vector;
    int i = 0, j = 0, pre = 0, diff, mi;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;
    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
    for (i = 0; i < sequence_number; i++)
        for (j = 0; j < N_insertions[i].size(); j++)
            len_sequences[i] += N_insertions[i][j].number;

    for (int k = 0; k < sequence_number; k++)
    {
        i = 0; j = 0; g_num = 0;
        for (int kk = 0; kk < sequence_number; kk++) ii[kk] = 0;
        std::vector<Insertion2>().swap(i_all_insertions);
        while (i < insertions[k].size() && j < N_insertions[k].size())
        {
            if (insertions[k][i].index == N_insertions[k][j].index)
            {
                g_num += insertions[k][i].number;
                if (N_insertions[k][j].number > insertions[k][i].number)
                {
                    i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index ,insertions[k][i].number, 0 }));  //gap 个 N
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
                }
                else
                {
                    i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index ,N_insertions[k][j].number ,insertions[k][i].number - N_insertions[k][j].number }));
                }
                i++;
                j++;
            }
            else if (insertions[k][i].index > N_insertions[k][j].index) //不等插新
            {
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
                j++;
            }
            else
            {
                g_num += insertions[k][i].number;
                i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
                i++;
            }
        }
        while (j < N_insertions[k].size())
        {
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size())
        {
            g_num += insertions[k][i].number;
            i_all_insertions.emplace_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.emplace_back(i_all_insertions);
    }

    for (i = 0; i < N_insertions[i].size(); i++)
    {
        std::vector<Insertion>().swap(N_insertions[i]);
        std::vector<Insertion>().swap(insertions[i]);
    }
    std::vector<std::vector<Insertion>>().swap(N_insertions);
    std::vector<std::vector<Insertion>>().swap(insertions);

    int min_size, min_i;
    std::vector<Insertion> multi;
    std::vector<Insertion> tmp;
    for (int k = 0; k < sequence_number; k++)
        if (k == 0 || more_insertions[k].size() < min_size)
        {
            min_size = more_insertions[k].size();
            min_i = k;
        }
    i = 0; j = 0;
    if (min_i == 0) diff = 1;
    else diff = 0;
    while (i < more_insertions[diff].size() && j < more_insertions[min_i].size())
    {
        if (more_insertions[diff][i].index == more_insertions[min_i][j].index)
        {
            if (more_insertions[diff][i].gap_num == 0 || more_insertions[min_i][j].gap_num == 0);
            else if (more_insertions[diff][i].gap_num > more_insertions[min_i][j].gap_num)
                multi.emplace_back(Insertion({ more_insertions[diff][i].index ,more_insertions[min_i][j].gap_num }));
            else multi.emplace_back(Insertion({ more_insertions[diff][i].index ,more_insertions[diff][i].gap_num }));
            i++; j++;
        }
        else if (more_insertions[diff][i].index > more_insertions[min_i][j].index) j++;
        else i++;
    }

    if (min_i == 0)
        for (int k = 2; k < sequence_number; k++)
        {
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.emplace_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.emplace_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }
    else
        for (int k = 1; k < sequence_number; k++)
        {
            if (k == min_i) continue;
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.emplace_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.emplace_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }

    int more = 0, all_size = 0, ti = 0, k = 0;
    for (i = 0; i < sequence_number; i++)
    {
        if (i == 0)
        {
            more = sequences[i].size();
            for (j = 0; j < all_insertions[i].size(); j++)
                more += (all_insertions[i][j].n_num + all_insertions[i][j].gap_num);
            all_size = more;
            mi = 0;
            for (j = 0; j < more_insertions[i].size(); j++)
            {
                all_size += more_insertions[i][j].n_num;
                if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
                    all_size -= multi[mi++].number;
                all_size += more_insertions[i][j].gap_num;
            }
            tmp_vector.resize(more);
        }

        ti = 0;
        k = 0;
        for (j = 0; j < all_insertions[i].size(); j++)
        {
            while (k < all_insertions[i][j].index)
                tmp_vector[ti++] = sequences[i][k++];
            for (int p = 0; p < all_insertions[i][j].n_num; p++)
                tmp_vector[ti++] = '\5';
            for (int p = 0; p < all_insertions[i][j].gap_num; p++)
                tmp_vector[ti++] = '\7';
        }
        while (k < sequences[i].size())tmp_vector[ti++] = sequences[i][k++];

        sequences[i].resize(all_size);
        mi = 0;
        k = 0;
        ti = 0;
        for (j = 0; j < more_insertions[i].size(); j++)
        {
            while (ti < more_insertions[i][j].index)
                sequences[i][k++] = tmp_vector[ti++];
            for (int p = 0; p < more_insertions[i][j].n_num; p++)
                sequences[i][k++] = '\5';
            if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            {
                for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                    sequences[i][k++] = '\7';
                mi++;
            }
            else
                for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                    sequences[i][k++] = '\7';
        }while (ti < more) sequences[i][k++] = tmp_vector[ti++];
        os << "> " << name[i]<< "\n";
        for (k = 0; k < sequences[i].size(); k++) os << chars[sequences[i][k]];
        os << "\n";
        std::vector<unsigned char>().swap(sequences[i]);
    }
    std::vector<unsigned char>().swap(sequences[0]);
    std::vector<unsigned char>().swap(tmp_vector);
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;
}

void utils::insert_and_write(std::ostream& os, std::istream& is, const std::vector<std::vector<Insertion>>& insertions)
{
    const size_t sequence_number = insertions.size();

    std::string each_line;
    std::string each_sequence;
    std::string each_sequence_aligned;
    for (unsigned count = 0, length, flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            if (flag)
            {
                if (count == 0)
                {
                    length = each_sequence.size();
                    for (auto insertion : insertions[0])
                        length += insertion.number;

                    each_sequence_aligned.reserve(length);
                    each_sequence.reserve(length);
                }

                utils::Insertion::insert_gaps(each_sequence.cbegin(), each_sequence.cend(),
                                              insertions[count].cbegin(), insertions[count].cend(), std::back_inserter(each_sequence_aligned), '-');

                if (arguments::output_matrix)
                    os << each_sequence_aligned;
                else
                    Fasta::cut_and_write(os, each_sequence_aligned);
                os << '\n';

                each_sequence.clear();
                each_sequence_aligned.clear();
                ++count;
            }

            if (arguments::output_matrix == false)
                os << each_line << '\n';
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
#if defined(__unix__) || defined(__unix) || defined(unix)
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
#endif
        }
    }

    utils::Insertion::insert_gaps(each_sequence.cbegin(), each_sequence.cend(),
                                  insertions.back().cbegin(), insertions.back().cend(), std::back_inserter(each_sequence_aligned), '-');

    if (arguments::output_matrix)
        os << each_sequence_aligned;
    else
        Fasta::cut_and_write(os, each_sequence_aligned);
}

void utils::write_to_fasta(std::ostream& os, std::istream& is, std::vector<std::vector<Insertion>>& insertions, size_t& II)
{
    std::string each_line;
    std::string each_sequence;
    std::string pint_str;
    std::string name;
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;
        if (each_line[0] == '>')
        {
            if (flag)
            {
                utils::write_to_str(pint_str, each_sequence, insertions[II++]);
                os << name << "\n"<< pint_str<<"\n";
                //insertions.erase(insertions.begin());
                //for (int k = 0; k < pint_str.size(); k++) os << pint_str[k];
                each_sequence.clear();
            }
            name = each_line;
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }
    utils::write_to_str(pint_str, each_sequence, insertions[II++]);
    os << name << "\n" << pint_str << "\n";
    //insertions.erase(insertions.begin());
    return;
}

void utils::write_to_str(std::string& ans, std::string& each_sequence, std::vector<Insertion>& insertions)
{
    size_t ti = 0;
    size_t k = 0;
    if (arguments::ALL_LEN != 0)
    {
        if (ans.size() != arguments::ALL_LEN)
            ans.resize(arguments::ALL_LEN);
    }
    else
    {
        arguments::ALL_LEN = each_sequence.size();
        for (size_t j = 0; j < insertions.size(); j++)
            arguments::ALL_LEN += insertions[j].number;
        ans.resize(arguments::ALL_LEN);
    }
    /*for (size_t j = 0; j < insertions.size(); j++)
    {
        while (k < insertions[j].index)
            ans[ti++] = each_sequence[k++];
        for (int p = 0; p < insertions[j].number; p++)
            ans[ti++] = '-';
    }
    while (k < each_sequence.size())ans[ti++] = each_sequence[k++];*/
    for (const auto& insertion : insertions)
    {
        std::copy(each_sequence.begin() + k, each_sequence.begin() + insertion.index, ans.begin() + ti);
        ti += insertion.index - k;

        std::fill(ans.begin() + ti, ans.begin() + ti + insertion.number, '-');
        ti += insertion.number;

        k = insertion.index;
    }

    std::copy(each_sequence.begin() + k, each_sequence.end(), ans.begin() + ti);
    return;
}

unsigned char* utils::copy_DNA(const std::vector<unsigned char>& sequence, unsigned char* A, size_t a_begin, size_t a_end)
{
    int i = 0;
    while (a_begin < a_end)
    {
        A[i++] = sequence[a_begin++];
    }
    return A;
}
