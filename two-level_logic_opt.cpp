#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>

#include <cmath>
#include <omp.h>

using namespace std;

int num_vars, num_products;

void ReadFile(string filename, vector<string> &sop)
{
    ifstream file;
    file.open(filename);

    file >> num_vars >> num_products;

    set<string> temp;
    string str;
    for (int i = 0; i < num_products; i++) {
        file >> str;
        temp.insert(str);
    }

    sop.assign(temp.begin(), temp.end());
    num_products = sop.size();

    file.close();
}

void ExpandToMinterms(vector<string> &sop, vector<vector<string>> &groups)
{
    vector<set<string>> temp(num_vars + 1);

   for (int i = 0; i < num_products; i++) {
        int dc_cnt = 0;
        for (int j = 0; j < num_vars; j++)
            if (sop[i][j] == '-')
                dc_cnt++;
        int dc_pow = pow(2, dc_cnt);
        for (int j = 0; j < dc_pow; j++) {
            string str = sop[i];
            int bit = j;
            int num_ones = 0;
            for (int k = 0; k < num_vars; k++) {
                if (str[k] == '-') {
                    str[k] = bit % 2 ? '1' : '0';
                    bit /= 2;
                }
                if (str[k] == '1')
                    num_ones++;
            }
            temp[num_ones].insert(str);
        }
    }

    for (int i = 0; i <= num_vars; i++)
        groups[i].assign(temp[i].begin(), temp[i].end());
}

int Merge(string str1, string str2)
{
    int idx = -1;
    int only_once = 0;
    for (int i = 0; i < num_vars; i++) {
        if (str1[i] != str2[i]) {
            if (str1[i] == '-' || str2[i] == '-') {
                idx = -1;
                break;
            }
            else {
                if (only_once == 0) {
                    only_once = 1;
                    idx = i;
                }
                else {
                    idx = -1;
                    break;
                }
            }
        }
    }
    return idx;
}

void PrimeImplicantGeneration(vector<string> &sop, vector<vector<string>> &groups)
{
    vector<string>().swap(sop);

    for (int cycle = num_vars; cycle > 0; cycle--) {
        vector<vector<string>> groups_temp(cycle);

        vector<vector<int>> products_used;
        for (int i = 0; i <= cycle; i++) {
            vector<int> temp(groups[i].size(), 0);
            products_used.emplace_back(temp);
        }

        set<string> temp[cycle];
        #pragma omp parallel num_threads(8)
        {
            #pragma omp for schedule(dynamic, 1)
            for (int group = 0; group < cycle; group++) {
                int len1 = groups[group].size();
                int len2 = groups[group + 1].size();
                for (int i = 0; i < len1; i++) {
                    for (int j = 0; j < len2; j++) {
                        int dc_pos = Merge(groups[group][i], groups[group + 1][j]);
                        if (dc_pos != -1) {
                            products_used[group][i] = 1;
                            products_used[group + 1][j] = 1;
                            string str = groups[group][i];
                            str[dc_pos] = '-';
                            temp[group].insert(str);
                        }
                    }
                }
            }
        }

        for (int i = 0; i <= cycle; i++) {
            int len = groups[i].size();
            for (int j = 0; j < len; j++) {
                if (products_used[i][j] == 0)
                    sop.emplace_back(groups[i][j]);
            }
        }

        for (int i = 0; i < cycle; i++)
            groups[i] = vector<string>(temp[i].begin(), temp[i].end());

        int done = 1;
        for (int i = 0; i < cycle; i++) {
            if (groups[i].size() > 0) {
                done = 0;
                break;
            }
        }
        if (done)
            break;
    }

    num_products = sop.size();
}

int MaxArray(int *arr, int len)
{
    int max = 0, idx = -1;
    for (int i = 0; i < len; i++)
        if (arr[i] > max) {
            max = arr[i];
            idx = i;
        }
    return idx;
}

void ColumnCovering(vector<string> &sop, vector<string> &sop_opt)
{
    set<int> temp;
    vector<vector<int>> sop_numbers(num_products);

   for (int i = 0; i < num_products; i++) {
        int dc_cnt = 0;
        for (int j = 0; j < num_vars; j++)
            if (sop[i][j] == '-')
                dc_cnt++;
        int dc_pow = pow(2, dc_cnt);
        for (int j = 0; j < dc_pow; j++) {
            string str = sop[i];
            int bit = j;
            int number = 0;
            for (int k = num_vars - 1; k >= 0; k--) {
                if (str[k] == '-') {
                    str[k] = bit % 2 ? '1' : '0';
                    bit /= 2;
                }
                if (str[k] == '1')
                    number += pow(2, num_vars - 1 - k);
            }
            sop_numbers[i].emplace_back(number);
            temp.insert(number);
        }
    }

    vector<int> numbers(temp.begin(), temp.end());
    int total_numbers = numbers.size();
    int **matrix = (int **)malloc(num_products * sizeof(int *));
    for (int i = 0; i < num_products; i++)
        matrix[i] = (int *)calloc(total_numbers, sizeof(int));

    for (int i = 0; i < num_products; i++) {
        int size = sop_numbers[i].size();
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < total_numbers; k++) {
                if (sop_numbers[i][j] == numbers[k]) {
                    matrix[i][k] = 1;
                    break;
                }
            }
        }
    }

    int *numbers_cnt = (int *)calloc(total_numbers, sizeof(int));
    for (int i = 0; i < total_numbers; i++)
        for (int j = 0; j < num_products; j++)
            numbers_cnt[i] += matrix[j][i];

    int *sop_numbers_cnt = (int *)malloc(num_products * sizeof(int));
    for (int i = 0; i < num_products; i++)
        sop_numbers_cnt[i] = sop_numbers[i].size();

    for (int i = 0; i < total_numbers; i++) {
        if (numbers_cnt[i] == 1) {
            for (int j = 0; j < num_products; j++) {
                if (matrix[j][i] == 1) {
                    sop_opt.emplace_back(sop[j]);
                    matrix[j][i] = 0;
                    sop_numbers_cnt[j] = 0;
                    for (int k = 0; k < total_numbers; k++) {
                        if (matrix[j][k] == 1) {
                            numbers_cnt[k] = 0;
                            matrix[j][k] = 0;
                            for (int l = 0; l < num_products; l++) {
                                if (matrix[l][k] == 1) {
                                    matrix[l][k] = 0;
                                    sop_numbers_cnt[l]--;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    free(numbers_cnt);

    int idx = MaxArray(sop_numbers_cnt, num_products);
    while (idx != -1) {
        sop_opt.emplace_back(sop[idx]);
        sop_numbers_cnt[idx] = 0;
        for (int i = 0; i < total_numbers; i++) {
            if (matrix[idx][i] == 1) {
                matrix[idx][i] = 0;
                for (int j = 0; j < num_products; j++) {
                    if (matrix[j][i] == 1) {
                        matrix[j][i] = 0;
                        sop_numbers_cnt[j]--;
                    }
                }
            }
        }
        idx = MaxArray(sop_numbers_cnt, num_products);
    }

    for (int i = 0; i < num_products; i++)
        free(matrix[i]);
    free(matrix);
    free(sop_numbers_cnt);
}

void Output(vector<string> &sop_opt, string filename)
{
    ofstream file;
    file.open(filename + ".out");

    num_products = sop_opt.size();
    int literals = 0;
    for (int i = 0; i < num_products; i++)
        for (int j = 0; j < num_vars; j++)
            if (sop_opt[i][j] != '-')
                literals++;

    file << literals << '\n' << num_products << '\n';
    for (int i = 0; i < num_products; i++)
        file << sop_opt[i] << '\n';

    file.close();
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        cout << "Usage: \n";
        cout << "./two-level_logic_opt <path/to/testcase>\n";
        exit(1);
    }

    vector<string> sop;
    ReadFile(argv[1], sop);

    vector<vector<string>> groups(num_vars + 1);
    ExpandToMinterms(sop, groups);

    PrimeImplicantGeneration(sop, groups);

    vector<string> sop_opt;
    ColumnCovering(sop, sop_opt);

    Output(sop_opt, argv[1]);

    return 0;
}
