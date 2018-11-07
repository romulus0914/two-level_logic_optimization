#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef struct string {
    int len;
    char *s;
} String;

char *filename;
int num_vars, num_products;
String *sop;
String **groups;

String StrNew(char *str)
{
    String string;
    string.len = num_vars + 1;
    string.s = (char *)malloc(string.len * sizeof(char));
    strcpy(string.s, str);
    return string;
}

void ReadFile()
{
    FILE *file = fopen(filename, "r");

    fscanf(file, "%d", &num_vars);
    fscanf(file, "%d", &num_products);

    ENTRY item;
    hcreate(num_products);

    sop = (String *)malloc(num_products * sizeof(String));

    int actual_num_products = 0;
    char temp[num_vars + 1];
    for (int i = 0; i < num_products; i++) {
        fscanf(file, "%s", temp);
        for (int j = 0; j < num_vars; j++)
            if (temp[j] == '-')
                temp[j] = '2';
        item.key = StrNew(temp).s;
        if (hsearch(item, FIND) == NULL) {
            item.data = temp;
            hsearch(item, ENTER);
            sop[actual_num_products++] = StrNew(temp);
        }
    }

    hdestroy();
    fclose(file);

    num_products = actual_num_products;
}

int PermutationC(int m, int n)
{
    long c = 1;
    for (int i = m; i > m - n; i--)
        c *= i;
    for (int i = n; i > 1; i--)
        c /= i;
    return c;
}

void ExpandToMinterms(int *groups_num_cnt)
{
    int *groups_num = (int *)calloc((num_vars + 1), sizeof(int));

    int one_cnt;
    int *dc_cnt = (int *)calloc(num_products, sizeof(int));
    for (int i = 0; i < num_products; i++) {
        one_cnt =  0;
        for (int j = 0; j < num_vars; j++) {
            if (sop[i].s[j] == '1')
                one_cnt++;
            else if (sop[i].s[j] == '2')
                dc_cnt[i]++;
        }
        for (int j = one_cnt; j <= one_cnt + dc_cnt[i]; j++)
            groups_num[j] += PermutationC(dc_cnt[i], j - one_cnt);
    }

    int max_groups_num = 0, total_groups_num = 0;
    for (int i = 0; i <= num_vars; i++) {
        total_groups_num += groups_num[i];
        if (groups_num[i] > max_groups_num)
            max_groups_num = groups_num[i];
    }

    groups = (String **)malloc((num_vars + 1) * sizeof(String *));
    for (int i = 0; i <= num_vars; i++)
        groups[i] = (String *)malloc(max_groups_num * sizeof(String));

    ENTRY item;
    hcreate(total_groups_num);

    for (int i = 0; i < num_products; i++) {
        int dc_pow = pow(2, dc_cnt[i]);
        for (int j = 0; j < dc_pow; j++) {
            String str = StrNew(sop[i].s);
            int bit = j;
            int num_ones = 0;
            for (int k = 0; k < num_vars; k++) {
                if (str.s[k] == '2') {
                    str.s[k] = bit % 2 ? '1' : '0';
                    bit /= 2;
                }
                if (str.s[k] == '1')
                    num_ones++;
            }
            item.key = str.s;
            if (hsearch(item, FIND) == NULL) {
                item.data = str.s;
                hsearch(item, ENTER);
                groups[num_ones][groups_num_cnt[num_ones]++] = str;
            }
        }
    }

    hdestroy();
    free(groups_num);
    free(dc_cnt);
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

int ColumnCovering(String *sop_opt)
{
    int *dc_cnt = (int *)calloc(num_products, sizeof(int));
    for (int i = 0; i < num_products; i++) {
        for (int j = 0; j < num_vars; j++) {
            if (sop[i].s[j] == '2') {
                dc_cnt[i]++;
                sop[i].s[j] = '-';
            }
        }
    }

    int total_numbers = 0;
    int *dc_pow = (int *)malloc(num_products * sizeof(int));
    int max_dc_pow = 0;
    for (int i = 0; i < num_products; i++) {
        dc_pow[i] = pow(2, dc_cnt[i]);
        if (dc_pow[i] > max_dc_pow)
            max_dc_pow = dc_pow[i];
        total_numbers += dc_pow[i];
    }

    int *numbers = (int *)malloc(total_numbers * sizeof(int));
    int **sop_numbers = (int **)malloc(num_products * sizeof(int *));
    for (int i = 0; i < num_products; i++)
        sop_numbers[i] = (int *)malloc(max_dc_pow * sizeof(int));
    int *sop_numbers_num = (int *)calloc(num_products, sizeof(int));

    ENTRY item;
    hcreate(total_numbers * 0.5);

    total_numbers = 0;
    for (int i = 0; i < num_products; i++) {
        for (int j = 0; j < dc_pow[i]; j++) {
            String str = StrNew(sop[i].s);
            int number = 0;
            int bit = j;
            for (int k = num_vars - 1; k >= 0; k--) {
                if (str.s[k] == '-') {
                    str.s[k] = bit % 2 ? '1' : '0';
                    bit /= 2;
                }
                if (str.s[k] == '1')
                    number += pow(2, num_vars - 1 - k);
            }
            sop_numbers[i][sop_numbers_num[i]++] = number;

            item.key = str.s;
            if (hsearch(item, FIND) == NULL) {
                item.data = str.s;
                hsearch(item, ENTER);
                numbers[total_numbers++] = number;
            }
        }
    }

    free(dc_cnt);
    free(dc_pow);
    hdestroy();

    int **matrix = (int **)malloc(num_products * sizeof(int *));
    for (int i = 0; i < num_products; i++)
        matrix[i] = (int *)calloc(total_numbers, sizeof(int));

    for (int i = 0; i < num_products; i++) {
        for (int j = 0; j < sop_numbers_num[i]; j++) {
            for (int k = 0; k < total_numbers; k++) {
                if (sop_numbers[i][j] == numbers[k]) {
                    matrix[i][k] = 1;
                    break;
                }
            }
        }
    }

    int num_products_opt = 0;

    int *numbers_cnt = (int *)calloc(total_numbers, sizeof(int));
    for (int i = 0; i < total_numbers; i++)
        for (int j = 0; j < num_products; j++)
            numbers_cnt[i] += matrix[j][i];

    for (int i = 0; i < total_numbers; i++) {
        if (numbers_cnt[i] == 1) {
            for (int j = 0; j < num_products; j++) {
                // which implicant contributes this 1
                if (matrix[j][i] == 1) {
                    sop_opt[num_products_opt++] = StrNew(sop[j].s);
                    matrix[j][i] = 0;
                    sop_numbers_num[j] = 0;
                    // set columns covered by this implicant to 0
                    for (int k = 0; k < total_numbers; k++) {
                        if (matrix[j][k] == 1) {
                            numbers_cnt[k] = 0;
                            matrix[j][k] = 0;
                            for (int l = 0; l < num_products; l++) {
                                if (matrix[l][k] == 1) {
                                    matrix[l][k] = 0;
                                    sop_numbers_num[l]--;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    free(numbers);
    free(numbers_cnt);

    int idx = MaxArray(sop_numbers_num, num_products);
    while (idx != -1) {
        sop_opt[num_products_opt++] = StrNew(sop[idx].s);
        sop_numbers_num[idx] = 0;
        for (int i = 0; i < total_numbers; i++) {
            if (matrix[idx][i] == 1) {
                matrix[idx][i] = 0;
                for (int j = 0; j < num_products; j++) {
                    if (matrix[j][i] == 1) {
                        matrix[j][i] = 0;
                        sop_numbers_num[j]--;
                    }
                }
            }
        }
        idx = MaxArray(sop_numbers_num, num_products);
    }

    for (int i = 0; i < num_products; i++)
        free(matrix[i]);
    free(matrix);
    for (int i = 0; i < num_products; i++)
        free(sop_numbers[i]);
    free(sop_numbers);
    free(sop_numbers_num);

    return num_products_opt;
}

void Output(String *sop_opt, int num_products_opt)
{
    int literals = 0;
    for (int i = 0; i < num_products_opt; i++)
        for (int j = 0; j < num_vars; j++)
            if (sop_opt[i].s[j] != '-')
                literals++;

    strcat(filename, ".out");
    FILE *out = fopen(filename, "w");
    fprintf(out, "%d\n%d\n", literals, num_products_opt);
    for (int i = 0; i < num_products_opt; i++)
        fprintf(out, "%s\n", sop_opt[i].s);
}

int Merge(char *str1, char *str2)
{
    int idx = -1;
    int only_once = 0;
    for (int i = 0; i < num_vars; i++) {
        if (str1[i] != str2[i]) {
            if (str1[i] == '2' || str2[i] == '2') {
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

void PrimeImplicantGeneration(int *groups_num)
{
    const int magical_num = 25;
    num_products = 0;

    for (int cycle = num_vars; cycle > 0; cycle--) {
        int max_group_num = 0;
        for (int i = 0; i <= cycle; i++)
            if (groups_num[i] > max_group_num)
                max_group_num = groups_num[i];

        String **groups_temp;
        groups_temp = (String **)malloc(cycle * sizeof(String *));
        for (int i = 0; i < cycle; i++)
            groups_temp[i] = (String *)malloc(max_group_num * magical_num * sizeof(String));

        int *groups_num_temp = (int *)calloc(cycle, sizeof(int));

        int **products_used = (int **)malloc((cycle + 1) * sizeof(int *));
        for (int i = 0; i < cycle + 1; i++)
            products_used[i] = (int *)calloc(max_group_num, sizeof(int));

        #pragma omp parallel num_threads(8)
        {
            #pragma omp for schedule(dynamic, 1)
            for (int group = 0; group < cycle; group++) {
                for (int i = 0; i < groups_num[group]; i++) {
                    for (int j = 0; j < groups_num[group + 1]; j++) {
                        int dc_pos = Merge(groups[group][i].s, groups[group + 1][j].s);
                        if (dc_pos != -1) {
                            products_used[group][i] = 1;
                            products_used[group + 1][j] = 1;

                            String str = StrNew(groups[group][i].s);
                            str.s[dc_pos] = '2';
                            groups_temp[group][groups_num_temp[group]++] = str;
                        }
                    }
                }
            }
        }

        // collect unused products as prime implicants
        for (int i = 0; i < cycle + 1; i++)
            for (int j = 0; j < groups_num[i]; j++)
                if (products_used[i][j] == 0)
                    sop[num_products++] = groups[i][j];

        // free memory and memcpy for next for loop
        for (int i = 0; i <= cycle; i++)
            free(groups[i]);
        free(groups);
        free(groups_num);

        groups = (String **)malloc(cycle * sizeof(String *));
        for (int i = 0; i < cycle; i++)
            groups[i] = (String *)malloc(max_group_num * magical_num * sizeof(String));

        groups_num = (int *)calloc(cycle, sizeof(int));
        for (int i = 0; i < cycle; i++)
            groups_num[i] = 0;

        ENTRY item;
        for (int i = 0; i < cycle; i++) {
            hcreate(groups_num_temp[i] * 1.25);
            for (int j = 0; j < groups_num_temp[i]; j++) {
                String str = StrNew(groups_temp[i][j].s);
                item.key = str.s;
                if (hsearch(item, FIND) == NULL) {
                    item.data = str.s;
                    hsearch(item, ENTER);
                    groups[i][groups_num[i]++] = str;
                }
            }
            hdestroy();
        }

        for (int i = 0; i < cycle; i++) {
            for (int j = 0; j < groups_num_temp[i]; j++)
                free(groups_temp[i][j].s);
            free(groups_temp[i]);
        }
        free(groups_temp);
        free(groups_num_temp);
        for (int i = 0; i < cycle + 1; i++)
            free(products_used[i]);
        free(products_used);

        // early stop
        int done = 1;
        for (int i = 0; i < cycle; i++)
            if (groups_num[i] > 0) {
                done = 0;
                break;
            }
        if (done) {
            for (int i = 0; i < cycle; i++) {
                for (int j = 0; j < groups_num[i]; j++)
                    free(groups[i][j].s);
                free(groups[i]);
            }
            free(groups);
            free(groups_num);
            break;
        }
    }
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage:\n");
        printf("    ./two-level_logic_opt <path/to/testcase>\n");
        exit(1);
    }
    filename = argv[1];

    ReadFile();

    int *groups_num = (int *)calloc((num_vars + 1), sizeof(int));
    ExpandToMinterms(groups_num);

    PrimeImplicantGeneration(groups_num);

    String *sop_opt = (String *)malloc(num_products * sizeof(String));
    int num_products_opt = ColumnCovering(sop_opt);

    Output(sop_opt, num_products_opt);

    free(sop_opt);

    return 0;
}
