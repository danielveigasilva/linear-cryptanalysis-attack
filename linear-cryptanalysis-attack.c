#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define uchar unsigned char 

#define KEY_SIZE 4
#define FILE_NAME_BIAS_CSV "biastable.csv"

/*
#define N_BITS 4
#define SBOX_SIZE 16
uchar SBOX[SBOX_SIZE] = { 14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7 };
*/

#define N_BITS 8
#define SBOX_SIZE 256
uchar SBOX[SBOX_SIZE] = {84, 132, 218, 91, 208, 159, 31, 168, 41, 176, 1, 19, 243, 100, 196, 34, 149, 77, 81, 167, 137, 94, 22, 16, 169, 190, 6, 56, 83, 0, 104, 90, 175, 139, 222, 7, 207, 20, 98, 64, 86, 24, 110, 8, 50, 130, 46, 170, 143, 131, 166, 74, 14, 101, 203, 160, 85, 60, 99, 191, 150, 182, 205, 231, 252, 198, 212, 179, 206, 189, 120, 58, 23, 242, 108, 129, 12, 118, 188, 211, 135, 180, 26, 13, 106, 235, 111, 192, 238, 214, 165, 154, 122, 236, 72, 247, 197, 39, 240, 224, 146, 125, 68, 158, 15, 45, 220, 173, 107, 226, 254, 217, 66, 121, 174, 225, 215, 42, 55, 43, 186, 147, 201, 80, 141, 28, 234, 95, 221, 144, 113, 164, 71, 187, 5, 52, 70, 127, 183, 33, 219, 202, 2, 177, 178, 228, 44, 117, 59, 255, 136, 102, 65, 57, 155, 87, 114, 210, 193, 62, 124, 233, 246, 92, 133, 97, 163, 195, 105, 89, 10, 232, 151, 152, 47, 229, 76, 199, 213, 30, 115, 245, 126, 239, 119, 153, 248, 145, 142, 156, 38, 96, 237, 128, 140, 223, 9, 54, 194, 75, 79, 249, 88, 67, 63, 244, 148, 32, 78, 171, 200, 21, 17, 230, 250, 162, 253, 53, 61, 181, 69, 3, 112, 109, 35, 209, 4, 82, 204, 241, 157, 29, 161, 18, 73, 184, 93, 40, 37, 27, 103, 134, 11, 36, 48, 123, 116, 216, 25, 172, 49, 227, 51, 185, 251, 138};


typedef struct expression_type {
    int x;
    int y;
    int bias;
    int hamming_x;
    float prob;
    int *xn;
    int *yn;
} expression;

typedef struct expressions_type {
    expression * expressions;
    int n_expression;
} expressions;


uchar e_block( uchar w, uchar k ) {
    uchar x = w ^ k;
    uchar y = SBOX[x];
    return y;
}

uchar * e_text( uchar * W, uchar * K, int w_size, int key_size) {
    uchar * C = (uchar *) calloc(w_size, sizeof(uchar));
    for (int i = 0; i < w_size; i ++)
        C[i] = e_block(W[i], K[i % key_size]);
    return C;
}



uchar dot_product(uchar w, uchar m, int n_bits){
    uchar w_and_m = w & m;
    uchar w_dot_m = 0;
    for (int i = 0; i < n_bits; i++)
        w_dot_m = w_dot_m ^ ((w_and_m >> i) & 1);
    return w_dot_m;
}

int mod(int n){
    return n < 0 ? n * -1 : n;
}

uchar bit( uchar x, uchar i ) {
    return (x >> (N_BITS-i)) & 1;
}


int ** get_bias_table(uchar * SBOX, int sbox_size, int n_bits) {
    
    int ** bias_table = (int ** ) calloc(sbox_size, sizeof(int*));

    for( int a = 0; a < sbox_size; a++ ) {
        bias_table[a] = (int * ) calloc(sbox_size, sizeof(int));

        for( int b = 0; b < sbox_size; b++ ) {
            int e = 0;
            for( int x = 0; x < sbox_size; x++ )
                if ((dot_product(a, x, n_bits) ^ dot_product(b, SBOX[x], n_bits)) == 0)
                    e += 1;
            bias_table[a][b] = e - pow(2, n_bits - 1);
        }
    }

    return bias_table;
}

expressions get_expressions(int ** bias_table, int sbox_size, int n_bits){
    expressions all_expres;

    all_expres.n_expression = (sbox_size - 1) * (sbox_size - 1);
    all_expres.expressions = (expression *) calloc(all_expres.n_expression, sizeof(expression));
    int index_expression = 0;
    for( int a = 1; a < sbox_size; a++ ) {
        for( int b = 1; b < sbox_size; b++ ) {
                all_expres.expressions[index_expression].x = a;
                all_expres.expressions[index_expression].y = b;
                all_expres.expressions[index_expression].bias = bias_table[a][b];
                if (bias_table[a][b] < 0)
                    all_expres.expressions[index_expression].prob = 1 - (0.5 + (float)bias_table[a][b]/(float)sbox_size);
                else
                    all_expres.expressions[index_expression].prob = 0.5 + (float)bias_table[a][b]/(float)sbox_size;

                all_expres.expressions[index_expression].xn = (int *) calloc((n_bits + 1), sizeof(int));
                all_expres.expressions[index_expression].yn = (int *) calloc((n_bits + 1), sizeof(int));
                for (int bit = 1; bit < n_bits + 1; bit ++){
                    all_expres.expressions[index_expression].xn[bit] = ((a >> (n_bits - bit)) & 1);
                    all_expres.expressions[index_expression].yn[bit] = ((b >> (n_bits - bit)) & 1);
                    all_expres.expressions[index_expression].hamming_x += all_expres.expressions[index_expression].xn[bit];
                }
                index_expression ++;
        }
    }


    expressions best_expres_for_xi;
    best_expres_for_xi.n_expression = N_BITS + 1;
    best_expres_for_xi.expressions = (expression *) calloc(best_expres_for_xi.n_expression, sizeof(expression));

    for (int i = 0; i < n_bits + 1; i ++){
        best_expres_for_xi.expressions[i].hamming_x = n_bits;
        best_expres_for_xi.expressions[i].prob = 0;
        best_expres_for_xi.expressions[i].x = -1;
        best_expres_for_xi.expressions[i].bias = 0;
    }

    for (int index_expres = 0; index_expres < all_expres.n_expression; index_expres ++){
            for (int index_x = 1; index_x < n_bits + 1; index_x++){
                if (all_expres.expressions[index_expres].xn[index_x] == 1 &&
                    all_expres.expressions[index_expres].hamming_x == 1 &&
                    mod(all_expres.expressions[index_expres].bias) >= mod(best_expres_for_xi.expressions[index_x].bias)){
                        
                        best_expres_for_xi.expressions[index_x].bias = all_expres.expressions[index_expres].bias;
                        best_expres_for_xi.expressions[index_x].x = all_expres.expressions[index_expres].x;
                        best_expres_for_xi.expressions[index_x].y = all_expres.expressions[index_expres].y;
                        best_expres_for_xi.expressions[index_x].prob = all_expres.expressions[index_expres].prob;
                        best_expres_for_xi.expressions[index_x].hamming_x = all_expres.expressions[index_expres].hamming_x;
                        best_expres_for_xi.expressions[index_x].xn = all_expres.expressions[index_expres].xn;
                        best_expres_for_xi.expressions[index_x].yn = all_expres.expressions[index_expres].yn;
                } 
            }
    }
    return best_expres_for_xi;
}


void print_bias_table_terminal(int ** bias_table, int sbox_size, int n_bits){
    printf("\n\n * TABELA DE BIAS\n\n");
    for( int a = 1; a < sbox_size; a++ ) {
        printf("    ");
        for( int b = 1; b < sbox_size; b++ ){
            printf("% d ", bias_table[a][b]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void print_expressions_terminal(expressions expres, int n_bits){

    printf("\n\n * EQUAÇÕES ESCOLHIDAS\n\n");
    for (int i = 0 ; i < expres.n_expression ; i ++){

        if (expres.expressions[i].x < 0)
            continue;

        printf("   * Para X = %d (hamming: %d) e Y = %d temos um bias de %d, resultando em:\n\n", expres.expressions[i].x, expres.expressions[i].hamming_x, expres.expressions[i].y, expres.expressions[i].bias);

        printf("          Pr[");
        int first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    printf(" XOR ");
                else
                    first_flag = 0;
                printf("x%d", j);
            }
        }
        for (int j = 1; j < n_bits + 1; j ++)
            if (expres.expressions[i].yn[j] == 1)
                printf(" XOR y%d", j);
        if (expres.expressions[i].bias < 0)
            printf(" = 1");
        else
            printf(" = 0");
        
        printf("] = %.3f", expres.expressions[i].prob);
        printf("\n\n          ");

        first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    printf(" XOR ");
                else
                    first_flag = 0;
                printf("k%d", j);
            }
        }

        printf(" = ");

        first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    printf(" XOR ");
                else
                    first_flag = 0;
                printf("w%d", j);
            }
        }

        for (int j = 1; j < n_bits + 1; j ++)
            if (expres.expressions[i].yn[j] == 1)
                printf(" XOR y%d", j);
        if (expres.expressions[i].bias < 0)
            printf(" XOR 1");

        printf("\n\n          ");


        printf("C CODE: ");

        first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    printf(" ^ ");
                else
                    first_flag = 0;
                printf("bit(K,%d)", j);
            }
        }

        printf(" == ");

        first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    printf(" ^ ");
                else
                    first_flag = 0;
                printf("bit(W,%d)", j);
            }
        }

        for (int j = 1; j < n_bits + 1; j ++)
            if (expres.expressions[i].yn[j] == 1)
                printf(" ^ bit(Y,%d)", j);
        if (expres.expressions[i].bias < 0)
            printf(" ^ 1");

        printf("\n\n\n");
    }
}


void print_bias_table_csv(int ** bias_table, int sbox_size, int n_bits, char * fileName){

    FILE * file_csv = fopen(fileName, "w");

    for (int b = 1; b < sbox_size; b++)
        fprintf(file_csv, " ;%x", b);
    fprintf(file_csv, "\n");
    for( int a = 1; a < sbox_size; a++ ) {
        fprintf(file_csv, "%x;", a);
        for( int b = 1; b < sbox_size; b++ ){
            fprintf(file_csv, "%d", bias_table[a][b]);
            if (b != sbox_size - 1)
                fprintf(file_csv, ";");
        }
        fprintf(file_csv,"\n");
    }
}


void print_bias_table_latex(int ** bias_table, int sbox_size, int n_bits, FILE * file){
    
    int sbox_size_limit_to_show = 16;
    
    fprintf(file, "\n\n"); 
    fprintf(file, "\\section{Tabela de Bias}\n");
    fprintf(file, "A partir da análise da S-BOX disponibilizada geramos a Tabela de Bias abaixo.\n");
    fprintf(file,"\\begin{center}\n"); 
    fprintf(file,"Tabela de Bias para SBOX (fator de 1/%d)*\n", sbox_size);
    fprintf(file,"    \\begin{tabular} ");
    fprintf(file,"{"); 
    for (int i = 0; i < sbox_size && i < sbox_size_limit_to_show ; i++)
        fprintf(file,"|c");
    
    if (sbox_size > sbox_size_limit_to_show)
        fprintf(file,"}\n");    
    else
        fprintf(file,"|}\n");

    fprintf(file,"        \\hline\n       ");
    
    for (int b = 1; b < sbox_size && b < sbox_size_limit_to_show ; b++){
        if (sbox_size > sbox_size_limit_to_show && b == sbox_size_limit_to_show - 1)
            fprintf(file," & \\textbf{...}");    
        else
            fprintf(file," & \\textbf{%x}", b);
    }

    fprintf(file,"\\\\\n");
    fprintf(file,"        \\hline\n       ");

    for( int a = 1; a < sbox_size && a < sbox_size_limit_to_show; a++ ) {
        
        if (sbox_size > sbox_size_limit_to_show && a > sbox_size_limit_to_show - 3){
            fprintf(file,"\\textbf{.} & ");
            for( int b = 1; b < sbox_size && b < sbox_size_limit_to_show - 1; b++ ){
                fprintf(file, " . & ");
            }
            fprintf(file,"\\\\\n");
        }    
        else{
            fprintf(file,"\\textbf{%x} & ", a);
            for( int b = 1; b < sbox_size && b < sbox_size_limit_to_show; b++ ){
                fprintf(file,"%d", bias_table[a][b]);
                if (b != sbox_size - 1)
                    fprintf(file," & ");
                if (b == sbox_size_limit_to_show - 1 - 1 && sbox_size > sbox_size_limit_to_show){
                    fprintf(file," . . .");
                    break;
                }
            }
            fprintf(file,"\\\\\n");
            fprintf(file,"        \\hline\n       ");
        }
    }


    fprintf(file,"    \\end{tabular}\n");
    fprintf(file,"\\end{center}\n");
    fprintf(file,"* Tabela completa disponível no anexo %s", FILE_NAME_BIAS_CSV);
}

void print_expressions_latex(expressions expres, int n_bits, FILE * file){

    fprintf(file, "\n\n"); 
    fprintf(file, "\\section{Equações}\n");
    fprintf(file, "Em seguida localizamos os maiores bias em módulo e selecionamos as fórmulas relativas a valores de X com peso de hamming = 1, afim de isolarmos um único K. Com isto obtemos %d equações:\n\n", N_BITS);
    fprintf(file, "\\newblock\n\n");
    for (int i = 0 ; i < expres.n_expression ; i ++){    

        if (expres.expressions[i].x < 0)
            continue;

        fprintf(file,"  Para X = %d (hamming: %d) e Y = %d com um bias de %d temos:\n", expres.expressions[i].x, expres.expressions[i].hamming_x, expres.expressions[i].y, expres.expressions[i].bias);

        fprintf(file,"\\[ Pr[");
        int first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    fprintf(file," \\oplus ");
                else
                    first_flag = 0;
                fprintf(file,"X_%d", j);
            }
        }
        for (int j = 1; j < n_bits + 1; j ++)
            if (expres.expressions[i].yn[j] == 1)
                fprintf(file," \\oplus Y_%d", j);
        if (expres.expressions[i].bias < 0)
            fprintf(file," = 1");
        else
            fprintf(file," = 0");
        
        fprintf(file,"] = %.3f", expres.expressions[i].prob);
        fprintf(file,"\\]\n");

        fprintf(file,"\\[");

        first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    fprintf(file," \\oplus ");
                else
                    first_flag = 0;
                fprintf(file,"K_%d", j);
            }
        }

        fprintf(file," = ");

        first_flag = 1;
        for (int j = 1; j < n_bits + 1; j ++){
            if (expres.expressions[i].xn[j] == 1){
                if (!first_flag)
                    fprintf(file," \\oplus ");
                else
                    first_flag = 0;
                fprintf(file,"W_%d", j);
            }
        }

        for (int j = 1; j < n_bits + 1; j ++)
            if (expres.expressions[i].yn[j] == 1)
                fprintf(file," \\oplus Y_%d", j);
        if (expres.expressions[i].bias < 0)
            fprintf(file," \\oplus 1");

        fprintf(file,"\\]\n\n");
    }
}

void create_latex_file(int ** bias_table, int sbox_size, expressions expres, int n_bits, char * fileName, int inte, int * Kdes){
    FILE * file_latex = fopen(fileName, "w");

    fprintf(file_latex, "\\documentclass{article}");
    fprintf(file_latex, "\\title{Criptoanálise Linear}");
    fprintf(file_latex, "\\author{Daniel Veiga da Silva Antunes }");
    fprintf(file_latex, "\\date{Junho 2023}");
    fprintf(file_latex, "\\begin{document}");
    fprintf(file_latex, "\\maketitle");

    print_bias_table_latex(bias_table, sbox_size, n_bits, file_latex);
    print_expressions_latex(expres, n_bits, file_latex);

    fprintf(file_latex, "\\section{Resultados}\n\n");
    fprintf(file_latex, "Foram executadas %d interações para que o seguinte valor de chave fosse obtido:\n", inte);
    fprintf(file_latex, "\n\\[ K = \\{\n");
    for (int i = 1; i < KEY_SIZE + 1; i++){
        fprintf(file_latex, " %d", Kdes[i]);
        if (i != KEY_SIZE)
            fprintf(file_latex, ",");
    }
    fprintf(file_latex, "\\} \\]\n");

    fprintf(file_latex, "\\end{document}");
}



int main() {

    int ** bias_table = get_bias_table(SBOX, SBOX_SIZE, N_BITS);
    print_bias_table_terminal(bias_table, SBOX_SIZE, N_BITS);
    print_bias_table_csv(bias_table, SBOX_SIZE, N_BITS, FILE_NAME_BIAS_CSV);
  
    expressions expression_table = get_expressions(bias_table, SBOX_SIZE, N_BITS);
    print_expressions_terminal(expression_table, N_BITS);

    srand( time(NULL) );
    
    int N = 4000000;

    uchar K[KEY_SIZE] = {7, 6, 11, 2};
    
    uchar Wp[N];
    uchar * Wc;
    
    for( int i = 0; i < N; i++ )
        Wp[i] = rand() & (SBOX_SIZE - 1);
    Wc = e_text(Wp, K, N, KEY_SIZE);
    
    float prob_1_Kij[KEY_SIZE + 1][N_BITS + 1];

    for (int i = 0; i < KEY_SIZE + 1; i ++)
        for (int j = 0; j < N_BITS + 1; j ++)
            prob_1_Kij[i][j] = 0;

    for( int i = 0; i < N; i++ ) {
        uchar W = Wp[i];
        uchar Y = Wc[i];
        
        for (int e = 0; e < expression_table.n_expression; e ++){

            if (expression_table.expressions[e].x < 0)
                continue;
            else {
                for (int x = 1; x < N_BITS + 1; x ++){
                    if (expression_table.expressions[e].xn[x] == 1){
                        int result = bit(W, x);
                        for (int y = 1; y < N_BITS + 1; y ++)
                            if (expression_table.expressions[e].yn[y] == 1)
                                result ^= bit(Y, y);
                        if (expression_table.expressions[e].bias < 0)
                            result ^= 1;

                        prob_1_Kij[(i % KEY_SIZE) + 1][x] += result == 1? 1 : 0;
                    }
                }
            }
        }
    }

    for (int i = 1; i < KEY_SIZE + 1; i++)
        for (int j = 1; j < N_BITS + 1; j++)
                prob_1_Kij[i][j] = (float) (prob_1_Kij[i][j]) / (N/KEY_SIZE);

    
    printf("\n\n");
    printf(" * CHAVE ENCONTRADA\n\n");
    printf("    K = {");
    int Kdes[KEY_SIZE + 1];
    for (int i = 1; i < KEY_SIZE + 1; i++){
        int exp_2 = 0;
        Kdes[i] = 0;
        for (int j = N_BITS; j > 0; j--){
            int bit = prob_1_Kij[i][j] < 0.5 ? 0 : 1;
            if (bit == 1)
                Kdes[i] += pow(2, exp_2);
            exp_2 ++;
        
        }
        printf(" %d ", Kdes[i]);
        if (i != KEY_SIZE)
            printf(",");
        else
            printf("}");
    }
    printf("\n\n");

    create_latex_file(bias_table, SBOX_SIZE, expression_table, N_BITS, "relatorio.tex", N, Kdes);
    return 0;
}
