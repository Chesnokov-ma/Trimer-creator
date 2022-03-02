#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include <cmath>


using namespace std;

struct PAIR
{
    double x;
    double y;
    PAIR (double x, double y)   {this->x = x; this->y = y;}
};

// расчет центральной точки
PAIR center_point_find(double x1, double y1, double x2, double y2, double x3, double y3)
{
    PAIR out_PAIR(0, 0);

    double a = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    double b = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));
    double c = sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1));

    double p = (a+b+c)/2;
    double S = 0.5 * abs(((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)));
    double r = S/p;

    out_PAIR.x = ((x1 + x2 - 2*x3)*((y3 - y1)*(x2 + x3 - 2*x1) + x1*(y2 + y3 - 2*y1)) - x3*(x2 + x3 - 2*x1)*(y1 + y2 - 2*y3))/((y2 + y3 - 2*y1)*(x1 + x2 - 2*x3) - (x2 + x3 - 2*x1)*(y1 + y2 - 2*y3));
    out_PAIR.y = ((out_PAIR.x - x1)*(y2 + y3 - 2*y1)/(x2 + x3 - 2*x1)) + y1;

    return out_PAIR;
}

int main(int argc, char** argv)
{
    // --help сообщение и папки вывода
    string help_message = "\nCreate lattice script (Trimer)\nm, n, b - are required\n{optional}\n-f -> output in mfsys/ and dat/ folders (created manualy)\n-s -> create square lattice\n";
    help_message += "-d -> create E matrix (.dat)\n--add -> addition to range (float)\n--check -> check before output\n-h (--help) -> help message\n{examples}\n./trimmer_create -m 7 -n 7 -b 700\n";
    help_message += "./trimer_create -mn 19 -b 850 -s --check\ntrimmer_create -m 7 -n 7 -b 500 --add 1\n";

    string def_folder0 = "", def_folder1 = "";

    bool round_format = true;  // "округленный" формат решетки
    bool dat_create = false;   // создавать ли файл с матрицей энергий

    // параметры
    double b = .0;
    int m = 0, n = 0;

    float addition = .0f;   // увеличить радиус круга на addition
    bool check_bef_out = false;     // спросить, делать ли вывод


    // получение параметров из консоли при запуске
    if(argc < 2)    
    {
        cout << "No arguments\nCan't create lattice" << endl;
        cout << help_message;
        return 0;
    }
    else   // и их проверка
    {
        for (int i = 1; i < argc; i++)
        {
            if (strcmp(argv[i], "-m") == 0)
            {
                m = atoi(argv[i + 1]);
                if ((m < 1 || m > 100) || (i + 1 >= argc)) { cout << "Wrong format {m}.\n"; return -1; }
                else    i++;
            }
            else if (strcmp(argv[i], "-n") == 0)
            {
                n = atoi(argv[i + 1]);
                if ((n < 1 || n > 100) || (i + 1 >= argc)) { cout << "Wrong format {n}.\n"; return -1; }
                else    i++;
            }
            else if (strcmp(argv[i], "-b") == 0)
            {
                b = double(atoi(argv[i + 1]));
                if ((b < 500 || b > 10000) || (i + 1 >= argc)) { cout << "Wrong format {b}.\n"; return -1; }
                else    i++;
            }
            else if ((strcmp(argv[i], "-mn") == 0) || (strcmp(argv[i], "-nm") == 0))
            {
                m = atoi(argv[i + 1]);
                n = atoi(argv[i + 1]);
                if ((n < 1 || n > 100) || (i + 1 >= argc)) { cout << "Wrong format {mn}.\n"; return -1; }
                else    i++;
            }
            else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            {
                cout << help_message;
                return 0;
            }
            else if (strcmp(argv[i], "--add") == 0)
            {
                addition = float(atoi(argv[i + 1]));
                if ((addition == 0) || (i + 1 >= argc))  { cout << "Wrong format {--add}.\nIf --add == 0 don't use it.\n"; return -1; }
                else    i++;
            }
            else if (strcmp(argv[i], "-s") == 0) { round_format = false; }
            else if (strcmp(argv[i], "-f") == 0) { def_folder0 = "mfsys/"; def_folder1 = "dat/"; }
            else if (strcmp(argv[i], "-d") == 0) { dat_create = true; }
            else if (strcmp(argv[i], "--check") == 0)   { check_bef_out = true; }
            else { cout << "Wrong format {else}.\n"; cout << help_message; return -1; }   
        }
    }

    // если указаны не все параметры или парметры n и m не равны
    if (b == .0 || m == 0 || n == 0)    {cout << "Not all parameters are specified.\n"; return -1;}   
    if (!(m == n) && round_format == true)    {cout << "Parameters m, n must be equal (m = n) in round format\n{use -s create lattice}.\n"; return -1;}     
    // else cout << m << ", " << n << ", " << b << endl;


    int N = n * m * 3;    // размер решетки

    //----------------------------------------------------------------------------------------------------------------------------

    // расчет матрицы парных энергий

    const double a = 450./2;
    vector<double> vx, vy, vmx, vmy;
    
    double phi_1 = 60 * M_PI/180.;
    double phi_2 = 180 * M_PI/180.;
    double phi_3 = 300 * M_PI/180.;
    double l = 300.;

    vector<double> tr_vx {a * cos(phi_1), a * cos(phi_2), a * cos(phi_3)};
    vector<double> tr_vy {a * sin(phi_1), a * sin(phi_2), a * sin(phi_3)};
    vector<double> tr_vmx {l * cos(phi_1), l * cos(phi_2), l * cos(phi_3)};
    vector<double> tr_vmy {l * sin(phi_1), l * sin(phi_2), l * sin(phi_3)};

    double dx = b, dy = b * sin(phi_1), offset = b * cos(phi_1);

    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < m; ++j)
        {
            if (i % 2 == 0)
            {
                for (auto k: tr_vx)
                    vx.push_back(k + dx * j);
                for (auto k: tr_vy)
                    vy.push_back(k + dy * i);
                for (auto k: tr_vmx)
                    vmx.push_back(k);
                for (auto k: tr_vmy)
                    vmy.push_back(k);
            }
            else
            {
                for (auto k: tr_vx)
                    vx.push_back(offset + k + dx * j);
                for (auto k: tr_vy)
                    vy.push_back(k + dy * i);
                for (auto k: tr_vmx)
                    vmx.push_back(k);
                for (auto k: tr_vmy)
                    vmy.push_back(k);
            }
        }
    }

    vector <int> in_range_spins;    // спины в радиусе
    in_range_spins.clear();

    if (round_format)   // для круглых решеток
    {
        PAIR center_point(0, 0);
        float b_f = b + addition;

        if (b == 500)   b_f += 1.0;   // увеличение радиуса для низких b
        if (b == 545)   b_f += 2.0;

        if (m % 2 != 0)
        {
            int ct = (int)(N / 2) + 1;  // треугольник в середине
            center_point = center_point_find(vx[ct],vy[ct],vx[ct-1],vy[ct-1],vx[ct-2],vy[ct-2]);   // расчет центральной точки для обрезания спинов вне радиуса
            // cout << center_point.x << " " << center_point.y << endl;
            // cout << ct << endl;

            int range = b_f * (int)(m / 2);   // масштабирование радиуса

            for (int i = 0; i < N; i++)     // поиск входящих в радиус диполей
            {
                if (((center_point.x - vx[i]) * (center_point.x - vx[i]) + (center_point.y - vy[i]) * (center_point.y - vy[i])) <= (range * range))
                {
                    in_range_spins.push_back(i);
                }
            }
            
            N = in_range_spins.size();  // кол-во входящих в радиус диполей
        }
        else    // если система с четным числом столбцов
        {
            int ct = (int)(N / 2) + 2 + 6;  // треугольник в середине
            center_point = center_point_find(vx[ct],vy[ct],vx[ct-1],vy[ct-1],vx[ct-2],vy[ct-2]);   // расчет центральной точки для обрезания спинов
            // cout << center_point.x << " " << center_point.y << endl;
            // cout << ct << endl;

            int range = b_f * (int)(m * 0.5);   // масштабирование радиуса ( можно попроб менять кооф, на который умножается m )

            // UPD: все равно лучше делать через нечетные

            for (int i = 0; i < N; i++)     // поиск входящих в радиус диполей
            {
                if (((center_point.x - vx[i]) * (center_point.x - vx[i]) + (center_point.y - vy[i]) * (center_point.y - vy[i])) <= (range * range))
                {
                    in_range_spins.push_back(i);
                }
            }
            
            N = in_range_spins.size();  // кол-во входящих в радиус диполей
        }
    }

    // проверка перед выводом
    string tmp_str = "Y";
    cout << "N = " << N << endl;
    if (check_bef_out)
    {
        cout << "Create .mfsys file? (y/n): ";
        cin >> tmp_str;
        if (strcmp(tmp_str.c_str(), "n") == 0 || strcmp(tmp_str.c_str(), "N") == 0)
            return 0;
    }
    if (strcmp(tmp_str.c_str(), "y") != 0 && strcmp(tmp_str.c_str(), "Y") != 0)
        return 0;

    //.mfsys
    FILE* f;
    string filename = def_folder0 + "trimer_N" + to_string(N) + "_b" + to_string(int(b)) + ".mfsys";
    f = fopen(filename.c_str(), "w");

    fprintf(f, "[header]\n");
    fprintf(f, "version=2\n");
    fprintf(f, "dimensions=2\n");
    fprintf(f, "type=standart\n");
    fprintf(f, "size=%llu\n", N);
    fprintf(f, "state=");
    for (int i = 0; i < N; i++)
        fprintf(f, "0");
    fprintf(f, "\nminstate=");
    for (int i = 0; i < N; i++)
        fprintf(f, "0");
    fprintf(f, "\nminscale=");
    for (int i = 0; i < N; i++)
        fprintf(f, "0");
    fprintf(f, "\ninteractionrange=0\n");
    fprintf(f, "sizescale=1\n");
    fprintf(f, "magnetizationscale=1\n");
    fprintf(f, "[parts]\n");
    int ind = 0;
    if (round_format)
    {
        for (int i = 0; i < N; i++)
        {
            int k = in_range_spins[i];

            fprintf(f, "%i\t%.16e\t%.16e\t0\t%.16e\t%.16e\t0\t0", ind, vx[k], vy[k], vmx[k], vmy[k]);
            if (ind != N - 1)
                fprintf(f, "\n");

            ind++;
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(f, "%i\t%.16e\t%.16e\t0\t%.16e\t%.16e\t0\t0", ind, vx[i], vy[i], vmx[i], vmy[i]);
            if (ind != N - 1)
                fprintf(f, "\n");

            ind++;
        }
    }
    
    fclose(f);

    if (!dat_create)
    {
        cout << "File succesfully created {.mfsys}.\n";
        return 0;
    }

    auto* matrix = new double[N * N]();     // создаем матрицу

    int index = 0;
    double r, Xij, Yij;
    double Ei = .0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                Xij = vx[i] - vx[j];
                Yij = vy[i] - vy[j];
                r = sqrt((double)(Xij * Xij + Yij * Yij));

                Ei = ((vmx[i] * vmx[j] + vmy[i] * vmy[j]) / (r * r * r) - 3. * (vmx[i] * Xij + vmy[i] * Yij) * (vmx[j] * Xij + vmy[j] * Yij) / (r * r * r * r * r));

                matrix[index] = Ei;
            }
            else
            {
                matrix[index] = .0;
            }

            index++;
        }
    }

    // .dat
    FILE* f_matrix;     
    string filename_1 = def_folder1 + "trimer_N" + to_string(vx.size()) + "_b" + to_string(int(b)) + ".dat" ;
    f_matrix = fopen(filename_1.c_str(), "w");
    for (int i = 0; i < N * N; i++)
    {
        if ((i + 1) % N == 0)
        {
            fprintf(f_matrix, "%lf\n",  matrix[i]);
        }
        else
        {
            fprintf(f_matrix, "%lf ",  matrix[i]);
        }
    }
    fclose(f_matrix);

    cout << "Files succesfully created {.mfsys, .dat}.\n";

    return 0;
}