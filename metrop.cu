/*
    // ОДНОСПИНОВЫЙ МЕТРОПОЛИС С УЧЕТОМ ТЕРМОДИНАМИЧЕСКОЙ ВЕРОЯТНОСТИ //
    
    * принимает .mfsys, расчитывает матрицу энергий
    * алг метрополиса с учетом термодинамической вероятности

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>

#include <cuda.h>

#define threads_in_block 1024

static void HandleError(cudaError_t err, const char *file, int line)  // проверка на ошибку при опрерации с device-памятью на host
{
    if (err != cudaSuccess)
    {
        printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__ ))

using namespace std;

class Timer
{
    private:
        std::chrono::time_point<std::chrono::steady_clock> m_StartTime;     // время, (теориетически) не зависящее от cpu
        std::chrono::time_point<std::chrono::steady_clock> m_EndTime;
        bool m_bRunning = false;
    
    public:
        void start()
        {
            m_StartTime = std::chrono::steady_clock::now();
            m_bRunning = true;
        }

        void stop()
        {
            m_EndTime = std::chrono::steady_clock::now();
            m_bRunning = false;
        }
    
        double Milliseconds()
        {
            std::chrono::time_point<std::chrono::steady_clock> endTime;
            
            if(m_bRunning)
            {
                endTime = std::chrono::steady_clock::now();
            }
            else
            {
                endTime = m_EndTime;
            }
            
            return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
        }
    
        double Seconds()
        {
            return Milliseconds() / 1000.0;
        }
};

// создание матрицы энергий из .mfsys файла
__global__ void matrix_create(double* matrixEn, double* x, double* y, double* mx, double* my, int N)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    // printf("%i на блоке %i\n", index, blockIdx.x);

    double r, Xij, Yij;
    int i = index / N;
    int j = index % N;
    if (index < N * N)
    {
        if (i < j)
        {
            Xij = x[i] - x[j];
            Yij = y[i] - y[j];
            r = sqrt((double)(Xij * Xij + Yij * Yij));

            matrixEn[index] = ((mx[i] * mx[j] + my[i] * my[j]) / (r * r * r) - 3. * (mx[i] * Xij + my[i] * Yij) * (mx[j] * Xij + my[j] * Yij) / (r * r * r * r * r));
        }
    }
}

// расчет энергия и по строкам и расчет полной энергии системы
__global__ void calc_E_cuda(double* matrixEn, uint8_t* spins, double* E_line, double* E, int N)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int i = index / N;
    int j = index % N;

    // расчет энергии по строкам
    if (index < N * N)  
    {
        if (i < j)
        {
            int si = 1, sj = 1;
            double tmp = .0;

            if (spins[i] == 0) si = -1;     // значения в excel инвертированы
            if (spins[j] == 0) sj = -1;     // поменять * -1 в условиях, если надо

            tmp += matrixEn[index] * si * sj;   

            // printf("index: %d  value: %lf  line: %d  colomn: %d\n", index, matrixEn[index], i, j);

            atomicAdd(&E_line[i], tmp);
            atomicAdd(&E_line[j], tmp);
            atomicAdd(E, tmp);

            // atomicAdd(int* address, int val)
            // atomicAdd(E_line + i, tmp);
        }
    }
}


// чтение параметров из консоли (имя .mfsys файла и т.д.)
void console_param_read(int argc, char** argv, string& filename)
{
    if(argc >= 2)    
    {
        for (int i = 1; i < argc; i++)      // to_string() - неправильная конвертация ()
        {
            if (string(argv[i]).find(".mfsys") != string::npos)
                filename = string(argv[i]);
            else cout << "Неправильный формат файла!\n";
        }
    }
}

// чтение результата из .mfsys
bool input_from_mfsys(const string& filename,vector <double>& vx, vector <double>& vy, vector <double>& vmx, vector <double>& vmy)
{
    fstream file(filename);
    if (!file.is_open())
    {
        cout << "File not found!" << endl;
        return false;
    }
    else
    {
        string s;
        for (file >> s; s != "[parts]"; file >> s)
            continue;
        int i, state;
        double x, y, z, mx, my, mz;
        while (!file.eof()){
            if(!(file >> i >> x >> y >> z >> mx >> my >> mz >> state))
                break;
            vx.push_back(x);
            vy.push_back(y);
            vmx.push_back(mx);
            vmy.push_back(my);
        }
        return true;
    }
}

// вывод для проверки
void print_2dmatrix(double* matrixEn, const int& N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << matrixEn[N*i + j] << ", ";
        cout << endl;
    }
}

int main( int argc, char** argv )
{   
    srand(time (NULL));

    string filename = "points376_80.mfsys";
    console_param_read(argc, argv, filename);     // чтение параметрова из консоли

    vector <double> vx, vy, vmx, vmy;
    if (!input_from_mfsys(filename, vx, vy, vmx, vmy)) return -1;   // если файл не найден

    const int N = vx.size();    // размер системы 
    auto matrixEn = new double[N * N]();    // матрица энергий
    double* dev_matrixEn;
    HANDLE_ERROR(cudaMalloc((void**)&dev_matrixEn, sizeof(double) * N * N));
    HANDLE_ERROR(cudaMemcpy(dev_matrixEn, matrixEn, sizeof(double) * N * N, cudaMemcpyHostToDevice));
    
    double * dev_x, * dev_y, * dev_mx, * dev_my;
    HANDLE_ERROR(cudaMalloc((void**)&dev_x, vx.size() * sizeof(double)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_y, vx.size() * sizeof(double)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_mx, vx.size() * sizeof(double)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_my, vx.size() * sizeof(double)));
    HANDLE_ERROR(cudaMemcpy(dev_x, vx.data(), vx.size() * sizeof(double), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_y, vy.data(), vx.size() * sizeof(double), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_mx, vmx.data(), vx.size() * sizeof(double), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_my, vmy.data(), vx.size() * sizeof(double), cudaMemcpyHostToDevice));

    int blocks = (int)(N * N / threads_in_block) + 1;
    matrix_create <<<blocks,threads_in_block>>> (dev_matrixEn, dev_x, dev_y, dev_mx, dev_my, N);    // создание матрицы энергий

    HANDLE_ERROR(cudaMemcpy(matrixEn, dev_matrixEn, sizeof(double) * N * N, cudaMemcpyDeviceToHost));

    // print_2dmatrix(matrixEn, N);
    
    uint8_t spins[N];    // Направление суперспинов
    double E_line[N];    // Энергия в каждом ряду
    for (int i = 0; i < N; i++)
    {
        spins[i] = 0;   // все вниз
        E_line[i] = .0;
    }
    // spins[N-1] = 1;

    // столбцы с плюсовой и минусовой энергией
    vector <int> pos_E_ind;
    vector <int> neg_E_ind;


    uint8_t* dev_spins;
    double* dev_E_line;
    HANDLE_ERROR(cudaMalloc((void**)&dev_spins, sizeof(uint8_t) * N));
    HANDLE_ERROR(cudaMemcpy(dev_spins, spins, sizeof(uint8_t) * N, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc((void**)&dev_E_line, sizeof(double) * N));
    HANDLE_ERROR(cudaMemcpy(dev_E_line, E_line, sizeof(double) * N, cudaMemcpyHostToDevice));


    // auto mexp = new double [N]; 
    double Z = .0, E = 0, E2 = .0, E_aver = .0, E2_aver = .0;
    double * dev_mexp, * dev_Z, * dev_E, * dev_E2;

    HANDLE_ERROR(cudaMemcpy(dev_matrixEn, matrixEn, sizeof(double) * N * N, cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMalloc( (void**)&dev_mexp, sizeof(double) * N));
    HANDLE_ERROR(cudaMalloc( (void**)&dev_Z, sizeof(double)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_E, sizeof(double)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_E2, sizeof(double)));

    HANDLE_ERROR(cudaMemset(dev_E, 0, sizeof(double)));
    HANDLE_ERROR(cudaMemset(dev_E2, 0, sizeof(double)));
    HANDLE_ERROR(cudaMemset(dev_Z, 0, sizeof(double)));

    Timer timer;
    timer.start();

    // Перед матрополисом происходит расчет начального E и массивов, чтобы отталкиваться от него на первой итерации
    calc_E_cuda <<<blocks,threads_in_block>>>(dev_matrixEn, dev_spins, dev_E_line, dev_E, N);

    HANDLE_ERROR(cudaMemcpy( E_line, dev_E_line, sizeof(double)*N, cudaMemcpyDeviceToHost));    // скачивание на cpu
    HANDLE_ERROR(cudaMemcpy( &E, dev_E, sizeof(double), cudaMemcpyDeviceToHost));

    for (int i = 0; i < N; i++)
    {
        if (E_line[i] > 0)
        {
            pos_E_ind.push_back(i);
        }
        else if (E_line[i] < 0)
        {
            neg_E_ind.push_back(i);
        }
        else if (E_line[i] == 0)
        {
            cout << "Энергия строки равна 0! Вероятнa ошибка!" << endl;
            return - 1;
        }

        cout << i << ": " << E_line[i] << endl;
    }
    cout << endl << "E = " << E << endl;

    cout << "Строки с положительными энергиями: ";
    for (int tmp : pos_E_ind)   cout << tmp << " ";     cout << endl;
    cout << "Строки с отрицательными энергиями: ";
    for (int tmp : neg_E_ind)   cout << tmp << " ";     cout << endl;

    // начало метрополиса
    // for (double T = 0.001; T<4; T+=0.1)
    for (double T = 0.001; T<0.101; T+=0.1)
    {
        E_aver = .0; // средние энергии
        E2_aver = .0;

        // for (int MK = 0; MK < N * 100; MK++)
        for (int MK = 0; MK < 1; MK++)
        {
            E2 = .0;

            int rand_spin = rand() % N;
            if (spins[rand_spin] == 1)  spins[rand_spin] = 0;
            else if (spins[rand_spin] == 0)  spins[rand_spin] = 1;      // случайный спин переворачивается

            for (int i = 0; i < N; i++) cout << unsigned(spins[i]) << ", ";  cout << endl;

            for (int i = 0; i < N; i++)    E_line[i] = .0;

            HANDLE_ERROR(cudaMemcpy(dev_spins, spins, sizeof(uint8_t) * N, cudaMemcpyHostToDevice));
            HANDLE_ERROR(cudaMemcpy(dev_E_line, E_line, sizeof(double) * N, cudaMemcpyHostToDevice));
            HANDLE_ERROR(cudaMemset(dev_E, 0, sizeof(double)));

            // 1 расчет E
            calc_E_cuda <<<blocks,threads_in_block>>>(dev_matrixEn, dev_spins, dev_E_line, dev_E, N);

            HANDLE_ERROR(cudaMemcpy( E_line, dev_E_line, sizeof(double)*N, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaMemcpy( &E2, dev_E, sizeof(double), cudaMemcpyDeviceToHost));

            for (int i = 0; i < N; i++)
            {
                cout << i << ": " << E_line[i] << endl;
            }
            cout << endl << "E2 = " << E2 << endl;

            if (E > E2)
            {
                // переворот в зависимости от теромодинамической вероятности
            }
        }

    }

    timer.stop();
    cout << "Время выполнения: " << timer.Milliseconds() << " мсек" << endl;

    return 0;
}