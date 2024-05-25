#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <ctime>
#include <sstream>
#include <chrono>

#if defined(_OPENMP)
#include <omp.h>
#define OMP_TREADS 8
#endif

using namespace std;
using Lattice = vector<int>;

struct TRange
{
    float T1;
    float T2;
    float dT;
    int aSteps;
    int mSteps;
};

struct InputParams {
    int L;
    float J;
    float Js;
    float h;
    int copies;
    vector<TRange> TRanges;
    int JsCount;
    float dJs;
};

struct CalculationJob {
    vector<InputParams> inputParams;
    bool read(const std::string& fileName) {
        ifstream ifs(fileName, ifstream::in);
        if (!ifs) {
            cout << "can't open the parameters file" << endl;
            return false;
        }
        string name;
        int L;
        int copies;
        float J;
        float h;
        float JsStart;
        float JsEnd;
        int JsCount;
        ifs >> name >> L; cout << name << " = " << L << endl;
        ifs >> name >> copies; cout << name << " = " << copies << endl;
        ifs >> name >> J; cout << name << " = " << J << endl;
        ifs >> name >> h; cout << name << " = " << h << endl;
        ifs >> name >> JsStart; cout << name << " = " << JsStart << endl;
        ifs >> name >> JsEnd; cout << name << " = " << JsEnd << endl;
        ifs >> name >> JsCount; cout << name << " = " << JsCount << endl;
        int TRangesCount;
        ifs >> name >> TRangesCount;
        vector<TRange> TRanges;
        for (int i = 0; i < TRangesCount; i++) {
            TRange r;
            ifs >> r.T1 >> r.T2 >> r.dT >> r.aSteps >> r.mSteps;
            TRanges.push_back(r);
        }

        for (int JsIdx = 0; JsIdx < JsCount; JsIdx++) {
            float Js = JsCount > 1 ? JsStart + JsIdx * (JsEnd - JsStart) / (JsCount - 1) : JsStart;
            float dJs = JsCount > 1 ? (JsEnd - JsStart) / (JsCount - 1) : 0;
            inputParams.push_back({ L,J,Js,h,copies,TRanges,JsCount,dJs  });
        }
        ifs.close();
        return true;
    }
};

struct CalculatedParams {
    float energy;
    float magneticMoment;
    float afmMoment; // AFM параметр порядка
    float afmLines;
};


int pbc(int x, int L) {
    //periodic boudery condition (ПГУ)
    if (x < 0)
        return x + L;
    return x % L;
}

//вычисление вычисление выходных параметров для заданного состояния
CalculatedParams computeState(vector<int>& lattice, const InputParams& ip) {
    int L = ip.L;
    float E = 0;
    float M = 0;
    float AFM = 0;
    float afmLines = 0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            int S = lattice[x * L + y];
            int Sr = lattice[pbc(x + 1, L) * L + y];//right
            int Sb = lattice[x * L + pbc(y + 1, L)]; //bottom
            int Srb = lattice[pbc(x - 1, L) * L + pbc(y + 1, L)];//rightBottom
            int Srt = lattice[pbc(x + 1, L) * L + pbc(y + 1, L)];//rightTop
            E += -ip.J * S * Sr;
            E += -ip.J * S * Sb;
            E += -ip.Js * S * Srt;
            E += -ip.Js * S * Srb;
            E += -ip.h * S;
            M += S;
            AFM += ((x + y) % 2 == 0) ? S : -S;
            afmLines += (x % 2 == 0) ? S : -S;
            afmLines += (y % 2 == 0) ? S : -S;
        }
    }
    return { E,M,AFM,afmLines };
}

void mcStep(Lattice& lattice, float T, const InputParams& ip, CalculatedParams& op, mt19937& generator, int x, int y) {
    int L = ip.L;
    uniform_real_distribution<float> realDist(0, 1); //сл. число от 0 до 1
    int S0 = lattice[x * L + y];
    int S1 = -S0; // новое значение спина на узле (x,y)
    //соседи
    int Sl = lattice[pbc(x - 1, L) * L + y]; //left
    int Sr = lattice[pbc(x + 1, L) * L + y];//right
    int St = lattice[x * L + pbc(y - 1, L)]; //top
    int Sb = lattice[x * L + pbc(y + 1, L)]; //bottom
    int Slt = lattice[pbc(x + 1, L) * L + pbc(y - 1, L)];
    int Slb = lattice[pbc(x - 1, L) * L + pbc(y - 1, L)];
    int Srt = lattice[pbc(x + 1, L) * L + pbc(y + 1, L)];
    int Srb = lattice[pbc(x - 1, L) * L + pbc(y + 1, L)];
    float dE = -ip.J * (Sl + Sr + St + Sb) * (S1 - S0) - ip.Js * (Slt + Slb + Srt + Srb) * (S1 - S0) - ip.h * (S1 - S0);
    if (realDist(generator) < expf(-dE / T)) {
        //принимаем состояние
        lattice[x * L + y] = S1;
        op.energy += dE;
        int dM = S1 - S0;
        op.magneticMoment += dM;
        op.afmMoment += ((x + y) % 2 == 0) ? dM : -dM;
        op.afmLines += (y % 2 == 0) ? dM : -dM;
        op.afmLines += (x % 2 == 0) ? dM : -dM;
    }
}
//https://en.wikipedia.org/wiki/Kahan_summation_algorithm
void addKahan(double& sum, double& c, double value) {
    auto y = value - c;    // Пока все хорошо: c - ноль.
    auto t = sum + y;         // Увы, sum велика, y мало, так что младшие разряды y потеряны.
    c = (t - sum) - y;   // (t - sum) восстанавливает старшие разряды y; вычитание y восстанавливает младшие разряды y
    sum = t;
}

template <typename T>
std::string toString(const T value, const int n = 4)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << value;
    return out.str();
}
//Для бинарного вывода данных в файл
template <typename type>
void writeBinary(ofstream& stream, type number) {
    stream.write(reinterpret_cast<char*>(&number), sizeof(type));
}

void macroMCStep(Lattice& lattice, float T, const InputParams& ip, CalculatedParams& op, mt19937& generator) {
    int L = ip.L;
    uniform_int_distribution<int> intDist(0, L - 1); //случайные узлы
    uniform_real_distribution<float> realDist(0, 1); //сл. число от 0 до 1
    for (int i = 0; i < L * L; i++) {
        // в среднем каждому узлу дали шанс поменять состояние
        int x = intDist(generator);
        int y = intDist(generator);
        mcStep(lattice, T, ip, op, generator, x, y);
    }
}

void compute(const InputParams& params, int tid) {
    int L = params.L;
    float J = params.J;
    float Js = params.Js;
    float h = params.h;
    float JsCount= params.JsCount;
    float dJs = params.dJs;
    int N = L * L;
    int K = params.copies; // количество реплик (копий системы)

    mt19937 generator{ random_device{}() };

    string prefix = "L" + toString(L) + "_J" + toString(J) + "_JS" +toString(Js) + "_h" + toString(h) + "_K" + toString(K);
    ofstream ofs(prefix + "_results.csv", ofstream::out); //поток вывода в файл
    ofs << fixed;
    ofstream ofsLattice(prefix + "_lattice.bin", ofstream::binary);
    writeBinary(ofsLattice, L);
    writeBinary(ofsLattice, J);
    writeBinary(ofsLattice, Js);
    vector<Lattice> lattices(K, Lattice(N, 1));
    ofs  << Js  << endl;

    cout << fixed;
    for (const auto& TRange : params.TRanges) {
        float T1 = TRange.T1;
        float T2 = TRange.T2;
        float dT = TRange.dT;
        int aSteps = TRange.aSteps; // шаги на узел
        int mSteps = TRange.mSteps;
        for (auto T = T1; T > T2 - dT / 2; T -= dT) {
            double avU = 0;
            double avC = 0;
            double avMagn = 0;
            double avMagnAFM = 0;
            double avMagnAFMLines = 0;
            double avX = 0;
            double avXafm = 0;
            double avXafmLines = 0;
            vector<double> U(K);
            vector<double> C(K);
            vector<double> Magn(K);
            vector<double> MagnAFM(K);
            vector<double> MagnAFMLines(K);
            vector<double> X(K);
            vector<double> Xafm(K);
            vector<double> XafmLines(K);
            for (int k = 0; k < K; k++)
            {
                double E = 0, E2 = 0, M = 0, M2 = 0, AFM = 0, AFM2 = 0,AFMLines = 0, AFMLines2 = 0;
                double Eerr = 0, E2err = 0, Merr = 0, M2err = 0, AFMerr = 0, AFM2err = 0, AFMLineserr = 0, AFMLines2err = 0;
                auto& lattice = lattices[k];
                // считаем полную энергию и маг. момент для начального состояния
                auto calcParams = computeState(lattices[k], params); // E,M, Maf....
                for (int step = 0; step < aSteps + mSteps; step++) {
                    macroMCStep(lattice, T, params, calcParams, generator);
                    if (step >= aSteps) {
                        //измеряем - накапливаем статистику
                        double e = calcParams.energy;
                        double m = calcParams.magneticMoment;
                        double afm = calcParams.afmMoment;
                        double afmLines = calcParams.afmLines;
                        addKahan(E, Eerr, e / mSteps);
                        addKahan(E2, E2err, e * e / mSteps);
                        addKahan(M, Merr, m / mSteps);
                        addKahan(M2, M2err, m * m / mSteps);
                        addKahan(AFM, AFMerr, afm / mSteps);
                        addKahan(AFM2, AFM2err, afm * afm / mSteps);
                        addKahan(AFMLines, AFMLineserr, afmLines / mSteps);
                        addKahan(AFMLines2, AFMLines2err, afmLines * afmLines / mSteps);
                    }
                }
                //считаем термодинамические средние
                U[k] = E / N;
                C[k] = (E2 - E * E) / N / T / T;
                Magn[k] = std::abs(M / N);
                MagnAFM[k] = std::abs(AFM / N);
                MagnAFMLines[k] = std::abs(AFMLines / N);
                X[k] = (M2 - M * M) / N / T;
                Xafm[k] = (AFM2 - AFM * AFM) / N / T;
                XafmLines[k] = (AFMLines2 - AFMLines * AFMLines) / N / T;
                for (int k = 0; k < K; k++) {
                    //вычисляем средние по репликам
                    avU += U[k] / K;
                    avC += C[k] / K;
                    avMagn += Magn[k] / K;
                    avMagnAFM += MagnAFM[k] / K;
                    avMagnAFMLines += MagnAFMLines[k] / K;
                    avX += X[k] / K;
                    avXafm += Xafm[k] / K;
                    avXafmLines += XafmLines[k] / K;
                }
            }
            if (tid == 0) {
                cout << "T = " << T << " ; E = " << avU << " ; M = " << avMagn << " ; Mafm = " << avMagnAFM << " ; MafmLines = " << avMagnAFMLines << " ; C = " << avC << " ; X = " << avX << " ; Xafm = " << avXafm << " ; XafmLines = " << avXafmLines << endl;
            }
            ofs << T << ";" << avU << ";" << avMagn << ";" << avMagnAFM << ";" << avMagnAFMLines << ";" << avC << ";" << avX << ";" << avXafm << ";" << avXafmLines << endl;
            ofs.flush();
            //сохраняем в файл T и состояние решетки (первую копию)
            writeBinary(ofsLattice, T);
            for (int x = 0; x < L; x++) {
                for (int y = 0; y < L; y++)
                    writeBinary(ofsLattice, lattices[0][x * L + y]);
            }
        }
    }
    ofs.close();
    ofsLattice.close();

}

int main(int argc, const char** argv)
{
    string paramsFile = "params_generalized.txt";
    if (argc > 1) {
        paramsFile = argv[1];
    }
    cout << paramsFile << endl;
    CalculationJob job;
    if (!job.read(paramsFile)) {
        cout << "can't open file " << paramsFile << endl;
        return 0;
    }
#if defined(_OPENMP)
    const int P = OMP_TREADS;
#else
    const int P = 1; // количество потоков
#endif

    int jobsCount = job.inputParams.size();
    auto start = std::chrono::high_resolution_clock::now();
#if defined(_OPENMP)
#pragma omp parallel num_threads(P)
#endif
    {
        int blocks = (jobsCount - 1) / P + 1;
        int tid = 0;
    #if defined(_OPENMP)
        tid = omp_get_thread_num();
    #endif

        for (int b = 0; b < blocks; b++) {
            int j = tid + b * P;
            if (j >= jobsCount)
                break;
            //cout << "computing h = " << job.inputParams[j].h << endl;
            compute(job.inputParams[j], tid);
        }

    }
    auto finish = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    cout << "time: " << milli.count() / 1000.0 << " s" << endl;
    return 0;
}
