#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <ctime>
#include <sstream>
#include <math.h>

using namespace std;

mt19937 generator{ random_device{}() };

struct TRange
{
    float T1;
    float T2;
    float dT;
    int aSteps;
    int mSteps;
};


struct ParamsNow
{
    int aSteps;
    int mSteps;
    float T;
    float Js;
};
//��������� ��������� �� �����
struct InputParams {
    int L;
    float J;
    float h;
    float Js;
    float JsEnd;
    float dJs;
    int Rep;
    vector<TRange> TRanges;
    void read(const std::string fileName) {
        ifstream ifs(fileName, ifstream::in);
        if (!ifs) {
            cout << "11111111" << endl;
            return;
        }
        string name;
        ifs >> name >> L; cout << name << " = " << L << endl;
        ifs >> name >> J; cout << name << " = " << J << endl;
        ifs >> name >> h; cout << name << " = " << h << endl;
        ifs >> name >> Js; cout << name << " = " << Js << endl;
        ifs >> name >> JsEnd; cout << name << " = " << JsEnd << endl;
        ifs >> name >> dJs; cout << name << " = " << dJs << endl;
        ifs >> name >> Rep; cout << name << " = " << Rep << endl;
        int count; //���-�� ����������
        ifs >> name >> count;
        for (int i = 0; i < count; i++) {
            TRange r;
            ifs >> r.T1 >> r.T2 >> r.dT >> r.aSteps >> r.mSteps;
            TRanges.push_back(r);
        }
    }
};

struct CalculatedParams {//��������� ������ ����������� ���������� ��������
    float energy;
    float magneticMoment;
    float afmMoment; // AFM �������� �������
    float AFMHorLine;
    float AFMVertLine;
};

struct Observables {//��������� ������� ����������� ����� ����������� �������� 1 �������
    double U;
    double C;
    double Magn;
    double MagnAFM;
    double MagnSumm;
    double X;
    double Xafm;
    double Xafmvl;
    double Xafmhl;
    vector<int> latticeNow; 
};

//periodic boudery condition (���)
int pbc(int x, int L) {
    if (x < 0)
        return x + L;
    return x % L;
}

//���������� �������� ���������� ��� ��������� ���������
CalculatedParams computeState(vector<int>& lattice, const InputParams& ip, const float Js) {//����� ��� ��������� ���������� �������
    int L = ip.L;
    float E = 0;
    float M = 0;
    float AFM = 0;
    float AFMVertLine = 0;
    float AFMHorLine = 0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            int S = lattice[x * L + y];
            int Sr = lattice[pbc(x + 1, L) * L + y];//right
            int Sb = lattice[x * L + pbc(y + 1, L)]; //bottom
            int Srb = lattice[pbc(x - 1, L) * L + pbc(y + 1, L)];//rightBottom
            int Srt = lattice[pbc(x + 1, L) * L + pbc(y + 1, L)];//rightTop
            E += -ip.J * S * Sr;
            E += -ip.J * S * Sb;
            E += -Js * S * Srt;
            E += -Js * S * Srb;
            E += -ip.h * S;
            M += S;
            AFM += ((x + y) % 2 == 0) ? S : -S;
            AFMVertLine += (y % 2 == 0) ? S : -S;
            AFMHorLine += (x % 2 == 0) ? S : -S;
        }
    }
    return { E,M,AFM,AFMVertLine,AFMHorLine };
}

//��� � ��������� ����������� ��������� 1 �����
void mcStep(vector<int>& lattice, float T, const InputParams& ip, const float Js, CalculatedParams& op, int x, int y) {
    int L = ip.L;
    uniform_real_distribution<float> realDist(0, 1); //��. ����� �� 0 �� 1
    int S0 = lattice[x * L + y];
    int S1 = -S0; // ����� �������� ����� �� ���� (x,y)
    //������
    int Sl = lattice[pbc(x - 1, L) * L + y]; //left
    int Sr = lattice[pbc(x + 1, L) * L + y];//right
    int St = lattice[x * L + pbc(y - 1, L)]; //top
    int Sb = lattice[x * L + pbc(y + 1, L)]; //bottom
    int Slt = lattice[pbc(x + 1, L) * L + pbc(y-1,L)]; 
    int Slb = lattice[pbc(x - 1, L) * L + pbc(y - 1, L)];
    int Srt = lattice[pbc(x + 1, L) * L + pbc(y + 1, L)]; 
    int Srb = lattice[pbc(x - 1, L) * L + pbc(y + 1, L)]; 
    float dE = -ip.J * (Sl + Sr + St + Sb) * (S1 - S0) - Js * (Slt + Slb + Srt + Srb) * (S1 - S0) - ip.h * (S1 - S0);
    if (realDist(generator) < expf(-dE / T)) {
        //��������� ���������
        lattice[x * L + y] = S1;
        op.energy += dE;
        int dM = S1 - S0;
        op.magneticMoment += dM;
        op.afmMoment += ((x + y) % 2 == 0) ? dM : -dM;
        op.AFMVertLine += (y % 2 == 0) ? dM : -dM;
        op.AFMHorLine += (x % 2 == 0) ? dM : -dM;
    }
}
//� ������ ������ �� ����������
void macroMCStepLegacy(vector<int>& lattice, float T, const InputParams& ip,const float Js, CalculatedParams& op) {
    int L = ip.L;
    uniform_int_distribution<int> intDist(0, L - 1); //��������� ����
    uniform_real_distribution<float> realDist(0, 1); //��. ����� �� 0 �� 1
    for (int i = 0; i < L * L; i++) {
        // � ������� ������� ���� ���� ���� �������� ���������
        int x = intDist(generator);
        int y = intDist(generator);
        mcStep(lattice, T, ip,Js, op, x, y);
    }
}

//������ ��������� �����������
void macroMCStep(vector<int>& lattice, float T, const InputParams& ip,const float Js, CalculatedParams& op, int eps) {
    int L = ip.L;
    for (int x = 0; x < L; x++) {
        for (int i = 0; i < L / 2; i++) {
            int y = pbc(2 * i - x + eps, L);
            mcStep(lattice, T, ip,Js, op, x, y);
        }
    }
}
//��� ���������� ������ ��� ������������ ������� ����� � ������
//https://en.wikipedia.org/wiki/Kahan_summation_algorithm
void addKahan(double& sum, double& c, double value) {
    auto y = value - c;    // ���� ��� ������: c - ����.
    auto t = sum + y;         // ���, sum ������, y ����, ��� ��� ������� ������� y ��������.
    c = (t - sum) - y;   // (t - sum) ��������������� ������� ������� y; ��������� y ��������������� ������� ������� y
    sum = t;
}
//��� �������������� ������ �������� ������
template <typename T>
std::string toString(const T value, const int n = 4)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << value;
    return out.str();
}
//��� ��������� ������ ������ � ����
template <typename type>
void writeBinary(ofstream& stream, type number) {
    stream.write(reinterpret_cast<char*>(&number), sizeof(type));
}

//������ ��� 1 �������
Observables computeObservables(const vector<int> latticeInit, const CalculatedParams cp, const ParamsNow& pn, const InputParams& ip) {
    double E = 0, E2 = 0, M = 0, M2 = 0, AFM = 0, AFM2 = 0, AFMVL = 0, AFMVL2 = 0, AFMHL = 0, AFMHL2 = 0;
    double Eerr = 0, E2err = 0, Merr = 0, M2err = 0, AFMerr = 0, AFM2err = 0, AFMVLerr = 0, AFMVL2err = 0, AFMHLerr = 0, AFMHL2err = 0;
    vector<int> lattice = latticeInit;
    CalculatedParams calcParams = cp;
    int L = ip.L;
    int N = L * L;
    float T = pn.T;
    float Js = pn.Js;
    int aSteps = pn.aSteps;
    int mSteps = pn.mSteps;

    for (int step = 0; step < aSteps + mSteps; step++) {
        // macroMCStepLegacy(lattice, T, params, calcParams);
        macroMCStep(lattice, T, ip, Js, calcParams, 0);
        macroMCStep(lattice, T, ip, Js, calcParams, 1);
        if (step >= aSteps) {
            //�������� - ����������� ����������
            double e = calcParams.energy;
            double m = calcParams.magneticMoment;
            double afm = calcParams.afmMoment;
            double afmvl = calcParams.AFMHorLine;
            double afmhl = calcParams.AFMVertLine;
            addKahan(E, Eerr, e / mSteps);
            addKahan(E2, E2err, e * e / mSteps);
            addKahan(M, Merr, m / mSteps);
            addKahan(M2, M2err, m * m / mSteps);
            addKahan(AFM, AFMerr, afm / mSteps);
            addKahan(AFM2, AFM2err, afm * afm / mSteps);
            addKahan(AFMVL, AFMVLerr, afmvl / mSteps);
            addKahan(AFMVL2, AFMVL2err, afmvl * afmvl / mSteps);
            addKahan(AFMHL, AFMHLerr, afmhl / mSteps);
            addKahan(AFMHL2, AFMHL2err, afmhl * afmhl / mSteps);
        }
    }
    //������� ����������������� �������
    double U = E / N;
    double C = (E2 - E * E) / N / T / T;
    double Magn = std::abs(M / N);
    double MagnAFM = std::abs(AFM / N);
    double MagnAFMVL = std::abs(AFMVL / N);
    double MagnAFMHL = std::abs(AFMHL / N);
    double MagnSumm = MagnAFMVL + MagnAFMHL;
    double X = (M2 - M * M) / N / T;
    double Xafm = (AFM2 - AFM * AFM) / N / T;
    double Xafmvl = (AFMHL2 - AFMHL * AFMHL) / N / T;
    double Xafmhl = (AFMVL2 - AFMVL * AFMVL) / N / T;
    return {U,C,Magn,MagnAFM ,MagnSumm,X,Xafm,Xafmvl,Xafmhl,lattice };
}


int main(int argc, char** argv)//��������� ��� ������� ����� �������
{
    //��������� ������
    InputParams params;
    ParamsNow paramsnow;
    string paramsFile = "params.txt";
    if (argc > 1) {
        paramsFile = argv[1];
    }
    params.read(paramsFile);//������������� ���������� �� �����
    int L = params.L;
    float J = params.J;
    float h = params.h;
    int N = L * L;
    int countJs = round((params.JsEnd- params.Js)/params.dJs+1);//���������� �������� ��������� ������ 2 �������
    int Rep = params.Rep;

    string prefix = "L" + toString(L) + "_J" + toString(J) + "_Js" + toString(params.Js) + "_JsEnd" + toString(params.JsEnd) + "_h" + toString(h);
    ofstream ofs(prefix + "_results.csv", ofstream::out); //����� ������ � ����
    ofs << fixed;
    ofstream ofsLattice(prefix + "_lattice.bin", ofstream::binary);
    writeBinary(ofsLattice, L);
    writeBinary(ofsLattice, J);
    writeBinary(ofsLattice, h);
    writeBinary(ofsLattice, countJs);//���������� ��������� � ����
    ofs << countJs << ";" << params.Js<< ";"<< params.dJs<<endl;

    auto start = clock();

    for (int Jstep = 1; Jstep <= countJs;Jstep++) {//���� �� ��������� ������ 2 �������
        float Js = params.Js + (Jstep - 1) * params.dJs;
        vector<int> lattice(N, 1);//�������� �������//����������� ��������� ������� ��� ������ �����������
        writeBinary(ofsLattice, Js);//������ �������� ������ 2 �������
        ofs << Js << endl;
        paramsnow.Js = Js;
       
 
        cout << fixed;
        for (const auto& TRange : params.TRanges) {//���� ��� ������ �������� ���������� ���� ���-�� �����
            float T1 = TRange.T1;
            float T2 = TRange.T2;
            float dT = TRange.dT;
            paramsnow.aSteps = TRange.aSteps;
            paramsnow.mSteps = TRange.mSteps;
            // ������� ������ ������� � ���. ������ ��� ���������� ���������
            CalculatedParams calcParams = computeState(lattice, params,Js); // E,M, Maf....

            //uniform_int_distribution<int> intDist(0, L - 1); //��������� �������� �� ��� ���
            for (paramsnow.T = T1; paramsnow.T > T2 - dT / 2; paramsnow.T -= dT) {
                double U=0, C=0, Magn=0, MagnAFM=0, MagnSumm=0, X=0, Xafm=0, Xafmvl=0, Xafmhl=0;
                vector<int> latticeNow;
                

                for (int i = 1; i <= Rep; i++) {//���������� �� ��������
                    Observables obs = computeObservables(lattice, calcParams, paramsnow, params);
                    U += obs.U / Rep;
                    C += obs.C / Rep;
                    Magn += obs.Magn / Rep;
                    MagnAFM += obs.MagnAFM / Rep;
                    MagnSumm += obs.MagnSumm / Rep;
                    X += obs.X / Rep;
                    Xafm += obs.Xafm / Rep;
                    Xafmvl += obs.Xafmvl / Rep;
                    Xafmhl += obs.Xafmhl / Rep;
                    latticeNow = obs.latticeNow;
                }

                cout << "T = " << paramsnow.T << " ; E = " << U << " ; M = " << Magn << " ; Mafm = " << MagnAFM << " ; MagnSumm = " << MagnSumm << " ; C = " << C << " ; X = " << X << " ; Xafm = " << Xafm << " ; Xafmvl = " << Xafmvl << "; Xafmhl = " << Xafmhl << endl;
                ofs << paramsnow.T << ";" << U << ";" << Magn << ";" << MagnAFM << ";" << MagnSumm << ";" << C << ";" << X << ";" << Xafm << ";" << Xafmvl << ";" << Xafmhl << endl;
                ofs.flush();
                //��������� � ���� T � ��������� ������� (������ �����)
                writeBinary(ofsLattice, paramsnow.T);
                for (int x = 0; x < L; x++) {
                    for (int y = 0; y < L; y++)
                        writeBinary(ofsLattice, latticeNow[x * L + y]);
                }
            }
        }
    }
    ofs.close();
    ofsLattice.close();
    cout << "time: " << (clock() - start) / 1e3 << " s" << endl;
    return 0;
}
