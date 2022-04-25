#include <iostream>
#include "cmath"

using namespace std;

#define e 1e-6
#define v2

void blad(double *x, double *y) {
    double temp[4];
    double blad;

    for (int i = 0; i < 4; ++i) {
        temp[i] = std::fabs(x[i] - y[i]);
    }

    blad = temp[0];

    for (int i = 1; i < 4; ++i) {
        if (temp[i] > blad)
            blad = temp[i];
    }
    cout << "////BLAD////" << endl;
    cout << blad << endl;
}

void wypelnijMacierz(double **macierz,int  N){
#ifdef v1
    double pom[4][4] =
                            {{1.0,-20.0,30.0,-4.0},
                            {2.0,-40.0,-6.0,50.0},
                            {9.0,-180.0,11.0,-12.0},
                            {-16.0,15.0,-140.0,13.0}};
#endif

#ifdef v2
    double E = e;
  double pom[4][4] =
          {{1.0 + E, 1.0, 1.0, 1.0},
       {1.0, 1.0 + E, 1.0, 1.0},
       {1.0, 1.0, 1.0 + E, 1.0},
       {1.0, 1.0, 1.0, 1.0 + E}};
#endif
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            macierz[i][j] = pom[i][j];
        }
    }
}

void printMacierz(double** macierz, int N, int* index){
    for(int i = 0; i < N; i++){
        for(int j = 0; j<N; j++){
            cout << macierz[index[i]][j] << "\t";
        }
        cout << endl;
    }
}

void printWektor(double *d){
    for(int i = 0; i<4;i++)
        cout << d[i] << endl;
}

void printWektor(double *d, int *index){
    for (int i = 0; i < 4; i++)
        cout << d[index[i]] << endl;
}

/**
 * funkcja wybierajaca elemant podstawowy przy dekompozycji Gaussa
 * @param macierz macierz na ktorej funkcja operuje
 * @param n indeks od ktorego funkcja zaczyna szukac
 * @param N wielkosc macierzy
 * @param index wektor indeksow
 * @return zwracany jest wybrany element, najwyzszy z wszystkich liczb znajdujacych sie ponizej startowej wartosci, w tej samej kolumnie
 */
int wyborElemPodst(double **macierz, int n, int N, int* index){
    int ret;
    for(int i = n; i<N-1; i++){
        if(fabs(macierz[index[i]][n]) < fabs(macierz[index[i+1]][n])){
            ret = index[i+1];
        }
        else{
            ret = index[i];
        }
    }
    return ret;
}

/**
 * funkcja wykonuje dekompozycje LU na macierzy przy pomocy dekompozycji Gaussa,
 * L i U znajdują się w jednej macierzy - glowna przekatna i trojkat ponad nia to macierz U,
 * natomiast reszta to macierz L
 * @param macierz macierz na ktorej wykonywana jest dekompozycja
 * @param index wektor indeksow potrzebna do wyboru elementow podstawowych
 * @param N wielkosc macierzy
 */
void gauss(double **macierz, int *index, int N){
    int row;
    double v;
    for(int k = 0; k < N-1; k++){
        if(macierz[index[k]][index[k]] == 0.0){
            row = wyborElemPodst(macierz,index[k], N, index);
            index[row] = index[k];
            index[k] = row;
        }
        for(int i = k+1; i<N;i++){
            v=macierz[index[i]][k];
            for (int j = k+1; j < N; ++j) {
                macierz[index[i]][j] = macierz[index[i]][j] - macierz[index[k]][j] * (v/macierz[index[k]][k]);
                macierz[index[i]][k] = v/macierz[index[k]][k];
            }
        }
    }
}

/**
 * funkcja wykonuje rownanie Ly = b aby obliczyc y
 * @param L macierz dolna uzyskana z dekompozycji
 * @param b wektor
 * @param index wektor indeksow
 * @param N wielkosc macierzy
 */
void macierzL(double **L, double *b, int *index, int N){
    double sum = 0.0;
    for (int i = 0; i <= N-1; i++) {
        for (int j = 0; j < i; j++) {
            sum = sum+L[index[i]][j]*b[index[j]];
        }
        b[index[i]] = (b[index[i]]-sum)/1.0;
        sum=0.0;
    }
}

/**
 * funkcja wykonuje rownanie Ux = b aby obliczyc x
 * @param U macierz gorna uzyskana z dekompozycji
 * @param b wektor
 * @param index wektor indeksow
 * @param N wielkosc macierzy
 */
void macierzU(double **U, double *b, int *index, int N) {
    double sum = 0.0;
    for (int i = 3; i >= 0; i--) {
        for (int j = i + 1; j <= N-1; j++) {
            sum = sum + U[index[i]][j] * b[index[j]];
        }
        b[index[i]] = (b[index[i]] - sum) / U[index[i]][i];
        sum = 0.0;
    }
}

int main() {
    int N = 4;
    double **macierz = new double *[N];
    for (int i = 0; i < N; ++i) {
        macierz[i] = new double[N];
    }
#ifdef v1
    double b[4] = {35.0, 104.0, -366.0, -354.0};
#endif
#ifdef v2
    double E = e;
    double b[4] = {6.0 + E, 6.0 + 2.0 * E, 6.0 + 2.0 * E, 6.0 + E};
    double b1, b2, b3, b4;
    b1 = 1;
    b2 = 2;
    b3 = 2;
    b4 = 1;
    double wynik[4] = {b1, b2, b3, b4};
#endif
    int index[4] = {0, 1, 2, 3};
    wypelnijMacierz(macierz, N);
    printMacierz(macierz, N, index);
    gauss(macierz, index, N);
    cout << "/////////MACIERZ PO GAUSSIE://////////" << endl;
    printMacierz(macierz, N, index);
    macierzL(macierz, b, index, N);
    macierzU(macierz, b, index, N);
    cout << "///////WYNIK://////" << endl;
    printWektor(b, index);
#ifdef v2
    blad(b, wynik);
#endif
    for (int i = 0; i < N; ++i) {
        delete [] macierz[i];
    }
    delete [] macierz;
    return 0;
}