//---------------------------------------------------------------------
//                                                                    |
//          ALGORITMO PARA RESOLVER LA ECUACIÓN DE SCHRÖDINGER        |
//                              Objetivos:                            |
//                                                                    |
//   1. Resolver la ecuación de Schrödinger unidimensional para un    |
//      potencial cuadrado.                                           |
//   2. Comprobar que se conserva la norma.                           | 
//                                                                    |
//---------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <fstream>     // Para trabajar con ficheros
#include <complex>     // Para trabajar con números complejos

#define N 1000         // Tamaño del retículo espacial
#define PI 3.141593

using namespace std;

// FUNCIONES QUE SE VAN A UTILIZAR
complex<double> CalculaPhi(int j, double k0_tilde);
double CalculaNorma2(complex<double> z);

/*---------------------------------------------------------------------
|                           FUNCIÓN PRINCIPAL                         |
---------------------------------------------------------------------*/
int main()
{
    // ------------------ DECLARACIÓN DE VARIABLES --------------------

    int ciclos;                 // Número de oscilaciones completas de la función de onda
    int iteraciones;            // Número de iteraciones que se realizan del algoritmo
    int j,n,k;                  // Contadores

    double lambda;              // Cte de proporcionalidad para la energía del fotón incidente
    double s_tilde;             // s/h^2=1/4k0_tilde^2
    double k0_tilde;            // k0*h, con k0*N*h=2*PI*ciclos
    double V[N+1];              // Potencial, es lambda*k0_tilde^2 si j\in[2N/5,3N/5]
    double Aplus, Aminus;       // Para calcular alpha, Aplus=1, Aminus=1
    double norma_phi;           // Norma de la función de onda
    double norma_txt;
    double densidad;            // Densidad de probabilidad: |phi|^2

    complex<double> phi[N+1];   // Función de onda
    complex<double> xi[N+1];    // phi[j,n+1]=xi[j,n]-phi[j,n]
    complex<double> alpha[N];   // xi[j+1]=alpha[j]*xi[j]+beta[j]
    complex<double> beta[N];    // xi[j+1]=alpha[j]*xi[j]+beta[j]
    complex<double> gamma[N];   // gamma[j]=1/(A0[j]+Aplus*alpha[j])
    complex<double> b[N];       // b[j]=2i*phi[j]/s_tilde
    complex<double> A0[N];      // Para calcular alpha, A0[j]=-2+2i/s_tilde-V[j]

    ofstream fich_densidad;     // Fichero para guardar la densidad de probabilidad
    ofstream fich_norma;        // Fichero para guardar la norma de phi (debe ser 1)
    ofstream fich_potencial;    // Fichero para guardar la forma del potencial


    // ------------------------ INICIALIZACIÓN ------------------------

    fich_densidad.open("densidad.txt");
    fich_norma.open("norma.txt");
    fich_potencial.open("potencial.txt");

    lambda = 5;
    ciclos = 50;                // Restringido a 1,...,N/4
    k0_tilde = 2*PI*ciclos/N;
    s_tilde = 1/(4*k0_tilde*k0_tilde);

    // Inicialización del potencial
    for(j=0; j<N-1; j++)
    {
        if (j>(2*N/5)&&j<(3*N/5)) V[j] = lambda*k0_tilde*k0_tilde;
        else V[j] = 0.0;
        fich_potencial << j << " " << V[j] << endl;
    }
    
    // Función de onda inicial y su norma (con condiciones de contorno)
    phi[0] = (0.0, 0.0);
    phi[N] = (0.0, 0.0);
    norma_phi = 0.0;
    for(j=1; j<N; j++)
    {
        phi[j] = CalculaPhi(j, k0_tilde);
        norma_phi = norma_phi + CalculaNorma2(phi[j]);
    }
    norma_phi = sqrt(norma_phi);

    // Normalizar la función de onda inicial
    for(j=1; j<N; j++) phi[j]=phi[j]/norma_phi;

    norma_txt=0.0;
    // Escribir la densidad de probabilidad en "densidad.txt"
    for(j=0; j<=N; j++)
    {
        densidad = CalculaNorma2(phi[j]);
        fich_densidad << j << " " << densidad << endl;
        norma_txt = norma_txt + densidad;
    }
    fich_densidad << endl;

    // Escribir la norma inicial en "norma.txt" (ahora debe ser 1)
    fich_norma << 0 << " " << sqrt(norma_txt) << endl;

    fich_densidad.flush();
    fich_norma.flush();


    // ----------------------- CÁLCULO DE ALPHA -----------------------

    // Inicializar alpha
    alpha[0] = (0.0, 0.0);
    alpha[N] = (0.0, 0.0);

    // Calcular las A
    Aplus = 1.0;
    Aminus = 1.0;
    for(j=1; j<N; j++) A0[j] = complex<double>(-2.0-V[j], 2.0/s_tilde);

    // Calcular alpha y gamma empezando en N-1
    gamma[N-1] = 1.0/A0[N-1];
    for(j=N-2; j>=0; j--)
    {
        alpha[j] = -Aminus*gamma[j+1];
        gamma[j] = 1.0/(A0[j]+alpha[j]);
    }


    // -------------------------- ALGORITMO ---------------------------

    beta[N]=0.0;
    iteraciones = 500;
    for(n=1;n<=iteraciones;n++){
        
        // CÁLCULO DE B: b[j]=4i*phi[j]/s_tilde
        for(j=1;j<N;j++) b[j] = complex<double>(0,4)*phi[j]/s_tilde;

        // CÁLCULO DE BETA: beta[j-1]=gamma[j]*(b[j]-Aplus*beta[j])
        for(j=N-2;j>=0;j--) beta[j] = gamma[j+1]*(b[j+1]-beta[j+1]);

        // CÁLCULO DE XI: xi[j+1]=alpha[j]*xi[j]+beta[j]
        xi[0] = 0;
        xi[N] = 0;
        for(j=1;j<N;j++) xi[j] = alpha[j-1]*xi[j-1]+beta[j-1];

        // CÁLCULO DE PHI: phi[j,n+1]=xi[j,n]-phi[j,n]
        norma_phi = 0;
        for(k=1;k<N;k++) phi[k] = xi[k]-phi[k];

        // GUARDAR LA DENSIDAD Y NORMA EN CADA ITERACIÓN
        norma_txt = 0;
        for(k=0;k<=N;k++){
            densidad = CalculaNorma2(phi[k]);
            fich_densidad << k << " " << densidad << endl;
            norma_txt = norma_txt+densidad;
        }
        fich_densidad<<endl;
        fich_norma << n << " "<<sqrt(norma_txt)<< endl;

        fich_densidad.flush();
        fich_norma.flush();
    }

    fich_densidad.close();
    fich_norma.close();
    fich_potencial.close();
    return 0;
}


/*---------------------------------------------------------------------
|                            OTRAS FUNCIONES                          |
---------------------------------------------------------------------*/

// Función CalculaPhi: devuelve el valor de la componente de Phi
complex<double> CalculaPhi(int j, double k0_tilde)
{
    double aux;
    complex<double> componente_phi;

    aux = exp(-8.0*(4*j-N)*(4*j-N)/(1.0*N*N));
    componente_phi = complex<double>(cos(k0_tilde*j),sin(k0_tilde*j))*aux;    

    return componente_phi;
}

// Función CalculaNorma2: devuelve el valor de la norma al cuadrado de un
// número complejo double
double CalculaNorma2(complex<double> z)
{
    double norma;
    norma = real(z)*real(z) + imag(z)*imag(z);

    return norma;
}