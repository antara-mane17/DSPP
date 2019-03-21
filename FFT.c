
/*
 * Program to perform Fast Fourier Transform (DFT) using Radix 2 DITFFT Algorithm
 *
 * Author: Srijal Poojari
 * Contact: srijal97@gmail.com
 * Description: Provides functions to calculate FFT of any signal of length N (radix 2)
 *
 */

#include<stdio.h>
#include<math.h>


//----------------------------------------------------------------------------------//

struct complex_number {   // structure to store and represent complex numbers
    double real;
    double imag;
};

struct complex_number twiddle_factor(int signal_length, int degree){  // calculates and returns the twiddle factor
    double theta = -2 * M_PI * degree / signal_length;

    // We know, twiddle_factor = exp(-j * 2*pi * n*k / N)
    // exp(theta) = cos(theta) + i*sin(theta)

    struct complex_number Wn = {cos(theta), sin(theta)};

    return Wn;
};

struct complex_number multiply_complex(struct complex_number a, struct complex_number b){  // multiplies 2 complex numbers
    struct complex_number result = {(a.real * b.real - a.imag * b.imag),
                                    (a.real * b.imag + a.imag * b.real)};

    return result;
};

struct complex_number add_complex(struct complex_number a, struct complex_number b){  // adds 2 complex numbers
    struct complex_number result = {(a.real + b.real),
                                    (a.imag + b.imag)};

    return result;
};

struct complex_number sub_complex(struct complex_number a, struct complex_number b){  // subtracts 2 complex numbers
    struct complex_number result = {(a.real - b.real),
                                    (a.imag - b.imag)};

    return result;
};

//----------------------------------------------------------------------------------//

int radix2_greater_than_equal_to(int n){
    // Calculates and returns a radix 2 number greater than or equal to n.

    int result = 1;

    while( result < n) {
        result = result << 1;  // Shift bits (increase power of 2) till desired number is achieved.
    }

    return result;
}

void fft(int N, struct complex_number a[], struct complex_number A[]) {
    // Calculate Fast Fourier Transform of signal a[] of length N. Result is A[].

    int k, n;

    if (N == 2){  // FFT flowchart (butterfly diagram) for 2 point signal.
        A[0] = add_complex(a[0], a[1]);
        A[1] = sub_complex(a[0], a[1]);
        return;
    }
    else {  // If N > 2, divide into further parts and first find its FFT. Uses recursion.

        // X[k] =  G[k]   +  W^(k) * H[k] //
        // N pt   N/2 pt           N/2 pt //
        // ------------------------------ //
        // where, G[k] = DFT[x(2*n)]      //
        //        H[k] = DFT[x(2*n + 1)]  //

        struct complex_number h[N/2];
        struct complex_number H[N/2];

        struct complex_number g[N/2];
        struct complex_number G[N/2];

        for (n = 0; n < N/2; n++){
            g[n] = a[2*n];
            h[n] = a[2*n + 1];
        }

        fft(N/2, g, G);
        fft(N/2, h, H);

        for (k = 0; k < N; k++){
            A[k] = add_complex(G[k % (N/2)], multiply_complex(twiddle_factor(N, k), H[k % (N/2)]));
        }

        return;
    }
}

void ifft(int N, struct complex_number A[], struct complex_number a[]) {

    // x[n] = (1/N) * FFT(X*[k])*
    int i;

    struct complex_number A_conj[N];  // Conjugate of A[k]
    struct complex_number a_conj[N];  // Conjugate of a[n]

    for(i = 0; i < N; i++) {  // Calculate conjugate
        A_conj[i].real = A[i].real;
        A_conj[i].imag = -1 * A[i].imag;
    }

    fft(N, A_conj, a_conj);  // find FFT of conjugate sequence

    for(i = 0; i < N; i++) {  // Again find conjugate and divide term by N.
        a[i].real = a_conj[i].real / N;
        a[i].imag = -1 * a_conj[i].imag / N;
    }
}

//----------------------------------------------------------------------------------//

int main()
{

    int selector, N;

    printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
    scanf("%d", &selector);

    printf( "Enter the length of x[n] i.e. N = ");
    scanf("%d", &N);

    // Convert N to radix 2 number for radix 2 FFT Algorithm.
    int N_rad2 = radix2_greater_than_equal_to(N);

    struct complex_number x[N_rad2];
    struct complex_number X[N_rad2];

    int i;

    for (i = 0; i < N_rad2; i++){
        // Initialize with zeros as we accept only N (not N_rad2) inputs from the user.
        x[i].real = 0;
        x[i].imag = 0;
    }

    printf( "Enter the values of x[n] : \n");

    if (selector == 1){ // If yes (real valued), accept only real values for ease of use.
        for(i = 0; i < N; i++) {
            scanf("%lf", &x[i].real);
            x[i].imag = 0;
        }
    }
    else{
        for(i = 0; i < N; i++) {
            printf("Real: ");
            scanf("%lf", &x[i].real);

            printf("Imaginary: ");
            scanf("%lf", &x[i].imag);
        }
    }

    //----------------FFT-------------------//

    fft(N_rad2, x, X);

    printf("FFT result, X[k]: \n");
    for(i = 0; i < N_rad2; i++) {
        printf("X[%d] = %lf + j%lf \n", i, X[i].real, X[i].imag);
    }

    //---------------IFFT-------------------//

    ifft(N_rad2, X, x);

    printf("IFFT result, x[n]: \n");
    for(i = 0; i < N_rad2; i++) {
        printf("x[%d] = %lf + j%lf \n", i, x[i].real, x[i].imag);
    }

    return 0;

}

