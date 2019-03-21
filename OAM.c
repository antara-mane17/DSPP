
/*
 * Program to perform Linear Filtering of long data sequence using the Overlap Add Method (OAM)
 *
 * Author: Srijal Poojari
 * Contact: srijal97@gmail.com
 * Description: Provides functions to calculate the linear convolution of 
 *              a long data sequence, of length N, using OAM.
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

void fast_circular_convolve(int N, struct complex_number x[], struct complex_number h[], struct complex_number y[]) {
    // Find circular convolution of 2 signals of length N. Result will be of length N.

    int N_rad2 = radix2_greater_than_equal_to(N);  // N for radix 2 fft algorithm

    struct complex_number x_new[N_rad2];
    struct complex_number h_new[N_rad2];
    struct complex_number y_new[N_rad2];

    int i;

    for(i = 0; i < N_rad2; i++) {  // Zero padding
        if(i < N) {
            x_new[i].real = x[i].real;
            x_new[i].imag = x[i].imag;

            h_new[i].real = h[i].real;
            h_new[i].imag = h[i].imag;
        }
        else {
            x_new[i].real = 0;
            x_new[i].imag = 0;

            h_new[i].real = 0;
            h_new[i].imag = 0;
        }
    }
    
    // We know, 
    // circular_convolution(x(n), h(n)) ===> IFFT( X[k] * H[K] )
    // where, X[k] = FFT[x(n)], H[k] = FFT[h(n)]

    struct complex_number X[N_rad2];
    struct complex_number H[N_rad2];
    struct complex_number Y[N_rad2];

    fft(N_rad2, x_new, X);
    fft(N_rad2, h_new, H);

    for(i = 0; i < N_rad2; i++) {
        Y[i] = multiply_complex(X[i], H[i]);
    }

    ifft(N_rad2, Y, y_new);

    for(i = 0; i < N; i++) {  // copy result of length N to y[n]
        y[i].real = y_new[i].real;
        y[i].imag = y_new[i].imag;
    }

    return;
}

void fast_linear_convolve(int len_x, int len_h, struct complex_number x[], struct complex_number h[], struct complex_number y[]) {
    // Find linear convolution of 2 signals of lengths len_x and len_h. Result will be of length (len_x + len_h - 1).

    int N = len_x + len_h - 1;  // length of result

    struct complex_number x_new[N];
    struct complex_number h_new[N];

    int i;
    
    // linear_convolution == circular_convolution, when N is same.
    
    // Therefore, to perform linear_convolution, we find circular_convolution
    // of the 2 signals by zero padding both signals to length N.

    for(i = 0; i < N; i++) {  // Zero padding
        if(i < len_x) {
            x_new[i].real = x[i].real;
            x_new[i].imag = x[i].imag;
        }
        else {
            x_new[i].real = 0;
            x_new[i].imag = 0;
        }

        if(i < len_h) {
            h_new[i].real = h[i].real;
            h_new[i].imag = h[i].imag;
        }
        else {
            h_new[i].real = 0;
            h_new[i].imag = 0;
        }
    }

    fast_circular_convolve(N, x_new, h_new, y);

    return;
}

//----------------------------------------------------------------------------------//

int main()
{

    int i, j, selector;
    int signal_length, M, N, L;
    
    // signal_length = length(x(n))
    //             M = length(h(n))
    //             L = length(decomposed x(n))
    //             N = L + M - 1; or length(decomposed y(n))

    printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
    scanf("%d", &selector);

    printf( "Length of x[n] = ");
    scanf("%d", &signal_length);

    struct complex_number x[signal_length];

    printf( "Enter the values of x[n] : \n");

    if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
        for(i = 0; i < signal_length; i++) {
            scanf("%lf", &x[i].real);
            x[i].imag = 0;
        }
    }
    else{
        for(i = 0; i < signal_length; i++) {
            printf("Real: ");
            scanf("%lf", &x[i].real);

            printf("Imaginary: ");
            scanf("%lf", &x[i].imag);
        }
    }

    printf( "Is h[n] a real valued signal? (1: Yes, 0: No): ");
    scanf("%d", &selector);

    printf( "Length of h[n] = ");
    scanf("%d", &M);

    struct complex_number h[M];

    printf( "Enter the values of h[n] : \n");

    if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
        for(i = 0; i < M; i++) {
            scanf("%lf", &h[i].real);
            h[i].imag = 0;
        }
    }
    else{
        for(i = 0; i < M; i++) {
            printf("Real: ");
            scanf("%lf", &h[i].real);

            printf("Imaginary: ");
            scanf("%lf", &h[i].imag);
        }
    }
    
    // Now to find value of N (radix 2) such that L < M in the relation N = L + M - 1

    N = 1; // Initial Assumption
    do {
        N = N << 1;    // get radix-2 value of N
        L = N - M + 1; // Length of decomposed x[n]
    } while( L < M );
    
    int num_of_decompositions = (signal_length / L) + 1;  // Number of decompositions required

    printf("Number of decompositions: %d \n", num_of_decompositions);

    int result_length = num_of_decompositions * L + M - 1;

    struct complex_number y[result_length];

    for(i = 0; i < result_length; i++){  // Intially, result is all 0
        y[i].real = 0;
        y[i].imag = 0;
    }
    
    //--------------- Perform OAM ------------------//

    for(i = 0; i < num_of_decompositions; i++){
        struct complex_number x_decom[L];  // decomposed x[n]
        struct complex_number y_decom[L+M-1];  // decomposed y[n]

        for(j = 0; j < L; j++) {  // Copy 'L' signal values from x[n] to x_decom[n]
            if ((L*i + j) < signal_length) {
                x_decom[j] = x[L*i + j]; 
            }
            else {
                x_decom[j].real = 0;
                x_decom[j].imag = 0;
            }
        }
        
        // Calculate y_decom[], length(y_decom[]) = L + M - 1

        fast_linear_convolve(L, M, x_decom, h, y_decom);

        for(j = 0; j < (L + M - 1); j++) {  // Overlap and add y_decom to y.  
            y[i*L + j] = add_complex(y[i*L + j], y_decom[j]);
        }
    }

    printf("OAM result, y[n]: \n");
    for(i = 0; i < result_length; i++) {
        printf("y[%d] = %lf \n", i, y[i].real);

    }

    return 0;

}

