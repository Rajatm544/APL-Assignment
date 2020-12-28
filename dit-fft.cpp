#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
using namespace std;

// converting binary String to decimal
int binToDec(string binary)
{
    string num = binary;
    int decVal = 0;

    int base = 1;
    int len = num.length();
    for (int i = len - 1; i >= 0; i--)
    {
        if (num[i] == '1')
            decVal += base;

        base = base * 2;
    }
    return decVal;
}

// converting decimal to binary in mirrored fashion, to emulate binary format of bit reversed order
int bitReverse(int n, int size)
{
    int num = int(log2(size));
    int binary[num];
    int i = 0;

    // converting decimal to binary in mirrored fashion
    // if n is zero, then add log2(size) number of zeroes as the binary form of n
    if (n == 0)
    {
        for (int i = 0; i < num; i++)
            binary[i] = 0;
    }
    while (n > 0)
    {
        binary[i] = n % 2;
        n = n / 2;
        i++;
    }

    // converting int[] to String
    string str;
    for (int j = 0; j < num; j++)
        str.push_back('0' + binary[j]);

    // convert binary String to a Decimal value
    return binToDec(str);
}

// Define a datatype 'dcomp' to declare complex numbers
typedef complex<double> dcomp;

int main()
{

    // input a choice to see what the user wants to compute
    // in case of IDFT, the only 2 differences are that the negetive power of the twiddle factor is multiplied to the subtracted number in the Cooley-Tukey algorithm
    // And the other change is that the sequence is divided by N at the end
    int choice;
    cout << "Enter 1 to compute DFT, Enter 2 to compute IDFT: ";
    cin >> choice;

    // Ensure choice is either 1 or 2
    while (choice != 1 && choice != 2)
    {
        cout << choice << " is not a valid option, please enter 1 or 2: ";
        cin >> choice;
    }

    int dftOrIdft;
    if (choice == 1)
        dftOrIdft = -1; // -1 multiplied to the power of twiddle factor is calculation of DFT
    else
        dftOrIdft = 1;

    int n; //length of the input sequence

    if (choice == 1)
        cout << "Enter the length of the DFT sequence in powers of 2: ";
    else
        cout << "Enter the length of the IDFT sequence in powers of 2: ";

    cin >> n;

    // Ensure n is a power of 2
    while (ceil(log2(n)) != floor(log2(n)))
    {
        cout << n << " is not a power of 2 \n";
        cout << "Enter the length of the input sequence in powers of 2: ";
        cin >> n;
    }

    // Display note to help input complex number
    cout << "\n**Note: You can input a complex number in the form (real,imaginary)\n";
    cout << "For example, to input 2+3j, you can type (2, 3)\n\n";

    // create an array of inputs
    dcomp arr[n];

    cout << "Enter " << n << " data points: ";
    for (int i = 0; i < n; i++)
    {
        cin >> arr[i];
    }

    // display the given input as x(n)
    cout << "\nThe given " << n << "-point input sequence is as follows: \n";

    if (choice == 1)
        cout << "x(n) = {";
    else
        cout << "X(k) = {";

    for (int i = 0; i < n; i++)
    {
        if (i == n - 1)
        {
            if (imag(arr[i]) < 0)
                cout << real(arr[i]) << " - " << -1 * imag(arr[i]) << "j}\n\n";
            else if (imag(arr[i]) > 0)
                cout << real(arr[i]) << " + " << imag(arr[i]) << "j}\n\n";
            else
                cout << real(arr[i]) << "}\n\n";
        }
        else
        {
            if (imag(arr[i]) < 0)
                cout << real(arr[i]) << " - " << -1 * imag(arr[i]) << "j, ";
            else if (imag(arr[i]) > 0)
                cout << real(arr[i]) << " + " << imag(arr[i]) << "j, ";
            else
                cout << real(arr[i]) << ", ";
        }
    }

    // verify if input is correct
    int isVerified = 0;
    cout << "Enter 1 to verify the input sequence, else press 0 to change any value: ";
    cin >> isVerified;

    while (isVerified == 0)
    {
        int changedIndex = 0;
        dcomp newVal;
        cout << "Enter index of value to be changed: ";
        cin >> changedIndex;

        // Ensure that index value is valid
        while (changedIndex > n - 1)
        {
            cout << changedIndex << " is not a valid index.";
            cout << "\nPlease enter a valid index: ";
            cin >> changedIndex;
        }

        cout << "Enter new value for index " << changedIndex << ": ";

        cin >> newVal;
        arr[changedIndex] = newVal;

        // display the given input
        cout << "\nThe new " << n << " point input sequence is as follows: \n";
        if (choice == 1)
            cout << "x(n) = {";
        else if (choice == 2)
            cout << " X(k) = {";

        for (int i = 0; i < n; i++)
        {
            if (i == n - 1)
            {
                if (imag(arr[i]) < 0)
                    cout << real(arr[i]) << " - " << -1 * imag(arr[i]) << "j}\n\n";
                else if (imag(arr[i]) > 0)
                    cout << real(arr[i]) << " + " << imag(arr[i]) << "j}\n\n";
                else
                    cout << real(arr[i]) << "}\n\n";
            }
            else
            {
                if (imag(arr[i]) < 0)
                    cout << real(arr[i]) << " - " << -1 * imag(arr[i]) << "j, ";
                else if (imag(arr[i]) > 0)
                    cout << real(arr[i]) << " + " << imag(arr[i]) << "j, ";
                else
                    cout << real(arr[i]) << ", ";
            }
        }

        cout << "Press 1 to verify, press 0 to change some input value: ";
        cin >> isVerified;
    }

    if (choice == 1)
        cout << "\nThe " << n << "-point DFT sequence, computed using DIF-FFT approach is as follows: \n";
    else
        cout << "\nThe " << n << "-point IDFT sequence, computed using DIF-IFFT approach is as follows: \n";

    // Obtain the index for bit reverse array
    int bitReversedIndex[n];
    for (int i = 0; i < n; i++)
    {
        bitReversedIndex[i] = bitReverse(i, n);
    }

    // create an array with inputs in bit reversed order
    dcomp bitRevInput[n];
    for (int i = 0; i < n; i++)
    {
        bitRevInput[i] = arr[bitReversedIndex[i]];
    }

    // Define values for Pi and log2(n)
    double pi = 2 * asin(1); // asin(1) return sin inverse of 1 in radians i.e pi/2
    int N = n;
    int logN = log2(N);

    // complex number to store twiddle factor
    dcomp wm;

    // To print complex numbers without scientific notations
    cout << fixed << setprecision(4);

    // Another array of N complex numbers to store values at the end of each stage
    dcomp tempArr[N];

    // ---------------------------------
    // Cooley-Tukey algorithm for DIF-FFT and DIF-IFFT
    // ---------------------------------

    // Number of stages is equal to log2(N)
    for (int i = 0; i < logN; i++)
    {
        // In each stage, there are N operations of addition of complex numbers
        for (int j = 0; j < N; j++)
        {
            int lenBy2PowI = (N / pow(2, i + 1));
            // In order to check if complex numbers need to be added or subtracted

            // If (j mod (N/2^i)) < value of N / (2 ^ (i + 1)), then just add the 2 correct complex numbers
            if ((j % int(N / pow(2, i))) < lenBy2PowI)
            {
                tempArr[j] = arr[j] + arr[j + lenBy2PowI];
            }
            else
            {
                // calculate the value of twiddle factor as wm = e^(-j*2 pi * k / N) but in polar form
                // k value will be same as (j / lenBy2PowI), N value will be (length of input / 2^i)

                // if power is multiplied by -1 => DFT
                // if power is multiplied by 1 => IDFT
                wm = polar(1.0, (dftOrIdft * 2) * M_PI * ((j % lenBy2PowI) / (N / pow(2, i))));

                // keep the intermediate value of subtraction of required numbers ready for multiplication with twiddle factor
                dcomp intermediate = arr[j - lenBy2PowI] - arr[j];

                // Store the result of this operation to the temp array
                tempArr[j] = intermediate * wm;
            }
        }

        // transfer tempArr to srr after each stage, so that tempArr can be overwritten as per next stage's result
        for (int k = 0; k < N; k++)
        {
            arr[k] = tempArr[k];
        }
    }

    // In case of IDFT, after the last stage of outputs, the result is divided by the length of the sequence
    dcomp length; //dcomp type because division of complex/int is not valid

    if (choice == 1)
        length = 1;
    else
        length = N;

    // Since output of DIF-FFT is in bit reversed order, store the result in correct order to arr itself
    for (int i = 0; i < N; i++)
    {
        arr[i] = tempArr[bitReversedIndex[i]] / length;
    }

    if (choice == 1)
    {
        // Output the DFT sequence
        for (int i = 0; i < N; i++)
        {
            if (real(arr[i]) == imag(arr[i]) && real(arr[i]) == 0)
                cout << "X[" << i << "] = 0\n";
            else if (imag(arr[i]) == 0)
                cout << "X[" << i << "] = " << real(arr[i]) << "\n";
            else if (imag(arr[i]) > 0)
                cout << "X[" << i << "] = " << real(arr[i]) << " + " << imag(arr[i]) << "j\n";
            else
                cout << "X[" << i << "] = " << real(arr[i]) << " - " << -1 * imag(arr[i]) << "j\n";
        }
    }
    else
    {
        // Output the DFT sequence
        for (int i = 0; i < N; i++)
        {
            if (real(arr[i]) == imag(arr[i]) && real(arr[i]) == 0)
                cout << "x(" << i << ") = 0\n";
            else if (imag(arr[i]) == 0)
                cout << "x(" << i << ") = " << real(arr[i]) << "\n";
            else if (imag(arr[i]) > 0)
                cout << "x(" << i << ") = " << real(arr[i]) << " + " << imag(arr[i]) << "j\n";
            else
                cout << "x(" << i << ") = " << real(arr[i]) << " - " << -1 * imag(arr[i]) << "j\n";
        }
    }
}
