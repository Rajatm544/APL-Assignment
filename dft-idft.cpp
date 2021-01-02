#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;

// shorthand datatype to store floating point complex numbers
typedef complex<float> fcomp;

class DFT
{
public:
    int N;
    double pi;
    int TwiddleFactorSign;

    vector<fcomp> inputs;         // an array of complex inputs
    vector<int> bitReversedIndex; // an array of indexes in the bit reveresed order
    vector<fcomp> output;         // an array to store the final output in the sequential order
    vector<fcomp> tempArr;        // an array of N complex numbers to store values at the end of each stage

    int bitReverse(int n, int size); // member function to obtain index in bit reversed order

    // constructor to assign default values of member data
    DFT()
    {
        TwiddleFactorSign = -1; // This sign is that of the power of the twiddle factor, which decides if computation is for DFT/IDFT
        N = 0;                  // initialise the length of the DFT sequence to default value
        pi = 2 * asin(1);       // asin(1) return sin inverse of 1 in radians i.e pi/2
    }

    // Get the value for the length of the input sequence
    void getLength()
    {
        cout << "Enter the length of the input sequence in powers of 2: ";
        cin >> N;

        // If a string is input instead of int, it defaults to 0, hence exit() the program to avoid errors
        if (N == 0)
        {
            cout << "\nInvalid input, please re-run the program\n";
            exit(0);
        }

        // Ensure N is a power of 2
        while (ceil(log2(N)) != floor(log2(N)))
        {
            cout << N << " is not a power of 2 \n";
            cout << "Enter the length of the input sequence in powers of 2: ";
            cin >> N;
        }

        // initilise the temArr vector to hold dummy values for N size, so that we don't get a segmentation fault
        for (int i = 0; i < N; i++)
            tempArr.push_back(0);
    }

    // Get the input sequence of complex numbers, and store as class member data
    void getInput()
    {
        // Display note to help input complex number
        cout << "\n**Note: You can input a complex number in the form (real,imaginary)\n";
        cout << "For example, to input 2+3j, you can type (2, 3)\n";
        cout << "Example Input: 1 0 (2, 3) 1\n\n";

        fcomp value;
        cout << "Enter " << N << " data points: ";
        for (int i = 0; i < N; i++)
        {
            cin >> value;
            inputs.push_back(value);
        }
    }

    // Display the input sequence after the user inputs it, to verify the input values
    void displayInput()
    {
        for (int i = 0; i < N; i++)
        {
            // For the last value, display a closing bracket and add a new line
            if (i == N - 1)
            {
                if (imag(inputs[i]) < 0)
                    cout << real(inputs[i]) << " - " << -1 * imag(inputs[i]) << "j}\n\n";
                else if (imag(inputs[i]) > 0)
                    cout << real(inputs[i]) << " + " << imag(inputs[i]) << "j}\n\n";
                else
                    cout << real(inputs[i]) << "}\n\n";
            }
            else
            {
                if (imag(inputs[i]) < 0)
                    cout << real(inputs[i]) << " - " << -1 * imag(inputs[i]) << "j, ";
                else if (imag(inputs[i]) > 0)
                    cout << real(inputs[i]) << " + " << imag(inputs[i]) << "j, ";
                else
                    cout << real(inputs[i]) << ", ";
            }
        }
    }

    // If the user wants to change any of the input value
    virtual void verifyInputs()
    {
        int changedIndex = 0;
        fcomp newVal;
        cout << "Enter index of value to be changed: ";
        cin >> changedIndex;

        // If a string is input insted of the int, then exit()
        if (changedIndex == 0)
        {
            cout << "\nInvalid input, please re-run the program\n";
            exit(0);
        }

        // Ensure that index value is valid
        while (changedIndex > N - 1)
        {
            cout << changedIndex << " is not a valid index.";
            cout << "\nPlease enter a valid index: ";
            cin >> changedIndex;
        }

        cout << "Enter new value for index " << changedIndex << ": ";

        cin >> newVal;
        inputs.at(changedIndex) = newVal;

        // display the newly given input sequence
        cout << "x(n) = {";
        this->displayInput();
    }

    // Create an array of integers to store the but reversed indexes
    void getBitReversedIndex()
    {
        for (int i = 0; i < N; i++)
        {
            bitReversedIndex.push_back(bitReverse(i, N));
        }
    }

    // calculate the DFT values for each input based on the Cooley-Tukey DIF-FFT method
    void calculateDFTOrIDFT()
    {
        // complex number to store twiddle factor
        fcomp wm;
        int logN = log2(N);

        // ---------------------------------
        // Cooley-Tukey algorithm for DIF-FFT
        // ---------------------------------

        // Number of stages is equal to log2(N)
        for (int i = 0; i < logN; i++)
        {
            // In each stage, there are N operations of addition of complex numbers
            for (int j = 0; j < N; j++)
            {
                int lenBy2PowI = int(N / pow(2, i + 1));

                // In order to check if complex numbers need to be added or subtracted

                // If (j mod (N/2^i)) < value of N / (2 ^ (i + 1)), then just add the 2 correct complex numbers
                if ((j % int(N / pow(2, i))) < lenBy2PowI)
                {
                    tempArr.at(j) = inputs.at(j) + inputs.at(j + lenBy2PowI);
                }
                else
                {
                    // calculate the value of twiddle factor as wm = e^(-j*2 pi * k / N) but in polar form
                    // k value will be same as (j / lenBy2PowI), N value will be (length of input / 2^i)

                    wm = polar(1.0, (TwiddleFactorSign * 2) * M_PI * ((j % lenBy2PowI) / (N / pow(2, i))));

                    // keep the intermediate value of subtraction of required numbers ready for multiplication with twiddle factor
                    fcomp intermediate = inputs.at(j - lenBy2PowI) - inputs.at(j);

                    // Store the result of this operation to the temp array
                    tempArr.at(j) = intermediate * wm;
                }
            }
            // transfer tempArr to input vector after each stage, so that tempArr can be overwritten as per next stage's result
            for (int k = 0; k < N; k++)
            {
                inputs.at(k) = tempArr.at(k);
            }
        }
    }

    // Store the computed DFT values in sequenctial order, as DIF-FFT gives it in bit reversed order
    virtual void setOutput()
    {
        for (int i = 0; i < N; i++)
            output.push_back(tempArr.at(bitReversedIndex.at(i)));
    }

    // Display the computed DFT values
    virtual void displayOutput()
    {
        for (int i = 0; i < N; i++)
        {
            if (real(output.at(i)) == imag(output.at(i)) && real(output.at(i)) == 0)
                cout << "X[" << i << "] = 0\n";
            else if (real(output.at(i)) == 0)
                cout << "X[" << i << "] = " << imag(output.at(i)) << "j\n";
            else if (imag(output.at(i)) == 0)
                cout << "X[" << i << "] = " << real(output.at(i)) << "\n";
            else if (imag(output.at(i)) > 0)
                cout << "X[" << i << "] = " << real(output.at(i)) << " + " << imag(output.at(i)) << "j\n";
            else
                cout << "X[" << i << "] = " << real(output.at(i)) << " - " << -1 * imag(output.at(i)) << "j\n";
        }
    }
};

// converting binary String to corresponding decimal number
int binaryToDecimal(string binary)
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
int DFT::bitReverse(int n, int size)
{
    int num = int(log2(size));
    int binary[num];
    int i = 0;

    // converting decimal to binary in mirrored fashion
    if (n == 0)
    {
        // if n is zero, then add log2(size) number of zeroes as the binary form of n
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

    // convert binary String to a corresponding decimal value
    return binaryToDecimal(str);
}

// Single inheritance of the DFT class
class IDFT : public DFT
{
public:
    int TwiddleFactorSign;

    // constructor to override the inherited constructor
    IDFT()
    {
        TwiddleFactorSign = 1; // This is the sign of the power of the twiddle factor, +1 indicates computation of IDFT
        N = 0;
        pi = 2 * asin(1); // asin(1) return sin inverse of 1 in radians i.e pi/2
    }

    // Redefine the verifyInputs() member function, to display X(k) instead of x(n) when needed
    void verifyInputs()
    {
        int changedIndex = 0;
        fcomp newVal;
        cout << "Enter index of value to be changed: ";
        cin >> changedIndex;

        // If a string is input insted of the int, then exit()
        if (changedIndex == 0)
        {
            cout << "\nInvalid input, please re-run the program\n";
            exit(0);
        }

        // Ensure that index value is valid
        while (changedIndex > N - 1)
        {
            cout << changedIndex << " is not a valid index.";
            cout << "\nPlease enter a valid index: ";
            cin >> changedIndex;
        }

        cout << "Enter new value for index " << changedIndex << ": ";

        cin >> newVal;
        inputs.at(changedIndex) = newVal;

        // display the newly given input sequence
        cout << "X(k) = {";
        this->displayInput();
    }

    // Store the computed DFT values in sequenctial order, as DIF-FFT gives it in bit reversed order
    void setOutput()
    {
        // In case of IDFT, after the last stage of outputs, the result is divided by the length of the sequence
        fcomp length = N; //fcomp type because division of complex/int is not valid
        for (int i = 0; i < N; i++)
        {
            output.push_back(tempArr.at(bitReversedIndex.at(i)) / length);
        }
    }

    // Display the computed DFT values
    void displayOutput()
    {
        for (int i = 0; i < N; i++)
        {
            if (real(output.at(i)) == imag(output.at(i)) && real(output.at(i)) == 0)
                cout << "x(" << i << ") = 0\n";
            else if (real(output.at(i)) == 0)
                cout << "x[" << i << "] = " << imag(output.at(i)) << "j\n";
            else if (imag(output.at(i)) == 0)
                cout << "x(" << i << ") = " << real(output.at(i)) << "\n";
            else if (imag(output.at(i)) > 0)
                cout << "x(" << i << ") = " << real(output.at(i)) << " + " << imag(output.at(i)) << "j\n";
            else
                cout << "x(" << i << ") = " << real(output.at(i)) << " - " << -1 * imag(output.at(i)) << "j\n";
        }
    }
};

int main()
{
    // input a choice to see what the user wants to compute
    // in case of IDFT, the only 2 differences are that the negetive power of the twiddle factor is multiplied to the subtracted number in the Cooley-Tukey algorithm
    // And the other change is that the sequence is divided by N at the end

    int choice;
    cout << "Enter 1 to compute DFT, Enter 2 to compute IDFT: ";
    cin >> choice;

    // Ensure choice is either 1 or 2
    while ((choice != 1 && choice != 2))
    {
        if (choice == 0)
        {
            cout << "\nInvalid input, please re-run the program\n";
            exit(0);
        }

        cout << choice << " is not a valid option, please enter 1 or 2: ";
        cin >> choice;
    }

    if (choice == 1)
    {
        DFT obj;

        obj.getLength();
        obj.getInput();

        // display the given input as x(n)
        cout << "\nThe given " << obj.N << "-point input sequence is as follows: \n";
        cout << "x(n) = {";
        obj.displayInput();

        // verify if input is correct
        int isVerified;
        cout << "Enter 1 to verify the input sequence, enter 2 to change any value: ";
        cin >> isVerified;

        // Ensure that the input is either 1 or 2
        while (isVerified != 1 && isVerified != 2)
        {
            // if user inputs a string, then the value defaults to 0 as int
            if (isVerified == 0)
            {
                cout << "\nInvalid input, please re-run the program\n";
                exit(0);
            }
            cout << isVerified << " is not a valid input. Please input 1 or 2: ";
            cin >> isVerified;
        }

        // Keep asking if the user wants to change any input value, unless 1 is input to verify the input sequence
        while (isVerified == 2)
        {
            obj.verifyInputs();
            cout << "Press 1 to verify, press 2 to change some input value: ";
            cin >> isVerified;
        }

        // To print complex numbers without scientific notations
        cout << fixed << setprecision(4);
        cout << "\nThe " << obj.N << "-point DFT sequence, computed using DIF-FFT approach is as follows: \n";

        obj.getBitReversedIndex();
        obj.calculateDFTOrIDFT();
        obj.setOutput();
        obj.displayOutput();
    }
    else
    {
        IDFT obj;

        obj.getLength();
        obj.getInput();

        // display the given input as x(n)
        cout << "\nThe given " << obj.N << "-point input sequence is as follows: \n";
        cout << "X(k) = {";
        obj.displayInput();

        // verify if input is correct
        int isVerified = 0;
        cout << "Enter 1 to verify the input sequence, enter 2 to change any of the input values: ";
        cin >> isVerified;

        // Ensure that the user enters either 1 or 2
        while (isVerified != 1 && isVerified != 2)
        {
            // Exit the program if user enter a string or a 0
            if (isVerified == 0)
            {
                cout << "\nInvalid input, please re-run the program\n";
                exit(0);
            }

            cout << isVerified << " is not a valid input. Please input 1 or 2: ";
            cin >> isVerified;
        }

        while (isVerified == 2)
        {
            obj.verifyInputs();
            cout << "Press 1 to verify, press 2 to change some input value: ";
            cin >> isVerified;
        }

        // To print complex numbers without scientific notations
        cout << fixed << setprecision(4);
        cout << "\nThe " << obj.N << "-point IDFT sequence, computed using DIF-IFFT approach is as follows: \n";

        obj.getBitReversedIndex();
        obj.calculateDFTOrIDFT();
        obj.setOutput();
        obj.displayOutput();
    }
}