#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;

// Define a datatype 'dcomp' to declare complex numbers
typedef complex<double> dcomp;

class DFT
{
public:
    int N;
    int TwiddleFactorSign;
    vector<dcomp> inputs;
    vector<int> bitReversedIndex;
    vector<dcomp> bitRevInput;

    vector<dcomp> output;
    // Another array of N complex numbers to store values at the end of each stage
    vector<dcomp> tempArr;
    // Define values for Pi and log2(n)
    double pi;

    int bitReverse(int n, int size);
    int binaryToDecimal(string binary);

    DFT()
    {
        TwiddleFactorSign = -1;
        N = 0;
        // Define values for Pi and log2(n)
        pi = 2 * asin(1); // asin(1) return sin inverse of 1 in radians i.e pi/2
    }

    void getLength()
    {
        cout << "Enter the length of the input sequence in powers of 2: ";
        cin >> N;

        // Ensure N is a power of 2
        while (ceil(log2(N)) != floor(log2(N)))
        {
            cout << N << " is not a power of 2 \n";
            cout << "Enter the length of the input sequence in powers of 2: ";
            cin >> N;
        }
    }

    void getInput()
    {
        dcomp value;
        cout << "Enter " << N << " data points: ";
        for (int i = 0; i < N; i++)
        {
            cin >> value;
            inputs.push_back(value);
        }
    }

    void displayInput()
    {
        for (int i = 0; i < N; i++)
        {
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

    virtual void verifyInputs()
    {
        int changedIndex = 0;
        dcomp newVal;
        cout << "Enter index of value to be changed: ";
        cin >> changedIndex;

        // Ensure that index value is valid
        while (changedIndex > N - 1)
        {
            cout << changedIndex << " is not a valid index.";
            cout << "\nPlease enter a valid index: ";
            cin >> changedIndex;
        }

        cout << "Enter new value for index " << changedIndex << ": ";

        cin >> newVal;
        inputs.insert(inputs.begin() + changedIndex, newVal);

        // display the given input
        cout << " X(k) = {";
        this->displayInput();
    }

    // Obtain the index for bit reverse array
    void getBitReversedIndex()
    {
        for (int i = 0; i < N; i++)
        {
            bitReversedIndex.push_back(bitReverse(i, N));
        }
    }

    // create an array with inputs in bit reversed order
    void setBitReversedInput()
    {
        for (int i = 0; i < N; i++)
        {
            bitRevInput.push_back(inputs.at(bitReversedIndex.at(i)));
        }
    }

    // calculate the DFT values for each input based on the Cooley-Tukey DIF-FFT method
    void calculateDFT()
    {
        int logN = log2(N);
        // initilise the temporary array to hold dummy values for N size, so that we don't get a segmentation fault
        for (int i = 0; i < N; i++)
        {
            tempArr.push_back(0);
        }

        // ---------------------------------
        // Cooley-Tukey algorithm for DIF-FFT
        // ---------------------------------

        // Number of stages is equal to log2(N)
        for (int i = 0; i < logN; i++)
        { // complex number to store twiddle factor
            dcomp wm;
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
                    dcomp intermediate = inputs.at(j - lenBy2PowI) - inputs.at(j);

                    // Store the result of this operation to the temp array
                    tempArr.at(j) = intermediate * wm;
                }
            }

            // transfer tempArr to srr after each stage, so that tempArr can be overwritten as per next stage's result
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
        {
            output.push_back(tempArr.at(bitReversedIndex.at(i)));
        }
    }

    // Display the computed DFT values
    virtual void displayOutput()
    {
        for (int i = 0; i < N; i++)
        {
            if (real(output[i]) == imag(output[i]) && real(output[i]) == 0)
                cout << "X[" << i << "] = 0\n";
            else if (imag(output[i]) == 0)
                cout << "X[" << i << "] = " << real(output[i]) << "\n";
            else if (imag(output[i]) > 0)
                cout << "X[" << i << "] = " << real(output[i]) << " + " << imag(output[i]) << "j\n";
            else
                cout << "X[" << i << "] = " << real(output[i]) << " - " << -1 * imag(output[i]) << "j\n";
        }
    }
};

// converting binary String to decimal
int DFT::binaryToDecimal(string binary)
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

// converting decimal to binary
int DFT::bitReverse(int n, int size)
{
    int num = int(log2(size));
    int binary[num];
    int i = 0;

    // converting decimal to binary in mirrored fashion
    if (n == 0)
        binary[0] = 0;
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
    return binaryToDecimal(str);
}

class IDFT : public DFT
{
public:
    int TwiddleFactorSign;

    IDFT()
    {
        TwiddleFactorSign = 1;
        N = 0;
        // Define values for Pi and log2(n)
        pi = 2 * asin(1); // asin(1) return sin inverse of 1 in radians i.e pi/2
    }

    void verifyInputs()
    {
        int changedIndex = 0;
        dcomp newVal;
        cout << "Enter index of value to be changed: ";
        cin >> changedIndex;

        // Ensure that index value is valid
        while (changedIndex > N - 1)
        {
            cout << changedIndex << " is not a valid index.";
            cout << "\nPlease enter a valid index: ";
            cin >> changedIndex;
        }

        cout << "Enter new value for index " << changedIndex << ": ";

        cin >> newVal;
        inputs.insert(inputs.begin() + changedIndex, newVal);

        // display the given input
        cout << " x(n) = {";
        this->displayInput();
    }

    // Store the computed DFT values in sequenctial order, as DIF-FFT gives it in bit reversed order
    void setOutput()
    {
        // In case of IDFT, after the last stage of outputs, the result is divided by the length of the sequence
        dcomp length = N; //dcomp type because division of complex/int is not valid
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
            if (real(output[i]) == imag(output[i]) && real(output[i]) == 0)
                cout << "x(" << i << ") = 0\n";
            else if (imag(output[i]) == 0)
                cout << "x(" << i << ") = " << real(output[i]) << "\n";
            else if (imag(output[i]) > 0)
                cout << "x(" << i << ") = " << real(output[i]) << " + " << imag(output[i]) << "j\n";
            else
                cout << "x(" << i << ") = " << real(output[i]) << " - " << -1 * imag(output[i]) << "j\n";
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
    while (choice != 1 && choice != 2)
    {
        cout << choice << " is not a valid option, please enter 1 or 2: ";
        cin >> choice;
    }

    if (choice == 1)
    {
        DFT obj;

        obj.getLength();

        // Display note to help input complex number
        cout << "\n**Note: You can input a complex number in the form (real,imaginary)\n";
        cout << "For example, to input 2+3j, you can type (2, 3)\n\n";

        obj.getInput();

        // display the given input as x(n)
        cout << "\nThe given " << obj.N << "-point input sequence is as follows: \n";
        cout << "x(n) = {";
        obj.displayInput();

        // verify if input is correct
        int isVerified = 0;
        cout << "Enter 1 to verify the input sequence, else press 0 to change any value: ";
        cin >> isVerified;

        while (isVerified == 0)
        {
            obj.verifyInputs();
            cout << "Press 1 to verify, press 0 to change some input value: ";
            cin >> isVerified;
        }

        // To print complex numbers without scientific notations
        cout << fixed << setprecision(4);
        cout << "\nThe " << obj.N << "-point DFT sequence, computed using DIF-FFT approach is as follows: \n";
        obj.getBitReversedIndex();
        obj.setBitReversedInput();
        obj.calculateDFT();
        obj.setOutput();
        obj.displayOutput();
    }
    else
    {

        IDFT obj;
        obj.getLength();

        // Display note to help input complex number
        cout << "\n**Note: You can input a complex number in the form (real,imaginary)\n";
        cout << "For example, to input 2+3j, you can type (2, 3)\n\n";

        obj.getInput();

        // display the given input as x(n)
        cout << "\nThe given " << obj.N << "-point input sequence is as follows: \n";
        cout << "X(k) = {";
        obj.displayInput();

        // verify if input is correct
        int isVerified = 0;
        cout << "Enter 1 to verify the input sequence, else press 0 to change any value: ";
        cin >> isVerified;

        while (isVerified == 0)
        {
            obj.verifyInputs();
            cout << "Press 1 to verify, press 0 to change some input value: ";
            cin >> isVerified;
        }

        // To print complex numbers without scientific notations
        cout << fixed << setprecision(4);
        cout << "\nThe " << obj.N << "-point IDFT sequence, computed using DIF-IFFT approach is as follows: \n";

        obj.getBitReversedIndex();
        obj.setBitReversedInput();
        obj.calculateDFT();
        obj.setOutput();
        obj.displayOutput();
    }
}