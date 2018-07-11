#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <cctype>
#include <complex>
#include <iomanip>
#include <vector>
#include <iostream>
#include <chrono>

using namespace std;

template<typename DataType>
vector<complex<DataType>> getCoefficients(size_t N);

template<typename DataType>
vector<complex<DataType>> myDFT(vector<complex<DataType>> inputSignal);

template<typename DataType>
vector<complex<DataType>> myFFT(vector<complex<DataType>> inputSignal);


template<typename DataType>
bool testSignal(vector<complex<DataType>>& inputSignal, vector<complex<DataType>>& fftSignal);

template<typename DataType>
vector<complex<DataType>> myFFTOptimized(vector<complex<DataType>>& inputSignal);

template<typename DataType>
void printSignal(vector<complex<DataType>> inputSignal);

int main()
{
	std::cout << "hello world " << std::endl;

	vector<complex<float>> inputSignalFloat;
	vector<complex<double>> inputSignalDouble;

	for (size_t i = 0; i < 1024; ++i)
	{
		float realPartFloat = static_cast<float>(std::rand() % 100);
		float imPartFloat = static_cast<float>(std::rand() % 100);

		double realPartDouble = static_cast<double>(std::rand() % 100);
		double imPartDouble = static_cast<double>(std::rand() % 100);

		inputSignalFloat.push_back(complex<float>(realPartFloat, imPartFloat));
		inputSignalDouble.push_back(complex<double>(realPartDouble, imPartDouble));
	}

	//std::cout << "inputsingal " << std::endl;
	//printSignal(inputSignalFloat);
	//std::cout << "DFT " << std::endl;
	//printSignal(myDFT(inputSignalFloat));
	//
	//std::cout << "FFT-Float " << std::endl;
	//printSignal(myFFT(inputSignalFloat));
	//
	//std::cout << "FFTOpt-Float " << std::endl;
	//printSignal(myFFTOptimized(inputSignalFloat));



	auto normalStart = std::chrono::high_resolution_clock::now();
	auto fftOutSignal = myFFT(inputSignalFloat);
	auto normalEnd = std::chrono::high_resolution_clock::now();
	std::uint64_t normalMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(normalEnd - normalStart).count();
	std::cout << "normal time " << normalMilliSeconds << std::endl;
	testSignal(inputSignalFloat, fftOutSignal);

	auto optStart = std::chrono::high_resolution_clock::now();
	auto fftOptOutSignal = myFFTOptimized(inputSignalFloat);
	auto optEnd = std::chrono::high_resolution_clock::now();
	std::uint64_t optMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(optEnd - optStart).count();
	std::cout << "opt time " << optMilliSeconds << std::endl;
	testSignal(inputSignalFloat, fftOptOutSignal);

	std::cin.get();
	return 0;
}

template<typename DataType>
bool testSignal(vector<complex<DataType>>& inputSignal, vector<complex<DataType>>& fftSignal)
{
	auto compareSignal = myDFT(inputSignal);

	const DataType DELTA = DataType(1e-3);// 
	for (std::size_t i = 0; i < compareSignal.size(); ++i)
	{
		DataType diffImag = std::abs(fftSignal[i].imag() - compareSignal[i].imag());

		// accurracy is 0.1%
		DataType diffImagAllowed = std::abs(compareSignal[i].imag())*DELTA;
		if (diffImag > diffImagAllowed) 
		{
			std::cout << "imag part is wrong at index " << i << " fftSignal " << 
				fftSignal[i].imag() << " != " << compareSignal[i].imag() << " difference " << diffImag << " allowed " << diffImagAllowed << std::endl;

			return false;
		}

		DataType diffReal = std::abs(fftSignal[i].real() - compareSignal[i].real());
		DataType diffRealAllowed = std::abs(compareSignal[i].real())*DELTA;
		if (diffReal > diffRealAllowed)
		{
			std::cout << "real part is wrong at index " << i << " fftSignal " <<
				fftSignal[i].real() << " != " << compareSignal[i].real() << " difference " << diffReal << " allowed " << diffRealAllowed << std::endl;

			return false;
		}
	}

	std::cout << "OK" << std::endl;
	return true;
}


template<typename DataType>
void printSignal(vector<complex<DataType>> inputSignal)
{
	std::cout << std::fixed << std::setprecision(5);
	for (size_t i = 0; i < inputSignal.size(); ++i)
	{
		std::cout << inputSignal[i].real();
		
		if (inputSignal[i].imag() >= 0.0f)
		{
			std::cout << "+";
		}

		std::cout << inputSignal[i].imag() << "i ";

		if (i +1 % 10 == 0)
		{
			std::cout << std::endl;
		}
	}

	std::cout << std::endl;
}

template<typename DataType>
vector<complex<DataType>> getCoefficients(size_t N)
{
	// -2*pi/N * j, j is from 0 to N-1
	vector<complex<DataType>> coefficients;
	coefficients.resize(N*N);
	using namespace std::complex_literals;

	DataType angularStep = DataType(-2.0) * static_cast<DataType>(M_PI) / static_cast<DataType>(N);

	DataType angularArgument;
	for (std::size_t j = 0; j < coefficients.size(); ++j)
	{
		angularArgument = angularStep * static_cast<DataType>(j);
		coefficients[j] = std::exp(angularArgument * std::complex<DataType>(DataType(0), DataType(1)));
	}

	return coefficients;
}

template<typename DataType>
vector<complex<DataType>> myDFT(vector<complex<DataType>> inputSignal)
{
	size_t N = inputSignal.size();
	vector<complex<DataType>> coefficients = getCoefficients<DataType>(N);

	vector<complex<DataType>> outputSignal;
	outputSignal.resize(N);

	for (size_t k = 0; k < N; ++k)
	{
		for (size_t j = 0; j < N; ++j)
		{
			outputSignal[k] += coefficients[k*j] * inputSignal[j];
		}
	}

	return outputSignal;
}


template<typename DataType>
vector<complex<DataType>> myFFT(vector<complex<DataType>> inputSignal)
{
	size_t N = inputSignal.size();
	vector<complex<DataType>> outSignal;
	outSignal.resize(N);
	if (N == 1) 
	{	
		outSignal[0] = inputSignal[0];
	}
	else
	{
		vector<complex<DataType>> coefficients = getCoefficients<DataType>(N);

		vector<complex<DataType>> evenSignal;
		evenSignal.resize(N / 2);

		vector<complex<DataType>> oddSignal;
		oddSignal.resize(N / 2);

		for (std::size_t i = 0; i < N / 2; i++)
		{
			evenSignal[i] = inputSignal[i] + inputSignal[N / 2 + i];
			oddSignal[i] = coefficients[i] * (inputSignal[i] - inputSignal[N / 2 + i]);
		}

		auto outEvenSignal = myFFT(evenSignal);
		auto outOddSignal = myFFT(oddSignal);

		for (size_t i = 0; i < N / 2; i++)
		{
			outSignal[2*i] = outEvenSignal[i];
			outSignal[2*i + 1] = outOddSignal[i];
		}
	}

	return outSignal;
}


// Utility function for reversing the bits
// of given index x
size_t bitReverse(size_t x, size_t log2n)
{
	size_t n = 0;
	for (size_t i = 0; i < log2n; i++)
	{
		n <<= 1;
		n |= (x & 1);
		x >>= 1;
	}
	return n;
}


template<typename DataType>
vector<complex<DataType>> myFFTOptimized(vector<complex<DataType>>& inputSignal)
{
	size_t N = inputSignal.size();
	vector<complex<DataType>> reInput;
	reInput.resize(N);

	// first step get given input signal into correct order
	// -> bit reversal
	size_t log2N = static_cast<size_t>(std::log2(N));
	for (size_t i = 0; i < N; ++i)
	{
		size_t reversedIndex = bitReverse(i, log2N);
		reInput[i] = inputSignal[reversedIndex];
	}

	size_t stageN = 1;
	complex<DataType> coefficient;

	// some variables outside the loop
	size_t keven;
	size_t kodd;
	complex<DataType> coeffOddValue;
	complex<DataType> evenValue;
	complex<DataType> oddValue;
	size_t butterflyOffset;
	size_t numButterflys;
	for (size_t stage = 0; stage < log2N; stage++)
	{
		using namespace std::complex_literals;
		DataType m = std::pow(2, static_cast<DataType>(stage + 1));
		complex<DataType> unitRoot = std::exp(DataType(-2.0) * static_cast<DataType>(M_PI) * complex<DataType>(0, 1) / m);

		// how many butterflies we have to calculate for next stage?
		numButterflys = N / (stageN * 2);
		for (size_t i = 0; i < numButterflys; ++i)
		{
			butterflyOffset = i*stageN * 2;
			coefficient = complex<DataType>(DataType(1.0), DataType(0.0));

			// we take fewer but bigger butterflies for higher stages
			for (std::size_t k = 0; k < stageN; k++)
			{
				keven = k + butterflyOffset;
				kodd = k + stageN + butterflyOffset;
				coeffOddValue = coefficient * reInput[kodd];

				// even is "left" and odd is "right" 
				evenValue = reInput[keven] + coeffOddValue;
				oddValue = reInput[keven] - coeffOddValue;
				reInput[keven] = evenValue;
				reInput[kodd] = oddValue;
				coefficient *= unitRoot;
			}
		}

		stageN *= 2;
	}

	return reInput;
}

template<typename DataType>
void four1(double* data, unsigned long nn)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	// reverse-binary reindexing
	n = nn << 1;
	j = 1;
	for (i = 1; i<n; i += 2) {
		if (j>i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j>m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};

	// here begins the Danielson-Lanczos section
	mmax = 2;
	while (n>mmax) {
		istep = mmax << 1;
		theta = -(2 * M_PI / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr*data[j - 1] - wi*data[j];
				tempi = wr * data[j] + wi*data[j - 1];

				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wtemp = wr;
			wr += wr*wpr - wi*wpi;
			wi += wi*wpr + wtemp*wpi;
		}
		mmax = istep;
	}
}