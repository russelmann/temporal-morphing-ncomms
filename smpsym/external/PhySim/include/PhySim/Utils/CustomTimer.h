#pragma once

#ifdef WIN32
#include <windows.h>
#else 
#include <sys/time.h>
#endif

#include <stdlib.h>

#include <PhySim/CommonIncludes.h>

#undef min
#undef max

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class CustomTimer
	{

	public:

		CustomTimer()
		{ 
			this->m_window = 10;
			this->m_identifier = string("");
			this->m_description = string("");
			this->initialize();
		}

		// Default constructor (identifier and description)
		CustomTimer(int window, const string& iden = "", const string& desc = "")
		{
			assert(window >= 1);
			this->m_window = window;
			this->m_identifier = iden;
			this->m_description = desc;
			this->initialize();
		}

		~CustomTimer() {}

		// Reinitialize
		void initialize()
		{
//#pragma omp critical
			{
				this->m_measurements.clear();

#ifdef WIN32
				QueryPerformanceFrequency(&m_frequency);
				m_startCount.QuadPart = 0;
				m_endCount.QuadPart = 0;
#else
				m_startCount.tv_sec = m_startCount.tv_usec = 0;
				m_endCount.tv_sec = m_endCount.tv_usec = 0;
#endif
				m_stopped = 0;
				m_staTimeInMicroSec = 0;
				m_endTimeInMicroSec = 0;
				m_cumTimeInMicroSec = 0;
			}
		}

		// Start timer
		void restart()
		{
			this->start();
		}

		// Do everything
		void stopStoreLog()
		{
			this->stop();
			double time = this->getElapsedTimeInMilliSec();
			this->store(time); // Store the time in millis
			this->log();
		}

		// Store value to vector 
		void store(double value)
		{
			this->m_measurements.push_back(value);
			if ((int) this->m_measurements.size() > this->m_window)
				this->m_measurements.erase(this->m_measurements.begin());
		}

		// Compute mean value
		double computeMean()
		{
			int n = (int) this->m_measurements.size();

			double timeSum = 0.0;
			for (int i = 0; i < n; ++i)
				timeSum += this->m_measurements[i];

			return timeSum / (double)n;
		}

		// Standard deviation
		double computeSD()
		{
			int n = (int) this->m_measurements.size();
			double mean = (double) this->computeMean();

			double diffSum = 0.0;
			for (int i = 0; i < n; ++i)
			{
				double diff = this->m_measurements[i] - mean;
				diffSum += diff*diff; // SD: E[(X - XBar)^2]
			}

			return diffSum / (double)n;
		}

		// Log state
		void log()
		{
//#pragma omp critical
			{
				double sd = this->computeSD();
				double mean = this->computeMean();
				logTime("%s \t %.9f \t %.9f\n", this->m_identifier.c_str(), mean, sd);
			}
		}

		const string& getIdentifier() const { this->m_identifier; }
		const string& getDescription() const { this->m_description; }

		void start()
		{
//#pragma omp critical
			{
				m_cumTimeInMicroSec = 0; //Reset 
				m_stopped = 0; // Reset stop flag
#ifdef WIN32
				QueryPerformanceCounter(&m_startCount);
#else
				gettimeofday(&m_startCount, NULL);
#endif
			}
		}

		void stop()
		{
//#pragma omp critical
			{
				m_stopped = 1; // Set stop flag
#ifdef WIN32
				QueryPerformanceCounter(&m_endCount);
#else
				gettimeofday(&m_endCount, NULL);
#endif
			}
		}

		void pause()
		{
//#pragma omp critical
			{
				this->stop();
				m_cumTimeInMicroSec += getElapsedTimeInMicroSec();
			}
		}

		void resume()
		{
//#pragma omp critical
			{
				double accumulated = this->m_cumTimeInMicroSec;
				this->start(); // This set accumulated time to 0
				this->m_cumTimeInMicroSec = accumulated; // Reset
			}
		}

		// Seconds
		double getElapsedTime()
		{
			return this->getElapsedTimeInSec();
		}

		// Seconds
		double getElapsedTimeInSec()
		{
			return this->getElapsedTimeInMicroSec() * 0.000001;
		}

		// Milli-seconds
		double getElapsedTimeInMilliSec()
		{
			return this->getElapsedTimeInMicroSec() * 0.001;
		}

		// Micro-seconds
		double getElapsedTimeInMicroSec()
		{
//#pragma omp critical
			{
#ifdef WIN32
				if (!m_stopped)
					QueryPerformanceCounter(&m_endCount);

				m_staTimeInMicroSec = m_startCount.QuadPart * (1000000.0 / m_frequency.QuadPart);
				m_endTimeInMicroSec = m_endCount.QuadPart * (1000000.0 / m_frequency.QuadPart);
#else
				if (!m_stopped)
					gettimeofday(&m_endCount, NULL);

				m_staTimeInMicroSec = (m_startCount.tv_sec * 1000000.0) + m_startCount.tv_usec;
				m_endTimeInMicroSec = (m_endCount.tv_sec * 1000000.0) + m_endCount.tv_usec;
#endif
			}

			return (m_endTimeInMicroSec - m_staTimeInMicroSec) + m_cumTimeInMicroSec;
		}

		int m_window;
		string m_identifier;
		string m_description;
		vector<double> m_measurements;


		double m_staTimeInMicroSec;	// Starting time in micro-seconds
		double m_endTimeInMicroSec;	// Ending time in micro-second
		double m_cumTimeInMicroSec;	// Storing cumulative time

		int  m_stopped; // Stop flag
#ifdef WIN32
		// Ticks per second
		LARGE_INTEGER m_frequency;
		LARGE_INTEGER m_startCount;
		LARGE_INTEGER m_endCount;
#else
		timeval m_startCount;
		timeval m_endCount;
#endif

	};
}