// This file is part of smpsym, an inverse design and simulation tool for
// research paper Guseinov R. et al "Programming temporal morphing of
// self-actuated shells"
//
// Copyright (C) 2019 Ruslan Guseinov <guseynov.ruslan@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef SINGLE_TASK_THREAD_H
#define SINGLE_TASK_THREAD_H

#include <future>

namespace smpup
{

	// Stoppable async thread handler for class with unique async task
	// *** Usage ***
	// public:
	//   void single_async_task() : task to run in an async thread
	//   bool run_task() : attempt to run the task (returns true on successful launch, false if already running)
	//   bool is_running_task() : check if the task is in progress
	//   void request_task_stop() : notify single_thread_task that a stop has been requested
	//   void wait_for_task() : wait for task to finish
	//   void shut_down_task() : request async task to stop and wait for it to finish
	// protected:
	//   bool task_stop_requested() : call from single_thread_task to check if stop has been requested
	// *** Recommendations ***
	// Always call shut_down_task() from child class destructor before touching data related to the task
	class SingleAsyncTask
	{
	private:
		const std::chrono::seconds zero_time{ 0 };
		mutable std::future<void> is_running_flag;
		mutable std::unique_ptr<std::promise<void>> stop_signal; // stop signal assiciated with current launch
		mutable std::future<void> future_stop_signal;

		// Member function from child class to run in a single async thread
		virtual void single_async_task(void* user_params = nullptr) { throw "single_async_task non-const undefined"; }
		virtual void single_async_task(void* user_params = nullptr) const { throw "single_async_task const undefined"; }

		bool run_task_(void* user_params, bool is_const) const
		{
			if (is_running_task()) return false;
			stop_signal = std::make_unique<std::promise<void>>();
			future_stop_signal = stop_signal->get_future();
			is_running_flag = std::async(std::launch::async, [this, user_params, is_const]() {
				if (is_const) single_async_task(user_params);
				else const_cast<SingleAsyncTask&>(*this).single_async_task(user_params);
			});
			return true;
		}

	protected:

		// Check if the thread is currently stopping; for use in single_thread_task()
		bool task_stop_requested() const
		{
			auto stopping_status = future_stop_signal.wait_for(zero_time);
			if (stopping_status == std::future_status::timeout) return false;
			return true;
		}

	public:

		SingleAsyncTask()
		{
			is_running_flag = std::promise<void>().get_future(); // to keep this future always valid
		}

		// Run the task; returns true if task has been launched
		bool run_task(void* user_params = nullptr) const { return run_task_(user_params, true); }

		// Run the task; returns true if task has been launched
		bool run_task(void* user_params = nullptr) { return run_task_(user_params, false); }

		// Check if the thread is currently running
		bool is_running_task() const
		{
			return is_running_flag.wait_for(zero_time) == std::future_status::timeout;
		}

		// Request async task to stop; must be handeled by single_thread_task()
		void request_task_stop() const
		{
			if (!is_running_task()) return;
			stop_signal->set_value();
		}

		// Wait for task to finish
		inline void wait_for_task() const { is_running_flag.wait(); }

		// Request async task to stop and wait for it to finish
		inline void shut_down_task() const
		{
			request_task_stop();
			wait_for_task();
		}

		~SingleAsyncTask()
		{
			shut_down_task();
		}
	};

}

#endif
