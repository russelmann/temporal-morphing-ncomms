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

#include "utils/utils.h"
#ifndef _WIN32
#include <unistd.h>
#endif

void run_process(std::string cmd)
{
#ifdef _WIN32
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;

	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));

	if (CreateProcessA(NULL, const_cast<char *>(cmd.c_str()),
		NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi))
	{
		WaitForSingleObject(pi.hProcess, INFINITE);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
	}
#else
	system(cmd.c_str());
#endif
}

void process_sleep(int milliseconds)
{
#ifdef _WIN32
	Sleep(milliseconds);
#else
	usleep(1000 * milliseconds);
#endif
}

void os_reveal_folder(std::string path)
{
#ifdef _WIN32
	run_process("explorer \"" + path +"\"");
#else
	run_process("finder \"" + path + "\"");
#endif
}

