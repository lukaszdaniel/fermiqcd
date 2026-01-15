/////////////////////////////////////////////////////////////////
/// @file mdp_timer.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains timing functions including functions to get cpu usage
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_TIMER_
#define MDP_TIMER_

#include <fstream>
#include <string>
#include <sstream>
#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#else
#include <unistd.h>
#endif
#include <sys/time.h>
#include "mdp_global_vars.h"
#include "mdp_macros.h"

namespace MDP
{
  double walltime()
  {
#ifdef NO_POSIX
    return (double)clock() / CLOCKS_PER_SEC;
#else
    double mic, time;
    double mega = 0.000001;
    struct timeval tp;
#ifdef HAVE_NO_TIMEZONE
    struct timezone
    {
      int tz_minuteswest; // minutes of Greenwich
      int tz_dsttime;     // type of dst correction
    } tzp = {0, 0};
#else
    struct timezone tzp = {0, 0};
#endif
    static mdp_int base_sec = 0;
    static mdp_int base_usec = 0;

    gettimeofday(&tp, &tzp);

    if (base_sec == 0)
    {
      base_sec = tp.tv_sec;
      base_usec = tp.tv_usec;
    }

    time = (double)(tp.tv_sec - base_sec);
    mic = (double)(tp.tv_usec - base_usec);
    time = (time + mic * mega);
    return (time);
#endif
  }

  std::string getname()
  {
#ifdef NO_POSIX
    return std::string("localhost");
#else
#ifdef _WIN32
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData); // Winsock initialisation
#endif

    static char tmp[1024];
    gethostname(tmp, sizeof(tmp));

#ifdef _WIN32
    WSACleanup(); // optional in case Winsock is no longer needed
#endif

    return std::string(tmp);
#endif
  }

  void getcpuusage(double &user, double &total)
  {
#ifdef NO_POSIX
    user = total = 0.0;
#else
    static mdp_int curr[4] = {0, 0, 0, 0};
    static mdp_int prev[4] = {0, 0, 0, 0};

    for (int i = 0; i < 4; ++i)
      prev[i] = curr[i];

    std::ifstream file("/proc/stat");
    if (!file)
    {
      user = total = 0.0;
      return;
    }

    std::string cpu;
    file >> cpu;

    if (cpu != "cpu")
      return;

    if (!(file >> curr[0] >> curr[1] >> curr[2] >> curr[3]))
    {
      error("Error while reading cpu stats");
    }

    double usage[4];
    double sum = 0.0;

    for (int i = 0; i < 4; ++i)
    {
      usage[i] = 1.0 * (curr[i] - prev[i]);
      sum += usage[i];
    }

    for (int i = 0; i < 4; ++i)
      usage[i] /= sum;

    user = 100.0 * usage[0];                          // user usage
    total = 100.0 * (usage[0] + usage[1] + usage[2]); // total cpu usage
#endif
  }
} // namespace MDP

#endif /* MDP_TIMER_ */
