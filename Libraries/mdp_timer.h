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

#include <string>

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
    static char tmp[1024];
    gethostname(tmp, 1024);
    return std::string(tmp);
#endif
  }

  void getcpuusage(double &user, double &total)
  {
#ifdef NO_POSIX
    user = total = 0;
#else
    static mdp_int t[4], s[4];
    double sum, usage[4];
    FILE *fp;
    s[0] = t[0];
    s[1] = t[1];
    s[2] = t[2];
    s[3] = t[3];
    fp = fopen("/proc/stat", "r");
    if (!fp)
      return;
    if (fscanf(fp, "cpu  %i%i%i%i", &t[0], &t[1], &t[2], &t[3]) != 4)
    {
      error("Error while reading cpu stats");
    }
    fclose(fp);
    usage[0] = (t[0] - s[0]);
    usage[1] = (t[1] - s[1]);
    usage[2] = (t[2] - s[2]);
    usage[3] = (t[3] - s[3]);
    sum = usage[0] + usage[1] + usage[2] + usage[3];
    usage[0] /= sum;
    usage[1] /= sum;
    usage[2] /= sum;
    usage[3] /= sum;
    user = 100.0 * (usage[0]);                        // user usage
    total = 100.0 * (usage[0] + usage[1] + usage[2]); // cpu usage
#endif
  }
} // namespace MDP

#endif /* MDP_TIMER_ */
