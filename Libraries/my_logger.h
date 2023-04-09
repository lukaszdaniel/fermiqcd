#ifndef MY_LOGGER_
#define MY_LOGGER_

#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

namespace MDP
{
  FILE *stdlog = NULL;
  static time_t time_start;
  static time_t time_finish;

  static clock_t CPU_start;
  static clock_t CPU_finish;

  static bool time_ok = false;
  static bool CPU_ok = false;

  int openlog(const char log[])
  {
    if (strncmp(log, "stderr", 7) == 0)
      stdlog = stderr;
    else if (strncmp(log, "stdout", 7) == 0)
      stdlog = stdout;
    else if ((stdlog = fopen(log, "w")) == NULL)
    {
      fprintf(stderr, "cannot open log file %s\n", log);
      exit(1);
    }
    return 0;
  }

  int closelog()
  {
    return fclose(stdlog);
  }

  void footer()
  {
    time_t current_time;

    char time_string[128];

    if ((current_time = time(NULL)) != -1)
      strftime(time_string, 64, "%d.%m.%y %H:%M:%S",
               localtime(&current_time));
    else
    {
      time_ok = false;
      strncpy(time_string, " time not available\n", 128);
    }

    time_finish = current_time;

    if ((CPU_finish = clock()) == -1)
    {
      CPU_ok = false;
    }

    if (CPU_ok)
    {
      fprintf(stdlog, "# CPU time %12li   sec\n",
              (CPU_finish - CPU_start) / CLOCKS_PER_SEC);
    }

    if (time_ok)
    {
      fprintf(stdlog, "# user time %13.1f sec\n",
              difftime(time_finish, time_start));
    }
    fprintf(stdlog, "# time stop  %s\n", time_string);

    fflush(stdlog);
  }

  void header(int argc, char *argv[], const char id[])
  {
    time_t current_time;
    char time_string[128];
    char host_name[128];

    if (stdlog == NULL)
    {
      fprintf(stderr, "stdlog not opened\n");
      exit(1);
    }

    if ((current_time = time(NULL)) != -1)
    {
      strftime(time_string, 64, "%d.%m.%y %H:%M:%S",
               localtime(&current_time));
      time_ok = true;
    }
    else
      strncpy(time_string, " time not available\n", 128);

    time_start = current_time;

    if ((CPU_start = clock()) != -1)
    {
      CPU_ok = true;
    }

    fprintf(stdlog, "# program    %s   %s\n", argv[0], id);
    fprintf(stdlog, "# time start %s\n", time_string);

    if (gethostname(host_name, 128) == 0)
    {
      fprintf(stdlog, "# host %s\n", host_name);
    }

    for (int i = 0; i < argc; i++)
      fprintf(stdlog, "# argv %2d %s\n", i, argv[i]);

    fprintf(stdlog, "# header end\n");

    if (atexit(footer))
    {
      fprintf(stdlog, "could not register function footer()\n");
    }

    fflush(stdlog);
  }
} // namespace MDP

#endif /* MY_LOGGER_ */
