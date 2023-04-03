// BEGIN FILE: mdp_all.h
// C headers
#include "sys/types.h"
#ifndef _WIN64
#include "sys/socket.h"
#endif
#include "sys/time.h"
#include <ctime>
#include "netinet/in.h"
#include "arpa/inet.h"
#include <cerrno>
#include "fcntl.h"
#include "netdb.h"
#include "signal.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "sys/stat.h"
#include "sys/uio.h"
#include "unistd.h"
#include "sys/wait.h"
#include "sys/un.h"
#include "sys/select.h"
#include "poll.h"
#include "strings.h"
#include "pthread.h"

// C++ headers and STL headers
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <map>

#ifndef HAVE_INET_NTOP
#define inet_ntop(a, b) inet_ntoa(b)
#define inet_pton(a, b, c) inet_aton(b, c)
#endif

void exit_message(int en, std::string message)
{
  std::cerr << "FROM PROCESS PID: " << getpid() << std::endl;
  std::cerr << "CHILD OF PROCESS PID: " << getppid() << std::endl;
  std::cerr << "FATAL ERROR: " << message << std::endl;
  std::cerr << "EXITING WITH ERROR NUMBER: " << en << std::endl;
  exit(en);
}
