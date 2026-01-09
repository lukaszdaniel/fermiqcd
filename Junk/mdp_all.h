#ifndef MDP_ALL_
#define MDP_ALL_

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX 1
#endif

#include <winsock2.h>
#include <ws2tcpip.h>
#include <windows.h>
#include <process.h> // _getpid
#include <io.h>

#else

#include <unistd.h>
#include <fcntl.h>
#include <csignal>
#include <pthread.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <sys/wait.h>
#include <sys/file.h>
#include <sys/un.h>
#include <sys/select.h>

#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <poll.h>
#endif

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

// C++ headers and STL headers
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <map>

#ifndef HAVE_INET_NTOP
#define inet_ptoa(a, b) inet_ntoa(b)
#define inet_pton(a, b, c) inet_aton(b, c)
#endif

namespace MDP
{
  void exit_message(int en, const std::string &message)
  {
    std::cerr << "FROM PROCESS PID: " << getpid() << std::endl;
    std::cerr << "CHILD OF PROCESS PID: " << getppid() << std::endl;
    std::cerr << "FATAL ERROR: " << message << std::endl;
    std::cerr << "EXITING WITH ERROR NUMBER: " << en << std::endl;
    std::exit(en);
  }
} // namespace MDP

#endif /* MDP_ALL_ */
