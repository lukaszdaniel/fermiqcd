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

#ifdef _WIN32
#define DATA_CAST (char *)
#define CONST_DATA_CAST (const char *)
#else
#define DATA_CAST
#define CONST_DATA_CAST
#endif

namespace MDP
{
  void exit_message(int en, const std::string &message)
  {
#ifdef _WIN32
    std::cerr << "FROM PROCESS PID: " << _getpid() << std::endl;
    std::cerr << "CHILD OF PROCESS PID: (unsupported on Windows)" << std::endl;
#else
    std::cerr << "FROM PROCESS PID: " << getpid() << std::endl;
    std::cerr << "CHILD OF PROCESS PID: " << getppid() << std::endl;
#endif

    std::cerr << "FATAL ERROR: " << message << std::endl;
    std::cerr << "EXITING WITH ERROR NUMBER: " << en << std::endl;
    std::exit(en);
  }

  class InternetAddress
  {
  private:
    struct sockaddr_in m_address;
    std::string m_ipaddress;
    int m_port;

  public:
    InternetAddress(const std::string &hostname = "127.0.0.1", int port = 0)
    {
      struct addrinfo hints;
      struct addrinfo *result = nullptr;

      std::memset(&hints, 0, sizeof(hints));
      hints.ai_family = AF_INET;       // IPv4 (use AF_UNSPEC for IPv6)
      hints.ai_socktype = SOCK_STREAM; // TCP

      std::string port_str = std::to_string(port);

      int ret = getaddrinfo(hostname.c_str(), port_str.c_str(), &hints, &result);
      if (ret != 0)
        throw std::string(gai_strerror(ret));

      // take the first correct address
      struct sockaddr_in *addr =
          reinterpret_cast<struct sockaddr_in *>(result->ai_addr);

      m_address = *addr;
      m_port = port;

      char ipstr[INET_ADDRSTRLEN];
      if (!inet_ntop(AF_INET, &addr->sin_addr, ipstr, sizeof(ipstr)))
      {
        freeaddrinfo(result);
        exit_message(1, "invalid IP address");
      }

      m_ipaddress = ipstr;

      freeaddrinfo(result);
    }

    struct sockaddr_in address() const
    {
      return m_address;
    }

    int getPort()
    {
      return ntohs(m_address.sin_port);
    }

    std::string getIPAddress()
    {
      return inet_ptoa(AF_INET, m_address.sin_addr);
    }

    int Connect(int sfd, int timeout = 0)
    {
      if (timeout == 0)
      {
#ifdef _WIN32
        return connect(sfd, (struct sockaddr *)&m_address,
                       (int)sizeof(m_address));
#else
        return connect(sfd, (struct sockaddr *)&m_address,
                       (socklen_t)sizeof(m_address));
#endif
      }

#ifdef _WIN32
      u_long mode = 1; // non-blocking
      ioctlsocket(sfd, FIONBIO, &mode);

      int ret = connect(sfd, (struct sockaddr *)&m_address,
                        (int)sizeof(m_address));

      if (ret == SOCKET_ERROR)
      {
        int err = WSAGetLastError();
        if (err != WSAEWOULDBLOCK && err != WSAEINPROGRESS)
          return ECONNREFUSED;

        WSAPOLLFD fds;
        fds.fd = sfd;
        fds.events = POLLWRNORM;
        fds.revents = 0;

        ret = WSAPoll(&fds, 1, timeout);
        if (ret <= 0)
          return ETIMEDOUT;
      }

      // back to blocking
      mode = 0;
      ioctlsocket(sfd, FIONBIO, &mode);

      int valopt = 0;
      int lon = sizeof(valopt);
      if (getsockopt(sfd, SOL_SOCKET, SO_ERROR,
                     (char *)&valopt, &lon) < 0 ||
          valopt != 0)
        return ECONNREFUSED;

      return 0;

#else
      int arg = fcntl(sfd, F_GETFL, NULL);
      fcntl(sfd, F_SETFL, arg | O_NONBLOCK);

      pollfd fds;
      fds.fd = sfd;
      fds.events = POLLWRNORM;
      fds.revents = 0;

      int ret = connect(sfd, (struct sockaddr *)&m_address,
                        (socklen_t)sizeof(m_address));

      if (ret < 0)
      {
        if (errno != EINPROGRESS)
          return ECONNREFUSED;

        ret = poll(&fds, 1, timeout);
        if (ret <= 0)
          return ETIMEDOUT;
      }

      // back to blocking
      fcntl(sfd, F_SETFL, arg);

      int valopt = 0;
      socklen_t lon = sizeof(valopt);
      if (getsockopt(sfd, SOL_SOCKET, SO_ERROR,
                     (void *)&valopt, &lon) < 0 ||
          valopt != 0)
        return ECONNREFUSED;

      return 0;
#endif
    }

    int Accept(int fd)
    {
      socklen_t ssize = sizeof(struct sockaddr_in);
      return accept(fd, (struct sockaddr *)&m_address, &ssize);
    }

    int sendTo(int sfd, void *data, int size, int options = 0)
    {
      return sendto(sfd, CONST_DATA_CAST data, size, options,
                    (struct sockaddr *)&m_address,
                    (socklen_t)sizeof(m_address));
    }

    int recvFrom(int sfd, void *data, int size, int options = 0)
    {
      socklen_t address_size = sizeof(m_address);
      return recvfrom(sfd, DATA_CAST data, size, options,
                      (struct sockaddr *)&m_address,
                      (socklen_t *)&address_size);
    }
  };

  int newTcpSocket(int flags = 0)
  {
    return socket(AF_INET, SOCK_STREAM, flags);
  }

  int newUdpSocket(int flags = 0)
  {
    return socket(AF_INET, SOCK_DGRAM, flags);
  }

  int newTcpClientSocket(const std::string &ipaddress, int port, int sleep_time = 10)
  {
    InternetAddress peer = InternetAddress(ipaddress, port);
    int sfd = newTcpSocket();
    while (1)
    {
      int r = peer.Connect(sfd);
      if (r == 0)
        return sfd;
      if (r != ETIMEDOUT)
      {
        close(sfd);
        return -1;
      }
#ifdef _WIN32
      Sleep(sleep_time * 1000); // Windows Sleep in ms
#else
      sleep(sleep_time); // POSIX sleep in s
#endif
    }
  }

  int Bind(int sfd, int port)
  {
    int server_fd = sfd;
    struct sockaddr_in address;
    memset(&address, 0, sizeof(address));
    address.sin_family = AF_INET;
    address.sin_port = htons(port);
    address.sin_addr.s_addr = htonl(INADDR_ANY);
    return ::bind(server_fd, (struct sockaddr *)&address, (socklen_t)sizeof(address));
  }

  int setSocketKeepAlive(int sfd, int on = 1)
  {
    return setsockopt(sfd, SOL_SOCKET, SO_KEEPALIVE, CONST_DATA_CAST &on, sizeof(on));
  }

  int setFileLock(int sfd)
  {
#ifdef _WIN32
    HANDLE h = (HANDLE)_get_osfhandle(sfd);
    if (h == INVALID_HANDLE_VALUE)
      return -1;

    OVERLAPPED ov = {};
    return LockFileEx(
               h,
               LOCKFILE_EXCLUSIVE_LOCK,
               0,
               MAXDWORD,
               MAXDWORD,
               &ov)
               ? 0
               : -1;
#else
    return flock(sfd, LOCK_EX);
#endif
  }

  int setFileLockShared(int sfd)
  {
#ifdef _WIN32
    HANDLE h = (HANDLE)_get_osfhandle(sfd);
    if (h == INVALID_HANDLE_VALUE)
      return -1;

    OVERLAPPED ov = {};
    return LockFileEx(
               h,
               0, // shared
               0,
               MAXDWORD,
               MAXDWORD,
               &ov)
               ? 0
               : -1;
#else
    return flock(sfd, LOCK_SH);
#endif
  }

  int setFileUnlock(int sfd)
  {
#ifdef _WIN32
    HANDLE h = (HANDLE)_get_osfhandle(sfd);
    if (h == INVALID_HANDLE_VALUE)
      return -1;

    OVERLAPPED ov = {};
    return UnlockFileEx(
               h,
               0,
               MAXDWORD,
               MAXDWORD,
               &ov)
               ? 0
               : -1;
#else
    return flock(sfd, LOCK_UN);
#endif
  }

  int setSocketNonblock(int sfd, int on = 1)
  {
#ifdef _WIN32
    u_long mode = on ? 1 : 0;
    return ioctlsocket(sfd, FIONBIO, &mode);
#else
    int flags = fcntl(sfd, F_GETFL, 0);
    if (flags < 0)
      return flags;
    if (on)
      flags |= O_NONBLOCK;
    else
      flags &= ~O_NONBLOCK;
    return fcntl(sfd, F_SETFL, flags);
#endif
  }

  int setSocketAsync(int sfd, int owner = 0, int on = 1)
  {
#ifdef _WIN32
    // Windows doesn't support FASYNC/F_SETOWN
    return setSocketNonblock(sfd + 0 * owner, on);
#else
    if (owner == 0)
      owner = getpid();
    int flags = fcntl(sfd, F_GETFL, 0);
    if (flags < 0)
      return flags;
    if (on)
      flags |= O_NONBLOCK | FASYNC;
    else
      flags &= ~(O_NONBLOCK | FASYNC);
    if (fcntl(sfd, F_SETFL, flags) < 0)
      return -1;
    if (fcntl(sfd, F_SETOWN, owner) < 0)
      return -1;
    return 1;
#endif
  }

  int setSocketReusable(int sfd, int on = 1)
  {
    return setsockopt(sfd, SOL_SOCKET, SO_REUSEADDR, CONST_DATA_CAST &on, sizeof(on));
  }

  int setSocketSendBroadcast(int sfd, int on = 1)
  {
    return setsockopt(sfd, SOL_SOCKET, SO_BROADCAST, CONST_DATA_CAST &on, sizeof(on));
  }

  int setSocketRecvBroadcast(int sfd, int on = 1)
  {
    return setSocketReusable(sfd, on);
  }

  int setSocketSendMulticast(int sfd, int ttl = 1)
  {
    return setsockopt(sfd, IPPROTO_IP, IP_MULTICAST_TTL, CONST_DATA_CAST &ttl, sizeof(ttl));
  }

  int setSocketRecvMulticast(int sfd, const std::string &from)
  {
    if (setSocketReusable(sfd, 1) < 0)
      return -1;

    struct ip_mreq mreq;
    mreq.imr_multiaddr.s_addr = InternetAddress(from, 0).address().sin_addr.s_addr;
    mreq.imr_interface.s_addr = htonl(INADDR_ANY);
    return setsockopt(sfd, IPPROTO_IP, IP_ADD_MEMBERSHIP, CONST_DATA_CAST &mreq, sizeof(mreq));
  }

  int forkTwice()
  {
#ifndef _WIN32
    pid_t pid;
    int status;

    pid = fork();
    if (pid == 0)
    {
      int pid2 = fork();
      switch (pid2)
      {
      case 0:
        return 0; // subprocess
      case -1:
        _exit(errno); // assumes all errnos are <256
      default:
        _exit(0); // parent of first fork
      }
    }
    else if (pid > 0)
    {
      if (waitpid(pid, &status, 0) < 0)
        return -1;
      else
        return 1; // parent
    }
    else
    {
      return -1; // fork error
    }
#endif
    return -1;
  }

  class Thread
  {
  private:
    pthread_t _thread_number;
    pthread_attr_t _thread_attributes;

  public:
    Thread();
    void threadStart();
    void threadStop();
    void threadSetJoinable();
    void threadJoin();
    virtual void run() = 0;
  };

  void *thread_function(void *p)
  {
    Thread *pt = (Thread *)p;
    pt->run();
    return pt;
  }

  Thread::Thread()
  {
    pthread_attr_init(&_thread_attributes);
  }

  void Thread::threadStart()
  {
    pthread_create(&_thread_number, &_thread_attributes, thread_function, (void *)this);
  }

  void Thread::threadStop()
  {
    pthread_exit(NULL);
  }

  void Thread::threadSetJoinable()
  {
    pthread_attr_setdetachstate(&_thread_attributes, PTHREAD_CREATE_JOINABLE);
  }

  void Thread::threadJoin()
  {
    pthread_join(_thread_number, NULL);
  }

  void _handler(int);

  class SignalHandler
  {
  public:
    virtual void handler(int s) = 0; // metoda wirtualna

    void catch_signal(int signalnum);

  private:
#ifndef _WIN32
    struct sigaction action;
#endif

    static void _handler(int s);
  };

  static std::map<int, SignalHandler *> _signal_handlers;

  void SignalHandler::_handler(int s)
  {
    auto it = _signal_handlers.find(s);
    if (it != _signal_handlers.end() && it->second)
    {
      it->second->handler(s);
    }
  }

  void SignalHandler::catch_signal(int signalnum)
  {
#ifndef _WIN32
    if (sigemptyset(&action.sa_mask) != 0)
      throw std::runtime_error("sigemptyset failed");

    if (sigaddset(&action.sa_mask, signalnum) != 0)
      throw std::runtime_error("sigaddset failed");

    action.sa_handler = _handler;
    action.sa_flags = 0;

    if (sigaction(signalnum, &action, nullptr) != 0)
      throw std::runtime_error("sigaction failed");
#else
    // Windows doesn't have sigaction, using signal()
    if (signal(signalnum, _handler) == SIG_ERR)
      throw std::runtime_error("signal failed");
#endif

    _signal_handlers[signalnum] = this;
  }

  void _handler(int s)
  {
    _signal_handlers[s]->handler(s);
  }

#ifndef _WIN32
  //
  // send the signum signal to the process pid
  //
  void signalSend(pid_t pid, int signum)
  {
    if (kill(pid, signum) < 0)
      throw std::runtime_error("kill failed");
  }

  //
  // block the signal signum
  //
  sigset_t signalBlock(int signum)
  {
    sigset_t set, oset;
    if (sigemptyset(&set) != 0)
      throw std::runtime_error("sigemptyset failed");
    if (sigaddset(&set, signum) != 0)
      throw std::runtime_error("sigaddset failed");
    if (sigprocmask(SIG_BLOCK, &set, &oset) != 0)
      throw std::runtime_error("sigprocmask failed");
    return oset;
  }

  //
  // check if signal signum is pending
  //
  bool signalPending(int signum)
  {
    sigset_t set;
    if (sigpending(&set) != 0)
      throw std::runtime_error("sigpending failed");
    return sigismember(&set, signum);
  }

  //
  // unblock the signal signum
  //
  sigset_t signalUnblock(int signum)
  {
    sigset_t set, oset;
    if (sigemptyset(&set) != 0)
      throw std::runtime_error("sigemptyset failed");
    if (sigaddset(&set, signum) != 0)
      throw std::runtime_error("sigaddset failed");
    if (sigprocmask(SIG_UNBLOCK, &set, &oset) != 0)
      throw std::runtime_error("sigprocmask failed");
    return oset;
  }

  //
  // restore signal blocks
  //
  sigset_t signalRestoreBlocks(sigset_t set)
  {
    sigset_t oset;
    if (sigprocmask(SIG_SETMASK, &set, &oset) != 0)
      throw std::runtime_error("sigprocmask failed");
    return oset;
  }
#endif

#undef DATA_CAST
#undef CONST_DATA_CAST
} // namespace MDP

#endif /* MDP_ALL_ */
