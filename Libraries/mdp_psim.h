/////////////////////////////////////////////////////////////////
/// @file mdp_psim.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_psim (the parallel simulator)
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PSIM_
#define MDP_PSIM_

#ifndef NO_POSIX

#include <unistd.h> // for getpid()
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <iostream>
#include <string>
#include <vector>
#include <map>
// #include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN64
#include <sys/socket.h>
#include <sys/errno.h>
#endif
#include <sys/file.h> // for flock()
#include <fcntl.h>
#include "mdp_global_vars.h"

namespace MDP
{
#ifndef FERMIQCD
  using mdp_int =  int; // this in case this lib is used without fermiqcd
#endif

  /// @brief Parallel SIMulator used by class mdp_communicator
  ///
  /// Attention: under MDP and/or FermiQCD this is already
  /// Instantiated inside class mdp_communicator.
  ///
  /// Example:
  /// @verbatim
  /// int main(int argc, char** argv) {
  ///    mdp_psim node(argc,argv);
  ///    int a=3, b=0;
  ///    if(node.id()==0) node.send(1,a);
  ///    if(node.id()==1) { node.recv(0,b); cout << b << endl;
  ///    return 0;
  /// }
  /// @endverbatim
  /// Compile with
  /// @verbatim
  ///    g++ [filename] -o a.out
  /// @endverbatim
  /// and run with
  /// @verbatim
  ///    ./a.out -PSIM_NPROCS=2
  /// @endverbatim
  /// Output should be 3.
  class mdp_psim
  {
  private:
    // Typed Constants
    static constexpr int PROCESS_COUNT_MIN = 1;        // minimum number of processes
    static constexpr int PROCESS_COUNT_MAX = 128;      // maximum number of processes
    static constexpr int CONN_LIST_IGNORE = -1;        // connections that cannot occur
    static constexpr int PROCESS_PARENT = 0;           // The parent process ID number
    static constexpr int COMM_RECV = 0;                // socket array indicator for reading
    static constexpr int COMM_SEND = 1;                // socket array indicator for writing
    static constexpr int COMM_TIMEOUT_DEFAULT = 86400; // 1 day default

    // common enum values for logging routines
    enum enumBegEnd
    {
      LOG_BEGIN,
      LOG_END
    };

    enum enumSendRecv
    {
      LOG_SR_SEND,
      LOG_SR_RECV
    };

    enum enumSendRecvStep
    {
      LOG_SR_START,
      LOG_SR_FAIL,
      LOG_SR_SUCCESS
    };

    // Set this to true if you want verbose testing description

    // Class variables
    int _verbatim;            // 0 for no output, 1 for some, 2 more
    int _processCount;        // Holds the number of processes
    std::string _logFileName; // filename of the log file
    bool _doLogging;          // do logging or not?
    FILE *_logfileFD;         // file descriptor for the logging file
    int _processID;           // process ID of "this" process

    /** @brief 2D array to hold all of the sockets
     *
     * First entry is the table of processes
     * Second entry is the socket array indicator for reading/writing
     */
    int (*_socketFD)[2];

    /** @brief defaults to COMM_TIMEOUT_DEFAULT
     */
    int _commTimeout;

    /** @brief Hash Map to hold out of sequence (send/receive) data
     */
    std::map<std::string, std::vector<char>> *_hash;

    /** @brief Initialise processes
     *
     * @note Used by the constructor ONLY
     */
    void psim_begin(int processCount, std::string logFileName, int verbatim)
    {
      _processCount = processCount;
      _logFileName = logFileName;
      _verbatim = verbatim;

      open_log();
      if ((processCount < PROCESS_COUNT_MIN) || (processCount > PROCESS_COUNT_MAX))
      {
        log("PSIM ERROR: Invalid number of processes");
        throw std::string("PSIM ERROR: Invalid number of processes");
      }

      initialize(processCount);
      create_sockets();
      fork_processes();
      close_sockets();

      char buffer[256];
      snprintf(buffer, 256, "process %i of %i created with pid=%i",
               _processID, _processCount, getpid());
      log(buffer, 1);
    }

    /** @brief Used by the destructor ONLY
     */
    void psim_end()
    {
      for (int source = 0; source < _processCount; source++)
      {
        for (int dest = 0; dest < _processCount; dest++)
        {
          if (dest == _processID)
            close(_socketFD[_processCount * source + dest][COMM_SEND]);
          if (source == _processID)
            close(_socketFD[_processCount * source + dest][COMM_RECV]);
        }
      }

      if (_socketFD != NULL)
        delete[] _socketFD;
      _socketFD = NULL;

      // only delete the _hash if it still exists
      if (_hash != NULL)
      {
        delete[] _hash;
        _hash = NULL;
      }

      char buffer[256];
      snprintf(buffer, 256, "process %i terminating", _processID);
      log(buffer, 1);

      close_log();
    }

    /** @brief Used by the constructor, this method sets up values and
     * some of the needed resources.
     */
    void initialize(int processCount)
    {
      _processCount = processCount;
      _processID = PROCESS_PARENT;
      _commTimeout = COMM_TIMEOUT_DEFAULT;

      _hash = new std::map<std::string, std::vector<char>>[_processCount];
      if (!_hash)
      {
        log("PSIM ERROR: failure to allocate hash");
        throw std::string("PSIM ERROR: failure to allocate hash");
      }

      _socketFD = new int[_processCount * _processCount][2];
      if (!_socketFD)
      {
        log("PSIM ERROR: failed to create socket array");
        throw std::string("PSIM ERROR: failed to create socket array");
      }
    }

    // *******************************************************************
    // ***         Private Method: create_sockets                      ***
    // ***                                                             ***
    // ***  Opens all of the sockets necessary for communication and   ***
    // ***  inserts them into an array.                                ***
    // ***                                                             ***
    // *******************************************************************

    void create_sockets()
    {
      char filename[512];
      for (int source = 0; source < _processCount; source++)
      {
        for (int dest = 0; dest < _processCount; dest++)
        {

          snprintf(filename, 512, ".fifo.%i.%i", source, dest);
          while (mkfifo(filename, 0666) < 0)
          {
            if (errno == EEXIST)
              unlink(filename);
            else
              throw std::string("PSIM ERROR: unable to mkfifo ") + filename;
          }
          _socketFD[_processCount * source + dest][COMM_RECV] = open(filename, O_RDONLY | O_NONBLOCK);
          _socketFD[_processCount * source + dest][COMM_SEND] = open(filename, O_WRONLY);
          int flags = fcntl(_socketFD[_processCount * source + dest][COMM_RECV], F_GETFL, 0);
          fcntl(_socketFD[_processCount * source + dest][COMM_RECV], F_SETFL, flags & ~O_NONBLOCK);
          if (_socketFD[_processCount * source + dest][COMM_SEND] <= 0 ||
              _socketFD[_processCount * source + dest][COMM_RECV] <= 0)
          {
            throw std::string("PSIM ERROR: unable to open fifo");
          }
        }
      }
      char buffer[256];
      for (int source = 0; source < _processCount; source++)
        for (int dest = 0; dest < _processCount; dest++)
        {
          snprintf(buffer, 256, "_socketFD[%i*%i+%i]={%i,%i}",
                   source, _processCount, dest,
                   _socketFD[_processCount * source + dest][COMM_SEND],
                   _socketFD[_processCount * source + dest][COMM_RECV]);
          log(buffer);
        }
    }
    void close_sockets()
    {
      for (int source = 0; source < _processCount; source++)
        for (int dest = 0; dest < _processCount; dest++)
        {
          if (dest != _processID)
            close(_socketFD[_processCount * source + dest][COMM_RECV]);
          if (source != _processID)
            close(_socketFD[_processCount * source + dest][COMM_SEND]);
        }
    }

    /** @brief Fork processes
     */
    void fork_processes()
    {
      _processID = 0;
      for (int i = 1; i < _processCount; i++)
      {
        int pid = fork();

        if (pid == -1)
        {
          log("PSIM ERROR: fork");
          throw("PSIM ERROR: fork");
        }
        else if (pid == 0)
        { // Child Process
          _processID = i;
          break;
        }
      }
    }

    /** @brief Verifies that the destination process ID is valid. This
     * is done before data is sent or received.
     */
    void check_process_id(int processID)
    {

      if ((processID == _processID) ||
          (processID < PROCESS_PARENT) ||
          (processID >= _processCount))
      {

        char buffer[256];
        snprintf(buffer, 256, "PSIM ERROR: process %i referencing %i.",
                 _processID, processID);
        log(buffer);
        throw std::string(buffer);
      }
    }

    /** @brief This method initializes the process log and sets it up for
     * appending messages.
     */
    void open_log()
    {
      _doLogging = false;
      if (_logFileName.length() == 0)
      {
        return;
      }

      // open and reset file
      if ((_logfileFD = fopen(_logFileName.c_str(), "w")) == NULL)
      {
        log("PSIM ERROR: unable to create logfile");
        throw std::string("PSIM ERROR: unable to create logfile");
      }
      // close the log file
      close_log();

      // reopen the log file in append mode
      if ((_logfileFD = fopen(_logFileName.c_str(), "a")) == NULL)
      {
        log("PSIM ERROR: unable to open logfile");
        throw std::string("PSIM ERROR: unable to open logfile");
      }
      _doLogging = true;
    }

    /** @brief Closes the log file.
     */
    void close_log()
    {
      if (_doLogging)
        fclose(_logfileFD);
    }

    /** @brief Centralizes the repetitive task of logging the steps
     * during send and receive.
     */
    void logSendRecv(int sourcedestProcessID,
                     std::string tag,
                     enumSendRecv method,
                     enumSendRecvStep step)
    {

      char buffer[256];
      const char cmdSendRecv[2][8] = {"send", "recv"};
      const char stepSendRecv[3][12] = {"starting...", "failed!", "success."};

      snprintf(buffer, 256, "%i %s(%i,%s) %s",
               _processID, cmdSendRecv[method],
               sourcedestProcessID, tag.c_str(), stepSendRecv[step]);
      log(buffer);
    }

    /** @brief get_source_index
     */
    int get_source_index(int source)
    {
      check_process_id(source);
      return _processCount * source + _processID;
    }

    /** @brief detDestIndex
     */
    int get_dest_index(int dest)
    {
      check_process_id(dest);
      return _processCount * _processID + dest;
    }

    /** @brief Handles the sending of binary data.
     */
    void send_buffer(int destProcessID,
                     const void *pdataToSend, mdp_int dataSize)
    {
      static int counter = 0;
      char filename[512];
      counter++;
      int destIndex = get_dest_index(destProcessID);
      snprintf(filename, 512, ".fifo.%i.%i.%i", _processID, destProcessID, counter);
      int fd = open(filename, O_WRONLY | O_CREAT, 0700);
      if (write(fd, pdataToSend, dataSize) != dataSize)
      {
        log("PSIM ERROR: failure to write to socket");
        throw std::string("PSIM ERROR: failure to write to socket");
      }
      if (write(_socketFD[destIndex][COMM_SEND], &counter, sizeof(counter)) != sizeof(counter))
      {
        log("PSIM ERROR: failure to write to socket");
        throw std::string("PSIM ERROR: failure to write to socket");
      }
      close(fd);
    }

    /** @brief Sends a data tag and a vector of chars (as binary data).
     */
    void send_binary(int destProcessID,
                     const std::string &tag,
                     const std::vector<char> &data)
    {

      size_t tagSize = tag.size();
      size_t dataSize = data.size();
      send_buffer(destProcessID, &tagSize, sizeof(tagSize));
      send_buffer(destProcessID, tag.c_str(), tagSize);
      send_buffer(destProcessID, &dataSize, sizeof(dataSize));
      send_buffer(destProcessID, &data[0], dataSize);
    }

    /** @brief Handles the receiving of binary data through the sockets.
     */
    void recv_buffer(int sourceProcessID,
                     void *pdataToReceive, mdp_int dataSize)
    {
      char filename[512];
      int counter;
      int sourceIndex = get_source_index(sourceProcessID);
      if (read(_socketFD[sourceIndex][COMM_RECV], &counter, sizeof(counter)) != sizeof(counter))
      {
        log("PSIM ERROR: timeout error in reading from socket");
        throw std::string("PSIM ERROR: timeout error in reading from socket");
      }
      snprintf(filename, 512, ".fifo.%i.%i.%i", sourceProcessID, _processID, counter);
      int fd = open(filename, O_RDONLY);
      if (read(fd, (char *)pdataToReceive, dataSize) != dataSize)
      {
        log("PSIM ERROR: timeout error in reading from socket");
        throw std::string("PSIM ERROR: timeout error in reading from socket");
      }
      else
      {
        unlink(filename);
      }
      close(fd);
    }

    /** @brief Receives data utilizing a data tag to make sure that the
     * data coming in is what was expected.
     */
    void recv_binary(int sourceProcessID,
                     const std::string &tag,
                     std::vector<char> &data)
    {

      static std::vector<char> dataTmp;
      static std::string tagReceived;
      int size;

      while (true)
      {
        const auto &itr = _hash[sourceProcessID].find(tag);

        if (itr != _hash[sourceProcessID].end())
        { // Found?
          data = _hash[sourceProcessID][tag];
          _hash[sourceProcessID].erase(itr);
          break;
        }
        else
        {
          recv_buffer(sourceProcessID, &size, sizeof(size));
          std::vector<char> buffer(size + 1);
          recv_buffer(sourceProcessID, buffer.data(), size);
          buffer[size] = '\0';
          tagReceived = buffer.data();

          if (tagReceived == tag)
          {
            recv_buffer(sourceProcessID, &size, sizeof(size));
            data.resize(size);
            recv_buffer(sourceProcessID, &data[0], size);
            break;
          }
          else
          {
            recv_buffer(sourceProcessID, &size, sizeof(size));
            dataTmp.resize(size);
            recv_buffer(sourceProcessID, &dataTmp[0], size);
            _hash[sourceProcessID][tagReceived] = dataTmp;
          }
        }
      }
    }

    // *******************************************************************
    // ***         Private Method: doBegEndLog                         ***
    // ***                                                             ***
    // ***  Centralized log method used by some of the public methods  ***
    // ***  to send a standardized message to the log at the beginning ***
    // ***  and the end of the routine.                                ***
    // ***                                                             ***
    // *******************************************************************

    void doBegEndLog(std::string method, enumBegEnd begEnd)
    {
      char buffer[256];
      char *be;

      if (begEnd == LOG_BEGIN)
        be = (char *)"BEGIN";
      else
        be = (char *)"END";

      snprintf(buffer, 256, "%i %s [%s]", _processID, be, method.c_str());
      log(buffer);
    }

  public:
    /** @brief Provide the number of processes to create and the name of
     * the logfile if desired and "" if no logfile is needed.
     */
    mdp_psim(int processCount, std::string logFileName = ".psim.log", int verbatim = 0)
    {
      psim_begin(processCount, logFileName, verbatim);
    }

    mdp_psim(int argc, char **argv)
    {
      int processCount = parse_argv_nprocs(argc, argv);
      std::string logFileName = parse_argv_logfile(argc, argv);
      int verbatim = parse_argv_verbatim(argc, argv);
      psim_begin(processCount, logFileName, verbatim);
    }

    /** @brief Deallocates space that was created within the process,
     * releases sockets, closes the log, etc.
     */
    virtual ~mdp_psim()
    {
      psim_end();
    }

    // *******************************************************************
    // ***         Public Method: log                                  ***
    // ***                                                             ***
    // ***  Accepts a string and appends the message to the common     ***
    // ***  log file.  Note: locking is not necessary because of the   ***
    // ***  definition of append.  It does not matter how many        ***
    // ***  processes share file pointers, writing will always occur   ***
    // ***  at the end of the file.                                    ***
    // ***                                                             ***
    // *******************************************************************

    void log(std::string message, int level = 2)
    {
      if (_doLogging)
      {
        int fd = fileno(_logfileFD);
        flock(fd, LOCK_EX);
        fwrite("PSIM LOG: ", 10, 1, _logfileFD);
        fwrite(message.c_str(), message.length(), 1, _logfileFD);
        fwrite("\n", 1, 1, _logfileFD);
        // Clear out the file buffer
        fflush(_logfileFD);
        flock(fd, LOCK_UN);
      }
      if (_verbatim >= level)
      {
        std::cout << "PSIM LOG: " << message << std::endl;
        std::cout.flush();
      }
    }

    /** @brief Returns an integer identifying which process is currently
     * executing.
     */
    int id()
    {
      return _processID;
    }

    /** @brief Returns an integer identifying the current number of
     * active processes.
     */
    int nprocs()
    {
      return _processCount;
    }

    // *******************************************************************
    // ***         Public Method: setCommTimeout                       ***
    // ***                                                             ***
    // ***  Sets the number of seconds that a process will wait to     ***
    // ***  receive data from another process before throwing an       ***
    // ***  exception.                                                 ***
    // ***                                                             ***
    // *******************************************************************

    void setCommTimeout(mdp_uint commTimeout)
    {
      _commTimeout = commTimeout;
    }

    // *******************************************************************
    // ***         Public Method: send                                 ***
    // ***                                                             ***
    // ***  This asynchronous method sends the data referenced by      ***
    // ***  "dataToSend" to "destProcessID".  The size of the data     ***
    // ***  is obtained by looking at the type "T".                    ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    void send(int destProcessID, std::string dataTag, T &dataToSend)
    {
      logSendRecv(destProcessID, dataTag, LOG_SR_SEND, LOG_SR_START);
      std::vector<char> data(sizeof(T));
      for (unsigned int k = 0; k < sizeof(T); k++)
        data[k] = ((char *)&dataToSend)[k];
      send_binary(destProcessID, dataTag, data);
      // std::cout << _processID << "->" << destProcessID << " " << dataTag << " " << dataToSend << std::endl;
      logSendRecv(destProcessID, dataTag, LOG_SR_SEND, LOG_SR_SUCCESS);
    }

    // *******************************************************************
    // ***         Public Method: send                                 ***
    // ***                                                             ***
    // ***  This asynchronous method sends the data at location        ***
    // ***  "pdataToSend" to "destProcessID".  The size of the data    ***
    // ***  being sent is provided in the integer: "dataSize".         ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    void send(int destProcessID, std::string dataTag,
              T *pdataToSend, mdp_int dataSize)
    {
      logSendRecv(destProcessID, dataTag, LOG_SR_SEND, LOG_SR_START);
      std::vector<char> data(sizeof(T) * dataSize);
      for (size_t k = 0; k < data.size(); k++)
        data[k] = ((char *)pdataToSend)[k];
      send_binary(destProcessID, dataTag, data);
      logSendRecv(destProcessID, dataTag, LOG_SR_SEND, LOG_SR_SUCCESS);
    }

    // *******************************************************************
    // ***         Public Method: recv                                 ***
    // ***                                                             ***
    // ***  This synchronous "blocking" method retrieves the data sent ***
    // ***  to "destProcessID".  The size of the data being sent is    ***
    // ***  provided in the integer: "dataSize".                       ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    void recv(int sourceProcessID, std::string dataTag, T &dataToReceive)
    {
      logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_START);
      std::vector<char> data;
      recv_binary(sourceProcessID, dataTag, data);
      if (data.size() != sizeof(T))
      {
        log("PSIM ERROR: recv invalid data)");
        throw std::string("PSIM ERROR: recv invalid data)");
      };
      for (unsigned int k = 0; k < sizeof(T); k++)
        ((char *)&dataToReceive)[k] = data[k];
      // std::cout << _processID << "<-" << sourceProcessID << " " << dataTag << " " << dataToReceive << std::endl;
      logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_SUCCESS);
    }

    // *******************************************************************
    // ***         Public Method: recv                                 ***
    // ***                                                             ***
    // ***  This synchronous "blocking" method retrieves the data sent ***
    // ***  to "destProcessID".  "dataSize" bytes are copied to        ***
    // ***  location "pdataToReceive".                                 ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    void recv(int sourceProcessID, std::string dataTag,
              T *pdataToReceive, mdp_int dataSize)
    {
      logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_START);
      std::vector<char> data;
      recv_binary(sourceProcessID, dataTag, data);
      if (data.size() != sizeof(T) * dataSize)
      {
        log("PSIM ERROR: recv invalid data size");
        throw std::string("PSIM ERROR: recv invalid data size");
      }
      for (size_t k = 0; k < data.size(); k++)
        ((char *)pdataToReceive)[k] = data[k];
      logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_SUCCESS);
    }

    // *******************************************************************
    // ***         Public Method: broadcast                            ***
    // ***                                                             ***
    // ***  Allows broadcasting data to all processes.                 ***
    // ***  sourceProcessID sends data (data) to all of the other      ***
    // ***  processes who receive the data through the same variable). ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    void broadcast(int sourceProcessID, T &data)
    {
      static std::string tag = "BROADCAST:0";
      doBegEndLog(tag, LOG_BEGIN);
      if (_processID == sourceProcessID)
      {
        for (int i = 0; i < _processCount; i++)
        {
          if (i != sourceProcessID)
          {
            send(i, tag, data);
          }
        }
      }
      else
      {
        recv(sourceProcessID, tag, data);
      }
      if (tag == "BROADCAST:0")
        tag = "BROADCAST:1";
      else
        tag = "BROADCAST:0";
      doBegEndLog(tag, LOG_END);
    }

    template <class T>
    void broadcast(int sourceProcessID, T *data, int dataSize)
    {
      static std::string tag = "BROADCASTV:0";
      doBegEndLog(tag, LOG_BEGIN);
      if (_processID == sourceProcessID)
      {
        for (int i = 0; i < _processCount; i++)
        {
          if (i != sourceProcessID)
            send(i, tag, data, dataSize);
        }
      }
      else
      {
        recv(sourceProcessID, tag, data, dataSize);
      }
      if (tag == "BROADCASTV:0")
        tag = "BROADCASTV:1";
      else
        tag = "BROADCASTV:0";
      doBegEndLog(tag, LOG_END);
    }

    // *******************************************************************
    // ***         Public Method: collect                              ***
    // ***                                                             ***
    // *** All parallel processes construct a list of the data passed  ***
    // *** by each process.  The list is broadcasted and returned by   ***
    // *** each processor.  This method is used to implement global    ***
    // *** sum and some other global operations.                       ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    std::vector<T> collect(int dest, T &data)
    {
      static std::string tag = "COLLECT";
      std::vector<T> dataList;
      T dataToReceive;
      dataList.resize(_processCount);
      doBegEndLog(tag, LOG_BEGIN);

      if (_processID != dest)
      {
        send(dest, tag, data);
      }
      else
      {
        dataList[dest] = data;

        for (int i = 0; i < _processCount; i++)
        {
          if (i != dest)
          {
            recv(i, tag, dataToReceive);
            dataList[i] = dataToReceive;
          }
        }
      }
      if (tag == "COLLECT:0")
        tag = "COLLECT:1";
      else
        tag = "COLLECT:0";
      doBegEndLog(tag, LOG_END);
      return dataList;
    }

    // *******************************************************************
    // ***         Public Method: combine                              ***
    // ***                                                             ***
    // *** All parallel processes construct a list of the data passed  ***
    // *** by each process.  The list is broadcasted and returned by   ***
    // *** each processor.  This method is used to implement global    ***
    // *** sum and some other global operations.                       ***
    // ***                                                             ***
    // *******************************************************************

    template <class T>
    std::vector<T> combine(T &data)
    {
      std::vector<T> dataList = collect(PROCESS_PARENT, data);
      std::cout << id() << " size=" << dataList.size() << std::endl;

      broadcast(PROCESS_PARENT, &dataList[0], dataList.size());
      std::cout << id() << " list=" << dataList[0] << dataList[1] << dataList[2] << std::endl;
      return dataList;
    }

    /** @brief Initiates a blocking point so that the processes pause
     * until ALL processes have reached the barrier.
     */
    void barrier()
    {
      int dummy;
      broadcast(PROCESS_PARENT, dummy);
      collect(PROCESS_PARENT, dummy);
    }

    /** @brief All parallel processes sum their data in parallel.  The sum
     * is returned.
     */
    template <class T>
    T add(T &item)
    {
      T total = 0;
      std::vector<T> dataList = collect(PROCESS_PARENT, item);
      if (_processID == PROCESS_PARENT)
        for (size_t i = 0; i < dataList.size(); i++)
        {
          total += dataList[i];
        }
      broadcast(PROCESS_PARENT, total);
      return total;
    }

    // *******************************************************************
    // ***         Public Method: auxiliary functions                  ***
    // *** Examples:                                                   ***
    // *** a.out -PSIM_NPROCS=2           (parallel processes)         ***
    // *** a.out -PSIM_LOGFILE=./test.log (log into test.log)          ***
    // *** a.out -PSIM_VERBATIM=1         (show all communications)    ***
    // ***                                                             ***
    // *******************************************************************

    static int parse_argv_nprocs(int argc, char **argv)
    {
      int n = 1;
      for (int i = 1; i < argc; i++)
        if (strncmp(argv[i], "-PSIM_NPROCS=", 13) == 0)
        {
          sscanf(argv[i] + 13, "%i", &n);
          break;
        }
      return n;
    }

    static std::string parse_argv_logfile(int argc, char **argv)
    {
      for (int i = 1; i < argc; i++)
        if (strncmp(argv[i], "-PSIM_LOGFILE=", 14) == 0)
        {
          return std::string(argv[i] + 14);
        }
      return std::string("");
    }

    static int parse_argv_verbatim(int argc, char **argv)
    {
      int n = 1;
      for (int i = 1; i < argc; i++)
        if (strncmp(argv[i], "-PSIM_VERBATIM=", 15) == 0)
        {
          sscanf(argv[i] + 15, "%i", &n);
          break;
        }
      return n;
    }
  };
} // namespace MDP

#endif // NO_POSIX

#endif /* MDP_PSIM_ */
