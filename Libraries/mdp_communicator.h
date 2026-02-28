/////////////////////////////////////////////////////////////////
/// @file mdp_communicator.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_array
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_COMMUNICATOR_
#define MDP_COMMUNICATOR_

#include <ctime>
#include <memory>
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "mdp_log.h"
#include "mdp_psim.h"
#include "mdp_timer.h"
#include "mdp_matrix.h"

namespace MDP
{
#ifdef PARALLEL
  using mdp_request = MPI_Request;
#else // not PARALLEL
  using mdp_request = int;
  using MPI_Comm = int;
#endif

  /// @brief DO NOT INSTANTIATE use object mdp instead
  ///
  /// Example:
  /// @verbatim
  /// int main(int argc, char**argv) {
  ///    mdp.open_wormholes(argc,argv);
  ///    // your code here
  ///    mdp << 3.14 << "\n";  // only process 0 prints
  ///    mdp.close_wormholes();
  ///    return 0;
  /// }
  /// @endverbatim
  class mdp_communicator : public mdp_log
  {
  private:
#ifndef PARALLEL
#ifndef NO_POSIX
    std::unique_ptr<mdp_psim> nodes;
#endif
#else
    MPI_Comm communicator;
#endif
    double mytime; // total time
    bool wormholes_open;
    int my_id;
    int my_nproc;

#ifndef PARALLEL
    double MPI_Wtime()
    {
      return walltime();
    }
#endif

  public:
    double comm_time; // time spent in communications

    mdp_communicator()
    {
      wormholes_open = false;
    }

    virtual ~mdp_communicator()
    {
    }

    template <class T>
    void put(T &obj, int destination)
    {
#ifdef PARALLEL
      mdp_request r;
      MPI_Isend(&obj, sizeof(T) / sizeof(char), MPI_CHAR, destination,
                tag(me(), destination), communicator, &r);
      wait(r);
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->send(destination, "", obj);
#endif
#endif
    }

    template <class T>
    void put(T &obj, int destination, [[maybe_unused]] mdp_request &r)
    {
#ifdef PARALLEL
      MPI_Isend(&obj, sizeof(T) / sizeof(char), MPI_CHAR, destination,
                tag(me(), destination), communicator, &r);
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->send(destination, "", obj);
#endif
#endif
    }

    template <class T>
    void put(T *objptr, mdp_int length, int destination)
    {
#ifdef PARALLEL
      mdp_request r;
      MPI_Isend(objptr, length * sizeof(T) / sizeof(char), MPI_CHAR, destination,
                tag(me(), destination), communicator, &r);
      wait(r);
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->send(destination, "", objptr, length);
#endif
#endif
    }

    template <class T>
    void put(T *objptr, mdp_int length, int destination, [[maybe_unused]] mdp_request &r)
    {
#ifdef PARALLEL
      if (MPI_Isend(objptr, length * sizeof(T) / sizeof(char), MPI_CHAR, destination,
                    tag(me(), destination), communicator, &r) != MPI_SUCCESS)
        error("Failure to send");
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->send(destination, "", objptr, length);
#endif
#endif
    }

    template <class T>
    void get(T &obj, int source)
    {
#ifdef PARALLEL
      MPI_Status status;
      MPI_Recv(&obj, sizeof(T) / sizeof(char), MPI_CHAR, source,
               tag(source, me()), communicator, &status);
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->recv(source, "", obj);
#endif
#endif
    }

    template <class T>
    void get(T *objptr, mdp_int length, int source)
    {
#ifdef PARALLEL
      MPI_Status status;
      if (MPI_Recv(objptr, length * sizeof(T) / sizeof(char), MPI_CHAR, source,
                   tag(source, me()), communicator, &status) != MPI_SUCCESS)
        error("Failure to receive");
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->recv(source, "", objptr, length);
#endif
#endif
    }

    // ////////////////////////////////
    // Compatibility functions
    // obj1 is input, obj2 is output
#ifdef PARALLEL
    void add(float &obj1, float &obj2)
    {
      MPI_Allreduce(&obj1, &obj2, 1, MPI_FLOAT, MPI_SUM, communicator);
    }

    void add(float *obj1, float *obj2, mdp_int length)
    {
      MPI_Allreduce(obj1, obj2, length, MPI_FLOAT, MPI_SUM, communicator);
    }

    void add(double &obj1, double &obj2)
    {
      MPI_Allreduce(&obj1, &obj2, 1, MPI_DOUBLE, MPI_SUM, communicator);
    }
#else // not PARALLEL
    template <typename T>
    void add(T &obj1, T &obj2)
    {
#ifndef NO_POSIX
      obj2 = nodes->add(obj1);
#endif
    }
#endif

#ifdef PARALLEL
    void add(double *obj1, double *obj2, mdp_int length)
    {
      MPI_Allreduce(obj1, obj2, length, MPI_DOUBLE, MPI_SUM, communicator);
    }
#else // not PARALLEL
    template <typename T>
    void add(T *obj1, T *obj2, mdp_int length)
    {
#ifndef NO_POSIX
      for (mdp_int i = 0; i < length; i++)
      {
        obj2[i] = nodes->add(obj1[i]);
      }
#endif
    }
#endif
    // end compatibility functions

    // /////////////////////////////////

#ifdef PARALLEL
    void add(mdp_int &obj1)
    {
      mdp_int obj2 = 0;
      MPI_Allreduce(&obj1, &obj2, 1, MPI_LONG, MPI_SUM, communicator);
      obj1 = obj2;
    }

    void add(float &obj1)
    {
      float obj2 = 0;
      MPI_Allreduce(&obj1, &obj2, 1, MPI_FLOAT, MPI_SUM, communicator);
      obj1 = obj2;
    }

    void add(double &obj1)
    {
      double obj2 = 0;
      MPI_Allreduce(&obj1, &obj2, 1, MPI_DOUBLE, MPI_SUM, communicator);
      obj1 = obj2;
    }

    void add(mdp_complex &obj1)
    {
      mdp_complex obj2 = 0;
#ifndef USE_DOUBLE_PRECISION
      MPI_Allreduce(&obj1, &obj2, 2, MPI_FLOAT, MPI_SUM, communicator);
#else
      MPI_Allreduce(&obj1, &obj2, 2, MPI_DOUBLE, MPI_SUM, communicator);
#endif
      obj1 = obj2;
    }
#else // not PARALLEL
    template <typename T>
    void add(T &obj1)
    {
#ifndef NO_POSIX
      obj1 = nodes->add(obj1);
#endif
    }
#endif

    void add(mdp_matrix &a)
    {
#ifdef PARALLEL
      add(a.address(), a.rows() * a.cols());
#else
#ifndef NO_POSIX
      for (unsigned int i = 0; i < a.size(); i++)
        a.address()[i] = nodes->add(a.address()[i]);
#endif
#endif
    }

#ifdef PARALLEL
    void add(mdp_int *obj1, mdp_int length)
    {
      auto obj2 = std::make_unique<mdp_int[]>(length);
      for (mdp_int i = 0; i < length; i++)
        obj2[i] = 0;
      MPI_Allreduce(obj1, obj2.get(), length, MPI_LONG, MPI_SUM, communicator);
      for (mdp_int i = 0; i < length; i++)
        obj1[i] = obj2[i];
    }

    void add(float *obj1, mdp_int length)
    {
      auto obj2 = std::make_unique<float[]>(length);
      for (mdp_int i = 0; i < length; i++)
        obj2[i] = 0;
      MPI_Allreduce(obj1, obj2.get(), length, MPI_FLOAT, MPI_SUM, communicator);
      for (mdp_int i = 0; i < length; i++)
        obj1[i] = obj2[i];
    }

    void add(double *obj1, mdp_int length)
    {
      auto obj2 = std::make_unique<double[]>(length);
      for (mdp_int i = 0; i < length; i++)
        obj2[i] = 0;
      MPI_Allreduce(obj1, obj2.get(), length, MPI_DOUBLE, MPI_SUM, communicator);
      for (mdp_int i = 0; i < length; i++)
        obj1[i] = obj2[i];
    }

    void add(mdp_complex *obj1, mdp_int length)
    {
      auto obj2 = std::make_unique<mdp_complex[]>(length);
      for (mdp_int i = 0; i < length; i++)
        obj2[i] = 0;
#ifndef USE_DOUBLE_PRECISION
      MPI_Allreduce(obj1, obj2.get(), 2 * length, MPI_FLOAT, MPI_SUM, communicator);
#else
      MPI_Allreduce(obj1, obj2.get(), 2 * length, MPI_DOUBLE, MPI_SUM, communicator);
#endif
      for (mdp_int i = 0; i < length; i++)
        obj1[i] = obj2[i];
    }
#else // not PARALLEL
    template <typename T>
    void add(T *obj1, mdp_int length)
    {
#ifndef NO_POSIX
      for (mdp_int i = 0; i < length; i++)
        obj1[i] = nodes->add(obj1[i]);
#endif
    }
#endif

    void add(mdp_matrix *a, mdp_int length)
    {
#ifdef PARALLEL
      for (mdp_int i = 0; i < length; i++)
        add(a[i]);
#else // not PARALLEL
#ifndef NO_POSIX
      for (mdp_int j = 0; j < length; j++)
        for (mdp_uint i = 0; i < a[j].size(); i++)
          a[j].address()[i] = nodes->add(a[j].address()[i]);
#endif
#endif
    }

    template <class T>
    void add(std::vector<T> &a)
    {
      add(&a[0], a.size());
    }

    template <class T>
    void broadcast(T &obj, int p)
    {
#ifdef PARALLEL
      MPI_Bcast(&obj, sizeof(T) / sizeof(char), MPI_CHAR, p, communicator);
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->broadcast(p, obj);
#endif
#endif
    }

    template <class T>
    void broadcast(T *obj, mdp_int length, int p)
    {
#ifdef PARALLEL
      MPI_Bcast(obj, length * sizeof(T) / sizeof(char), MPI_CHAR, p, communicator);
#else // not PARALLEL
#ifndef NO_POSIX
      nodes->broadcast(p, obj, length);
#endif
#endif
    }

    void wait([[maybe_unused]] mdp_request &r)
    {
#ifdef PARALLEL
      MPI_Status status;
      MPI_Wait(&r, &status);
#endif
    }

    void wait([[maybe_unused]] mdp_request *r, [[maybe_unused]] int length)
    {
#ifdef PARALLEL
      MPI_Status status;
      MPI_Waitall(length, r, &status);
#endif
    }

    /** @brief Returns the unique id of this process
     */
    int me()
    {
      return my_id;
    }

    /** @brief Returns the total number of parallel processes for this job
     */
    int nproc()
    {
      return my_nproc;
    }

    void barrier()
    {
#ifdef PARALLEL
      MPI_Barrier(communicator);
#endif
    }

    int tag(int i, int j)
    {
      return i * nproc() + j;
    }

    void reset_time()
    {
      mytime = MPI_Wtime();
      comm_time = 0;
    }

    /** @brief returns the time in seconds since call to mdp_communicator::open_wormholes
     */
    double time()
    {
      return MPI_Wtime() - mytime;
    }

    /** @brief starts communications
     *
     * parses command line argument for MPI or PSIM parameters
     */
    void open_wormholes(int argc, char **argv
#ifdef PARALLEL
                        ,
                        MPI_Comm communicator_ = MPI_COMM_WORLD
#endif
    )
    {
      if (wormholes_open)
        return;

#ifdef PARALLEL
      communicator = communicator_;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(communicator, &(my_id));
      MPI_Comm_size(communicator, &(my_nproc));
      print = true;
      for (int i = 0; i < nproc(); i++)
      {
        if (me() == i)
          (*this) << "MPI PROCESS " << me()
                  << " of " << nproc()
                  << " STARTING ON " << getname() << "\n";
        barrier();
      }
#else // not PARALLEL
#ifndef NO_POSIX
      nodes = std::make_unique<mdp_psim>(argc, argv);
      my_id = nodes->id();
      my_nproc = nodes->nprocs();
#else
      my_id = 0;
      my_nproc = 1;
#endif
#endif
      setPrinting(me() == 0);

      begin_function("PROGRAM");
      begin_function("open_wormholes");
      (*this) << "<head>\n"
                 "*************************************************\n"
                 "* Starting [Matrix Distributed Processing]      *\n"
                 "* created by Massimo Di Pierro                  *\n"
                 "* copyrighted by www.metacryption.com           *\n"
                 "*************************************************\n"
                 "</head>\n";

      reset_time();
      double a, b;
      getcpuusage(a, b);
      wormholes_open = true;
      end_function("open_wormholes");
    }

    /** @brief prints statistics about parallel processes
     */
    void print_stats()
    {
#ifndef NO_POSIX
      auto a = std::make_unique<double[]>(nproc());
      auto b = std::make_unique<double[]>(nproc());
      auto c = std::make_unique<double[]>(nproc());
      for (int i = 0; i < nproc(); i++)
        a[i] = b[i] = c[i] = 0;
      getcpuusage(a[me()], b[me()]);
      c[me()] = 100.0 * comm_time / (time());
      add(a.get(), nproc());
      add(b.get(), nproc());
      add(c.get(), nproc());
      char buffer[256];
      for (int i = 0; i < nproc(); i++)
      {
        snprintf(buffer, 256,
                 "* Process %i stats: CPU=%.2f%% PROCESS=%.2f%% COMM=%.2f%%\n",
                 i, a[i], b[i], c[i]);
        (*this) << buffer;
      }
      (*this) << "* (above numbers make no sense under windows)\n";
#endif
    }

    /** @brief closes parallel communications
     */
    void close_wormholes()
    {
      begin_function("close_wormholes");
      (*this) << "<foot>\n"
                 "*************************************************\n"
                 "* Ending [Matrix Distributed Processing]        *\n";
      print_stats();
      (*this) << "*************************************************\n"
                 "</foot>\n";
      (*this) << "PROCESS " << me() << " ENDING AFTER " << time() << " sec.\n";
      wormholes_open = false;
#ifdef PARALLEL
      MPI_Finalize();
#endif
      end_function("close_wormholes");
      end_function("PROGRAM");
#ifndef PARALLEL
#ifndef NO_POSIX
      if (nodes)
        nodes.reset();
#endif
#endif
    }

    /** @brief forces the process to exit
     */
    void abort()
    {
#ifdef PARALLEL
      MPI_Abort(communicator, 1);
#else // not PARALLEL
      (*this) << "calling abort...";
      exit(-1);
#endif
    }
  };

  /// the only communicator object
  mdp_communicator mdp;

  void _mpi_error_message(const std::string &message, const std::string &file = "unknown", int line = 0)
  {
    mdp.error_message(message, file, line);
  }

  /// Logs in xml the start of a function with message fun
  void begin_function(const std::string &fun)
  {
    mdp.begin_function(fun);
  }

  /// Logs in xml the end of a function with message fun
  void end_function(const std::string &fun)
  {
    mdp.end_function(fun);
  }
} // namespace MDP

#endif /* MDP_COMMUNICATOR_ */
