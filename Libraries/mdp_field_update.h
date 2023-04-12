/////////////////////////////////////////////////////////////////
/// @file mdp_field_update.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains mdp_field::update()
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_FIELD_UPDATE_
#define MDP_FIELD_UPDATE_

namespace MDP
{
  /** @brief The only communication function for a field object
   *
   * To be invoked every time field variables are assigned and
   * need to be synchronized between the parallel processes
   */
  template <class T>
  void mdp_field<T>::update(int np, int d, int ncomp)
  {
    T *dynamic_buffer = nullptr;
    T *where_to = nullptr;
    mpi.comm_time -= mpi.time();
    mdp_request request;
    mdp_int start_to_send = 0;
    mdp_int process = 0;
    mdp_int length = 0;
    int ni, nf;

    if (d == -1)
    {
      d = 0;
      ncomp = m_field_components;
    }

    if ((ncomp == m_field_components) && (d != 0))
      error("update(): packet is too big");

    if (np < 2)
    {
      ni = nf = np;
    }
    else
    {
      ni = 0;
      nf = 1;
    }

    for (mdp_int dp = 1; dp < Nproc; dp++)
    {
      process = (ME + dp) % Nproc;

      if (np < 2)
      {
        length = lattice().len_to_send0(process, np);
      }
      else
      {
        length = lattice().len_to_send0(process, 0) + lattice().len_to_send0(process, 1);
      }

      if (np == 1)
      {
        start_to_send = lattice().len_to_send0(process, 0);
      }
      else
      {
        start_to_send = 0;
      }

      if (length > 0)
      {
        dynamic_buffer = new T[length * ncomp];
        for (mdp_int idx = 0; idx < length; idx++)
          for (mdp_int k = 0; k < ncomp; k++)
          {
            dynamic_buffer[idx * ncomp + k] =
                *(m_data.get() + lattice().to_send0(process, start_to_send + idx) * m_field_components + d * ncomp + k);
          }
        mpi.put(dynamic_buffer, length * ncomp, process, request);
        std::cout.flush();
      }
      else
      {
        dynamic_buffer = nullptr;
      }

      process = (ME - dp + Nproc) % Nproc;
      length = lattice().stop0(process, nf) - lattice().start0(process, ni);
      if (length > 0)
      {
        if (ncomp == m_field_components)
        {
          where_to = m_data.get() + lattice().start0(process, ni) * m_field_components;
          mpi.get(where_to, length * m_field_components, process);
          where_to = nullptr;
        }
        else
        {
          where_to = new T[length * ncomp];
          mpi.get(where_to, length * ncomp, process);
          for (mdp_int idx = 0; idx < length; idx++)
            for (mdp_int k = 0; k < ncomp; k++)
            {
              *(m_data.get() + (lattice().start0(process, ni) + idx) * m_field_components +
                d * ncomp + k) = where_to[idx * ncomp + k];
            }
          delete[] where_to;
        }
      }

      process = (ME + dp) % Nproc;

      if (dynamic_buffer != nullptr)
      {
        mpi.wait(request);
        delete[] dynamic_buffer;
      }
    }
    mpi.comm_time += mpi.time();
  }
} // namespace MDP

#endif /* MDP_FIELD_UPDATE_ */
