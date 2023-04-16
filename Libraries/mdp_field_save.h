/////////////////////////////////////////////////////////////////
/// @file mdp_field_save.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains file IO operations for class mdp_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_FIELD_SAVE_
#define MDP_FIELD_SAVE_

namespace MDP
{
  /// Auxiliary function
  bool mdp_default_user_write(FILE *fp,
                              void *p,
                              mdp_int psize,
                              mdp_int header_size,
                              mdp_int position,
                              const mdp_lattice &lattice)
  {
    if (fseek(fp, position * psize + header_size, SEEK_SET) ||
        fwrite(p, psize, 1, fp) != 1)
      return false;
    return true;
  }

  /// Best way to save a field
  template <class T>
  bool mdp_field<T>::save(std::string filename,
                          int processIO,
                          mdp_int max_buffer_size,
                          bool save_header,
                          mdp_int skip_bytes,
                          bool (*user_write)(FILE *, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &))
  {

    filename = next_to_latest_file(filename);

    mdp_int header_size = 0;
    mdp_int psize = m_field_components * sizeof(T);
    mdp_int idx_gl, nvol_gl = lattice().global_volume();
    double mytime = mpi.time();
    m_header.reset();
    if (ME == processIO)
    {
      mdp_int *buffer_size = new mdp_int[Nproc];
      mdp_int *buffer_ptr = new mdp_int[Nproc];
      mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
      T *short_buffer = new T[m_field_components];
      int process;
      for (process = 0; process < Nproc; process++)
        buffer_ptr[process] = 0;
      std::cout << "Saving file " << filename
                << " from process " << processIO
                << " (buffer = " << max_buffer_size << " sites)\n";
      fflush(stdout);
      FILE *fp = fopen(filename.c_str(), "wb+");
      if (fp == 0)
        error("Unable to open file");

      m_header.set_time();

      if (save_header)
      {
        header_size = sizeof(mdp_field_file_header);
        if (fseek(fp, skip_bytes, SEEK_SET) ||
            fwrite(&m_header, header_size, 1, fp) != 1)
          error("Unable to write file header");
      }

      skip_bytes += header_size;
      bool exception = false;
      fseek(fp, skip_bytes, SEEK_SET);
      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        process = where_global(idx_gl);
        if ((process != NOWHERE) && (process != processIO))
        {
          if (buffer_ptr[process] == 0)
          {
            mpi.get(buffer_size[process], process);
            mpi.get(&(large_buffer(process, 0, 0)),
                    buffer_size[process] * m_field_components, process);
          }
          for (mdp_uint k = 0; k < m_field_components; k++)
            short_buffer[k] = large_buffer(process, buffer_ptr[process], k);
          buffer_ptr[process]++;
          if (buffer_ptr[process] == buffer_size[process])
            buffer_ptr[process] = 0;
        }
        if (process == processIO)
        {
          for (mdp_uint k = 0; k < m_field_components; k++)
            short_buffer[k] = *(m_data.get() + lattice().local(idx_gl) * m_field_components + k);
        }
        if (process != NOWHERE)
        {
          if (user_write)
          {
            if (!user_write(fp, short_buffer,
                            m_field_components * sizeof(T),
                            skip_bytes,
                            idx_gl, lattice()))
              error("propably out ofdisk space");
          }
          else
          {
            if (exception && fseek(fp, idx_gl * psize + skip_bytes, SEEK_SET))
            {
              error("probably out of disk space");
            }
            if (fwrite(short_buffer, psize, 1, fp) != 1)
            {
              error("probably out of disk space");
            }
          }
        }
        else
        {
          exception = true;
        }
      }
      delete[] buffer_size;
      delete[] buffer_ptr;
      delete[] short_buffer;
      fclose(fp);
    }
    else
    {
      int process;
      mdp_int buffer_size = 0, idx, idx_gl;
      mdp_int *local_index = new mdp_int[max_buffer_size];
      mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);
      mdp_request request;
      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        process = where_global(idx_gl);
        if (process == ME)
        {
          local_index[buffer_size] = lattice().local(idx_gl);
          buffer_size++;
        }
        if ((buffer_size == max_buffer_size) ||
            ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
        {
          for (idx = 0; idx < buffer_size; idx++)
            for (mdp_uint k = 0; k < m_field_components; k++)
              local_buffer(idx, k) = *(m_data.get() + local_index[idx] * m_field_components + k);
          mpi.put(buffer_size, processIO, request);
          mpi.wait(request);
          mpi.put(&(local_buffer(0, 0)), buffer_size * m_field_components,
                  processIO, request);
          mpi.wait(request);
          buffer_size = 0;
        }
      }
      delete[] local_index;
    }
    if (ME == 0 && mdp_shutup == false)
    {
      printf("... Saving time: %f (sec)\n", mpi.time() - mytime);
      fflush(stdout);
    }
    return true;
  }
} // namespace MDP

#endif /* MDP_FIELD_SAVE_ */
