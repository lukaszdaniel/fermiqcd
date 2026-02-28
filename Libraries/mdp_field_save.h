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

#include <cstdio>
#include "mdp_global_vars.h"
#include "mdp_utils.h"
#include "mdp_array.h"
#include "mdp_lattice.h"
#include "mdp_field.h"

namespace MDP
{
  /// Auxiliary function
  bool mdp_default_user_write(FILE *fp,
                              void *p,
                              mdp_int psize,
                              mdp_int header_size,
                              mdp_int position,
                              [[maybe_unused]] const mdp_lattice &lattice)
  {
    if (fseek(fp, position * psize + header_size, SEEK_SET) ||
        fwrite(p, psize, 1, fp) != 1)
      return false;
    return true;
  }

  bool mdp_default_user_write(std::ofstream &file,
                              void *p,
                              mdp_int psize,
                              mdp_int header_size,
                              mdp_int position,
                              [[maybe_unused]] const mdp_lattice &lattice)
  {
    // Offset in bytes
    std::streamoff offset = static_cast<std::streamoff>(position) * psize + header_size;

    // Move pointer to the correct position
    file.seekp(offset, std::ios::beg);
    if (!file)
      return false;

    // Write data
    file.write(reinterpret_cast<char *>(p), psize);
    if (!file)
      return false;

    return true;
  }

  /// Best way to save a field
  template <class T>
  bool mdp_field<T>::save_old(std::string filename,
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
    double mytime = mdp.time();

    m_header.reset();

    if (isSubProcess(processIO))
    {
      auto buffer_size = std::make_unique<mdp_int[]>(Nproc);
      auto buffer_ptr = std::make_unique<mdp_int[]>(Nproc);

      mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
      auto short_buffer = std::make_unique<T[]>(m_field_components);

      for (int process = 0; process < Nproc; process++)
      {
        buffer_size[process] = 0;
        buffer_ptr[process] = 0;
      }

      std::cout << "Saving file " << filename
                << " from process " << processIO
                << " (buffer = " << max_buffer_size << " sites)\n";

      FILE *fp = fopen(filename.c_str(), "wb+");
      if (fp == nullptr)
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
        int process = where_global(idx_gl);

        if ((process != NOWHERE) && (process != processIO))
        {
          if (buffer_ptr[process] == 0)
          {
            mdp.get(buffer_size[process], process);
            mdp.get(&(large_buffer(process, 0, 0)),
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
          T *src = m_data.get() + lattice().local(idx_gl) * m_field_components;

          for (mdp_uint k = 0; k < m_field_components; ++k)
            short_buffer[k] = src[k];
        }

        if (process != NOWHERE)
        {
          if (user_write)
          {
            if (!user_write(fp, short_buffer.get(),
                            m_field_components * sizeof(T),
                            skip_bytes,
                            idx_gl, lattice()))
              error("probably out of disk space");
          }
          else
          {
            if (exception && fseek(fp, idx_gl * psize + skip_bytes, SEEK_SET))
            {
              error("probably out of disk space");
            }

            if (fwrite(short_buffer.get(), psize, 1, fp) != 1)
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

      fclose(fp);
    }
    else
    {
      mdp_int buffer_size = 0;
      auto local_index = std::make_unique<mdp_int[]>(max_buffer_size);
      mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);
      mdp_request request;

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);

        if (process == ME)
          local_index[buffer_size++] = lattice().local(idx_gl);

        if ((buffer_size == max_buffer_size) ||
            ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
        {
          for (mdp_int i = 0; i < buffer_size; i++)
            for (mdp_uint k = 0; k < m_field_components; k++)
              local_buffer(i, k) = m_data.get()[local_index[i] * m_field_components + k];

          mdp.put(buffer_size, processIO, request);
          mdp.wait(request);

          mdp.put(&(local_buffer(0, 0)), buffer_size * m_field_components,
                  processIO, request);
          mdp.wait(request);

          buffer_size = 0;
        }
      }
    }

    if (isMainProcess() && !mdp_shutup)
    {
      std::cout << "... Saving time: " << (mdp.time() - mytime) << " (sec)\n";
    }

    return true;
  }

  template <class T>
  bool mdp_field<T>::save(std::string filename,
                          int processIO,
                          mdp_int max_buffer_size,
                          bool save_header,
                          mdp_int skip_bytes,
                          bool (*user_write)(std::ofstream &, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &))
  {
    filename = next_to_latest_file(filename);

    mdp_int header_size = 0;
    mdp_int psize = m_field_components * sizeof(T);
    mdp_int idx_gl, nvol_gl = lattice().global_volume();
    double mytime = mdp.time();

    m_header.reset();

    if (isSubProcess(processIO))
    {
      auto buffer_size = std::make_unique<mdp_int[]>(Nproc);
      auto buffer_ptr = std::make_unique<mdp_int[]>(Nproc);

      mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
      auto short_buffer = std::make_unique<T[]>(m_field_components);

      for (int process = 0; process < Nproc; process++)
      {
        buffer_size[process] = 0;
        buffer_ptr[process] = 0;
      }

      std::cout << "Saving file " << filename
                << " from process " << processIO
                << " (buffer = " << max_buffer_size << " sites)\n";

      std::ofstream fp(filename, std::ios::binary | std::ios::trunc);
      if (!fp)
        error("Unable to open file");

      m_header.set_time();

      if (save_header)
      {
        header_size = sizeof(mdp_field_file_header);

        fp.seekp(skip_bytes, std::ios::beg);
        fp.write(reinterpret_cast<const char *>(&m_header), header_size);

        if (!fp)
          error("Unable to write file header");
      }

      skip_bytes += header_size;
      bool exception = false;

      fp.seekp(skip_bytes, std::ios::beg);

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);

        if ((process != NOWHERE) && (process != processIO))
        {
          if (buffer_ptr[process] == 0)
          {
            mdp.get(buffer_size[process], process);
            mdp.get(&(large_buffer(process, 0, 0)),
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
          T *src = m_data.get() + lattice().local(idx_gl) * m_field_components;

          for (mdp_uint k = 0; k < m_field_components; ++k)
            short_buffer[k] = src[k];
        }

        if (process != NOWHERE)
        {
          if (user_write)
          {
            if (!user_write(fp, short_buffer.get(),
                            m_field_components * sizeof(T),
                            skip_bytes,
                            idx_gl, lattice()))
              error("probably out of disk space");
          }
          else
          {
            if (exception && fp.seekp(idx_gl * psize + skip_bytes, std::ios::beg))
            {
              error("probably out of disk space");
            }

            fp.write(reinterpret_cast<const char *>(short_buffer.get()), psize);

            if (!fp)
              error("probably out of disk space");
          }
        }
        else
        {
          exception = true;
        }
      }
    }
    else
    {
      mdp_int buffer_size = 0;
      auto local_index = std::make_unique<mdp_int[]>(max_buffer_size);
      mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);
      mdp_request request;

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);

        if (process == ME)
          local_index[buffer_size++] = lattice().local(idx_gl);

        if ((buffer_size == max_buffer_size) ||
            ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
        {
          for (mdp_int i = 0; i < buffer_size; i++)
            for (mdp_uint k = 0; k < m_field_components; k++)
              local_buffer(i, k) = m_data.get()[local_index[i] * m_field_components + k];

          mdp.put(buffer_size, processIO, request);
          mdp.wait(request);

          mdp.put(&(local_buffer(0, 0)), buffer_size * m_field_components,
                  processIO, request);
          mdp.wait(request);

          buffer_size = 0;
        }
      }
    }

    if (isMainProcess() && !mdp_shutup)
      std::cout << "... Saving time: " << (mdp.time() - mytime) << " (sec)\n";

    return true;
  }
} // namespace MDP

#endif /* MDP_FIELD_SAVE_ */
