/////////////////////////////////////////////////////////////////
/// @file mdp_field_load.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains file IO operations for class mdp_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_FIELD_LOAD_
#define MDP_FIELD_LOAD_

#include <cstdio>
#include "mdp_global_vars.h"
#include "mdp_utils.h"
#include "mdp_array.h"
#include "mdp_lattice.h"
#include "mdp_field.h"

namespace MDP
{
  /// Auxiliary function
  bool mdp_default_user_read(FILE *fp,
                             void *p,
                             mdp_int psize,
                             mdp_int header_size,
                             mdp_int position,
                             [[maybe_unused]] const mdp_lattice &lattice)
  {
    if (fseek(fp, (size_t)position * psize + header_size, SEEK_SET) ||
        fread(p, psize, 1, fp) != 1)
      return false;
    return true;
  }

  bool mdp_default_user_read(std::ifstream &file,
                             void *p,
                             mdp_int psize,
                             mdp_int header_size,
                             mdp_int position,
                             [[maybe_unused]] const mdp_lattice &lattice)
  {
    // Offset in bytes
    std::streamoff offset = static_cast<std::streamoff>(position) * psize + header_size;

    // Move pointer to the correct position
    file.seekg(offset, std::ios::beg);
    if (!file)
      return false;

    // Read data
    file.read(reinterpret_cast<char *>(p), psize);
    if (!file)
      return false;

    return true;
  }

  /// Best way to load a field
  template <class T>
  bool mdp_field<T>::load_old(std::string filename,
                          int processIO,
                          mdp_int max_buffer_size,
                          bool load_header,
                          mdp_int skip_bytes,
                          bool (*user_read)(FILE *, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &),
                          bool try_switch_endianess)
  {
    if (!file_exists(filename))
      throw std::string("file ") + filename + std::string(" does not exist");

    filename = latest_file(filename);
    if (filename == "?")
      return false;

    mdp_int header_size = 0;
    size_t idx_gl, nvol_gl = lattice().global_volume();
    size_t psize = m_field_components * sizeof(T);
    double mytime = mdp.time();
    bool reversed_header_endianess = false;

    struct stat statbuf;

    if (isSubProcess(processIO))
    {
      auto buffer_size = std::make_unique<mdp_int[]>(Nproc);
      mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
      auto short_buffer = std::make_unique<T[]>(m_field_components);
      mdp_request request;

      for (int process = 0; process < Nproc; process++)
        buffer_size[process] = 0;
      std::cout << "Loading file " << filename
                << " from process " << processIO
                << " (buffer = " << max_buffer_size << " sites)\n";
      fflush(stdout);
      stat(filename.c_str(), &statbuf);
      int total_size = statbuf.st_size;

      FILE *fp = fopen(filename.c_str(), "rb");
      if (fp == nullptr)
        error("Unable to open file");

      if (load_header)
      {
        mdp_field_file_header tmp_header;
        header_size = sizeof(mdp_field_file_header);

        if (fseek(fp, skip_bytes, SEEK_SET) ||
            fread(&tmp_header, header_size, 1, fp) != 1)
        {
          fprintf(stderr, "mdp_field.load(): Unable to load file header\n");
          return false;
        }

        reversed_header_endianess = mdp_field_file_header::switch_header_endianess(tmp_header);

        std::cout << "reverse: " << reversed_header_endianess << std::endl;

        if (tmp_header.endianess != m_header.endianess)
          fprintf(stderr, "Unrecognized endianess... trying to read anyway\n");

        // UGLY BUT FIXES INCOMPATIBIITY
        int actual_size = tmp_header.box[0];
        for (mdp_int d = 1; d < tmp_header.ndim; d++)
          actual_size *= tmp_header.box[d];

        tmp_header.sites = actual_size;
        header_size += total_size - (tmp_header.bytes_per_site * actual_size + header_size);

        if (tmp_header.ndim != m_header.ndim)
        {
          fprintf(stderr, "mdp_field.load(): wrong ndim\n");
          return false;
        }

        for (mdp_int d = 0; d < lattice().ndim(); d++)
          if (tmp_header.box[d] != m_header.box[d])
          {
            fprintf(stderr, "mdp_file.load(): wrong lattice size\n");
            return false;
          }
        if (tmp_header.bytes_per_site != m_header.bytes_per_site)
        {
          fprintf(stderr, "mdp_file.load(): wrong type of field (%i bytes per site?)\n", tmp_header.bytes_per_site);
          return false;
        }
        if (tmp_header.sites != m_header.sites)
        {
          fprintf(stderr, "mdp_field.load(): wrong number of sites\n");
          return false;
        }
        m_header = tmp_header;
      }

      skip_bytes += header_size;

      bool exception = false;
      fseek(fp, skip_bytes, SEEK_SET);

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);

        if (process != NOWHERE)
        {
          if (user_read)
          {
            if (!user_read(fp, short_buffer.get(),
                           m_field_components * sizeof(T),
                           skip_bytes,
                           idx_gl, lattice()))
              error("unexpected end of file");
          }
          else
          {
            if (exception && fseek(fp, idx_gl * psize + skip_bytes, SEEK_SET))
            {
              std::cout << "debug info: " << idx_gl * psize + skip_bytes << " " << psize << std::endl;
              error("unexpected end of file");
            }

            if (fread(short_buffer.get(), psize, 1, fp) != 1)
            {
              std::cout << "failure to read" << std::endl;
              error("unexpected end of file");
            }
          }
        }
        else
        {
          exception = true;
        }

        if ((process != NOWHERE) && (process != processIO))
        {
          for (mdp_uint k = 0; k < m_field_components; k++)
            large_buffer(process, buffer_size[process], k) = short_buffer[k];

          buffer_size[process]++;

          if (buffer_size[process] == max_buffer_size)
          {
            mdp.put(&(large_buffer(process, 0, 0)),
                    max_buffer_size * m_field_components, process, request);
            mdp.wait(request);
            buffer_size[process] = 0;
          }

          if (idx_gl == nvol_gl - 1)
          {
            for (int process = 0; process < Nproc; process++)
              if ((process != ME) &&
                  (buffer_size[process] != max_buffer_size) &&
                  (buffer_size[process] > 0))
              {
                mdp.put(&(large_buffer(process, 0, 0)),
                        buffer_size[process] * m_field_components,
                        process, request);
                mdp.wait(request);
              }
          }
        }

        if (process == processIO)
        {
          T *dst = m_data.get() + lattice().local(idx_gl) * m_field_components;

          for (mdp_uint k = 0; k < m_field_components; k++)
            dst[k] = short_buffer[k];
        }
      }

      fclose(fp);
    }
    else
    {
      mdp_int buffer_size = 0;
      auto local_index = std::make_unique<mdp_int[]>(max_buffer_size);
      mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);
        if (process == ME)
        {
          local_index[buffer_size++] = lattice().local(idx_gl);
        }

        if ((buffer_size == max_buffer_size) ||
            ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
        {
          mdp.get(&(local_buffer(0, 0)), buffer_size * m_field_components, processIO);

          for (mdp_int i = 0; i < buffer_size; i++)
            for (mdp_uint k = 0; k < m_field_components; k++)
              m_data.get()[local_index[i] * m_field_components + k] = local_buffer(i, k);

          buffer_size = 0;
        }
      }
    }

    update();
    mdp.broadcast(reversed_header_endianess, processIO);

    if (try_switch_endianess && reversed_header_endianess)
    {
      mdp << "switching endianess...\n";
#ifdef USE_DOUBLE_PRECISION
      switch_endianess_8bytes();
#else
      switch_endianess_4bytes();
#endif
    }

    if (isMainProcess() && !mdp_shutup)
    {
      printf("... Loading time: %f (sec)\n", mdp.time() - mytime);
      fflush(stdout);
    }
    return true;
  }

  template <class T>
  bool mdp_field<T>::load(std::string filename,
                          int processIO,
                          mdp_int max_buffer_size,
                          bool load_header,
                          mdp_int skip_bytes,
                          bool (*user_read)(std::ifstream &, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &),
                          bool try_switch_endianess)
  {
    if (!file_exists(filename))
      throw std::string("file ") + filename + std::string(" does not exist");

    filename = latest_file(filename);
    if (filename == "?")
      return false;

    mdp_int header_size = 0;
    size_t idx_gl, nvol_gl = lattice().global_volume();
    size_t psize = m_field_components * sizeof(T);
    double mytime = mdp.time();
    bool reversed_header_endianess = false;

    struct stat statbuf;

    if (isSubProcess(processIO))
    {
      auto buffer_size = std::make_unique<mdp_int[]>(Nproc);
      mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
      auto short_buffer = std::make_unique<T[]>(m_field_components);
      mdp_request request;

      for (int process = 0; process < Nproc; process++)
        buffer_size[process] = 0;
      std::cout << "Loading file " << filename
                << " from process " << processIO
                << " (buffer = " << max_buffer_size << " sites)\n";

      stat(filename.c_str(), &statbuf);
      int total_size = statbuf.st_size;

      std::ifstream fp(filename, std::ios::binary);
      if (!fp)
        error("Unable to open file");

      if (load_header)
      {
        mdp_field_file_header tmp_header;
        header_size = sizeof(mdp_field_file_header);

        fp.seekg(skip_bytes, std::ios::beg);
        if (!fp.read(reinterpret_cast<char *>(&tmp_header), header_size))
        {
          fprintf(stderr, "mdp_field.load(): Unable to load file header\n");
          return false;
        }

        reversed_header_endianess = mdp_field_file_header::switch_header_endianess(tmp_header);

        std::cout << "reverse: " << reversed_header_endianess << std::endl;

        if (tmp_header.endianess != m_header.endianess)
          fprintf(stderr, "Unrecognized endianess... trying to read anyway\n");

        // UGLY BUT FIXES INCOMPATIBIITY
        int actual_size = tmp_header.box[0];
        for (mdp_int d = 1; d < tmp_header.ndim; d++)
          actual_size *= tmp_header.box[d];

        tmp_header.sites = actual_size;
        header_size += total_size - (tmp_header.bytes_per_site * actual_size + header_size);

        if (tmp_header.ndim != m_header.ndim)
        {
          fprintf(stderr, "mdp_field.load(): wrong ndim\n");
          return false;
        }

        for (mdp_int d = 0; d < lattice().ndim(); d++)
          if (tmp_header.box[d] != m_header.box[d])
          {
            fprintf(stderr, "mdp_file.load(): wrong lattice size\n");
            return false;
          }
        if (tmp_header.bytes_per_site != m_header.bytes_per_site)
        {
          fprintf(stderr, "mdp_file.load(): wrong type of field (%i bytes per site?)\n", tmp_header.bytes_per_site);
          return false;
        }
        if (tmp_header.sites != m_header.sites)
        {
          fprintf(stderr, "mdp_field.load(): wrong number of sites\n");
          return false;
        }
        m_header = tmp_header;
      }

      skip_bytes += header_size;

      bool exception = false;
      fp.seekg(skip_bytes, std::ios::beg);

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);

        if (process != NOWHERE)
        {
          if (user_read)
          {
            if (!user_read(fp, short_buffer.get(),
                           m_field_components * sizeof(T),
                           skip_bytes,
                           idx_gl, lattice()))
              error("unexpected end of file");
          }
          else
          {
            if (exception && fp.seekg(idx_gl * psize + skip_bytes, std::ios::beg))
            {
              std::cout << "debug info: " << idx_gl * psize + skip_bytes << " " << psize << std::endl;
              error("unexpected end of file");
            }

            if (!fp.read(reinterpret_cast<char *>(short_buffer.get()), psize))
            {
              std::cout << "failure to read" << std::endl;
              error("unexpected end of file");
            }
          }
        }
        else
        {
          exception = true;
        }

        if ((process != NOWHERE) && (process != processIO))
        {
          for (mdp_uint k = 0; k < m_field_components; k++)
            large_buffer(process, buffer_size[process], k) = short_buffer[k];

          buffer_size[process]++;

          if (buffer_size[process] == max_buffer_size)
          {
            mdp.put(&(large_buffer(process, 0, 0)),
                    max_buffer_size * m_field_components, process, request);
            mdp.wait(request);
            buffer_size[process] = 0;
          }

          if (idx_gl == nvol_gl - 1)
          {
            for (int process = 0; process < Nproc; process++)
              if ((process != ME) &&
                  (buffer_size[process] != max_buffer_size) &&
                  (buffer_size[process] > 0))
              {
                mdp.put(&(large_buffer(process, 0, 0)),
                        buffer_size[process] * m_field_components,
                        process, request);
                mdp.wait(request);
              }
          }
        }

        if (process == processIO)
        {
          T *dst = m_data.get() + lattice().local(idx_gl) * m_field_components;

          for (mdp_uint k = 0; k < m_field_components; k++)
            dst[k] = short_buffer[k];
        }
      }
    }
    else
    {
      mdp_int buffer_size = 0;
      auto local_index = std::make_unique<mdp_int[]>(max_buffer_size);
      mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);

      for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        int process = where_global(idx_gl);
        if (process == ME)
        {
          local_index[buffer_size++] = lattice().local(idx_gl);
        }

        if ((buffer_size == max_buffer_size) ||
            ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
        {
          mdp.get(&(local_buffer(0, 0)), buffer_size * m_field_components, processIO);

          for (mdp_int i = 0; i < buffer_size; i++)
            for (mdp_uint k = 0; k < m_field_components; k++)
              m_data.get()[local_index[i] * m_field_components + k] = local_buffer(i, k);

          buffer_size = 0;
        }
      }
    }

    update();
    mdp.broadcast(reversed_header_endianess, processIO);

    if (try_switch_endianess && reversed_header_endianess)
    {
      mdp << "switching endianess...\n";
#ifdef USE_DOUBLE_PRECISION
      switch_endianess_8bytes();
#else
      switch_endianess_4bytes();
#endif
    }

    if (isMainProcess() && !mdp_shutup)
      printf("... Loading time: %f (sec)\n", mdp.time() - mytime);

    return true;
  }
} // namespace MDP

#endif /* MDP_FIELD_LOAD_ */
