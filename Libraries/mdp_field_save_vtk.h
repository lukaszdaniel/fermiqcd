/////////////////////////////////////////////////////////////////
/// @file mdp_field_save_vtk.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains file IO operations for class mdp_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_FIELD_SAVE_VTK_
#define MDP_FIELD_SAVE_VTK_

#include "mdp_field.h"
#include "mdp_utils.h"
#include "mdp_array.h"

namespace MDP
{
  template <class T>
  bool mdp_field<T>::save_vtk(std::string filename,
                              int t,
                              int component,
                              int processIO,
                              bool ASCII)
  {

    filename = next_to_latest_file(filename);
    std::string filename_tmp = filename + ".tmp";

    int max_buffer_size = 1024;
    int timeslice;
    mdp_int nvol_gl = lattice().global_volume();
    double mytime = mdp.time();
    m_header.reset();
    if (lattice().n_dimensions() < 3 || lattice().n_dimensions() > 4)
      error("mdp_field::save_vtk only works for ndim=3 and 4");
    if (isSubProcess(processIO))
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
      FILE *fp = nullptr;

      float fval;
      char tmp[1024];
      char header[1024];

      int space_volume = lattice().size() / lattice().size(0);
      int LZ = lattice().size(1);
      int LY = lattice().size(2);
      int LX = lattice().size(3);
      if (lattice().n_dimensions() == 3)
      {
        space_volume = lattice().size();
        LZ = lattice().size(0);
        LY = lattice().size(1);
        LX = lattice().size(2);
        t = -1;
      }

      snprintf(header, 1024,
               "# vtk DataFile Version 2.0\n"
               "%s\n"
               "%s\n"
               "DATASET STRUCTURED_POINTS\n"
               "DIMENSIONS %i %i %i\n"
               "ORIGIN     0   0   0\n"
               "SPACING    1   1   1\n"
               "POINT_DATA %i",
               filename.c_str(),
               ((ASCII) ? "ASCII" : "BINARY"),
               LX, LY, LZ, LX * LY * LZ);

      fp = fopen(filename_tmp.c_str(), "wb");
      fwrite(header, sizeof(char), strlen(header), fp);

      for (mdp_int idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
      {
        process = where_global(idx_gl);
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
          for (mdp_uint k = 0; k < m_field_components; k++)
            short_buffer[k] = *(m_data.get() + lattice().local(idx_gl) * m_field_components + k);
        }
        if (process != NOWHERE)
        {
          timeslice = idx_gl / space_volume;
          if (idx_gl % space_volume == 0 && (t < 0 || timeslice == t))
          {
            snprintf(header, 1024, "\nSCALARS scalars_t%i float\nLOOKUP_TABLE default\n", timeslice);
            fwrite(header, sizeof(char), strlen(header), fp);
          }
          if (t < 0 || timeslice == t || lattice().n_dimensions() == 3)
          {
            for (mdp_uint fc = 0; fc < m_field_components; fc++)
              if (component == -1 || fc == mdp_uint(component))
              {
                fval = (float)short_buffer[fc];
                if (!ASCII)
                {
                  switch_endianess_byte4(fval);
                  if (fwrite(&fval, sizeof(fval), 1, fp) != 1)
                    error("probably out of disk space");
                }
                else
                {
                  snprintf(tmp, 1024, "%e\n", fval);
                  if (fwrite(tmp, strlen(tmp), 1, fp) != 1)
                    error("probably out of disk space");
                }
              }
          }
        }
      }
      if (fp)
      {
        fclose(fp);
        remove(filename.c_str());
        rename(filename_tmp.c_str(), filename.c_str());
      }
      delete[] buffer_size;
      delete[] buffer_ptr;
      delete[] short_buffer;
    }
    else
    {
      int process;
      mdp_int buffer_size = 0;
      mdp_int *local_index = new mdp_int[max_buffer_size];
      mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);
      mdp_request request;
      for (mdp_int idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
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
          for (mdp_int idx = 0; idx < buffer_size; idx++)
            for (mdp_uint k = 0; k < m_field_components; k++)
              local_buffer(idx, k) = *(m_data.get() + local_index[idx] * m_field_components + k);
          mdp.put(buffer_size, processIO, request);
          mdp.wait(request);
          mdp.put(&(local_buffer(0, 0)), buffer_size * m_field_components,
                  processIO, request);
          mdp.wait(request);
          buffer_size = 0;
        }
      }
      delete[] local_index;
    }
    if (isMainProcess() && !mdp_shutup)
    {
      printf("... Saving time: %f (sec)\n", mdp.time() - mytime);
      fflush(stdout);
    }
    return true;
  }

  mdp_field<mdp_real> &cumulate_field(mdp_field<mdp_real> &field, std::string filename)
  {
    static std::map<std::string, int> counter;
    static std::map<std::string, mdp_field<mdp_real> *> fields;
    mdp_site p(field.lattice());

    if (counter.find(filename) == counter.end())
    {
      fields[filename] = new mdp_field<mdp_real>(field.lattice(), field.size_per_site());
      forallsites(p)
      {
        for (int i = 0; i < field.size_per_site(); i++)
          (*fields[filename])(p, i) = field(p, i);
      }
      counter[filename] = 1;
    }
    else
    {
      int k = counter[filename];
      forallsites(p)
      {
        for (int i = 0; i < field.size_per_site(); i++)
          (*fields[filename])(p, i) = (k * (*fields[filename])(p, i) + field(p, i)) / (k + 1);
      }
      counter[filename]++;
    }
    return (*fields[filename]);
  }

#if 0
  void save_vtk(mdp_field<mdp_real> &field,
                std::string filename,
                int t = -1,
                int component = -1,
                int processIO = 0,
                bool ASCII = false)
  {
    static std::map<std::string, int> counter;
    static std::map<std::string, mdp_field<mdp_real> *> fields;
    mdp_site p(field.lattice());
    int k = 0;
    if (counter.find(filename) == counter.end())
    {
      fields[filename] = new mdp_field<mdp_real>(field.lattice(), field.size_per_site());
      forallsites(p)
      {
        for (int i = 0; i < field.size_per_site(); i++)
          (*fields[filename])(p, i) = field(p, i);
      }
      counter[filename] = 1;
    }
    else
    {
      k = counter[filename];
      forallsites(p)
      {
        for (int i = 0; i < field.size_per_site(); i++)
          (*fields[filename])(p, i) = (k * (*fields[filename])(p, i) + field(p, i)) / (k + 1);
      }
      counter[filename]++;
    }
    if (filename[filename.size() - 1] != '*')
    {
      filename = filename.substr(0, filename.size() - 4);
      field.save_vtk(filename + ".vtk", t, component, processIO, ASCII);
      fields[filename]->save_vtk(filename + "_average.vtk", t, component, processIO, ASCII);
    }
    else
    {
      std::string filename2 = filename.substr(0, filename.size() - 1);
      field.save_vtk(filename2 + "_" + std::to_string(k) + ".vtk", t, component, processIO, ASCII);
      fields[filename]->save_vtk(filename2 + "_average_" + std::to_string(k) + ".vtk", t, component, processIO, ASCII);
    }
  }
#endif
} // namespace MDP

#endif /* MDP_FIELD_SAVE_VTK_ */
