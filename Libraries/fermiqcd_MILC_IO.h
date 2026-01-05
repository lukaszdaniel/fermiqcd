/////////////////////////////////////////////////////////////////
/// @file fermiqcd_MILC_IO.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Functions to read a MILC gauge configuration without conversion
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_MILC_IO_
#define FERMIQCD_MILC_IO_

#include <memory>
#include <string>
#include <cstdio>
#include <fstream>
#include "mdp_global_vars.h"
#include "mdp_site.h"
#include "fermiqcd_gauge_field.h"

namespace MDP
{
  struct NoEndianSwitch
  {
    static inline void apply(float &) noexcept {}
  };

  struct EndianSwitch
  {
    static inline void apply(float &v) noexcept
    {
      switch_endianess_byte4(v);
    }
  };

  template <typename EndianPolicy>
  bool milc_read_as_float(std::ifstream &fp,
                          void *data,
                          mdp_int psize,
                          mdp_int header_size,
                          mdp_int position,
                          const mdp_lattice &lattice)
  {
    auto *out = static_cast<double *>(data);

    mdp_site x(lattice);
    x.set_global(position);
    position = (((x(0) * lattice.size(3) + x(3)) * lattice.size(2) + x(2)) * lattice.size(1) + x(1));

    constexpr std::size_t float_size = sizeof(float);
#ifdef USE_DOUBLE_PRECISION
    const std::size_t bytes_to_read = psize / 2;
#else
    const std::size_t bytes_to_read = psize;
#endif
    const std::size_t float_count = bytes_to_read / float_size;

    // zero-copy: we read directly into the float buffer
    auto *buffer = static_cast<float *>(data);

    fp.seekg(static_cast<std::streamoff>(position * bytes_to_read + header_size), std::ios::beg);
    if (!fp.read(reinterpret_cast<char *>(buffer), static_cast<std::streamsize>(bytes_to_read)))
    {
      return false;
    }

    // in-place conversion from end (to not overwrite data)
    for (std::size_t i = float_count; i-- > 0;)
    {
      EndianPolicy::apply(buffer[i]);
      out[i] = static_cast<double>(buffer[i]);
    }

    return true;
  }

  inline bool milc_read_as_float_noswitch(std::ifstream &fp,
                                          void *data,
                                          mdp_int psize,
                                          mdp_int header_size,
                                          mdp_int position,
                                          const mdp_lattice &lattice)
  {
    return milc_read_as_float<NoEndianSwitch>(fp, data, psize, header_size, position, lattice);
  }

  inline bool milc_read_as_float_switch(std::ifstream &fp,
                                        void *data,
                                        mdp_int psize,
                                        mdp_int header_size,
                                        mdp_int position,
                                        const mdp_lattice &lattice)
  {
    return milc_read_as_float<EndianSwitch>(fp, data, psize, header_size, position, lattice);
  }

  bool load_milc(gauge_field &U, const std::string &filename,
                 mdp_int max_buffer_size = 128, int processIO = 0)
  {
    struct
    {
      mdp_uint magic_number; /* Identifies file format */
      mdp_int dims[4];       /* Full lattice dimensions */
      char time_stamp[64];   /* Date and time stamp - used to
  check consistency between the
  ASCII header file and the
  lattice file */
      mdp_int header_bytes;  /* NOT WRITTEN TO THE FILE but
  helpful for finding the data */
      mdp_int what_is_this;
      mdp_int order; /* 0 means no coordinate list is
  attached and the values are in
  coordinate serial order.
  Nonzero means that a
  coordinate list is attached,
  specifying the order of values */
      // mdp_int boo;
    } milc_header;

    bool endian_swap = false;
    mdp_uint size = sizeof(milc_header);

    std::ifstream fp(filename, std::ios::binary);
    if (!fp)
      return false;

    // Read header
    if (fp.read(reinterpret_cast<char *>(&milc_header), size).gcount() != size)
    {
      return false;
    }

    // Check magic number for endianess
    if (milc_header.magic_number == 0x874e0000)
    {
      switch_endianess_byte4(milc_header.dims[0]);
      switch_endianess_byte4(milc_header.dims[1]);
      switch_endianess_byte4(milc_header.dims[2]);
      switch_endianess_byte4(milc_header.dims[3]);
      endian_swap = true;
    }

    if (U.lattice().ndim() != 4 ||
        milc_header.dims[0] != U.lattice().size(1) ||
        milc_header.dims[1] != U.lattice().size(2) ||
        milc_header.dims[2] != U.lattice().size(3) ||
        milc_header.dims[3] != U.lattice().size(0))
    {
      return false;
    }

    if (endian_swap)
      return U.load(filename, processIO, max_buffer_size, false, size,
                    milc_read_as_float_switch, false);
    else
      return U.load(filename, processIO, max_buffer_size, false, size,
                    milc_read_as_float_noswitch, false);
  }
} // namespace MDP

#endif /* FERMIQCD_MILC_IO_ */
