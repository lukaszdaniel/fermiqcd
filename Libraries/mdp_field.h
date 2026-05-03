/////////////////////////////////////////////////////////////////
/// @file mdp_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_FIELD_
#define MDP_FIELD_

#include <memory>
#include <cstring>
#include <fstream>
#include "mdp_global_vars.h"
#include "mdp_complex.h"
#include "mdp_endianess_converter.h"
#include "mdp_lattice.h"
#include "mdp_site.h"

namespace MDP
{
  template <typename Src, typename Dst>
  bool mdp_write_convert(std::ofstream &file,
                         void *data,
                         mdp_int psize,
                         mdp_int header_size,
                         mdp_int position,
                         [[maybe_unused]] const mdp_lattice &lattice)
  {
    auto *src = static_cast<Src *>(data);

    const size_t src_count = psize / sizeof(Src);
    const size_t bytes_to_write = src_count * sizeof(Dst);

    auto buffer = std::make_unique<Dst[]>(src_count);

    for (mdp_uint i = 0; i < src_count; ++i)
      buffer[i] = static_cast<Dst>(src[i]);

    const std::streamoff offset = static_cast<std::streamoff>(position) * bytes_to_write + header_size;

    file.seekp(offset, std::ios::beg);
    if (!file)
      return false;

    file.write(reinterpret_cast<const char *>(buffer.get()), bytes_to_write);
    if (!file)
      return false;

    return true;
  }

  inline bool mdp_write_double_as_float(std::ofstream &file,
                                        void *data,
                                        mdp_int psize,
                                        mdp_int header_size,
                                        mdp_int position,
                                        const mdp_lattice &lattice)
  {
    return mdp_write_convert<double, float>(file, data, psize, header_size, position, lattice);
  }

  inline bool mdp_write_float_as_double(std::ofstream &file,
                                        void *data,
                                        mdp_int psize,
                                        mdp_int header_size,
                                        mdp_int position,
                                        const mdp_lattice &lattice)
  {
    return mdp_write_convert<float, double>(file, data, psize, header_size, position, lattice);
  }

  template <typename Dst, typename Src>
  bool mdp_read_convert(std::ifstream &file,
                        void *data,
                        mdp_int psize,
                        mdp_int header_size,
                        mdp_int position,
                        [[maybe_unused]] const mdp_lattice &lattice)
  {
    auto *dst = static_cast<Dst *>(data);

    const size_t dst_count = psize / sizeof(Dst);
    const size_t bytes_to_read = dst_count * sizeof(Src);

    auto buffer = std::make_unique<Src[]>(dst_count);

    const std::streamoff offset = static_cast<std::streamoff>(position) * bytes_to_read + header_size;

    file.seekg(offset, std::ios::beg);
    if (!file)
      return false;

    file.read(reinterpret_cast<char *>(buffer.get()), bytes_to_read);
    if (!file)
      return false;

    for (mdp_uint i = 0; i < dst_count; ++i)
      dst[i] = static_cast<Dst>(buffer[i]);

    return true;
  }

  inline bool mdp_read_double_as_float(std::ifstream &file,
                                       void *data,
                                       mdp_int psize,
                                       mdp_int header_size,
                                       mdp_int position,
                                       const mdp_lattice &lattice)
  {
    return mdp_read_convert<double, float>(file, data, psize, header_size, position, lattice);
  }

  inline bool mdp_read_float_as_double(std::ifstream &file,
                                       void *data,
                                       mdp_int psize,
                                       mdp_int header_size,
                                       mdp_int position,
                                       const mdp_lattice &lattice)
  {
    return mdp_read_convert<float, double>(file, data, psize, header_size, position, lattice);
  }

  /** @brief header for field file IO
   *
   * Used to store the binary file header that precedes the data
   * when storing an object of class mdp_field<> in a file
   */
  class mdp_field_file_header
  {
  public:
    char file_id[60];
    char program_version[60];
    char creation_date[60];
    mdp_uint endianess;
    mdp_suint ndim;
    mdp_uint box[10];
    mdp_uint bytes_per_site;
    mdp_uint sites;

    mdp_field_file_header()
    {
      reset();
    }

    void reset()
    {
      strcpy(file_id, "File Type: MDP FIELD\n");
      strcpy(program_version, mdp_program_name);
      endianess = mdp_local_endianess;
      // ndim = 0;
      // bytes_per_site = 0;
      // sites = 0;
      // for (int i = 0; i < 10; i++)
      //   box[i] = 1;
      // set_time();
    }

    void set_time()
    {
      time_t time_and_date;
      time(&time_and_date);
      strcpy(creation_date, ctime(&time_and_date));
      for (size_t i = strlen(creation_date); i < sizeof(creation_date); i++)
        creation_date[i] = '\0';
    }

    /** @brief tries to switch the endianess of the numerical members of the header
     */
    static bool switch_header_endianess(mdp_field_file_header &header)
    {
      if (header.endianess == mdp_local_endianess)
        return false;

      switch_endianess(header.endianess);

      if (header.endianess == mdp_local_endianess)
      {
        switch_endianess(header.endianess);
        switch_endianess(header.ndim);
        switch_endianess(header.box[0]);
        switch_endianess(header.box[1]);
        switch_endianess(header.box[2]);
        switch_endianess(header.box[3]);
        switch_endianess(header.box[4]);
        switch_endianess(header.box[5]);
        switch_endianess(header.box[6]);
        switch_endianess(header.box[7]);
        switch_endianess(header.box[8]);
        switch_endianess(header.box[9]);
        switch_endianess(header.bytes_per_site);
        switch_endianess(header.sites);
        return true;
      }
      else
      {
        switch_endianess(header.endianess);
        return false;
      }
    }
  };

  /** @brief most generic field object
   *
   * mdp_field is a vector field of n components of type T.
   * scalar field is a special case where n = 1.
   *
   * Example:
   * @verbatim
   *    constexpr Box box = {10,10,10};
   *    mdp_lattice lattice(box);
   *    mdp_field<mdp_real> psi(lattice,10);
   *    mdp_site x(lattice);
   *    forallsites(x)
   *      for(int i=0; i<10; i++)
   *         psi(x,i)=0.0;
   *    psi.update(); // synchronization
   *    psi.save("myfield");
   *    psi.load("myfield");
   * @endverbatim
   */
  template <class T = mdp_complex>
  class mdp_field
  {
  private:
    void fill_header()
    {
      mdp_suint i;
      m_header.bytes_per_site = sizeof(T) * m_field_components;
      m_header.sites = lattice().size();
      m_header.ndim = lattice().n_dimensions();
      for (i = 0; i < lattice().n_dimensions(); i++)
        m_header.box[i] = lattice().size(i);
      for (; i < 10; i++)
        m_header.box[i] = 0;
    }

    /** @brief only used by mdp_field::load() and mdp_field::save()
     */
    mdp_uint where_global(mdp_uint global_idx) const
    {
      return lattice().where_global(global_idx);
    }

  protected:
    const mdp_lattice *m_lattice; /* this points to the lattice for this field  */
    std::unique_ptr<T[]> m_data;  /* this is to store the main field            */
    mdp_uint m_size;              /* this is the total number field components on the lattice */
    mdp_uint m_field_components;  /* this is the number of field components per site */

    /** @brief the field file header, contains data only if field was read from file
     */
    mdp_field_file_header m_header;

    /** @brief Allows dynamic allocation of a field that is not allocated
     *
     * @param a lattice where the field resides on
     * @param n size of the field vector. For n = 1 field is a scalar field
     */
    void allocate_field(const mdp_lattice &a, mdp_uint n = 0)
    {
      deallocate_field();
      if (n == 0)
        error("You cannot have a field of zero size!");
      m_field_components = n;
      m_size = a.enclosing_volume() * m_field_components;
      m_data = std::make_unique<T[]>(m_size);
      m_lattice = &a;
      fill_header();
    }

    void deallocate_memory()
    {
      m_data = nullptr;
      m_size = m_field_components = 0;
    }

  public:
    /** @brief declare empty field (zero size)
     */
    mdp_field() : m_lattice(nullptr), m_data(nullptr), m_size(0), m_field_components(0)
    {
    }

    /** @brief declares a field on lattice a and allocates a vector of n T at each site
     */
    mdp_field(const mdp_lattice &a, mdp_uint n = 1)
    {
      allocate_field(a, n);
    }

    /** @brief Copy constructor
     */
    mdp_field(const mdp_field &field)
    {
      allocate_field(field.lattice(), field.m_field_components);
      mdp_uint i_min = physical_local_start(EVENODD);
      mdp_uint i_max = physical_local_stop(EVENODD);
      for (mdp_uint i = i_min; i < i_max; i++)
        m_data[i] = field.m_data[i];
    }

    /** @brief checks if a field is allocated of has zero-size
     */
    bool allocated()
    {
      return (m_data != nullptr);
    }

    /** @brief dynamically deallocate field
     */
    void deallocate_field()
    {
      deallocate_memory();
    }

    void reset_field()
    {
      deallocate_memory();
    }

    virtual ~mdp_field()
    {
      deallocate_memory();
    }

    /** @brief returns component i of the vector of objects T stored at site x
     */
    T &operator()(mdp_site x, mdp_uint i = 0)
    {
#ifdef CHECK_ALL
      if (!(x.is_here()))
      {
        error("You are looking for a site that is not here");
      }
#endif
#ifdef CHECK_BOUNDARY
      if (i >= m_field_components)
      {
        error(std::format("field rows can be indexed up to {}", m_field_components - 1));
      }
#endif
      return m_data[x.local_index() * m_field_components + i];
    }

    /** @brief returns const component i of the vector of objects T stored at site x
     */
    const T &operator()(mdp_site x, mdp_uint i = 0) const
    {
#ifdef CHECK_ALL
      if (!(x.is_here()))
      {
        error("You are looking for a site that is not here");
      }
#endif
#ifdef CHECK_BOUNDARY
      if (i >= m_field_components)
      {
        error(std::format("field rows can be indexed up to {}", m_field_components - 1));
      }
#endif
      return m_data[x.local_index() * m_field_components + i];
    }

    /** @brief returns component i of the vector of objects T stored at site x
     */
    T &operator()(mdp_uint idx, mdp_uint i = 0)
    {
#ifdef CHECK_BOUNDARY
      if (i >= m_field_components)
      {
        error("field rows can be indexed up to " + (m_field_components - 1));
      }
#endif
      return m_data[idx * m_field_components + i];
    }

    /** @brief returns const component i of the vector of objects T stored at site x
     */
    const T &operator()(mdp_uint idx, mdp_uint i = 0) const
    {
#ifdef CHECK_BOUNDARY
      if (i >= m_field_components)
      {
        error("field rows can be indexed up to " + (m_field_components - 1));
      }
#endif
      return m_data[idx * m_field_components + i];
    }

    /** @brief returns the address of the vector of objects T stored at site x
     */
    T *operator[](mdp_site x) const
    {
      return address(x, 0);
    }

    T &operator[](mdp_int i) const
    {
      return m_data[i];
    }

    T *address(mdp_site x, mdp_int i = 0) const
    {
#ifdef CHECK_ALL
      if (!(x.is_here()))
      {
        error("You are looking for a site that is not here");
      }
#endif
      return m_data.get() + x.local_index() * m_field_components + i;
    }

    /** @brief shifts the entire fields in direction mu of i steps
     *
     * (i can be positive or negative)
     * note that if i=1, field(x-mu) is assigned to field(x)
     * function requires communication
     */
    void shift(int i, mdp_suint mu)
    {
      mdp_field tmp(lattice(), m_field_components);
      mdp_site x(lattice());

      while (i != 0)
      {
        update();
        if (i == +1)
        {
          forallsites(x)
          {
            for (mdp_uint k = 0; k < m_field_components; k++)
              tmp(x, k) = (*this)(x - mu, k); // mind here
          }
          i--;
        }
        else if (i == -1)
        {
          forallsites(x)
          {
            for (mdp_uint k = 0; k < m_field_components; k++)
              tmp(x, k) = (*this)(x + mu, k); // mind here
          }
          i++;
        }
        (*this) = tmp;
      }
    }

    mdp_field &operator=(const mdp_field &a)
    {
      if (this == &a)
        return *this;

      if (&lattice() != &a.lattice() ||
          m_size != a.m_size ||
          m_field_components != a.m_field_components)
        error("mdp_field: operator=() incompatible fields");

      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] = a.m_data[i];

      return *this;
    }

    void operator=(const T a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] = a;
    }

    void operator+=(const mdp_field &a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] += a.m_data[i];
    }

    void operator-=(const mdp_field &a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] -= a.m_data[i];
    }

    template <class T2>
    void operator*=(const T2 a)
    {
      if (a == T2(1.0)) [[unlikely]]
        return;

      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] *= a;
    }

    template <class T2>
    void operator/=(const T2 a)
    {
      if (a == T2(1.0)) [[unlikely]]
        return;

      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] /= a;
    }

    /** @brief returns by reference the lattice this field is defined on
     */
    const mdp_lattice &lattice() const
    {
      return *m_lattice;
    }

    /** Dimension of the field
     */
    mdp_suint ndim() const
    {
      return m_lattice->n_dimensions();
    }

    /** @brief returns the total memory in bytes occupied by the field
     */
    mdp_int field_size()
    {
      return lattice().size() * m_field_components * sizeof(T);
    }

    /** @brief returns the total space in bytes required to store the field
     */
    mdp_int file_size()
    {
      return sizeof(mdp_field_file_header) + field_size();
    }

    void switch_endianess_4bytes()
    {
      // This method might not work for complex<double> because it assumes that the data type T
      // can be safely cast to int32_t for endianess switching, which may not be true for complex<double>.
      mdp_int *p = nullptr;

      if (sizeof(T) * m_field_components % 4 != 0)
        error("Field not % 4");
      mdp_site x(lattice());
      forallsitesandcopies(x)
      {
        p = (mdp_int *)address(x);
        for (mdp_uint i = 0; i < sizeof(T) * m_field_components / 4; i++)
        {
          switch_endianess(*(p + i));
        }
      }
    }

    void switch_endianess_8bytes()
    {
      // I am not sure if this works for complex<double>
      int64_t *p;

      if (sizeof(T) * m_field_components % 8 != 0)
        error("Field not % 8");
      mdp_site x(lattice());
      forallsitesandcopies(x)
      {
        p = (int64_t *)address(x);
        for (mdp_uint i = 0; i < (sizeof(T) * m_field_components / 8); i++)
        {
          switch_endianess(*(p + i));
        }
      }
    }

    /** @brief lattice size in units of sizeof(T)
     */
    mdp_uint global_size() const
    {
      return m_field_components * lattice().global_volume();
    }

    mdp_uint physical_size() const
    {
      return m_size;
    }

    mdp_uint size_per_site() const
    {
      return m_field_components;
    }

    mdp_uint physical_local_start(mdp_parity i = EVENODD) const
    {
      if (i == EVENODD)
        i = EVEN;
      return m_field_components * lattice().start0(ME, i);
    }

    mdp_uint physical_local_stop(mdp_parity i = EVENODD) const
    {
      if (i == EVENODD)
        i = ODD;
      return m_field_components * lattice().stop0(ME, i);
    }

    T *physical_address(mdp_uint i = 0)
    {
      return m_data.get() + i;
    }

    /**
     * @brief Synchronizes field data between MPI processes (halo exchange).
     *
     * This method exchanges boundary (ghost/halo) data between neighboring
     * processes to ensure consistency of distributed field values across
     * the lattice decomposition.
     *
     * Each process sends selected components of its local field to other
     * processes that own neighboring sites and receives corresponding data
     * to update its ghost regions.
     *
     * @param np     Parity selector:
     *               - EVENODD (default): update both even and odd sites
     *               - EVEN or ODD: update only selected subset
     * @param d      Starting component index (used when updating only part
     *               of the field). If -1, all components are updated.
     * @param ncomp  Number of components to update starting from @p d.
     *
     * @note
     * - This function performs a nearest-neighbor communication (halo exchange),
     *   not a global reduction or averaging.
     * - Only boundary/ghost data are exchanged; interior sites are untouched.
     * - The communication pattern is defined by the lattice partitioning.
     *
     * @warning
     * - Must be called after modifying field values that are used by other
     *   processes (e.g., before stencil operations or nearest-neighbor access).
     * - Not required if only strictly local (non-neighbor-dependent) computations
     *   are performed.
     *
     * @usage
     * Typical usage pattern:
     * @code
     * // modify local field values
     * field(x) = ...;
     *
     * // synchronize ghost cells before using neighbor values
     * field.update();
     *
     * // now safe to use field at neighboring sites
     * @endcode
     *
     * @note
     * - Partial updates (using @p d and @p ncomp) allow reducing communication
     *   when only some components have changed.
     * - EVEN/ODD splitting can be used for performance optimization in
     *   checkerboard-style algorithms.
     */
    void update(mdp_parity np = EVENODD, int d = -1, mdp_uint ncomp = 1);

    /**
     * @brief Collective parallel read of a field from a binary file.
     *
     * This method loads field data from a single binary file using the process
     * identified by @p processIO and distributes the data across MPI processes.
     * The I/O process reads the file in global lattice order and sends buffered
     * chunks to the appropriate processes.
     *
     * Data are expected to be stored as contiguous blocks of
     * @c m_field_components elements of type @c T per lattice site,
     * in global lattice ordering (as produced by save()).
     *
     * @param filename              Name of the input file (resolved via latest_file()).
     * @param processIO             Rank of the process responsible for file I/O.
     * @param max_buffer_size       Maximum number of lattice sites transferred
     *                              per communication buffer.
     * @param load_header           If true, reads and validates an
     *                              @c mdp_field_file_header from the file.
     * @param skip_bytes            Number of bytes to skip before reading data
     *                              (useful for custom headers).
     * @param user_read             Optional callback for custom reading.
     *                              If provided, it is called for each lattice site:
     *                              @code
     *                              bool user_read(std::ifstream& fp,
     *                                             void* data,
     *                                             mdp_int nbytes,
     *                                             mdp_int offset,
     *                                             mdp_int global_index,
     *                                             const mdp_lattice& lattice);
     *                              @endcode
     *                              If null, raw binary data are read.
     * @param try_switch_endianess  If true, automatically converts data endianess
     *                              when the file header indicates mismatch.
     *
     * @return true on success, false if header validation fails.
     *
     * @throws std::format if the file does not exist.
     * @throws error on I/O failures or unexpected end of file.
     *
     * @note
     * - Only the @p processIO performs file reads.
     * - Other processes receive their data via buffered communication.
     * - If @p load_header is enabled, the file header is validated against the
     *   current lattice (dimensions, size, data layout).
     * - Endianess mismatch is detected from the header and optionally corrected
     *   after loading.
     * - Missing lattice sites (NOWHERE) trigger repositioning in the input stream.
     * - The method assumes the file was written with a compatible layout (see save()).
     * - After loading, update() is called and endianess information is broadcast.
     */
    bool load(std::string filename,
              mdp_uint processIO = 0,
              mdp_int max_buffer_size = 1024,
              bool load_header = true,
              mdp_int skip_bytes = 0,
              bool (*user_read)(std::ifstream &, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &) = nullptr,
              bool try_switch_endianess = true);

    /**
     * @brief Collective parallel write of the field to a binary file.
     *
     * This method gathers field data distributed across MPI processes and writes
     * them into a single binary file using the process identified by @p processIO.
     * Other processes send their local data in buffered chunks.
     *
     * Data are written in global lattice order as contiguous blocks of
     * @c m_field_components elements of type @c T per lattice site.
     *
     * The method supports writing an optional file header and allows full
     * customization of the output format via a user-provided callback.
     *
     * @param filename        Name of the output file (may be modified by
     *                        next_to_latest_file() to avoid overwriting).
     * @param processIO       Rank of the process responsible for file I/O.
     * @param max_buffer_size Maximum number of lattice sites sent in one buffer
     *                        from each process.
     * @param load_header     (Unused / legacy) retained for API compatibility.
     * @param skip_bytes      Number of bytes to skip at the beginning of the file
     *                        before writing data (useful for custom headers).
     * @param user_write      Optional callback for custom writing.
     *                        If provided, it is called for each lattice site:
     *                        @code
     *                        bool user_write(std::ofstream& fp,
     *                                        void* data,
     *                                        mdp_int nbytes,
     *                                        mdp_int offset,
     *                                        mdp_int global_index,
     *                                        const mdp_lattice& lattice);
     *                        @endcode
     *                        If null, raw binary data are written.
     *
     * @return true on success.
     *
     * @throws error if the file cannot be opened or a write operation fails.
     *
     * @note
     * - Only the @p processIO performs file operations.
     * - Other processes act as data providers via buffered communication.
     * - If @p save_header is enabled, an @c mdp_field_file_header is written
     *   at the beginning of the file (after @p skip_bytes).
     * - If @p user_write is provided, it overrides the default raw binary write.
     * - The method assumes a consistent global lattice ordering across processes.
     * - Missing lattice sites (NOWHERE) trigger repositioning in the output stream.
     */
    bool save(std::string filename,
              mdp_uint processIO = 0,
              mdp_int max_buffer_size = 1024,
              bool load_header = true,
              mdp_int skip_bytes = 0,
              bool (*user_write)(std::ofstream &, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &) = nullptr);

    /**
     * @brief Collective parallel write of the field to a VTK file (STRUCTURED_POINTS).
     *
     * This method gathers field data distributed across MPI processes and writes it
     * into a single VTK file using the process identified by @p processIO.
     * Other processes send their local data in buffered chunks to the I/O process.
     *
     * Supported only for 3D or 4D lattices:
     * - In 4D, data are interpreted as a sequence of 3D time slices.
     * - In 3D, the entire field is written as a single dataset (t_slice ignored).
     *
     * The output file is written via a temporary file and then atomically renamed
     * to avoid partial/corrupted results.
     *
     * @param filename   Base name of the output VTK file (a temporary ".tmp" file is used internally).
     * @param t_slice    Time slice to save (4D only).
     *                   -1 writes all time slices.
     * @param component  Field component to save.
     *                   -1 writes all components.
     * @param processIO  Rank of the process responsible for file I/O.
     *                   This process gathers data from all others.
     * @param ASCII      If true, writes ASCII VTK; otherwise writes binary (with endian conversion).
     *
     * @return true on success.
     *
     * @throws std::ios_base::failure if the output file cannot be opened.
     * @throws error if lattice dimensionality is not 3 or 4, or on I/O failure.
     *
     * @note
     * - Data are written in global lattice order.
     * - In parallel mode, non-I/O processes do not perform file operations.
     * - Binary output uses mdp_real precision and applies endian conversion.
     * - File name may be modified by next_to_latest_file() to avoid overwriting.
     */
    bool save_vtk(std::string filename,
                  int t_slice = -1,
                  int component = -1,
                  mdp_uint processIO = 0,
                  bool ASCII = false);

    bool save_as_float(const std::string &filename,
                       mdp_uint processIO = 0,
                       mdp_int max_buffer_size = 1024,
                       bool load_header = true,
                       mdp_int skip_bytes = 0)
    {
#ifdef USE_SINGLE_PRECISION
      save(filename, processIO, max_buffer_size, load_header, skip_bytes, nullptr);
#else
      m_header.bytes_per_site /= 2;
      save(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_write_double_as_float);
      m_header.bytes_per_site *= 2;
#endif
      return true;
    }

    bool load_as_float(const std::string &filename,
                       mdp_uint processIO = 0,
                       mdp_int max_buffer_size = 1024,
                       bool load_header = true,
                       mdp_int skip_bytes = 0)
    {

#ifdef USE_SINGLE_PRECISION
      load(filename, processIO, max_buffer_size, load_header, skip_bytes, nullptr, true);
#else
      m_header.bytes_per_site /= 2;
      load(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_read_double_as_float, true);
      m_header.bytes_per_site *= 2;
#endif
      return true;
    }

    bool load_as_double(const std::string &filename,
                        mdp_uint processIO = 0,
                        mdp_int max_buffer_size = 1024,
                        bool load_header = true,
                        mdp_int skip_bytes = 0)
    {
#ifdef USE_SINGLE_PRECISION
      m_header.bytes_per_site *= 2;
      load(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_read_float_as_double, true);
      m_header.bytes_per_site /= 2;
#else
      load(filename, processIO, max_buffer_size, load_header, skip_bytes, nullptr, true);
#endif
      return true;
    }

    bool save_as_double(const std::string &filename,
                        mdp_uint processIO = 0,
                        mdp_int max_buffer_size = 1024,
                        bool load_header = true,
                        mdp_int skip_bytes = 0)
    {
#ifdef USE_SINGLE_PRECISION
      m_header.bytes_per_site *= 2;
      save(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_write_float_as_double);
      m_header.bytes_per_site /= 2;
#else
      save(filename, processIO, max_buffer_size, load_header, skip_bytes, nullptr);
#endif
      return true;
    }
  };

  /** @brief Other usefull aliases.
   *
   * Parked here for now.
   */
  using mdp_complex_field = mdp_field<mdp_complex>;
  using mdp_complex_scalar_field = mdp_field<mdp_complex>;
  using mdp_complex_vector_field = mdp_field<mdp_complex>;

  using mdp_real_field = mdp_field<mdp_real>;
  using mdp_real_scalar_field = mdp_field<mdp_real>;
  using mdp_real_vector_field = mdp_field<mdp_real>;

  using mdp_int_field = mdp_field<mdp_int>;
  using mdp_int_scalar_field = mdp_field<mdp_int>;
  using mdp_int_vector_field = mdp_field<mdp_int>;

  mdp_real norm_square(mdp_complex_field &psi,
                       mdp_parity parity = EVENODD)
  {
    mdp_real n2 = 0;
    mdp_uint i_min = psi.physical_local_start(parity);
    mdp_uint i_max = psi.physical_local_stop(parity);

    for (mdp_uint i = i_min; i < i_max; i++)
      n2 += abs2(psi[i]);

    mdp.add(n2);
    return n2;
  }

  mdp_complex scalar_product(mdp_complex_field &psi,
                             mdp_complex_field &chi,
                             mdp_parity parity = EVENODD)
  {
    mdp_complex n2 = 0;
    mdp_uint i_min = psi.physical_local_start(parity);
    mdp_uint i_max = psi.physical_local_stop(parity);

    for (mdp_uint i = i_min; i < i_max; i++)
      n2 += conj(psi[i]) * chi[i];

    mdp.add(n2);

    return n2;
  }

  mdp_real real_scalar_product(mdp_complex_field &psi,
                               mdp_complex_field &chi,
                               mdp_parity parity = EVENODD)
  {

    mdp_real n2 = 0;
    mdp_uint i_min = psi.physical_local_start(parity);
    mdp_uint i_max = psi.physical_local_stop(parity);

    for (mdp_uint i = i_min; i < i_max; i++)
    {
      n2 +=
          real(chi[i]) * real(psi[i]) +
          imag(chi[i]) * imag(psi[i]);
    }

    mdp.add(n2);
    return n2;
  }

  mdp_real imag_scalar_product(mdp_complex_field &psi,
                               mdp_complex_field &chi,
                               mdp_parity parity = EVENODD)
  {
    mdp_real n2 = 0;
    mdp_uint i_min = psi.physical_local_start(parity);
    mdp_uint i_max = psi.physical_local_stop(parity);

    for (mdp_uint i = i_min; i < i_max; i++)
    {
      n2 +=
          real(psi[i]) * imag(chi[i]) +
          imag(psi[i]) * real(chi[i]);
    }
    mdp.add(n2);
    return n2;
  }

  void mdp_add_scaled_field(mdp_complex_field &psi,
                            mdp_real alpha,
                            mdp_complex_field &chi,
                            mdp_parity parity = EVENODD)
  {
    mdp_uint i_min = psi.physical_local_start(parity);
    mdp_uint i_max = psi.physical_local_stop(parity);

    for (mdp_uint i = i_min; i < i_max; i++)
      psi[i] += alpha * chi[i];
  }

  void mdp_add_scaled_field(mdp_complex_field &psi,
                            mdp_complex alpha,
                            mdp_complex_field &chi,
                            mdp_parity parity = EVENODD)
  {
    mdp_uint i_min = psi.physical_local_start(parity);
    mdp_uint i_max = psi.physical_local_stop(parity);

    // this needs optimization.
    for (mdp_uint i = i_min; i < i_max; i++)
      psi[i] += alpha * chi[i];
  }

  mdp_complex operator*(mdp_complex_field &psi,
                        mdp_complex_field &chi)
  {
    return scalar_product(psi, chi);
  }

  mdp_real relative_residue(mdp_complex_field &p,
                            mdp_complex_field &q,
                            mdp_parity parity = EVENODD)
  {
    mdp_real residue = 0, num = 0, den = 0;
    mdp_uint i_min = p.physical_local_start(parity);
    mdp_uint i_max = q.physical_local_stop(parity);

    // this needs optimization.
    for (mdp_uint i = i_min; i < i_max;)
    {
      num += abs2(p[i]);
      den += abs2(q[i]);
      if (++i % p.size_per_site() == 0)
      {
        residue += (den == 0) ? 1.0 : (num / den);
        num = den = 0;
      }
    }
    mdp.add(residue);
    return std::sqrt(residue / p.lattice().global_volume());
  }
} // namespace MDP

#endif /* MDP_FIELD_ */
