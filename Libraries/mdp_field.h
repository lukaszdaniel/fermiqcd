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

namespace MDP
{
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
    uint32_t endianess;
    int32_t ndim;
    int32_t box[10];
    int32_t bytes_per_site;
    int32_t sites;

    mdp_field_file_header()
    {
      reset();
    }

    void reset()
    {
      strcpy(file_id, "File Type: MDP FIELD\n");
      strcpy(program_version, mdp_program_name);
      endianess = mdp_local_endianess;
    }

    void set_time()
    {
      time_t time_and_date;
      time(&time_and_date);
      strcpy(creation_date, ctime(&time_and_date));
      for (size_t i = strlen(creation_date) + 1; i < sizeof(creation_date); i++)
        creation_date[i] = '\0';
    }

    /** @brief tries to switch the endianess of the numerical members of the header
     */
    friend bool switch_header_endianess(mdp_field_file_header &header)
    {
      if (header.endianess == mdp_local_endianess)
        return false;

      switch_endianess_byte4(header.endianess);

      if (header.endianess == mdp_local_endianess)
      {
        switch_endianess_byte4(header.endianess);
        switch_endianess_byte4(header.ndim);
        switch_endianess_byte4(header.box[0]);
        switch_endianess_byte4(header.box[1]);
        switch_endianess_byte4(header.box[2]);
        switch_endianess_byte4(header.box[3]);
        switch_endianess_byte4(header.box[4]);
        switch_endianess_byte4(header.box[5]);
        switch_endianess_byte4(header.box[6]);
        switch_endianess_byte4(header.box[7]);
        switch_endianess_byte4(header.box[8]);
        switch_endianess_byte4(header.box[9]);
        switch_endianess_byte4(header.bytes_per_site);
        switch_endianess_byte4(header.sites);
        return true;
      }
      else
      {
        switch_endianess_byte4(header.endianess);
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
   *    int box[]={10,10,10};
   *    mdp_lattice lattice(3,box);
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
  template <class T>
  class mdp_field
  {
  private:
    void fill_header()
    {
      int i;
      m_header.bytes_per_site = sizeof(T) * m_field_components;
      m_header.sites = lattice().size();
      m_header.ndim = lattice().n_dimensions();
      for (i = 0; i < lattice().n_dimensions(); i++)
        m_header.box[i] = lattice().size(i);
      for (; i < 10; i++)
        m_header.box[i] = 0;
    }

  protected:
    mdp_lattice *m_ptr;          /* this points to the lattice for this field  */
    std::unique_ptr<T[]> m_data; /* this is to store the main field            */
    mdp_uint m_size;             /* this is the size of the field in sizeof(T) */
    mdp_uint m_field_components; /* this is the number of field components per site */

    /** @brief the field file header, contains data only if field was read from file
     */
    mdp_field_file_header m_header;

    /** @brief Allows dynamic allocation of a field that is not allocated
     *
     * @param a lattice where the field resides on
     * @param n size of the field vector. For n = 1 field is a scalar field
     */
    void allocate_field(mdp_lattice &a, mdp_uint n = 0)
    {
      deallocate_field();
      if (n == 0)
        n = m_field_components;
      else
        m_field_components = n;
      if (m_field_components == 0)
        error("You cannot have a field of zero size!");
      m_size = a.enclosing_volume() * m_field_components;
      m_data = std::make_unique<T[]>(m_size);
      if (m_data == nullptr)
        error("OUT OF MEMORY !!!!");
      m_ptr = &a;
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
    mdp_field() : m_ptr(nullptr), m_data(nullptr), m_size(0), m_field_components(0)
    {
    }

    /** @brief declares a field on lattice a and allocates a vector of n T at each site
     */
    mdp_field(mdp_lattice &a, mdp_uint n = 1)
    {
      allocate_field(a, n);
    }

    /** @brief Copy constructor
     */
    mdp_field(const mdp_field &field)
    {
      allocate_field(field.lattice(), field.m_field_components);
      mdp_int i_min = physical_local_start(EVENODD);
      mdp_int i_max = physical_local_stop(EVENODD);
      for (mdp_int i = i_min; i < i_max; i++)
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
      return m_data[x.local_index() * m_field_components + i];
    }

    /** @brief returns component i of the vector of objects T stored at site x
     */
    T &operator()(mdp_uint idx, mdp_uint i = 0)
    {
      return m_data[idx * m_field_components + i];
    }

    /** @brief returns the address of the vector of objects T stored at site x
     */
    T *operator[](mdp_site x)
    {
      return address(x, 0);
    }

    T &operator[](mdp_int i)
    {
      return m_data[i];
    }

    T *address(mdp_site x, int i = 0) const
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
    void shift(int i, int mu)
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

    void operator=(const mdp_field &a)
    {
      if (&lattice() != &a.lattice() ||
          m_size != a.m_size ||
          m_field_components != a.m_field_components)
        error("mdp_field: operator=() incompatible fields");

      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] = a.m_data[i];
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
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] *= a;
    }

    template <class T2>
    void operator/=(const T2 a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] /= a;
    }

    /** @brief returns by reference the lattice this field is defined on
     */
    mdp_lattice &lattice() const
    {
      return *m_ptr;
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

    /** @brief only used by mdp_field::load() and mdp_field::save()
     */
    int where_global(mdp_int global_idx)
    {
      return lattice().where_global(global_idx);
    }

    void switch_endianess_4bytes()
    {
      // I am not sure if this works for complex<double>
      mdp_int *p = nullptr;

      if (sizeof(T) * m_field_components % 4 != 0)
        error("Field not % 4");
      mdp_site x(lattice());
      forallsitesandcopies(x)
      {
        p = (mdp_int *)address(x);
        for (mdp_uint i = 0; i < sizeof(T) * m_field_components / 4; i++)
        {
          switch_endianess_byte4(*(p + i));
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
          switch_endianess_byte8(*(p + i));
        }
      }
    }

    /** @brief lattice size in units of sizeof(T)
     */
    mdp_int global_size()
    {
      return m_field_components * lattice().global_volume();
    }

    mdp_int physical_size()
    {
      return m_size;
    }

    mdp_int size_per_site()
    {
      return m_field_components;
    }

    mdp_int physical_local_start(int i = 2)
    {
      if (i == 2)
        i = 0;
      return m_field_components * lattice().start0(ME, i);
    }

    mdp_int physical_local_stop(int i = 2)
    {
      if (i == 2)
        i = 1;
      return m_field_components * lattice().stop0(ME, i);
    }

    T *physical_address(mdp_int i = 0)
    {
      return m_data.get() + i;
    }

    /** @brief communication function for a field object
     *
     * The only communication function for a field object
     * to be invoked every time field variables are assigned and
     * need to be synchronized between the parallel processes
     * The most important communication function in MDP.
     * it must be called after each field variables are modified.
     * it restores the synchronization between parallel processes.
     */
    void update(int np = 2, int d = -1, mdp_uint ncomp = 1);

    /** @brief Best way to load a field
     *
     * IO functions: load(filename, processIO, buffersize)
     *               save(filename, processIO, buffersize)
     * filename should possibly include the path.
     * processIO is the process that physically perform the IO.
     * buffersize if the size of the buffer associated to the
     *   communication to each process. buffersize*Nproc
     *   should fit in the memory of processIO.
     *   By default buffersize=1024 and it works reasonably fast.
     */
    bool load(std::string filename,
              int processIO = 0,
              mdp_int max_buffer_size = 1024,
              bool load_header = true,
              mdp_int skip_bytes = 0,
              bool (*user_read)(FILE *, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &) = nullptr,
              bool try_switch_endianess = true);

    /** @brief Best way to save a field
     */
    bool save(std::string filename,
              int processIO = 0,
              mdp_int max_buffer_size = 1024,
              bool load_header = true,
              mdp_int skip_bytes = 0,
              bool (*user_write)(FILE *, void *, mdp_int, mdp_int, mdp_int, const mdp_lattice &) = nullptr);

    /** @brief Best way to save a field to a VTK file
     *
     * @param filename String giving name of file to create
     * @param t Int that specifies which timeslice to save.  Value of -1 saves all timeslices.
     * @param component Int specifying which component to save.  Value of -1 saves all components.
     * @param processIO Int naming the ID of processor that should perform the save
     * @param ASCII Bool flag to save data in ASCII when true, binary when false
     * @return Returns true on success, throws an error for unexpected conditions.
     */
    bool save_vtk(std::string filename,
                  int t = -1,
                  int component = -1,
                  int processIO = 0,
                  bool ASCII = false);

#ifdef INCLUDE_DEPRECATED_IO
    void load(char filename[],
              int processIO,
              mdp_int max_buffer_size,
              char *header,
              mdp_int header_size = 0,
              mdp_int (*sort_x)(mdp_lattice &, mdp_int) = nullptr,
              int auto_switch_endianess = true);

    void save(char filename[],
              int processIO,
              mdp_int max_buffer_size,
              char *header,
              mdp_int header_size = 0,
              mdp_int (*sort_x)(mdp_lattice &, mdp_int) = nullptr,
              const char *mode = "w");
#endif
  };

  /** @brief Other usefull aliases.
   *
   * Parked here for now.
   */
  using mdp_complex_scalar_field = mdp_field<mdp_complex>;
  using mdp_complex_vector_field = mdp_field<mdp_complex>;

  using mdp_real_field = mdp_field<mdp_real>;
  using mdp_real_scalar_field = mdp_field<mdp_real>;
  using mdp_real_vector_field = mdp_field<mdp_real>;

  using mdp_int_field = mdp_field<mdp_int>;
  using mdp_int_scalar_field = mdp_field<mdp_int>;
  using mdp_int_vector_field = mdp_field<mdp_int>;
} // namespace MDP

#endif /* MDP_FIELD_ */
