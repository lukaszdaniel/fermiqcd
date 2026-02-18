/////////////////////////////////////////////////////////////////
/// @file mdp_lattice.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_lattice
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_LATTICE_
#define MDP_LATTICE_

#define MDP_LATTICE

#include <memory>
#include <string>
#include <algorithm>
#include <fstream>
#include <array>
#include <cstddef>
#include "mdp_global_vars.h"
#include "mdp_prng.h"
#include "mdp_partitionings.h"
#include "mdp_topologies.h"

#define MAX_DIM 10

namespace MDP
{
  class mdp_site;

  struct Box
  {
    constexpr Box(std::initializer_list<int> init) noexcept
        : m_ndim(init.size())
    {
      mdp_int i = 0;
      for (int v : init)
      {
        m_box[i++] = v;
      }
    }

    template <mdp_int N>
    constexpr Box(const int (&arr)[N]) noexcept
        : m_ndim(N)
    {
      for (mdp_int i = 0; i < N && i < MAX_DIM; ++i)
      {
        m_box[i] = arr[i];
      }
    }

    constexpr mdp_int dim() const noexcept { return m_ndim; }

    constexpr mdp_int dim(mdp_int i) const noexcept { return m_box[i]; }

    constexpr mdp_int operator[](mdp_int i) const noexcept
    {
      return m_box[i];
    }

    constexpr mdp_int &operator[](mdp_int i) noexcept
    {
      return m_box[i];
    }

    constexpr mdp_int volume() const noexcept
    {
      mdp_int v = 1;
      for (mdp_int i = 0; i < m_ndim; ++i)
        v *= m_box[i];
      return v;
    }

    constexpr const mdp_int *get() const noexcept { return m_box.data(); }

  private:
    std::array<mdp_int, MAX_DIM> m_box{};
    mdp_int m_ndim{0};
  };

  /// @brief distributed lattice object
  ///
  /// Example:
  /// @verbatim
  ///    constexpr Box box = {3,3,3};
  ///    int seed=0, border_width=1;
  ///    mdp_lattice lattice(box,default_partitioning0,
  ///                        torus_topology,seed,border_width);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      std::cout << lattice.random(x).plain() << std::endl;
  /// @endverbatim
  class mdp_lattice
  {
  private:
    int m_ndim;                   /* number of dimensions       */
    int m_ndir;                   /* number of directions       */
    int m_next_next;              /* 1, 2 or 3 is the thickness of the boundary */
    int *m_nx;                    /* box containing the lattice */
    mdp_int m_local_volume;       /* local volume               */
    mdp_int m_global_volume;      /* global volume              */
    mdp_int m_internal_volume;    /* internal volume            */
    mdp_int *m_global_from_local; /* map local to global        */
    mdp_int *m_local_from_global; /* map global to local        */
#ifdef MDP_NO_LG
    mutable std::fstream m_lg_file; /* temporary file to store local_from_global map if not enough memory */
#endif
    mdp_int **m_up; /* move up in local index     */
    mdp_int **m_dw; /* move dw in local index     */
    int **m_co;     /* coordinate x in local idx  */
    int *m_wh;      /* in which process? is idx   */
    int *m_parity;  /* parity of local mdp_site idx   */
    mdp_int m_start[_NprocMax_][2];
    mdp_int m_stop[_NprocMax_][2];
    mdp_int m_len_to_send[_NprocMax_][2];
    mdp_int *m_to_send[_NprocMax_];
    bool m_local_random_generator;
    std::unique_ptr<mdp_prng[]> m_random_obj;
    mdp_int m_random_seed;
    using where_fn = int (*)(const int x[], const int ndim, const int nx[]);
    using neighbour_fn = void (*)(const int mu,
                                  int x_dw[],
                                  const int x[],
                                  int x_up[],
                                  const int ndim,
                                  const int nx[]);

    where_fn m_where;
    neighbour_fn m_neighbour;

    // helper function
    inline bool is_me(const int x[]) const noexcept
    {
      const int owner_id = m_where(x, m_ndim, m_nx);
      return owner_id == ME;
    }

    // helper function
    inline bool is_valid_process(const int x[]) const noexcept
    {
      const int owner_id = m_where(x, m_ndim, m_nx);
      return (owner_id < Nproc);
    }

    // helper function
    void init_basic_parameters(const Box &box,
                               int ndir_,
                               int (*where_)(const int[], const int, const int[]),
                               void (*neighbour_)(const int, int[], const int[], int[], const int, const int[]),
                               int next_next_,
                               bool local_random_)
    {
      deallocate_memory();

      // //////////////////////////////////////////////////////////////////
      // (*where)(x,nc) must return the processor rank of where x is stored
      // (*neghbour)(mu,x_dw,x,x_up,ndim, nx) must fill x_dw and x_up
      //            according with current position x and direction mu
      // //////////////////////////////////////////////////////////////////

      m_ndim = box.dim();
      m_ndir = ndir_;
      m_where = where_;
      m_neighbour = neighbour_;
      m_next_next = next_next_;
      m_local_random_generator = local_random_;

      m_local_volume = 0;
      m_internal_volume = 0;
      m_global_volume = box.volume();

      if (m_ndim != m_ndir)
        mdp << "Warning: The number of dimensions (ndim) does not match the number of directions (ndir). This may cause unexpected behavior in lattice operations.\n";

      mdp << "Lattice dimension: " << box[0];
      for (mdp_int mu = 1; mu < m_ndim; mu++)
        mdp << " x " << box[mu];
      mdp << "\n";

      ///////////////////////////////////////////////////////////////////
      // Dynamically allocate lattice size
      ///////////////////////////////////////////////////////////////////
      m_nx = new int[m_ndim];
      for (mdp_int mu = 0; mu < m_ndim; ++mu)
        m_nx[mu] = box[mu];
    }

    // helper function
#ifndef MDP_NO_LG
    void compute_local_sites(int *x, int *x_up, int *x_dw, mdp_int *local_mdp_sites)
#else
    void compute_local_sites(int *x, int *x_up, int *x_dw, std::fstream &lms_file)
#endif
    {
      int x_up_dw[MAX_DIM], x_up_up[MAX_DIM], x_up_up_dw[MAX_DIM], x_up_up_up[MAX_DIM];
      int x_dw_up[MAX_DIM], x_dw_dw[MAX_DIM], x_dw_dw_dw[MAX_DIM], x_dw_dw_up[MAX_DIM];

      // ///////////////////////////////////////////////////////////////////
      // Fill table local_mdp_sites with the global coordinate of mdp_sites in ME
      // m_local_volume counts the mdp_sites stored in local_mdp_sites
      // ///////////////////////////////////////////////////////////////////
      //
      // cases of interest...:
      //
      // m_next_next < 0 => only local mdp_sites NEW!
      // m_next_next = 0 => only x+mu and x-mu (wilson) NEW!
      // m_next_next = 1 => only x+-mu, x+-mu+-nu mu!=nu (clover)
      // m_next_next = 2 => only x+-mu, x+-mu+-nu
      // m_next_next = 3 => only x+-mu, x+-mu+-nu, x+-mu+-nu+-rho (asqtad)
      //
      // ///////////////////////////////////////////////////////////////////

      for (mdp_int mu = 0; mu < m_ndim; mu++)
        x[mu] = 0;

      do
      {
        const mdp_int global_idx = global_coordinate(x);

        if (is_me(x))
        {
#ifndef MDP_NO_LG
          local_mdp_sites[m_local_volume] = global_idx;
#else
          lms_file.seekp(m_local_volume * sizeof(mdp_int), std::ios::beg);
          lms_file.write(reinterpret_cast<const char *>(&global_idx), sizeof(mdp_int));

          if (!lms_file)
            error("mdp_lattice::allocate_lattice()\n"
                  "Unable to write to temporary file");
#endif
          ++m_local_volume;
        }
        else if (m_next_next >= 0)
        {
          bool is_boundary = false;
          for (mdp_int mu = 0; mu < m_ndir; mu++)
          {
            // calculate up and down points in mu direction given the neighbour topology
            (*m_neighbour)(mu, x_dw, x, x_up, m_ndim, m_nx);

            if (!is_valid_process(x_up) || !is_valid_process(x_dw))
              error("Incorrect partitioning");

            if (is_me(x_up) || is_me(x_dw))
              is_boundary = true;

            // ////////////////////////////////////////////////////
            // cases:
            // 1) mu+nu
            // 2) mu+nu and mu+mu
            // 3) mu+nu+rho and mu+mu+rho and mu+mu+mu
            // One may want to optimize case 3. It was thought for
            // improved staggered fermions.
            // ////////////////////////////////////////////////////
            for (mdp_int nu = 0; nu < m_ndir; nu++)
            {
              if ((nu != mu) || (m_next_next > 1))
              {
                (*m_neighbour)(nu, x_dw_dw, x_dw, x_dw_up, m_ndim, m_nx);
                (*m_neighbour)(nu, x_up_dw, x_up, x_up_up, m_ndim, m_nx);

                if (!is_valid_process(x_up_dw) ||
                    !is_valid_process(x_up_up) ||
                    !is_valid_process(x_dw_dw) ||
                    !is_valid_process(x_dw_up))
                  error("Incorrect partitioning");

                if (is_me(x_up_dw) ||
                    is_me(x_up_up) ||
                    is_me(x_dw_dw) ||
                    is_me(x_dw_up))
                  is_boundary = true;

                // mu-nu-rho terms (asqtad case)
                if (m_next_next == 3)
                {
                  for (mdp_int rho = 0; rho < m_ndir; rho++)
                  {
                    (*m_neighbour)(rho, x_dw_dw_dw, x_dw_dw, x_dw_dw_up, m_ndim, m_nx);
                    (*m_neighbour)(rho, x_up_up_dw, x_up_up, x_up_up_up, m_ndim, m_nx);

                    if (!is_valid_process(x_up_up_up) ||
                        !is_valid_process(x_up_up_dw) ||
                        !is_valid_process(x_dw_dw_up) ||
                        !is_valid_process(x_dw_dw_dw))
                      error("Incorrect partitioning");

                    if (is_me(x_up_up_up) ||
                        is_me(x_up_up_dw) ||
                        is_me(x_dw_dw_up) ||
                        is_me(x_dw_dw_dw))
                      is_boundary = true;
                  }
                }
              }
            }
          }

          if (is_boundary)
          {
#ifndef MDP_NO_LG
            local_mdp_sites[m_local_volume] = global_idx;
#else
            lms_file.seekp(m_local_volume * sizeof(mdp_int), std::ios::beg);
            lms_file.write(reinterpret_cast<const char *>(&global_idx), sizeof(mdp_int));

            if (!lms_file)
              error("mdp_lattice::allocate_lattice()\n"
                    "Unable to write to temporary file");
#endif
            ++m_local_volume;
          }
        }

        x[0]++;
        for (mdp_int mu = 0; mu < m_ndim - 1; mu++)
          if (x[mu] >= m_nx[mu])
          {
            x[mu] = 0;
            x[mu + 1]++;
          }
      } while (x[m_ndim - 1] < m_nx[m_ndim - 1]);
    }

    // helper function
    void allocate_basic_arrays()
    {
      // /////////////////////////////////////////////////////////////////
      // Dynamically allocate some other arrays
      // /////////////////////////////////////////////////////////////////

      // local to global mapping table
      m_global_from_local = new mdp_int[m_local_volume];

#ifndef MDP_NO_LG
      m_local_from_global = new mdp_int[m_global_volume];
      std::fill_n(m_local_from_global, m_global_volume, NOWHERE);
#else
      std::string lg_filename(6, '\0');
      std::generate(lg_filename.begin(), lg_filename.end(),
                    []
                    { return "abcdef123456"[rand() % 12]; });

      m_lg_file.open(lg_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);

      if (!m_lg_file)
        error("mdp_lattice::mdp_lattice()\n"
              "Unable to create temporary m_local_from_global file");

      for (mdp_int global_idx = 0; global_idx < m_global_volume; ++global_idx)
      {
        m_lg_file.seekp(global_idx * sizeof(mdp_int), std::ios::beg);
        m_lg_file.write(reinterpret_cast<const char *>(&NOWHERE), sizeof(mdp_int));

        if (!m_lg_file)
          error("mdp_lattice::allocate_lattice()\n"
                "Unable to write to temporary file");
      }
#endif

      m_wh = new mdp_int[m_local_volume];
      m_parity = new int[m_local_volume];

      for (int process = 0; process < Nproc; ++process)
      {
        m_start[process][0] = 0;
        m_start[process][1] = 0;
        m_stop[process][0] = 0;
        m_stop[process][1] = 0;
      }
    }

    // helper function
    void allocate_connectivity_arrays()
    {
      // neighbour and coordinates tables
      m_dw = new mdp_int *[m_local_volume];
      m_up = new mdp_int *[m_local_volume];
      m_co = new mdp_int *[m_local_volume];

      for (mdp_int idx = 0; idx < m_local_volume; ++idx)
      {
        m_dw[idx] = new mdp_int[m_ndir];
        m_up[idx] = new mdp_int[m_ndir];
        m_co[idx] = new mdp_int[m_ndir];
      }
    }

    // helper function
    void build_local_global_maps(int *x, mdp_int *local_mdp_sites)
    {
#ifdef MDP_NO_LG
      mdp_int lms_tmp = 0;
#endif
      for (int process = 0; process < Nproc; process++)
      {
        if (process > 0)
        {
          m_start[process][0] = m_stop[process][0] = m_stop[process - 1][1];
        }

        for (int parity = 0; parity < 2; ++parity)
        {
          if (parity > 0)
          {
            m_start[process][1] = m_stop[process][1] = m_stop[process][0];
          }

          for (mdp_int old_idx = 0; old_idx < m_local_volume; ++old_idx)
          {
#ifndef MDP_NO_LG
            translate_to_coordinates(local_mdp_sites[old_idx], x);
#else
            lms_file.seekg(old_idx * sizeof(mdp_int), std::ios::beg);
            lms_file.read(reinterpret_cast<char *>(&lms_tmp), sizeof(mdp_int));

            if (!lms_file)
              error("mdp_lattice::allocate_lattice()\n"
                    "Unable to read from temporary file");
            translate_to_coordinates(lms_tmp, x);
#endif
            int x_parity = compute_parity(x);
            if (((*m_where)(x, m_ndim, m_nx) == process) && (x_parity == parity))
            {
              mdp_int new_idx = m_stop[process][parity];
#ifndef MDP_NO_LG
              m_local_from_global[local_mdp_sites[old_idx]] = new_idx;
              m_global_from_local[new_idx] = local_mdp_sites[old_idx];
#else

              lms_file.seekg(old_idx * sizeof(mdp_int), std::ios::beg);
              lms_file.read(reinterpret_cast<char *>(&lms_tmp), sizeof(mdp_int));

              if (!lms_file)
                error("mdp_lattice::allocate_lattice()\n"
                      "Unable to read from temporary file");

              m_lg_file.seekp(lms_tmp * sizeof(mdp_int), std::ios::beg);
              m_lg_file.write(reinterpret_cast<const char *>(&new_idx), sizeof(mdp_int));

              if (!m_lg_file)
                error("mdp_lattice::allocate_lattice()\n"
                      "Unable to write to temporary file");
              m_global_from_local[new_idx] = lms_tmp;
#endif
              m_wh[new_idx] = process;
              m_parity[new_idx] = x_parity;
              m_stop[process][parity]++;
            }
          }
        }
      }
    }

    // helper function
    void build_neighbour_tables(int *x, int *x_dw, int *x_up)
    {
      for (mdp_int local_idx = 0; local_idx < m_local_volume; ++local_idx)
      {
        const mdp_int global = m_global_from_local[local_idx];

        translate_to_coordinates(global, x);

        for (mdp_int mu = 0; mu < m_ndim; ++mu)
          m_co[local_idx][mu] = x[mu];

        for (mdp_int mu = 0; mu < m_ndir; ++mu)
        {
          (*m_neighbour)(mu, x_dw, x, x_up, m_ndim, m_nx);

          const mdp_int g_up = global_coordinate(x_up);
          const mdp_int g_dw = global_coordinate(x_dw);

          if (m_wh[local_idx] == ME)
          {
            m_dw[local_idx][mu] = local(g_dw);
            m_up[local_idx][mu] = local(g_up);
          }
          else
          {
            if (local(g_dw) != NOWHERE)
              m_dw[local_idx][mu] = local(g_dw);
            else
              m_dw[local_idx][mu] = NOWHERE;

            if (local(g_up) != NOWHERE)
              m_up[local_idx][mu] = local(g_up);
            else
              m_up[local_idx][mu] = NOWHERE;
          }
        }
      }
    }

    /** @brief share information with other processors
     */
    void communicate_results_to_all_processes()
    {
      mdp_int buffer[2];
      mdp_int length;
      int process;
// #define LATTICE_DEBUG
#ifdef LATTICE_DEBUG
      int process2;
#endif
      mdp_request request;

      // sending length ///////////////////////////
#ifdef LATTICE_DEBUG // debugging code below
      if (Nproc % 2 == 1 || m_where != default_partitioning0)
      {
#endif
        for (int dp = 1; dp < Nproc; dp++)
        {
          process = (ME + dp) % Nproc;
          for (int np = 0; np < 2; np++)
          {
            buffer[np] = m_stop[process][np] - m_start[process][np];
          }
          mdp.put(buffer, 2, process, request);
          process = (ME - dp + Nproc) % Nproc;
          mdp.get(m_len_to_send[process], 2, process);
          mdp.wait(request);
          process = (ME + dp) % Nproc;
          length = m_stop[process][1] - m_start[process][0];
          std::unique_ptr<mdp_int[]> dynamic_buffer = std::make_unique<mdp_int[]>(length);
          for (int idx = 0; idx < length; idx++)
            dynamic_buffer[idx] = m_global_from_local[m_start[process][0] + idx];
          mdp.put(dynamic_buffer.get(), length, process, request);
          process = (ME - dp + Nproc) % Nproc;
          length = m_len_to_send[process][0] + m_len_to_send[process][1];
          m_to_send[process] = new mdp_int[length];
          mdp.get(m_to_send[process], length, process);
          for (int idx = 0; idx < length; idx++)
            m_to_send[process][idx] = local(m_to_send[process][idx]);
          mdp.wait(request);
        }
#ifdef LATTICE_DEBUG // debugging code below
      }
      else
      {
        for (int dp = 1; dp < Nproc; dp++)
        {
          for (int k = 0; k < 2; k++)
          {

            process = (ME + dp) % Nproc;
            process2 = (ME - dp + Nproc) % Nproc;

            if ((k + ME) % 2 == 0)
            {
              for (int np = 0; np < 2; np++)
                buffer[np] = m_stop[process][np] - m_start[process][np];
              mdp.put(buffer, 2, process, request);
              length = m_stop[process][1] - m_start[process][0];
              std::unique_ptr<mdp_int[]> dynamic_buffer = std::make_unique<mdp_int[]>(length);
              for (int idx = 0; idx < length; idx++)
                dynamic_buffer[idx] = m_global_from_local[m_start[process][0] + idx];
              mdp.put(dynamic_buffer.get(), length, process, request);
            }
            else
            {
              mdp.get(m_len_to_send[process2], 2, process2);
              length = m_len_to_send[process2][0] + m_len_to_send[process2][1];
              m_to_send[process2] = new mdp_int[length];
              mdp.get(m_to_send[process2], length, process);
              for (int idx = 0; idx < length; idx++)
                m_to_send[process2][idx] = local(m_to_send[process2][idx]);
            }
          }
        }
      }
#endif
    }

    /** @brief initialize random number generator for each local mdp_site
     */
    void initialize_random(mdp_int random_seed_)
    {
      if (mdp_random_seed_filename && random_seed_ == 0)
      {
        if (isMainProcess())
        {
          {
            std::ifstream fp(mdp_random_seed_filename, std::ios::binary);
            if (fp)
            {
              fp.read(reinterpret_cast<char *>(&random_seed_), sizeof(random_seed_));
              if (!fp)
              {
                random_seed_ = 0;
              }
            }
            mdp << "Reading from file '" << mdp_random_seed_filename << "' lattice().random_seed = " << random_seed_ << "\n";
          }

          {
            std::ofstream fp(mdp_random_seed_filename, std::ios::binary | std::ios::trunc);
            if (fp)
            {
              mdp_int tmp = random_seed_ + 1;
              fp.write(reinterpret_cast<char *>(&tmp), sizeof(tmp));
              mdp << "Writing to   file '" << mdp_random_seed_filename << "' lattice().random_seed = " << tmp << "\n";
            }
          }
        }
        mdp.broadcast(random_seed_, 0);
      }
      else
      {
        mdp << "Adopting random_seed = " << random_seed_ << "\n";
      }

      m_random_seed = random_seed_;
      if (m_local_random_generator)
      {
        mdp << "Using a local random generator\n";
        m_random_obj = std::make_unique<mdp_prng[]>(m_internal_volume);
        for (mdp_int idx = 0; idx < m_internal_volume; idx++)
          m_random_obj[idx].initialize(m_global_from_local[idx + m_start[ME][0]] + m_random_seed);
      }
    }

    /** @brief reallocate a lattice dynamically
     *
     * @param box size of the lattice
     * @param ndir_ number of directions
     * @param where pointer to a partitioning function
     * @param neighbour_ pointer to a topology function
     * @param random_seed_ seed to be used by the parallel prng
     * @param next_next_ size of the buffer between neighbor processes
     * @param local_random_ true is local random generator is required
     */
    void allocate_lattice(const Box &box,
                          int ndir_,
                          int (*where_)(const int[], const int, const int[]) = default_partitioning0,
                          void (*neighbour_)(const int, int[], const int[], int[], const int, const int[]) = torus_topology,
                          mdp_int random_seed_ = 0,
                          int next_next_ = 1,
                          bool local_random_ = true)
    {
      mdp.begin_function("allocate_lattice");

      mdp << "Initializing mdp_lattice...\n";
      init_basic_parameters(box, ndir_, where_, neighbour_, next_next_, local_random_);

#ifndef MDP_NO_LG
      std::unique_ptr<mdp_int[]> local_mdp_sites = std::make_unique<mdp_int[]>(m_global_volume);
#else
      std::string lms_filename(8, '\0');
      std::generate(lms_filename.begin(), lms_filename.end(), []()
                    { return "abcdef123456"[rand() % 12]; });
      std::fstream lms_file;
      lms_file.open(lms_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);

      if (!lms_file)
        error("mdp_lattice::mdp_lattice()\n"
              "Unable to create temporary lms file");
#endif

      int x[MAX_DIM];
      int x_up[MAX_DIM];
      int x_dw[MAX_DIM];
#ifndef MDP_NO_LG
      compute_local_sites(x, x_up, x_dw, local_mdp_sites.get());
#else
      compute_local_sites(x, x_up, x_dw, lms_file);
#endif

      allocate_basic_arrays();
      allocate_connectivity_arrays();
      build_local_global_maps(x, local_mdp_sites.get());
      build_neighbour_tables(x, x_dw, x_up);

      m_internal_volume = m_stop[ME][1] - m_start[ME][0];

      mdp << "Communicating...\n";
      communicate_results_to_all_processes();

      mdp << "Initializing random per mdp_site...\n";
      initialize_random(random_seed_);

      mdp << "Lattice created.\n";
      mdp.end_function("allocate_lattice");
    }

    /** @brief dynamically deallocate a lattice
     */
    void deallocate_memory()
    {
      if (m_local_volume == 0)
        return;

      delete[] m_nx;

      for (int process = 0; process < Nproc; process++)
      {
        if (process != ME)
        {
          if (m_len_to_send[process][0] + m_len_to_send[process][1] != 0)
            delete[] m_to_send[process];
        }
      }

      for (int new_idx = 0; new_idx < m_local_volume; ++new_idx)
      {
        delete[] m_dw[new_idx];
        delete[] m_up[new_idx];
        delete[] m_co[new_idx];
      }

      delete[] m_dw;
      delete[] m_up;
      delete[] m_co;
      delete[] m_wh;
      delete[] m_global_from_local;
#ifndef MDP_NO_LG
      delete[] m_local_from_global;
#else
      m_lg_file.close();
#endif
      delete[] m_parity;
      // delete[] m_random_obj;
    }

  public:
    /** Which process point x[] is in?
     *
     * @param x Point x[] to be inspected
     * @return Process ID
     */
    int where(const int x[]) const
    {
      return (*m_where)(x, m_ndim, m_nx);
    }

    mdp_int process(mdp_int local_idx) const
    {
      return m_wh[local_idx];
    }

    /** mu-th coordinate for point x[] based on local index
     *
     * @param local_idx Local index for point x[] to be inspected
     * @return Value of mu-th coordinate
     */
    mdp_int coordinate(mdp_int local_idx, mdp_int mu) const
    {
      return m_co[local_idx][mu];
    }

    /** @brief Calculate global coordinate
     *
     * Calculate ordinal (global) coordinate of
     * a n-dimensional point x[].
     */
    mdp_int global_coordinate(const int x[])
    {
      mdp_int global_idx = 0;
      for (mdp_int mu = 0; mu < m_ndim - 1; mu++)
        global_idx = (global_idx + x[mu]) * m_nx[mu + 1];
      return global_idx + x[m_ndim - 1];
    }

    /** @brief Translate global coordinate into x[] coordinates
     *
     * Given the ordinal coordinate, recover
     * x = (x0, x1, ... x9) coordinates.
     */
    void translate_to_coordinates(mdp_int global_idx, int x[])
    {
      for (mdp_int mu = m_ndim - 1; mu > 0; mu--)
      {
        x[mu] = global_idx % m_nx[mu];
        global_idx = (global_idx - x[mu]) / m_nx[mu];
      }
      x[0] = global_idx;
    }

    /** Which process point x[] is in based on its global index?
     *
     * @param global_idx Global index for point x[] to be inspected
     * @return Process ID
     */
    mdp_int where_global(mdp_int global_idx)
    {
      int x[MAX_DIM];
      translate_to_coordinates(global_idx, x);

      return (*m_where)(x, m_ndim, m_nx);
    }

    /** @brief Calculate parity of point x[]
     *
     * Check the parity of the sum of x coordinates.
     */
    int compute_parity(const int x[])
    {
      int p = 0;
      for (mdp_int mu = 0; mu < m_ndim; mu++)
        p = p + x[mu];
      return (p % 2);
    }

    // base constructor
    mdp_lattice() : m_local_volume(0), m_random_obj(nullptr)
    {
    }

    /** @brief declares a lattice object
     *
     * @param box size of the lattice
     * @param ndir_ number of directions
     * @param where_ pointer to a partitioning function
     * @param neighbour_ pointer to a topology function
     * @param random_seed_ seed to be used by the parallel prng
     * @param next_next_ size of the buffer between neighbour processes
     * @param local_random_ true is local random generator is required
     */
    mdp_lattice(const Box &box,
                int ndir_,
                int (*where_)(const int[], const int, const int[]) = default_partitioning0,
                void (*neighbour_)(const int, int[], const int[], int[], const int, const int[]) = torus_topology,
                mdp_int random_seed_ = 0,
                int next_next_ = 1,
                bool local_random_ = true) : mdp_lattice()
    {
      if (ndir_ < 0)
        ndir_ = box.dim();
      allocate_lattice(box, ndir_, where_, neighbour_,
                       random_seed_, next_next_, local_random_);
    }

    /** @brief declares a lattice object
     *
     * @param box size of the lattice
     * @param where_ pointer to a partitioning function
     * @param neighbour_ pointer to a topology function
     * @param random_seed_ seed to be used by the parallel prng
     * @param next_next_ size of the buffer between neighbour processes
     * @param local_random_ true is local random generator is required
     */
    mdp_lattice(const Box &box,
                int (*where_)(const int[], const int, const int[]) = default_partitioning0,
                void (*neighbour_)(const int, int[], const int[], int[], const int, const int[]) = torus_topology,
                mdp_int random_seed_ = 0,
                int next_next_ = 1,
                bool local_random_ = true) : mdp_lattice(box, box.dim(), where_, neighbour_, random_seed_, next_next_, local_random_)
    {
    }

    /** @brief deallocate all the dynamically allocated arrays
     */
    virtual ~mdp_lattice()
    {
      deallocate_memory();
    }

    mdp_prng &random(mdp_site);

    // //////////////////////////////////
    // functions for external access ...
    // to be used to access member variables
    // /////////////////////////////////

    /** @brief number of dimensions of the lattice
     */
    mdp_int n_dimensions() const
    {
      return m_ndim;
    }

    mdp_int ndim() const
    {
      return m_ndim;
    }

    /** @brief number of directions one can move on the lattice; usually same as ndim
     */
    mdp_int n_directions() const
    {
      return m_ndir;
    }

    /** @brief size of the lattice in direction mu
     */
    mdp_int size(const int mu) const
    {
      return m_nx[mu];
    }

    /** @brief Size of each dimension
     *
     * @return Pointer to the array of sizes
     */
    const mdp_int *dims() const
    {
      return m_nx;
    }

    /** @brief number of lattice sites stored locally by current process
     */
    mdp_int local_volume() const
    {
      return m_internal_volume;
    }

    /** @brief number of lattice sites stored locally by current process
     * together with the number of copies of the neigbouring sites
     */
    mdp_int enclosing_volume() const
    {
      return m_local_volume;
    }

    /** @brief total lattice volume
     */
    mdp_int global_volume() const
    {
      return m_global_volume;
    }

    /** @brief number of sites of the lattice
     */
    mdp_int size() const
    {
      return m_global_volume;
    }

    mdp_int boundary_thickness() const
    {
      return m_next_next;
    }

    mdp_int move_up(const mdp_int idx, const int mu) const
    {
      return m_up[idx][mu];
    }

    mdp_int move_down(const mdp_int idx, const int mu) const
    {
      return m_dw[idx][mu];
    }

#ifdef SSE2
    mdp_int **up()
    {
      return m_up;
    }

    mdp_int **down()
    {
      return m_dw;
    }
#endif

    /** @brief Local index based on global index
     */
    mdp_int local(mdp_int global_idx) const
    {
#ifndef MDP_NO_LG
      return m_local_from_global[global_idx];
#else
      mdp_int lg_tmp;
      m_lg_file.seekg(global_idx * sizeof(mdp_int), std::ios::beg);
      m_lg_file.read(reinterpret_cast<char *>(&lg_tmp), sizeof(mdp_int));

      if (!m_lg_file)
        error("mdp_lattice::allocate_lattice()\n"
              "Unable to read from temporary file");

      return lg_tmp;
#endif
    }

    /** @brief Global index based on local index
     */
    mdp_int global(mdp_int local_idx) const
    {
      return m_global_from_local[local_idx];
    }

    int site_parity(const mdp_int idx) const
    {
      return m_parity[idx];
    }

    mdp_int start_index(const int process_id, int parity = EVENODD) const
    {
      if (parity == EVENODD)
        parity = 0;
      return m_start[process_id][parity];
    }

    mdp_int stop_index(const int process_id, int parity = EVENODD) const
    {
      if (parity == EVENODD)
        parity = 1;
      return m_stop[process_id][parity];
    }

    mdp_int start0(mdp_int process_id, mdp_int parity) const
    {
      return m_start[process_id][parity];
    }

    mdp_int stop0(mdp_int process_id, mdp_int parity) const
    {
      return m_stop[process_id][parity];
    }

    mdp_int len_to_send0(mdp_int process_id, mdp_int parity) const
    {
      return m_len_to_send[process_id][parity];
    }

    mdp_int to_send0(mdp_int process_id, mdp_int local_idx) const
    {
      return m_to_send[process_id][local_idx];
    }
  };
} // namespace MDP

#endif /* MDP_LATTICE_ */
