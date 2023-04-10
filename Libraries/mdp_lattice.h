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

#include <string>

namespace MDP
{
  const mdp_int NOWHERE = INT_MAX;

  class mdp_site;

  /// @brief distributed lattice object
  ///
  /// Example:
  /// @verbatim
  ///    int box[]={3,3,3};
  ///    int seed=0, border_width=1;
  ///    mdp_lattice lattice(3,box,default_partitioning0,
  ///                        torus_topology,seed,border_width);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      cout << lattice.random(x).plain() << endl;
  /// @endverbatim
  class mdp_lattice
  {
  public:
    int ndim;        /* number of dimensions       */
    int ndir;        /* number of directions       */
    int next_next;   /* 1 or 2 is the thikness of the boudary */
    int *nx;         /* box containing the lattice */
    mdp_int nvol;    /* local volume               */
    mdp_int nvol_gl; /* global volume              */
    mdp_int nvol_in; /* internal volume            */
    mdp_int *gl;     /* map local to global        */
    mdp_int *lg;     /* map global to local        */
    FILE *lg_file;   /* temporary file to store lg if not enough memory */
    mdp_int **up;    /* move up in local index     */
    mdp_int **dw;    /* move dw in local index     */
    int **co;        /* coordinate x in local idx  */
    int *wh;         /* in which process? is idx   */
    int *parity;     /* parity of local mdp_site idx   */
    mdp_int start[_NprocMax_][2];
    mdp_int stop[_NprocMax_][2];
    mdp_int len_to_send[_NprocMax_][2];
    mdp_int *to_send[_NprocMax_];

  private:
    bool local_random_generator;
    mdp_prng *random_obj;
    mdp_int random_seed;
    // ////////////////////////////////////////////////////
    // share information with other processors
    // ////////////////////////////////////////////////////
    void communicate_results_to_all_processes()
    {
      mdp_int buffer[2];
      mdp_int *dynamic_buffer;
      mdp_int length;
      int dp, process, np, idx;
#if 0
    int process2;
#endif
      mdp_request request;

      // sending length ///////////////////////////
#if 0 // debugging code below
    if (Nproc % 2 == 1 || m_where != default_partitioning0)
    {
#endif
      for (dp = 1; dp < Nproc; dp++)
      {
        process = (ME + dp) % Nproc;
        for (np = 0; np < 2; np++)
          buffer[np] = stop[process][np] - start[process][np];
        mpi.put(buffer, 2, process, request);
        process = (ME - dp + Nproc) % Nproc;
        mpi.get(len_to_send[process], 2, process);
        mpi.wait(request);
        process = (ME + dp) % Nproc;
        length = stop[process][1] - start[process][0];
        dynamic_buffer = new mdp_int[length];
        for (idx = 0; idx < length; idx++)
          dynamic_buffer[idx] = gl[start[process][0] + idx];
        mpi.put(dynamic_buffer, length, process, request);
        process = (ME - dp + Nproc) % Nproc;
        length = len_to_send[process][0] + len_to_send[process][1];
        to_send[process] = new mdp_int[length];
        mpi.get(to_send[process], length, process);
        for (idx = 0; idx < length; idx++)
          to_send[process][idx] = local(to_send[process][idx]);
        mpi.wait(request);
        delete[] dynamic_buffer;
      }
#if 0 // debugging code below
    }
    else
    {
      for (dp = 1; dp < Nproc; dp++)
      {
        for (int k = 0; k < 2; k++)
        {

          process = (ME + dp) % Nproc;
          process2 = (ME - dp + Nproc) % Nproc;

          if ((k + ME) % 2 == 0)
          {
            for (np = 0; np < 2; np++)
              buffer[np] = stop[process][np] - start[process][np];
            mpi.put(buffer, 2, process, request);
            length = stop[process][1] - start[process][0];
            dynamic_buffer = new mdp_int[length];
            for (idx = 0; idx < length; idx++)
              dynamic_buffer[idx] = gl[start[process][0] + idx];
            mpi.put(dynamic_buffer, length, process, request);
          }
          else
          {
            mpi.get(len_to_send[process2], 2, process2);
            length = len_to_send[process2][0] + len_to_send[process2][1];
            to_send[process2] = new mdp_int[length];
            mpi.get(to_send[process2], length, process);
            for (idx = 0; idx < length; idx++)
              to_send[process2][idx] = local(to_send[process2][idx]);
          }
        }
        delete[] dynamic_buffer;
      }
    }
#endif
    }

    int (*m_where)(const int *, const int, const int *);
    void (*m_neighbour)(const int mu, int *x_dw, const int *x, int *x_up, const int ndim, const int *nx);

  public:
    auto &where() const
    {
      return m_where;
    }

    /** @brief Calculate global coordinate
     *
     * Calculate ordinal (global) coordinate of
     * a n-dimentional point x[].
     */
    mdp_int global_coordinate(int *x)
    {
      mdp_int global_idx = 0;
      int mu;
      for (mu = 0; mu < ndim - 1; mu++)
        global_idx = (global_idx + x[mu]) * nx[mu + 1];
      return global_idx + x[ndim - 1];
    }

    /** @brief Translate global coordinate into x[] coordinates
     *
     * Given the ordinal coordinate, recover
     * x = (x0, x1, ... x9) coordinates.
     */
    void global_coordinate(mdp_int global_idx, int *x)
    {
      for (int mu = ndim - 1; mu > 0; mu--)
      {
        x[mu] = global_idx % nx[mu];
        global_idx = (global_idx - x[mu]) / nx[mu];
      }
      x[0] = global_idx;
    }

    int compute_parity(int *x)
    {
      int mu = 0;
      int p = 0;
      for (mu = 0; mu < ndim; mu++)
        p = p + x[mu];
      return (p % 2);
    }

    mdp_lattice()
    {
      nvol = 0;
      m_random_obj = nullptr;
    }

    /// declares a lattice object
    /// @param ndim_ dimensions of the lattice
    /// @param nx_ size of the lattice
    /// @param where pointer to a partitioning function
    /// @param neighbour_ pointer to a topology function.
    /// @param random_seed_ seed to be used by the parallel prng
    /// @param next_next_ size of the buffer between neighbour processes
    /// @param local_random_ true is local random generator is required
    mdp_lattice(int ndim_,
                int nx_[],
                int (*where_)(const int *, const int, const int *) = default_partitioning0,
                void (*neighbour_)(const int, int *, const int *, int *, const int, const int *) = torus_topology,
                mdp_int random_seed_ = 0,
                int next_next_ = 1,
                bool local_random_ = true)
    {
      nvol = 0;
      m_random_obj = nullptr;
      allocate_lattice(ndim_, ndim_, nx_, where_, neighbour_,
                       random_seed_, next_next_, local_random_);
    }

    /// for weird stuff
    mdp_lattice(int ndim_,
                int ndir_,
                int nx_[],
                int (*where_)(const int *, const int, const int *) = default_partitioning0,
                void (*neighbour_)(const int, int *, const int *, int *, const int, const int *) = torus_topology,
                mdp_int random_seed_ = 0,
                int next_next_ = 1,
                bool local_random_ = true)
    {
      nvol = 0;
      m_random_obj = nullptr;
      allocate_lattice(ndim_, ndir_, nx_, where_, neighbour_,
                       random_seed_, next_next_, local_random_);
    }

    /// reallocate a lattice dynamically
    /// @param ndim_ dimensions of the lattice
    /// @param nx_ size of the lattice
    /// @param where pointer to a partitioning function
    /// @param neighbour_ pointer to a topology function.
    /// @param random_seed_ seed to be used by the parallel prng
    /// @param next_next_ size of the buffer between neighbour processes
    /// @param local_random_ true is local random generator is required
    void allocate_lattice(int ndim_,
                          int nx_[],
                          int (*where_)(const int *, const int, const int *) = default_partitioning0,
                          void (*neighbour_)(const int, int *, const int *, int *, const int, const int *) = torus_topology,
                          mdp_int random_seed_ = 0,
                          int next_next_ = 1,
                          bool local_random_ = true)
    {
      allocate_lattice(ndim_, ndim_, nx_, where_, neighbour_,
                       random_seed_, next_next_, local_random_);
    }

    /// for weird stuff
    void allocate_lattice(int ndim_,
                          int ndir_,
                          int nx_[],
                          int (*where_)(const int *, const int, const int *) = default_partitioning0,
                          void (*neighbour_)(const int, int *, const int *, int *, const int, const int *) = torus_topology,
                          mdp_int random_seed_ = 0,
                          int next_next_ = 1,
                          bool local_random_ = true)
    {
      mpi.begin_function("allocate_lattice");
      local_random_generator = local_random_;
      if (ndim_ != ndir_)
        mpi << "It is getting complicated: you have ndim!=ndir\n";
      deallocate_memory();
      // //////////////////////////////////////////////////////////////////
      // (*where)(x,nc) must return the processor rank of where x is stored
      // (*neghbour)(mu,x_dw,x,x_up,ndim, nx) must fill x_dw and x_up
      //            according with current position x and direction mu
      // //////////////////////////////////////////////////////////////////
      mpi << "Initializing a mdp_lattice...\n";

      int mu, nu, rho, process, np;
      bool is_boundary;
      mdp_int global_idx, old_idx, new_idx;
      ndim = ndim_;
      ndir = ndir_;
      m_where = where_;
      m_neighbour = neighbour_;
      next_next = next_next_;
      nvol = 0;
      nvol_gl = 1;

      mpi << "Lattice dimension: " << nx_[0];
      for (mu = 1; mu < ndim; mu++)
        mpi << " x " << nx_[mu];
      mpi << "\n";

      ///////////////////////////////////////////////////////////////////
      // Dynamically allocate some arrays
      ///////////////////////////////////////////////////////////////////
      nx = new int[ndim];

      for (mu = 0; mu < ndim; mu++)
        nvol_gl *= (nx[mu] = nx_[mu]);
      int x[10];
      int x_up[10], x_up_dw[10], x_up_up[10], x_up_up_dw[10], x_up_up_up[10];
      int x_dw[10], x_dw_up[10], x_dw_dw[10], x_dw_dw_dw[10], x_dw_dw_up[10];
#if !defined(MDP_NO_LG)
      mdp_int *local_mdp_sites = new mdp_int[nvol_gl];
#else
      mdp_int lms_tmp = 0;
      FILE *lms_file = tmpfile();
      if (lms_file == 0)
        error("mdp_lattice::mdp_lattice()\n"
              "Unable to create temporary lms file");
#endif
      // ///////////////////////////////////////////////////////////////////
      // Fill table local_mdp_sites with the global coordinate of mdp_sites in ME
      // nvol counts the mdp_sites stored in local_mdp_sites
      // ///////////////////////////////////////////////////////////////////
      //
      // cases of interest...:
      //
      // next_next < 0 => only local mdp_sites NEW!
      // next_next = 0 => only x+mu and x-mu (wilson) NEW!
      // next_next = 1 => only x+-mu, x+-mu+-nu mu!=nu (clover)
      // next_next = 2 => only x+-mu, x+-mu+-nu
      // next_next = 3 => only x+-mu, x+-mu+-nu, x+-mu+-nu+-rho (asqtad)
      //
      // ///////////////////////////////////////////////////////////////////

      for (mu = 0; mu < ndim; mu++)
        x[mu] = 0;

      do
      {
        global_idx = global_coordinate(x);
        if ((*m_where)(x, ndim, nx) == ME)
        {
#if !defined(MDP_NO_LG)
          local_mdp_sites[nvol] = global_idx;
#else
          if (fseek(lms_file, nvol * sizeof(mdp_int), SEEK_SET) != 0 ||
              fwrite(&global_idx, sizeof(mdp_int), 1, lms_file) != 1)
            error("mdp_lattice::allocate_lattice()\n"
                  "Unable to write to temporary file");
#endif
          nvol++;
        }
        else if (next_next >= 0)
        {
          is_boundary = false;
          for (mu = 0; mu < ndir; mu++)
          {
            (*m_neighbour)(mu, x_dw, x, x_up, ndim, nx);
            if (((*m_where)(x_up, ndim, nx) >= Nproc) ||
                ((*m_where)(x_dw, ndim, nx) >= Nproc))
              error("Incorrect patitioning");
            if (((*m_where)(x_up, ndim, nx) == ME) ||
                ((*m_where)(x_dw, ndim, nx) == ME))
              is_boundary = true;
            // ////////////////////////////////////////////////////
            // cases:
            // 1) mu+nu
            // 2) mu+nu and mu+mu
            // 3) mu+nu+rho and mu+mu+rho and mu+mu+mu
            // One may want to optimize case 3. Is was thought for
            // improved staggered fermions.
            // ////////////////////////////////////////////////////
            for (nu = 0; nu < ndir; nu++)
              if ((nu != mu) || (next_next > 1))
              {
                (*m_neighbour)(nu, x_dw_dw, x_dw, x_dw_up, ndim, nx);
                (*m_neighbour)(nu, x_up_dw, x_up, x_up_up, ndim, nx);
                if (((*m_where)(x_up_dw, ndim, nx) >= Nproc) ||
                    ((*m_where)(x_up_up, ndim, nx) >= Nproc) ||
                    ((*m_where)(x_dw_dw, ndim, nx) >= Nproc) ||
                    ((*m_where)(x_dw_up, ndim, nx) >= Nproc))
                  error("Incorrect patitioning");
                if (((*m_where)(x_up_dw, ndim, nx) == ME) ||
                    ((*m_where)(x_up_up, ndim, nx) == ME) ||
                    ((*m_where)(x_dw_dw, ndim, nx) == ME) ||
                    ((*m_where)(x_dw_up, ndim, nx) == ME))
                  is_boundary = true;
                if (next_next == 3)
                  for (rho = 0; rho < ndir; rho++)
                  {
                    (*m_neighbour)(rho, x_dw_dw_dw, x_dw_dw, x_dw_dw_up, ndim, nx);
                    (*m_neighbour)(rho, x_up_up_dw, x_up_up, x_up_up_up, ndim, nx);
                    if (((*m_where)(x_up_up_up, ndim, nx) >= Nproc) ||
                        ((*m_where)(x_up_up_dw, ndim, nx) >= Nproc) ||
                        ((*m_where)(x_dw_dw_up, ndim, nx) >= Nproc) ||
                        ((*m_where)(x_dw_dw_dw, ndim, nx) >= Nproc))
                      error("Incorrect patitioning");
                    if (((*m_where)(x_up_up_up, ndim, nx) == ME) ||
                        ((*m_where)(x_up_up_dw, ndim, nx) == ME) ||
                        ((*m_where)(x_dw_dw_up, ndim, nx) == ME) ||
                        ((*m_where)(x_dw_dw_dw, ndim, nx) == ME))
                      is_boundary = true;
                  }
              }
          }
          if (is_boundary == true)
          {
#if !defined(MDP_NO_LG)
            local_mdp_sites[nvol] = global_idx;
#else
            if (fseek(lms_file, nvol * sizeof(mdp_int), SEEK_SET) != 0 ||
                fwrite(&global_idx, sizeof(mdp_int), 1, lms_file) != 1)
              error("mdp_lattice::allocate_lattice()\n"
                    "Unable to write to temporary file");
#endif
            nvol++;
          }
        }
        x[0]++;
        for (mu = 0; mu < ndim - 1; mu++)
          if (x[mu] >= nx[mu])
          {
            x[mu] = 0;
            x[mu + 1]++;
          }
      } while (x[ndim - 1] < nx[ndim - 1]);

      // /////////////////////////////////////////////////////////////////
      // Dynamically allocate some other arrays
      // /////////////////////////////////////////////////////////////////

      dw = new mdp_int *[nvol];
      up = new mdp_int *[nvol];
      co = new int *[nvol];
      for (new_idx = 0; new_idx < nvol; new_idx++)
      {
        dw[new_idx] = new mdp_int[ndir];
        up[new_idx] = new mdp_int[ndir];
        co[new_idx] = new int[ndir];
      }
      gl = new mdp_int[nvol];
#if !defined(MDP_NO_LG)
      lg = new mdp_int[nvol_gl];
#else
      lg_file = tmpfile();
      if (lg_file == 0)
        error("mdp_lattice::mdp_lattice()\n"
              "Unable to create temporary lg file");
#endif
      wh = new int[nvol];
      for (global_idx = 0; global_idx < nvol_gl; global_idx++)
#if !defined(MDP_NO_LG)
        lg[global_idx] = NOWHERE;
#else
        if (fseek(lg_file, global_idx * sizeof(mdp_int), SEEK_SET) != 0 ||
            fwrite(&NOWHERE, sizeof(mdp_int), 1, lg_file) != 1)
          error("mdp_lattice::allocate_lattice()\n"
                "Unable to write to temporary file");
#endif
      parity = new int[nvol];
      // /////////////////////////////////////////////////////////////////
      start[0][0] = stop[0][0] = 0;
      for (process = 0; process < Nproc; process++)
      {
        if (process > 0)
          start[process][0] = stop[process][0] = stop[process - 1][1];
        for (np = 0; np < 2; np++)
        {
          if (np > 0)
            start[process][1] = stop[process][1] = stop[process][0];
          for (old_idx = 0; old_idx < nvol; old_idx++)
          {
#if !defined(MDP_NO_LG)
            global_coordinate(local_mdp_sites[old_idx], x);
#else
            if (fseek(lms_file, old_idx * sizeof(mdp_int), SEEK_SET) != 0 ||
                fread(&lms_tmp, sizeof(mdp_int), 1, lms_file) != 1)
              error("mdp_lattice::allocate_lattice()\n"
                    "Unable to read to temporary file");
            global_coordinate(lms_tmp, x);
#endif
            if (((*m_where)(x, ndim, nx) == process) && (compute_parity(x) == np))
            {
              new_idx = stop[process][np];
#if !defined(MDP_NO_LG)
              lg[local_mdp_sites[old_idx]] = new_idx;
              gl[new_idx] = local_mdp_sites[old_idx];
#else

              if (fseek(lms_file, old_idx * sizeof(mdp_int), SEEK_SET) != 0 ||
                  fread(&lms_tmp, sizeof(mdp_int), 1, lms_file) != 1)
                error("mdp_lattice::allocate_lattice()\n"
                      "Unable to read to temporary file");
              if (fseek(lg_file, lms_tmp * sizeof(mdp_int), SEEK_SET) != 0 ||
                  fwrite(&new_idx, sizeof(mdp_int), 1, lg_file) != 1)
                error("mdp_lattice::allocate_lattice()\n"
                      "Unable to write to temporary file");
              gl[new_idx] = lms_tmp;
#endif
              wh[new_idx] = process;
              parity[new_idx] = compute_parity(x);
              stop[process][np]++;
            }
          }
        }
      }
      // deallcate temporary array
#if !defined(MDP_NO_LG)
      delete[] local_mdp_sites;
#else
      fclose(lms_file);
#endif
      // /////////////////////////
      for (new_idx = 0; new_idx < nvol; new_idx++)
      {
        global_coordinate(gl[new_idx], x);
        for (mu = 0; mu < ndim; mu++)
          co[new_idx][mu] = x[mu];
        for (mu = 0; mu < ndir; mu++)
        {
          (*m_neighbour)(mu, x_dw, x, x_up, ndim, nx);
          if (wh[new_idx] == ME)
          {
            dw[new_idx][mu] = local(global_coordinate(x_dw));
            up[new_idx][mu] = local(global_coordinate(x_up));
          }
          else
          {
            if (local(global_coordinate(x_dw)) != NOWHERE)
              dw[new_idx][mu] = local(global_coordinate(x_dw));
            else
              dw[new_idx][mu] = NOWHERE;
            if (local(global_coordinate(x_up)) != NOWHERE)
              up[new_idx][mu] = local(global_coordinate(x_up));
            else
              up[new_idx][mu] = NOWHERE;
          }
        }
      }
      nvol_in = stop[ME][1] - start[ME][0];
      mpi << "Communicating...\n";
      communicate_results_to_all_processes();
      mpi << "Initializing random per mdp_site...\n";
      if (mdp_random_seed_filename && random_seed_ == 0)
      {
        if (ME == 0)
        {
          FILE *fp = fopen(mdp_random_seed_filename, "r");
          if (fp != nullptr)
          {
            if (fread(&random_seed_, sizeof(random_seed_), 1, fp) != 1)
              random_seed_ = 0;
            fclose(fp);
          };
          fp = fopen(mdp_random_seed_filename, "w");
          if (fp != nullptr)
          {
            random_seed_ += 1;
            fwrite(&random_seed_, sizeof(random_seed_), 1, fp);
            random_seed_ -= 1;
            fclose(fp);
          }
          mpi << "Reading from file " << mdp_random_seed_filename << " lattice().random_seed=" << random_seed_ << "\n";
          mpi << "Writing to   file " << mdp_random_seed_filename << " lattice().random_seed=" << random_seed_ + 1 << "\n";
        }
        mpi.broadcast(random_seed_, 0);
      }
      else
      {
        mpi << "Adopting random_seed=" << random_seed_ << "\n";
      }
      initialize_random(random_seed_);
      mpi << "Lattice created.\n";
      mpi.end_function("allocate_lattice");
    }

    /** @brief deallocate all the dynamically allocated arrays
     */
    virtual ~mdp_lattice()
    {
      deallocate_memory();
    }

    /** @brief dynamically deallocate a lattice
     */
    void deallocate_memory()
    {
      if (nvol == 0)
        return;

      int process, new_idx;
      delete[] nx;

      for (process = 0; process < Nproc; process++)
        if (process != ME)
        {
          if (len_to_send[process][0] + len_to_send[process][1] != 0)
            delete[] to_send[process];
        }

      for (new_idx = 0; new_idx < nvol; new_idx++)
      {
        delete[] dw[new_idx];
        delete[] up[new_idx];
        delete[] co[new_idx];
      }

      delete[] dw;
      delete[] up;
      delete[] co;
      delete[] wh;
      delete[] gl;
#ifndef MDP_NO_LG
      delete[] lg;
#else
      fclose(lg_file);
#endif
      delete[] parity;
      delete[] random_obj;
    }

    /** @brief initialize random number generator for each local mdp_site
     */
    void initialize_random(mdp_int random_seed_ = 0)
    {
      random_seed = random_seed_;
      if (local_random_generator)
      {
        mdp << "Using a local random generator\n";
        if (random_obj != 0)
          delete[] random_obj;
        random_obj = new mdp_prng[nvol_in];
        for (mdp_int idx = 0; idx < nvol_in; idx++)
          random_obj[idx].initialize(gl[idx + start[ME][0]] + random_seed);
      }
    }

    mdp_prng &random(mdp_site);
    // //////////////////////////////////
    // functions for external access ...
    // to be used to access member variables
    // /////////////////////////////////

    /** @brief number of dimensions of the lattice
     */
    int n_dimensions() const
    {
      return ndim;
    }

    /** @brief number of directions one can move on the lattice; usually same as ndim
     */
    int n_directions() const
    {
      return ndir;
    }

    /** @brief size of the lattice in direction mu
     */
    mdp_int size(const int mu) const
    {
      return nx[mu];
    }

    /** @brief Size of each dimension
     *
     * @return Pointer to the array of sizes
     */
    const mdp_int *dims() const
    {
      return nx;
    }

    /** @brief number of lattice sites stored locally by current process
     */
    mdp_int local_volume() const
    {
      return nvol_in;
    }

    /** @brief number of lattice sites stored locally by current process
     * together with the number of copies of the neigbouring sites
     */
    mdp_int enclosing_volume() const
    {
      return nvol;
    }

    /** @brief total lattice volume
     */
    mdp_int global_volume() const
    {
      return nvol_gl;
    }

    /** @brief number of sites of the lattice
     */
    mdp_int size() const
    {
      return nvol_gl;
    }

    mdp_int boundary_thickness() const
    {
      return next_next;
    }

    mdp_int move_up(const mdp_int idx, const int mu) const
    {
      return up[idx][mu];
    }

    mdp_int move_down(const mdp_int idx, const int mu) const
    {
      return dw[idx][mu];
    }

    mdp_int local(mdp_int idx) const
    {
#ifndef MDP_NO_LG
      return lg[idx];
#else
      mdp_int lg_tmp;
      if (fseek(lg_file, idx * sizeof(mdp_int), SEEK_SET) != 0 ||
          fread(&lg_tmp, sizeof(mdp_int), 1, lg_file) != 1)
        error("mdp_lattice::allocate_lattice()\n"
              "Unable to read to temporary file");
      return lg_tmp;
#endif
    }

    mdp_int global(mdp_int idx) const
    {
      return gl[idx];
    }

    int site_parity(const mdp_int idx) const
    {
      return parity[idx];
    }

    mdp_int start_index(const int process, int p = EVENODD) const
    {
      if (p == EVENODD)
        p = 0;
      return start[process][p];
    }

    mdp_int stop_index(const int process, int p = EVENODD) const
    {
      if (p == EVENODD)
        p = 1;
      return stop[process][p];
    }
  };
} // namespace MDP

#endif /* MDP_LATTICE_ */
