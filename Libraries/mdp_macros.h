/////////////////////////////////////////////////////////////////
/// @file mdp_macros.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains various mdp macros
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_MACROS_
#define MDP_MACROS_

namespace MDP
{
/// Loop on all local sites of this process
#define forallsites(x) \
     for (x.start(); x.is_in(); x.next())

/// Loop on all local sites of this process with given parity
/// If pofx is EVENODD=2 then loops on even and odd sites
#define forallsitesofparity(x, pofx)                                        \
     for (x.start(), x.set_local(x.lattice().start0(ME, (pofx % 2)));       \
          x.local_index() < x.lattice().stop0(ME, (pofx + (pofx % 2)) / 2); \
          x.next())

/// Loop on all sites stored by this process
#define forallsitesandcopies(x) \
     for (x.start(), x.set_local(0); x.local_index() < x.lattice().enclosing_volume(); x.next())

/// Loop on all sites stored by this process with given parity
// if pofx is EVENODD=2 then loops on even and odd sites
#define forallsitesandcopiesofparity(x, pofx)                                           \
     for (int __process = 0; __process < Nproc; __process++)                            \
          for (x.start(), x.set_local(x.lattice().start0(__process, (pofx % 2)));       \
               x.local_index() < x.lattice().stop0(__process, (pofx + (pofx % 2)) / 2); \
               x.next())

/// Returns the unique id of this process
#define ME mdp.me()

/// Check for Sub Process
#define isSubProcess(i) (mdp.me() == (i))

/// Check for Parent (Main) Process
#define isMainProcess() (isSubProcess(0))

/// Returns the total number of parallel processes for this job
#define Nproc mdp.nproc()

/// Reports a runtime error and the line that caused it
#define error(a) _mpi_error_message(a, __FILE__, __LINE__);
} // namespace MDP

#endif /* MDP_MACROS_ */
