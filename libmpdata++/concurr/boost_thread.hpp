/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/concurr/detail/concurr_common.hpp>

// TODO: make it work with clang as well!
#if !defined(_REENTRANT)
#  error _REENTRANT not defined, please use something like -pthread flag for gcc
#endif
#include <boost/thread.hpp>

#include <cstdlib> // std::getenv()

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      std::map<boost::thread::id, int> boost_thread_id;
    };

    template <
      class solver_t,
      bcond::bcond_e bcx,
      bcond::bcond_e bcy = bcond::null,
      bcond::bcond_e bcz = bcond::null
    >
    class boost_thread : public detail::concurr_common<solver_t, bcx, bcy, bcz>
    {
      using parent_t = detail::concurr_common<solver_t, bcx, bcy, bcz>;
 
      class mem_t : public parent_t::mem_t 
      {
        boost::barrier b;

        public:

        int rank()
        {
          return detail::boost_thread_id[boost::this_thread::get_id()];
        }

	static int size() 
	{
	  const char *env_var("OMP_NUM_THREADS");
	  return (std::getenv(env_var) != NULL)
	    ? std::atoi(std::getenv(env_var)) // TODO: check if convesion OK?
	    : boost::thread::hardware_concurrency();
	}


        // TODO:  could one constructor be enough if n_dims a template param?
        mem_t(int s0) : 
          b(size()),
          parent_t::mem_t(s0, size()) 
        {}; 

        mem_t(int s0, int s1) : 
          b(size()),
          parent_t::mem_t(s0, s1, size()) 
        {}; 

        mem_t(int s0, int s1, int s2) : 
          b(size()),
          parent_t::mem_t(s0, s1, s2, size()) 
        {}; 

	void barrier()
	{
	  b.wait();
	}
      };

      public:

      void solve(int nt)
      {
        boost::thread_group threads; // TODO: member field?
        for (int i = 0; i < this->algos.size(); ++i) 
        {  
          std::unique_ptr<boost::thread> thp;
          thp.reset(new boost::thread(
            &solver_t::solve, boost::ref(this->algos[i]), nt
          ));
          detail::boost_thread_id[thp->get_id()] = i;
          threads.add_thread(thp.release());
        }
        threads.join_all();
      }

// could it be just one ctor with int[solver_t::n_dims]?
      // 1D ctor
      boost_thread(
	const int s0,
	typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
        parent_t(s0, params, new mem_t(s0), mem_t::size())
      {}

      // 2D ctor
      boost_thread(
	const int s0,
	const int s1,
	typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
	parent_t(s0, s1, params, new mem_t(s0, s1), mem_t::size(), 1)
      {}

      // 3D ctor
      boost_thread(
	const int s0,
	const int s1,
	const int s2,
	typename solver_t::params_t &params = typename solver_t::params_t()
      ) :
	parent_t(s0, s1, s2, params, new mem_t(s0, s1, s2), mem_t::size(), 1, 1)
      {}
    };
  }; // namespace concurr
}; // namespace libmpdataxx
