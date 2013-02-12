/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "concurr_common.hpp"

#if !defined(_REENTRANT)
#  error _REENTRANT not defined, pleas use something like -pthread flag for gcc
#endif
#include <boost/thread.hpp>

namespace advoocat
{
  namespace concurr
  {
    template <
      class solver_t,
      bcond::bcond_e bcx,
      bcond::bcond_e bcy = bcond::null,
      bcond::bcond_e bcz = bcond::null
    >
    class boost_thread : public detail::concurr_common<solver_t, bcx, bcy, bcz>
    {
      using parent_t = detail::concurr_common<solver_t, bcx, bcy, bcz>;
 
      std::unique_ptr<boost::barrier> b;

      int size() 
      {
        return boost::thread::hardware_concurrency();
      }

      void init()
      {
        b.reset(new boost::barrier(size()));
      }

      void barrier()
      {
        b->wait();
      }

      public:

      void advance(int nt)
      {
        boost::thread_group threads; // member field?
// TODO!
//        for (int i = 0; i < this->algos.size(); ++i) threads.add_thread(new boost::thread(
//          &solver_t::solve, this->algos[i], nt 
//        ));
        threads.join_all();
      }

      // 1D ctor
      boost_thread(
	const int s0,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
        parent_t(s0, params, size())
      {
        init();
      }

      // 2D ctor
      boost_thread(
	const int s0,
	const int s1,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) : 
	parent_t(s0, s1, params, size(), 1)
      {
        init();
      }

      // 3D ctor
      boost_thread(
	const int s0,
	const int s1,
	const int s2,
	const typename solver_t::params_t &params = typename solver_t::params_t()
      ) :
	parent_t(s0, s1, s2, params, size(), 1, 1)
      {
        init();
      }
    };
  }; // namespace concurr
}; // namespace advoocat
