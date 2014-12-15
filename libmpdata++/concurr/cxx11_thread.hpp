// C++11 Thread-based shared-memory concurrency for libmpdata++
//
// author[s]: Sylwester Arabas
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include <libmpdata++/concurr/detail/concurr_common.hpp>

#include <thread>
#include <mutex>
#include <condition_variable>

#include <boost/ptr_container/ptr_vector.hpp>

#include <cstdlib> // std::getenv()

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
      std::map<std::thread::id, int> cxx11_thread_id;

      // based on boost barrier's code
      class barrier
      {
	std::mutex m_mutex;
	std::condition_variable m_cond;
	std::size_t m_generation, m_count;
        const std::size_t m_threshold;

	public:

	explicit barrier(const std::size_t count) : 
          m_count(count), 
          m_threshold(count),
          m_generation(0) 
        { }

	bool wait()
	{
          std::unique_lock<std::mutex> lock(m_mutex);
          unsigned int gen = m_generation;

          if (--m_count == 0)
          {
            m_generation++;
            m_count = m_threshold;
            m_cond.notify_all();
            return true;
          }

          while (gen == m_generation)
            m_cond.wait(lock);
          return false;
	}
      };
    };

    template <
      class solver_t,
      bcond::bcond_e bcxl,
      bcond::bcond_e bcxr,
      bcond::bcond_e bcyl = bcond::null,
      bcond::bcond_e bcyr = bcond::null,
      bcond::bcond_e bczl = bcond::null,
      bcond::bcond_e bczr = bcond::null
    >
    class cxx11_thread : public detail::concurr_common<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
    {
      using parent_t = detail::concurr_common<solver_t, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;
 
      class mem_t : public parent_t::mem_t 
      {
        detail::barrier b;

        public:

        int rank()
        {
          return detail::cxx11_thread_id[std::this_thread::get_id()];
        }

	static int size() 
	{
	  const char *env_var("OMP_NUM_THREADS");
	  return (std::getenv(env_var) != NULL)
	    ? std::atoi(std::getenv(env_var)) // TODO: check if conversion OK?
	    : std::thread::hardware_concurrency();
	}

        // ctor
        mem_t(const std::array<int, solver_t::n_dims> &grid_size) :
          b(size()),
          parent_t::mem_t(grid_size, size()) 
        {}; 

	void barrier()
	{
	  b.wait();
	}
      };

      public:

      void solve(int nt)
      {
        boost::ptr_vector<std::thread> threads(mem_t::size());
        for (int i = 0; i < this->algos.size(); ++i) 
        {  
          threads.push_back(new std::thread(
            &solver_t::solve, std::ref(this->algos[i]), nt
          ));
          detail::cxx11_thread_id[threads.back().get_id()] = i;
        }
        for (auto &th : threads) th.join();
      }

      // ctor
      cxx11_thread(const typename solver_t::rt_params_t &p) : 
        parent_t(p, new mem_t(p.grid_size), mem_t::size())
      {}

    };
  }; // namespace concurr
}; // namespace libmpdataxx
