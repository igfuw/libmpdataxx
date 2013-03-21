/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
// TODO: 2D assumed - shouldn't the name be suffixed with _2d?
// TODO: document divideby 
      template <class inhomo_solver_t, int u, int w, int divideby = -1> // TODO: -1 => some enum or smthng
      class solver_velocity_common : public inhomo_solver_t
      {
	protected:
	using parent_t = inhomo_solver_t;
        public:
	typedef typename parent_t::real_t real_t;
        protected:

	// member fields
        const rng_t im, jm;
	const real_t dx, dz;

	void hook_ante_loop(const int nt)
	{
          rng_t &i = this->i, &j = this->j;
          // allow extrapolation at the first time-step
          this->state(u, -1)(i, j) = this->state(u)(i, j); // TODO: loop over {u,w,divideby}
          this->state(w, -1)(i, j) = this->state(w)(i, j);
          if (divideby != -1) this->state(divideby, -1)(i, j) = this->state(divideby)(i, j);
          parent_t::hook_ante_loop(nt);
	}

// TODO: tmp_idx = -1 do jakiegos enuma czy consta 

	void extrp(int e) // extrapolate in time to t+1/2
	{            // psi[n-1] will not be used anymore, and it will be intentionally overwritten!
          rng_t &i = this->i, &j = this->j;
	  auto tmp = this->state(e, -1);

	  tmp(i,j) /= -2;
	  tmp(i,j) += 3./2 * this->state(e)(i,j);

	  this->xchng(e, -1);      // filling halos 
	}

	template<int d, class arr_t> 
	void intrp(
	  const arr_t psi,
	  const rng_t &i, 
	  const rng_t &j, 
	  const real_t &dx 
	)   
	{   
	  using idxperm::pi;
	  using namespace arakawa_c;
   
	  this->mem->C[d](pi<d>(i+h,j)) = this->dt / dx * .5 * (
            psi(pi<d>(i,    j)) + 
            psi(pi<d>(i + 1,j))
          );
	}  

	void hook_ante_step()
	{
	  extrp(u);      // extrapolate velocity field in time (t+1/2)
	  extrp(w);

          if (divideby != -1) // known at compile-time
          {
	    extrp(divideby);           
            // TODO: what if divideby == 0?
	    intrp<0>(this->state(u, -1) / this->state(divideby, -1), im, this->j^this->halo, dx);
	    intrp<1>(this->state(w, -1) / this->state(divideby, -1), jm, this->i^this->halo, dz);
          }
          else
          {
	    intrp<0>(this->state(u, -1), im, this->j^this->halo, dx);
	    intrp<1>(this->state(w, -1), jm, this->i^this->halo, dz);
          }

          // TODO: document why we compute courants before forcings
          parent_t::hook_ante_step(); // forcings
	}

	public:

	struct params_t : parent_t::params_t 
        { 
          real_t dx, dz;
        };

	// ctor
	solver_velocity_common(
	  typename parent_t::ctor_args_t args,
	  const params_t &p
	) : 
	  parent_t(args, p),
          im(args.i.first() - 1, args.i.last()),
          jm(args.j.first() - 1, args.j.last()),
          dx(p.dx), dz(p.dz)
	{
// TODO: assert dx / dz were set!
        } 
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
