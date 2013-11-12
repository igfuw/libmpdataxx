/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
// TODO: 2D assumed - shouldn't the name be suffixed with _2d?
      template <class inhomo_solver_t, int u, int w, int density = -1> // TODO: -1 => some enum or smthng
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
          parent_t::hook_ante_loop(nt);
          // write to stash current state 
          // to artificially allow extrapolation at the first time-step
          if (density == -1) 
          {
	    stash[0](this->ijk) = this->state(u)(this->ijk);
	    stash[1](this->ijk) = this->state(w)(this->ijk);
          }
          else
          {
	    stash[0](this->ijk) = this->state(u)(this->ijk) / this->state(density)(this->ijk); // TODO: what if density == 0?
	    stash[1](this->ijk) = this->state(w)(this->ijk) / this->state(density)(this->ijk); // TODO: what if density == 0?
          } 
	}

        template <int d>
	void extrp(int e) // extrapolate velocity field in time to t+1/2
	{                 // (write the result to stash since we don't need previous state any more)
	  stash[d](this->ijk) /= -2.;

          if (density == -1) 
            stash[d](this->ijk) += 3./2 * this->state(e)(this->ijk);
          else
            stash[d](this->ijk) += 3./2 * (this->state(e)(this->ijk) / this->state(density)(this->ijk)); // TODO: what if density == 0?

	  this->xchng(stash[d], this->i^this->halo, this->j^this->halo);      // filling halos 
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
    
          if (this->mem->G.numElements() == 0)
          {
	    this->mem->GC[d](pi<d>(i+h,j)) = this->dt / dx * .5 * (
	      psi(pi<d>(i,    j)) + 
	      psi(pi<d>(i + 1,j))
	    );
          } 
          else
          { 
            // TODO
            assert(false);
          }
	}  

	void hook_ante_step()
	{ //extrapolate velocity field in time (t+1/2)
	  extrp<0>(u);     
	  extrp<1>(w);
          //interpolate from velocity field to courant field (mpdata needs courant numbers from t+1/2)
	  intrp<0>(stash[0], im, this->j^this->halo, dx);
	  intrp<1>(stash[1], jm, this->i^this->halo, dz);

          this->mem->barrier();

          // filling the stash with data from current velocity field 
          // (so that in the next time step they can be used for extrapolation in time)
          if (density == -1)
          {
	    stash[0](this->ijk) = this->state(u)(this->ijk);
	    stash[1](this->ijk) = this->state(w)(this->ijk);
          }
          else
          {
	    stash[0](this->ijk) = this->state(u)(this->ijk) / this->state(density)(this->ijk);
	    stash[1](this->ijk) = this->state(w)(this->ijk) / this->state(density)(this->ijk);
          }

          // intentionally after stash !!!
          // (we have to stash data from the current time step before applying any forcings to it)
          parent_t::hook_ante_step(); // applying forcings
	}

        // member fields
        arrvec_t<typename parent_t::arr_t> &stash;

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
          dx(p.dx), dz(p.dz),
          stash(args.mem->tmp[__FILE__][0])
	{
// TODO: assert dx / dz were set!
        } 

	// TODO: merge the two allocs into one!

	// 1D version
	static void alloc(typename parent_t::mem_t *mem, const int nx)
	{
	  parent_t::alloc(mem, nx);
	  parent_t::alloc_tmp_sclr(mem, nx, __FILE__, parent_t::n_dims); // psi[n-1] secret stash for velocity extrapolation in time
	}

	// 2D version
	static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
	{
	  parent_t::alloc(mem, nx, ny);
	  parent_t::alloc_tmp_sclr(mem, nx, ny, __FILE__, parent_t::n_dims); // psi[n-1] secret stash for velocity extrpolation in time
	}

	// TODO: 3D version
      }; 
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
