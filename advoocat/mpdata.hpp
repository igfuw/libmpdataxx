/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

namespace mpdata 
{

// TODO
/*
  struct varsgn
  {
    template<class T> auto aon(const T &x) -> decltype(abs(x))
    {
      return abs(x);
    }
  };

  struct posdef
  {
    template<class T> T aon(const T &x)
    {
      return x;
    }
  };
*/

    template<class nom_t, class den_t>
    inline auto frac(
      const nom_t &nom, const den_t &den
    ) return_macro(
      where(den > 0, nom / den, 0)
    ) 

    template<int d>
    inline auto A(const arr_t &psi, 
      const rng_t &i, const rng_t &j
    ) return_macro(
      frac(
	  abs(psi(pi<d>(i+1, j))) 
	- abs(psi(pi<d>(i,   j))),
	// ----------------------
	  abs(psi(pi<d>(i+1, j))) 
	+ abs(psi(pi<d>(i,   j)))
      ) 
    ) 

    template<int d>
    inline auto B(const arr_t &psi, 
      const rng_t &i, const rng_t &j
    ) return_macro(
     frac(//aon=abs or not -> depending on the template paramater
	  abs(psi(pi<d>(i+1, j+1))) 
	+ abs(psi(pi<d>(i,   j+1))) 
	- abs(psi(pi<d>(i+1, j-1))) 
	- abs(psi(pi<d>(i,   j-1))),
	// ------------------------
	  abs(psi(pi<d>(i+1, j+1))) 
	+ abs(psi(pi<d>(i,   j+1))) 
	+ abs(psi(pi<d>(i+1, j-1))) 
	+ abs(psi(pi<d>(i,   j-1)))
      ) / 2
    )

    template<int d>
    inline auto C_bar(
      const arr_t &C, 
      const rng_t &i, 
      const rng_t &j
    ) return_macro(
      (
	C(pi<d>(i+1, j+h)) + 
	C(pi<d>(i,   j+h)) +
	C(pi<d>(i+1, j-h)) + 
	C(pi<d>(i,   j-h)) 
      ) / 4
    )

    template <int dim>
    inline auto antidiff_2D(
      const arr_t &psi, 
      const rng_t &i, const rng_t &j,
      const arrvec_t &C
    ) return_macro(
      abs(C[dim](pi<dim>(i+h, j))) 
      * (1 - abs(C[dim](pi<dim>(i+h, j)))) 
      * A<dim>(psi, i, j) 
      - C[dim](pi<dim>(i+h, j)) 
      * C_bar<dim>(C[dim-1], i, j)
      * B<dim>(psi, i, j)
    ) 

}; // namespace mpdata

template<int n_iters, class bcx_t, class bcy_t>
struct mpdata_2D : solver_2D<bcx_t, bcy_t>
{
  // member fields
  arrvec_t tmp[2];
  rng_t im, jm;

  // ctor
  mpdata_2D(int nx, int ny) : 
    solver_2D<bcx_t, bcy_t>(nx, ny, 1), 
    im(this->i.first() - 1, this->i.last()),
    jm(this->j.first() - 1, this->j.last())
  {
    int n_tmp = n_iters > 2 ? 2 : 1;
    for (int n = 0; n < n_tmp; ++n)
    {
      tmp[n].push_back(new arr_t(
        this->i^h, this->j^this->hlo));
      tmp[n].push_back(new arr_t(
        this->i^this->hlo, this->j^h));
    }
  }

  // method invoked by the solver
  void advop()
  {
    for (int step = 0; step < n_iters; ++step) 
    {
      if (step == 0) 
        donorcell::op_2D(this->psi, 
          this->n, this->C, this->i, this->j);
      else
      {
        this->cycle();
        this->bcx.fill_halos(this->psi[this->n]);
        this->bcy.fill_halos(this->psi[this->n]);

        // choosing input/output for antidiff C
        const arrvec_t 
          &C_unco = (step == 1) 
            ? this->C 
            : (step % 2) 
              ? tmp[1]  // odd steps
              : tmp[0], // even steps
          &C_corr = (step  % 2) 
            ? tmp[0]    // odd steps
            : tmp[1];   // even steps

        // calculating the antidiffusive C 
        C_corr[0](this->im+h, this->j) = 
          mpdata::antidiff_2D<0>(
            this->psi[this->n], 
            this->im, this->j, C_unco
          );
        this->bcy.fill_halos(C_corr[0]);

        C_corr[1](this->i, this->jm+h) = 
          mpdata::antidiff_2D<1>(
            this->psi[this->n], 
            this->jm, this->i, C_unco
        );
        this->bcx.fill_halos(C_corr[1]);

        // donor-cell step 
        donorcell::op_2D(this->psi, 
          this->n, C_corr, this->i, this->j);
      }
    }
  }
};
