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

    // 1D
    inline auto A(
      const arr_1d_t &psi, 
      const rng_t &i 
    ) return_macro(
      frac(
	  abs(psi(i+1)) 
	- abs(psi(i  )),
	// ----------------------
	  abs(psi(i+1)) 
	+ abs(psi(i  ))
      ) 
    ) 

    // 2D
    template<int d>
    inline auto A(
      const arr_2d_t &psi, 
      const rng_t &i, 
      const rng_t &j
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
    inline auto B(
      const arr_2d_t &psi, 
      const rng_t &i, 
      const rng_t &j
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
      const arr_2d_t &C, 
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

    inline auto antidiff(
      const arr_1d_t &psi, 
      const rng_t &i, 
      const arr_1d_t &C
    ) return_macro(
      abs(C(i+h)) 
      * (1 - abs(C(i+h))) 
      * A(psi, i) 
    ) 

    template <int dim>
    inline auto antidiff(
      const arr_2d_t &psi, 
      const rng_t &i, 
      const rng_t &j,
      const arrvec_t<arr_2d_t> &C
    ) return_macro(
      abs(C[dim](pi<dim>(i+h, j))) 
      * (1 - abs(C[dim](pi<dim>(i+h, j)))) 
      * A<dim>(psi, i, j) 
      - C[dim](pi<dim>(i+h, j)) 
      * C_bar<dim>(C[dim-1], i, j)
      * B<dim>(psi, i, j)
    ) 

}; // namespace mpdata

namespace solvers
{
  template<int n_iters, class bcx_t, int n_eqs = 1>
  class mpdata_1d : public solver_1d<bcx_t, n_eqs>
  {
    // member fields
    arrvec_t<arr_1d_t> tmp[2];
    rng_t im;

    protected:

    // method invoked by the solver
    void advop(int e)
    {
      for (int step = 0; step < n_iters; ++step) 
      {
        if (step == 0) 
          donorcell::op_1d(this->psi[e], this->n, this->C[0], this->i);
        else
        {
          this->cycle();
          this->bcx.fill_halos(this->psi[e][this->n]);

          // choosing input/output for antidiff C
          const arr_1d_t
            &C_unco = (step == 1) 
              ? this->C[0] 
              : (step % 2) 
                ? tmp[1][0]  // odd steps
                : tmp[0][0], // even steps
            &C_corr = (step  % 2) 
              ? tmp[0][0]    // odd steps
              : tmp[1][0];   // even steps

          // calculating the antidiffusive C 
          C_corr(this->im+h) = 
            mpdata::antidiff(
              this->psi[e][this->n], 
              this->im, C_unco[0]
            );

          // donor-cell step 
          donorcell::op_1d(this->psi[e], 
            this->n, C_corr[0], this->i);
        }
      }
    }

    public:

    // ctor
    mpdata_1d(int nx) : 
      solver_1d<bcx_t, n_eqs>(nx, /* halo = */1), 
      im(this->i.first() - 1, this->i.last())
    {
      int n_tmp = n_iters > 2 ? 2 : 1;
      for (int n = 0; n < n_tmp; ++n)
      {
        tmp[n].push_back(new arr_1d_t(
          this->i^h));
      }
    }

  };

  template<int n_iters, class bcx_t, class bcy_t, int n_eqs = 1>
  class mpdata_2d : public solver_2d<bcx_t, bcy_t, n_eqs>
  {
    // member fields
    arrvec_t<arr_2d_t> tmp[2];
    rng_t im, jm;

    protected:

    // method invoked by the solver
    void advop(int e)
    {
      for (int step = 0; step < n_iters; ++step) 
      {
        if (step == 0) 
          donorcell::op_2d(
            this->psi[e], 
            this->n, this->C, this->i, this->j);
        else
        {
          this->cycle();
          this->bcx.fill_halos(this->psi[e][this->n], this->j^this->halo);
          this->bcy.fill_halos(this->psi[e][this->n], this->i^this->halo);

          // choosing input/output for antidiff C
          const arrvec_t<arr_2d_t>
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
            mpdata::antidiff<0>(
              this->psi[e][this->n], 
              this->im, this->j, C_unco
            );
          this->bcy.fill_halos(C_corr[0], this->i^h);

          C_corr[1](this->i, this->jm+h) = 
            mpdata::antidiff<1>(
              this->psi[e][this->n], 
              this->jm, this->i, C_unco
          );
          this->bcx.fill_halos(C_corr[1], this->j^h);

          // donor-cell step 
          donorcell::op_2d(this->psi[e], 
            this->n, C_corr, this->i, this->j);
        }
      }
    }

    public:

    // ctor
    mpdata_2d(int nx, int ny) : 
      solver_2d<bcx_t, bcy_t, n_eqs>(nx, ny, 1), 
      im(this->i.first() - 1, this->i.last()),
      jm(this->j.first() - 1, this->j.last())
    {
      int n_tmp = n_iters > 2 ? 2 : 1;
      for (int n = 0; n < n_tmp; ++n)
      {
        tmp[n].push_back(new arr_2d_t(
          this->i^h, this->j^this->halo 
        ));
        tmp[n].push_back(new arr_2d_t(
          this->i^this->halo, this->j^h
        ));
      }
    }
  };
}; // namespace solvers
