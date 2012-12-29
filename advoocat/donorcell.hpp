/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

namespace donorcell
{
//listing11
  template<class T1, class T2, class T3> 
  inline auto F(
    const T1 &psi_l, const T2 &psi_r, const T3 &C
  ) return_macro(
    (
      (C + abs(C)) * psi_l + 
      (C - abs(C)) * psi_r
    ) / 2
  )

  inline auto donorcell( 
    const arr_1d_t &psi, 
    const arr_1d_t &C, 
    const rng_t &i
  ) return_macro(
    F(
      psi(i  ), 
      psi(i+1), 
        C(i+h)
    ) -
    F(
      psi(i-1), 
      psi(i  ), 
        C(i-h)
    )
  )

  template<int d>  
  inline auto donorcell( 
    const arr_2d_t &psi, 
    const arr_2d_t &C, 
    const rng_t &i, 
    const rng_t &j
  ) return_macro(
    F(
      psi(pi<d>(i,   j)), 
      psi(pi<d>(i+1, j)), 
        C(pi<d>(i+h, j))
    ) -
    F(
      psi(pi<d>(i-1, j)), 
      psi(pi<d>(i,   j)), 
        C(pi<d>(i-h, j))
    )
  )

  void op_1d(
    const arrvec_t<arr_1d_t> &psi, 
    const int n,
    const arr_1d_t &C, 
    const rng_t &i
  ) { 
    psi[n+1](i) = psi[n](i)
      - donorcell(psi[n], C, i);
  }

  void op_2d(
    const arrvec_t<arr_2d_t> &psi, const int n,
    const arrvec_t<arr_2d_t> &C, 
    const rng_t &i, const rng_t &j
  ) { 
    psi[n+1](i,j) = psi[n](i,j)
      - donorcell<0>(psi[n], C[0], i, j)
      - donorcell<1>(psi[n], C[1], j, i); 
  }
}; // namespace donorcell 

namespace solvers
{

template<class bcx_t, int n_eqs = 1>
struct donorcell_1d : solver_1d<bcx_t, n_eqs> 
{
  donorcell_1d(int nx) :
    solver_1d<bcx_t, n_eqs>(nx, /* halo = */ 1)
  {}  

  void advop(int e)
  {
    donorcell::op_1d(
      this->psi[e], this->n, this->C[0], this->i
    );
  }
};

template<class bcx_t, class bcy_t, int n_eqs = 1>
struct donorcell_2d : solver_2d<bcx_t, bcy_t, n_eqs> 
{
  donorcell_2d(int nx, int ny) :
    solver_2d<bcx_t, bcy_t, n_eqs>(nx, ny, /* halo = */ 1)
  {}  

  void advop(int e)
  {
    donorcell::op_2d(
      this->psi[e], this->n, this->C, this->i, this->j
    );
  }
};
}; // namespace solvers

