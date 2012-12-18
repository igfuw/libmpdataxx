/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/
/** @mainpage
 * List of examples (subset of list of tests):
 * - test_gnuplot-iostream.cpp
 */
//listing00
// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw
//listing01
typedef double real_t;
//listing02
#include <blitz/array.h>
using arr_t = blitz::Array<real_t, 2>;
using rng_t = blitz::Range;
using idx_t = blitz::RectDomain<2>;
//listing03
#define return_macro(expr) \
  -> decltype(safeToReturn(expr)) \
{ return safeToReturn(expr); } 
//listing04
#include <boost/ptr_container/ptr_vector.hpp>
struct arrvec_t : boost::ptr_vector<arr_t> {
  const arr_t &operator[](const int i) const {   
    return this->at(
      (i + this->size()) % this->size()
    ); 
  }
};
//listing05
struct hlf_t {} h;

inline rng_t operator+(
  const rng_t &i, const hlf_t &
) { 
  return i; 
} 

inline rng_t operator-(
  const rng_t &i, const hlf_t &
) { 
  return i-1; 
}
//listing06
template<class n_t>
inline rng_t operator^(
  const rng_t &r, const n_t &n
) { 
  return rng_t(
    (r - n).first(), 
    (r + n).last()
  ); 
} 
//listing07
template<int d> 
inline idx_t pi(const rng_t &i, const rng_t &j);

template<>
inline idx_t pi<0>(
  const rng_t &i, const rng_t &j
) {
  return idx_t({i,j});
};

template<>
inline idx_t pi<1>(
  const rng_t &j, const rng_t &i
) {
  return idx_t({i,j});
}; 
//listing08
template<class bcx_t, class bcy_t>
struct solver_2D
{
  // member fields
  arrvec_t psi, C;
  int n, hlo;
  rng_t i, j;
  bcx_t bcx;
  bcy_t bcy;

  // ctor
  solver_2D(int nx, int ny, int hlo) :
    hlo(hlo),
    n(0), 
    i(0, nx-1), 
    j(0, ny-1),  
    bcx(i, j, hlo), 
    bcy(j, i, hlo)
  {
    for (int l = 0; l < 2; ++l) 
      psi.push_back(new arr_t(i^hlo, j^hlo));
    C.push_back(new arr_t(i^h, j^hlo));
    C.push_back(new arr_t(i^hlo, j^h));
  }

  // accessor methods
  arr_t state() {
    return psi[n](i,j).reindex({0,0});
  }

  arr_t courant(int d) 
  { 
    return C[d]; 
  }

  // helper methods invoked by solve()
  virtual void advop() = 0;

  void cycle() 
  { 
    n = (n + 1) % 2 - 2; 
  }

  // integration logic
  void solve(const int nt) 
  {
    for (int t = 0; t < nt; ++t) 
    {
      bcx.fill_halos(psi[n]);
      bcy.fill_halos(psi[n]);
      advop();
      cycle();
    }
  }
};
//listing09
template<int d>
struct cyclic
{
  // member fields
  idx_t left_halo, rght_halo;
  idx_t left_edge, rght_edge;;

  // ctor
  cyclic(
    const rng_t &i, const rng_t &j, int hlo
  ) :
    left_halo(pi<d>(
      rng_t(i.first()-hlo, i.first()-1), j//^hlo // TODO Range::all()
    )),
    rght_edge(pi<d>(
      rng_t(i.last()-hlo+1, i.last()  ), j//^hlo // TODO Range::all()
    )),
    rght_halo(pi<d>(
      rng_t(i.last()+1, i.last()+hlo  ), j//^hlo // TODO Range::all()
    )),
    left_edge(pi<d>(
      rng_t(i.first(), i.first()+hlo-1), j//^hlo // TODO Range::all()
    ))
  {} 

  // method invoked by the solver
  void fill_halos(const arr_t &a)
  {
    a(left_halo) = a(rght_edge);     
    a(rght_halo) = a(left_edge);     
  }
};
//listing10
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
//listing12
  template<int d>  
  inline auto donorcell( 
    const arr_t &psi, const arr_t &C, 
    const rng_t &i, const rng_t &j
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
//listing13
  void op_2D(
    const arrvec_t &psi, const int n,
    const arrvec_t &C, 
    const rng_t &i, const rng_t &j
  ) { 
    psi[n+1](i,j) = psi[n](i,j)
      - donorcell<0>(psi[n], C[0], i, j)
      - donorcell<1>(psi[n], C[1], j, i); 
  }
//listing14
}; 
//listing15
template<class bcx_t, class bcy_t>
struct donorcell_2D : solver_2D<bcx_t, bcy_t> 
{
  donorcell_2D(int nx, int ny) :
    solver_2D<bcx_t, bcy_t>(nx, ny, 1)
  {}  

  void advop()
  {
    donorcell::op_2D(
      this->psi, this->n, this->C, 
      this->i, this->j
    );
  }
};
//listing16
namespace mpdata
{
//listing17
  template<class nom_t, class den_t>
  inline auto frac(
    const nom_t &nom, const den_t &den
  ) return_macro(
    where(den > 0, nom / den, 0)
  ) 
//listing18
  template<int d>
  inline auto A(const arr_t &psi, 
    const rng_t &i, const rng_t &j
  ) return_macro(
    frac(
      psi(pi<d>(i+1, j)) - psi(pi<d>(i,j)),
      psi(pi<d>(i+1, j)) + psi(pi<d>(i,j))
    ) 
  ) 
//listing19
  template<int d>
  inline auto B(const arr_t &psi, 
    const rng_t &i, const rng_t &j
  ) return_macro(
   frac(
      psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) -
      psi(pi<d>(i+1, j-1)) - psi(pi<d>(i, j-1)),
      psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) +
      psi(pi<d>(i+1, j-1)) + psi(pi<d>(i, j-1))
    ) / 2
  )
//listing20
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
//listing21
  template<int d>
  inline auto antidiff_2D(
    const arr_t &psi, 
    const rng_t &i, const rng_t &j,
    const arrvec_t &C
  ) return_macro(
    abs(C[d](pi<d>(i+h, j))) 
    * (1 - abs(C[d](pi<d>(i+h, j)))) 
    * A<d>(psi, i, j) 
    - C[d](pi<d>(i+h, j)) 
    * C_bar<d>(C[d-1], i, j)
    * B<d>(psi, i, j)
  ) 
//listing22
};
//listing23
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
//listing24
