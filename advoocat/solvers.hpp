/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

namespace solvers
{
  template <int n_eqs>
  class solver_common
  {
    protected: 

    std::vector<int> n;

    virtual void advop(int e) = 0;

    // helper methods invoked by solve()
    void adv()
    {
      for (int e = 0; e < n_eqs; ++e) advop(e);
    }

    void cycle(int e = -1) 
    { 
      if (e == -1) 
        for (e = 0; e < n_eqs; ++e) cycle(e);
      else
        n[e] = (n[e] + 1) % 2 - 2; 
    }

    public:
 
    solver_common() :
      n(n_eqs, 0) 
    {}
  };

  template<class bcx_t, int n_eqs>
  class solver_1d : public solver_common<n_eqs>
  {
    //typedef arr_1d_t arr_t;

    protected:

    bcx_t bcx;
 
    // member fields
    arrvec_t<arr_1d_t> C;

    // psi contains model state including halo region
    arrvec_t<arr_1d_t> psi[n_eqs];

    int halo;
    rng_t i;

    void xchng() 
    {
      for (int e = 0; e < n_eqs; ++e) 
        bcx.fill_halos( psi[e][ this->n[e] ] );
    }

    // ctor
    solver_1d(int nx, int halo) :
      halo(halo),
      i(0, nx-1), 
      bcx(rng_t(0, nx-1), halo)
    {
      for (int e = 0; e < n_eqs; ++e) // equations
        for (int l = 0; l < 2; ++l) // time levels
          psi[e].push_back(new arr_1d_t(i^halo));

      C.push_back(new arr_1d_t(i^h));
    }

    public:

    // integration logic
    void solve(const int nt) 
    {
      for (int t = 0; t < nt; ++t) 
      {
        xchng();
        this->adv();
        this->cycle();
      }
    }

    // accessor method for psi (hides the halo region)
    arr_1d_t state(int e = 0) 
    {
      return psi[e][ this->n[e] ](i).reindex({0});
    }

    // accessor method for the Courant number field
    arr_1d_t courant() 
    { 
      return C[0]; 
    }
  };

  template<class bcx_t, class bcy_t, int n_eqs>
  class solver_2d : public solver_common<n_eqs>
  {
    //typedef arr_2d_t arr_t;

    protected:
  
    bcx_t bcx;
    bcy_t bcy;

    arrvec_t<arr_2d_t> C, psi[n_eqs];
    int halo;
    std::vector<int> n;
    rng_t i, j;

    void xchng() 
    {
      for (int e = 0; e < n_eqs; ++e) 
      {
        bcx.fill_halos(psi[e][ this->n[e] ], j^halo);
        bcy.fill_halos(psi[e][ this->n[e] ], i^halo);
      }
    }

    // ctor
    solver_2d(int nx, int ny, int halo) :
      halo(halo),
      i(0, nx-1), 
      j(0, ny-1),  
      bcx(rng_t(0, nx-1), rng_t(0, ny-1), halo), 
      bcy(rng_t(0, ny-1), rng_t(0, nx-1), halo)
    {
      for (int e = 0; e < n_eqs; ++e) // equations
        for (int l = 0; l < 2; ++l) // time levels
          psi[e].push_back(new arr_2d_t(i^halo, j^halo));

      C.push_back(new arr_2d_t(i^h, j^halo));
      C.push_back(new arr_2d_t(i^halo, j^h));
    }

    public:

    // accessor methods
    arr_2d_t state(int e = 0) 
    {
      return psi[e][ this->n[e] ](i,j).reindex({0,0});
    }

    arr_2d_t courant(int d) 
    { 
      return C[d]; 
    }

    // integration logic
    void solve(const int nt) 
    {
      for (int t = 0; t < nt; ++t) 
      {
        xchng();
        this->adv();
        this->cycle();
      }
    }
  };
}; // namespace solvers
