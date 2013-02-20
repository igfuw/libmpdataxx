/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

template <class parent_t>
class bombel : public parent_t
{
  real_t Tht_amb;

  void forcings(real_t dt)  //explicit forcings (to be applied before the eliptic solver)
  {
    auto W   = this->psi(w);
    auto Tht = this->psi(tht);
    rng_t &i = this->i, &j = this->j;

    W(i,j) += (dt * si:: seconds) * formulae::g<real_t>() * si::seconds / si::metres * (Tht(i,j) - Tht_amb) / Tht_amb;
  }

  public:

  struct params_t : parent_t::params_t { real_t Tht_amb; };

  // ctor
  bombel( // TODO: ctor_arg_t
    typename parent_t::mem_t *mem, 
    typename parent_t::bc_p &bcxl,
    typename parent_t::bc_p &bcxr,
    typename parent_t::bc_p &bcyl,
    typename parent_t::bc_p &bcyr,
    const rng_t &i,
    const rng_t &j,
    const params_t &p
  ) :
    parent_t(mem, bcxl, bcxr, bcyl, bcyr, i, j, p),
    Tht_amb(p.Tht_amb)
  {}
};
