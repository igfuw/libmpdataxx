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
    auto W   = this->state(w);
    auto Tht = this->state(tht); // TODO: relieson global tht!!!!!
    rng_t &i = this->i, &j = this->j;

    W(i,j) += (dt * si:: seconds) * formulae::g<real_t>() * si::seconds / si::metres * (Tht(i,j) - Tht_amb) / Tht_amb;
  }

  public:

  struct params_t : parent_t::params_t { real_t Tht_amb; };

  // ctor
  bombel( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) :
    parent_t(args, p),
    Tht_amb(p.Tht_amb)
  {}
};
