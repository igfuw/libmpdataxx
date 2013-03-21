/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// TODO: includes!

template <class parent_t>
class bombel : public parent_t
{
  real_t Tht_amb;

  // explicit forcings (to be applied before the eliptic solver)
  void forcings(real_t dt)  
  {
    auto W   = this->state(w);
    auto Tht = this->state(tht); // TODO: relies on global tht!!!!!
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
