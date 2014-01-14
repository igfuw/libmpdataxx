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
  void update_forcings(arrvec_t<typename parent_t::arr_t> &rhs)  
  {
    parent_t::update_forcings(rhs); // incl. zeroing the rhs arrays

    auto dW   = rhs.at(w); // TODO: relies on global w !!!
    auto Tht = this->state(tht); // TODO: relies on global tht!!!!!
    auto &ijk = this->ijk;

    dW(ijk) += (formulae::g<real_t>() * si::seconds * si::seconds / si::metres) /*/ 300*/ * (Tht(ijk) - Tht_amb) / Tht_amb; // TODO: get rid of units here?
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t Tht_amb; 
  };

  // ctor
  bombel( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) :
    parent_t(args, p),
    Tht_amb(p.Tht_amb)
  {}
};
