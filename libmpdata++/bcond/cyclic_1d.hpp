/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/bcond/bcond.hpp>

namespace libmpdataxx
{
  namespace bcond
  {
    template <typename real_t>
    class cyclic_left_1d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // methods invoked by the solver
      void bcinit(const arr_t &a) {}

      void fill_halos_sclr(const arr_t &a)
      {
	a(this->left_halo_sclr) = a(this->rght_intr_sclr);
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av)
      {
	av[0](this->left_halo_vctr) = av[0](this->rght_intr_vctr);
      }
    };

    template <typename real_t>
    class cyclic_rght_1d : public bcond_t<real_t>
    {
      using parent_t = bcond_t<real_t>;
      using arr_t = blitz::Array<real_t, 1>;
      using parent_t::parent_t; // inheriting ctor

      public:

      // methods invoked by the solver
      void bcinit(const arr_t &a) {}

      void fill_halos_sclr(const arr_t &a)
      {
	a(this->rght_halo_sclr) = a(this->left_intr_sclr);
      }

      void fill_halos_vctr_alng(const arrvec_t<arr_t> &av)
      {
	av[0](this->rght_halo_vctr) = av[0](this->left_intr_vctr);
      }
      
    };
  }; // namespace bcond
}; // namespace libmpdataxx
