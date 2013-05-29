/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  */

#pragma once

namespace libmpdataxx
{
  namespace bcond
  {
    enum bcond_e { null, cyclic }; 

    template <typename real_t>
    class bcond_t
    {
      public:

      // 1D
      virtual void fill_halos_sclr(const blitz::Array<real_t, 1> &) 
      {
        assert(false && "bcond::fill_halos_sclr() called!");
      };
      virtual void fill_halos_vctr(const blitz::Array<real_t, 1> &) 
      {
        assert(false && "bcond::fill_halos_vctr() called!");
      };

      // 2D
      virtual void fill_halos_sclr(const blitz::Array<real_t, 2> &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_sclr() called!");
      };
      virtual void fill_halos_vctr(const blitz::Array<real_t, 2> &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_vctr() called!");
      };

      // 3D
      virtual void fill_halos_sclr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_sclr() called!");
      };
      virtual void fill_halos_vctr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) 
      {
        assert(false && "bcond::fill_halos_vctr() called!");
      };
    };

    template <typename real_t>
    class shared : public bcond_t<real_t> // TODO: move to a bcond_shared file and document!
    {
      public:

      virtual void fill_halos_sclr(const blitz::Array<real_t, 1> &) { };
      virtual void fill_halos_sclr(const blitz::Array<real_t, 2> &, const rng_t &) { };
      virtual void fill_halos_sclr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) { };

      virtual void fill_halos_vctr(const blitz::Array<real_t, 1> &) { };
      virtual void fill_halos_vctr(const blitz::Array<real_t, 2> &, const rng_t &) { };
      virtual void fill_halos_vctr(const blitz::Array<real_t, 3> &, const rng_t &, const rng_t &) { };
    };

  }; // namespace bcond
}; // namespace libmpdataxx
