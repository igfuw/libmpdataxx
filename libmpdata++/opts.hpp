/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <cfenv> // for fesetround() 
//#pragma STDC FENV_ACCESS ON (as of august 2014 seems not implemented in either gcc or clang - gcc bug 34678)

namespace libmpdataxx
{
  namespace opts
  {
    using opts_t = unsigned long;

    // http://stackoverflow.com/questions/523724/c-c-check-if-one-bit-is-set-in-i-e-int-variable
    constexpr opts_t bit(const opts_t &x)
    {
      return opts_t(1) << x;
    }

    constexpr bool isset(const opts_t &x, const opts_t &y)
    {
      return 0 != (x & y);
    }

    enum
    {
      fct = opts::bit(0), // flux-corrected transport
      abs = opts::bit(1), // use the abs() trick to handle variable-sign signal
      tot = opts::bit(2), // third-order accuracy terms
      eps = opts::bit(3), // use frac=nom/(den+eps) instead of frac=where(den!=0,nom/den,0) 
      npa = opts::bit(4), // use nprt=(x-abs(x))/2 instead of nprt=min(0,x), and analogous formulae for pprt
      iga = opts::bit(5), // infinite-gauge option
      nug = opts::bit(6), // non-unit G (default G = 1) - see Smolarkiewicz 2006 eq (25) and discussion below for info on G
      dfl = opts::bit(7), // devergent flows
      nkh = opts::bit(8)  // don't use Kahan summation algorithm in the donor-cell formulae
    };

  }; // namespace opts

  struct ct_params_default_t
  {
    enum { opts = opts::iga | opts::fct };
    enum { hint_norhs = 0 };
    enum { fp_round_mode = FE_TOWARDZERO };
    //enum { fp_bits = 64 };  // TODO
    struct ix {};
    static constexpr int hint_scale(const int &e) { return 0; } // base-2 logarithm
  };

}; // namespace libmpdataxx
