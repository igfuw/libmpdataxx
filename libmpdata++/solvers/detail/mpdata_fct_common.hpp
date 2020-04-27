/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      const int fct_min_halo = 2; // TODO move to fct::formulae? & document why 2

      template <typename ct_params_t, int minhalo>
      class mpdata_fct_common : public mpdata_osc<
        ct_params_t, detail::max(minhalo, fct_min_halo)
      >
      {
        using parent_t = mpdata_osc<
          ct_params_t, detail::max(minhalo, fct_min_halo)
        >;

        protected:

        // member fields
        typename parent_t::arr_t psi_min, psi_max, beta_up, beta_dn; 
        arrvec_t<typename parent_t::arr_t> GC_mono; 

        arrvec_t<typename parent_t::arr_t> &GC(int iter) 
        {
          if (iter > 0) return GC_mono;
          return parent_t::GC(iter);
        }

        void beta_barrier(const int &iter)
        {
          if (!opts::isset(ct_params_t::opts, opts::iga)) // this->flux would be overwritten by donor-cell
            this->mem->barrier();
        }

        public:

        // ctor
        mpdata_fct_common(
          typename parent_t::ctor_args_t args,
          const typename parent_t::rt_params_t &p
        ) : 
          parent_t(args, p),
          psi_min(args.mem->tmp[__FILE__][0][0]),
          psi_max(args.mem->tmp[__FILE__][0][1]),
          GC_mono(args.mem->tmp[__FILE__][1]),
          beta_up(args.mem->tmp[__FILE__][2][0]),
          beta_dn(args.mem->tmp[__FILE__][2][1])
        {}

        static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
          parent_t::alloc(mem, n_iters);
          parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // psi_min and psi_max
          parent_t::alloc_tmp_vctr(mem, __FILE__);    // GC_mono
          parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // beta_up, beta_dn
        }
      };

      // partial specialisations
      template<typename ct_params_t, int minhalo, class enableif = void> 
      class mpdata_fct
      {}; 
    } // namespace detail
  } // namespace solvers
} // namescpae libmpdataxx
