/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/concurr/detail/concurr_ref_common_hlpr.hpp>

namespace libmpdataxx
{
  namespace concurr
  {
    namespace detail
    {
    /*
      template <typename T, typename = int>
      struct has_solver_family : std::false_type { };

      template <typename T>
      struct has_solver_family <T, decltype((void) T::solver_family, 0)> : std::true_type { };
      */


      // helpers to find solvers that have the grid_refinement
      /*
      template< class, class = void >
      struct has_grid_refinement : std::false_type { };
      template< class T >
      struct has_grid_refinement<T, std::void_t<typename T::grid_refinement>> : std::true_type { };
      */

/*
      template<class solver_t>
      constexpr bool is_derived_from_fra()
      {
        if(!has_parent_t<solver_t>::value) return false;
        return std::is_same_v<typename solver_t::solver_family, solvers::mpdata_rhs_vip_prs_sgs_fra_family_tag> ? true : is_derived_from_fra<typename solver_t::parent_t>();
      }
      */

      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr,
        class enableif = void
      >
      class concurr_common : public concurr_ref_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {
        public:
        using parent_t = concurr_ref_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;

        // ctor
        concurr_common(
          const typename solver_t_::rt_params_t &p,
          typename solver_t_::mem_t *mem_p,
          const int &size
        ) :
          parent_t(p, mem_p, size)
        {
          // allocate per-thread structures
          this->init(p, mem_p->grid_size, size);
        }
      };
/*
      // specialization for solvers that do not have grid refinement
      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
        //typename std::enable_if_t<not std::is_same_v<solver_t_::solver_family, solvers::mpdata_rhs_vip_prs_sgs_fra_family_tag{}>>
        //typename std::enable_if_t<not is_derived_from_fra<solver_t_>()>
      >
      class concurr_common <
        solver_t_,
        bcxl, bcxr,
        bcyl, bcyr,
        bczl, bczr,
        //typename std::enable_if_t<true>
        //typename std::enable_if_t<not is_derived_from_fra<solver_t_>()>
        //typename std::enable_if_t<not is_derived_from_fra<solver_t_>()>
        typename std::enable_if_t<solver_t_::grid_refinement>
      > : public concurr_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {};
      */

      // specialzie for solvers with refined grid
      template<
        class solver_t_,
        bcond::bcond_e bcxl, bcond::bcond_e bcxr,
        bcond::bcond_e bcyl, bcond::bcond_e bcyr,
        bcond::bcond_e bczl, bcond::bcond_e bczr
        //typename std::enable_if_t<not std::is_same_v<solver_t_::solver_family, solvers::mpdata_rhs_vip_prs_sgs_fra_family_tag{}>>
        //typename std::enable_if_t<not is_derived_from_fra<solver_t_>()>
      >
      class concurr_common <
        solver_t_,
        bcxl, bcxr,
        bcyl, bcyr,
        bczl, bczr,
        //typename std::enable_if_t<true>
        //typename std::enable_if_t<not is_derived_from_fra<solver_t_>()>
        //typename std::enable_if_t<not is_derived_from_fra<solver_t_>()>
        typename std::enable_if_t<solver_t_::grid_refinement>
      > : public concurr_ref_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>
      {
        public:
        using parent_t = concurr_ref_common_hlpr<solver_t_, bcxl, bcxr, bcyl, bcyr, bczl, bczr>;

        // ctor
        concurr_common(
          const typename solver_t_::rt_params_t &p,
          typename solver_t_::mem_t *mem_p,
          const int &size
        ) :
          parent_t(p, mem_p, size)
        {
          // allocate per-thread structures
          this->init(p, mem_p->grid_size, size);
        }
      };
    } // namespace detail
  } // namespace concurr
} // namespace libmpdataxx
