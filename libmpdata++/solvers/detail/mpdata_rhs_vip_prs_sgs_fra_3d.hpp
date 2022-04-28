/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/solvers/detail/mpdata_rhs_vip_prs_sgs_fra_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      template<typename ct_params_t, int minhalo, class enableif = void>
      class mpdata_rhs_vip_prs_sgs_fra_dim
      {
//        static_assert(false);
      };

      template<typename ct_params_t, int minhalo>
      class mpdata_rhs_vip_prs_sgs_fra_dim<
        ct_params_t,
        minhalo,
        typename std::enable_if_t<ct_params_t::n_dims == 3>
      > : public detail::mpdata_rhs_vip_prs_sgs_fra_common<ct_params_t, minhalo>
      {
        using parent_t = detail::mpdata_rhs_vip_prs_sgs_fra_common<ct_params_t, minhalo>;
        using parent_t::parent_t; // inheriting constructors
        using real_t = typename ct_params_t::real_t;

        protected:
  
        // helper variables for grid refinement
        rng_t mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2;     // positions between already known values (to be filled during given iteration)
        rng_t ijk_r2r_0_h_with_halo, ijk_r2r_1_h, ijk_r2r_2_h; // all positions at resolution of given iteration
        int stride, hstride;                                   // stride and half stride

        // CDF of stretching parameter d; assuming CDF = A d^B + C under condition that CDF(0.5)=0 and CDF(1)=1
        static real_t CDF_of_d(const real_t &d)
        {
          const real_t B = -0.39091161; // comes from a fit to the E. Akinlabi data
          return (pow(d,B) - pow(0.5,B)) / (1. - pow(0.5,B));
        }
        // inverse function: d of a CDF
        static real_t d_of_CDF(const real_t &CDF)
        {
          const real_t B = -0.39091161; // comes from a fit to the E. Akinlabi data
          return pow((1.-pow(0.5,B))*CDF + pow(0.5,B), 1./B);
        }

        struct d_of_CDF_fctr
        {
          real_t operator()(const real_t &CDF) const
          {
            return CDF < 0 ? -d_of_CDF(-CDF) : d_of_CDF(CDF);
          }
          BZ_DECLARE_FUNCTOR(d_of_CDF_fctr);
        };


        // interpolation similar to mpdata_rhs_vip...
        template<int d, class arr_t>
        void intrp(
          arr_t arr,
          const rng_t &i,
          const rng_t &j,
          const rng_t &k,
          const int &dist
        )
        {
          using idxperm::pis;

          // debug output
          std::cerr << "ranges in interpolation" << std::endl;
          if(d==0)
            std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
          if(d==1)
            std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
          if(d==2)
            std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;

          arr(pis<d>(i, j, k)) = real_t(.5) * (
            arr(pis<d>(i - dist, j, k)) +
            arr(pis<d>(i + dist, j, k))
          );
        }

        // fractral reconstruction (calculates 2 mid-points based on 3)
        // notation following Akinlabi et al. "Fractal reconstruciton..." 2019
        template<int d, class arr_t>
        void rcnstrct(
          arr_t arr,
          const rng_t &i,
          const rng_t &j,
          const rng_t &k,
          const int &dist // number of indices between a reconstructed point and a known point
        )
        {
          // TODO: move some to formulas::
          using idxperm::pis;

          // second reconstructed position (j=2, between i=1 and i=2)
          const rng_t i_next = i + 2*dist;
          
          // debug output
          std::cerr << "ranges in rcnstrct" << std::endl;
          if(d==0)
            std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
          if(d==1)
            std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
          if(d==2)
            std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;
          // do interpolation, useful for testing if ranges are correct
          /*
          arr(pis<d>(i, j, k)) = real_t(.5) * (
            arr(pis<d>(i - dist, j, k)) +
            arr(pis<d>(i + dist, j, k))
            );
          arr(pis<d>(i_next, j, k)) = real_t(.5) * (
            arr(pis<d>(i_next - dist, j, k)) +
            arr(pis<d>(i_next + dist, j, k))
            );
          return;
          */
          // end of the interpolation test


          // stretching parameter, TODO: draw it from the distribution
//          const real_t d_1 = -pow(2, -1./3.),
//                       d_2 =  pow(2, -1./3.);


// solution following Scotti and Meneveau Physica D 1999
/*

          // helper references
          // NOTE: these are actually the resolved values only in the first iteration, in the subsequent iterations some of these are reconstructed values
          //       should we point to the resolved values in all iterations?? 
          const arr_t u_im1 = arr(pis<d>(i      - dist, j, k)), // resolved u_(i-1)
                      u_i   = arr(pis<d>(i      + dist, j, k)), // resolved u_(i)
                      u_ip1 = arr(pis<d>(i_next + dist, j, k)); // resolved u_(i+1)

          // c_j and f_j naming comes from the Akinlabi paper, TODO: rename to a_i b_i
          // NOTE: in multiple iterations of reconstruction, they only change due to random stretching parameters? 
          //       or the stretching parameter should be drawn from the distro once? then no need for c_j and f_j to have the same size as reconstructed fields...
          arr_t a_1 = this->c_j(pis<d>(i,      j, k)),
                a_2 = this->c_j(pis<d>(i_next, j, k)),
                b_1 = this->f_j(pis<d>(i,      j, k)),
                b_2 = this->f_j(pis<d>(i_next, j, k));

          // Eqs. (13) and (14) from Scotti Meneveau 1999 
          a_1 = u_i - u_im1 - d_1*(u_ip1 - u_im1);
          a_2 = u_ip1 - u_i - d_2*(u_ip1 - u_im1);
          b_1 = u_im1*(real_t(1) - d_1); 
          b_2 = u_i - d_2*u_im1; 

          // W(u_0(xi)) at xi=1/4, between u_i-1 and u_i
          // u_0(2*xi) = u_i-1 + 1/2 * (u_i+1 - u_i-1)
          // Eq. (6) and caption of fig. 1 Scotti 1999
          // NOTE: in multiple iterations, should xi be between the 3 resolved points? or should we consider the reconstructed points as new resolved points?
          arr(pis<d>(i,      j, k)) = d_1 * 
//                                      (u_im1 + 0.5 * (u_ip1 - u_im1)) + // u_0(2*xi)
                                      u_i + 
                                      0.5 * a_1 + b_1;                  // q_1(2*xi)
          // W(u_0(xi)) at xi=3/4, between u_i and u_i+1
          // u_0(2*xi-1) = u_i-1 + 1/2 * (u_i+1 - u_i-1)
          // NOTE: u_0(2*3/4-1) == u_0(2*1/4), so calculate it once in the previous step and store it...
          arr(pis<d>(i_next, j, k)) = d_2 * 
//                                      (u_im1 + 0.5 * (u_ip1 - u_im1)) + // u_0(2*xi-1)
                                      u_i + 
                                      0.5 * a_2 + b_2;                  // q_2(2*xi-1)
*/
          // Eq. (5) Scotti and Meneveau PRL 1997, a_1 and a_2 differ from Scotti Meneveau 1999 !
          // NOTE: in the 1999 paper, a_1 and a_2 are multipled by 2*xi and 2*xi-1, respectively
          //       in the 1997 paper, the are both multipled by xi
          /*
          a_1 = 2 * (u_im1 - u_i - d_1*(u_ip1 - u_im1));
          a_2 = 2 * (u_ip1 - u_i - d_2*(u_ip1 - u_im1));
          b_1 = u_im1*(real_t(1) - d_1); 
          b_2 = u_i - d_2*u_im1; 
          */

          //const real_t xi[3] = {0, 0.5, 1};


// solution following Akinlabi et al.


          const int x[3] = {0, 1, 2};

          // helper references
          const arr_t u_0 = arr(pis<d>(i      - dist, j, k)),
                      u_1 = arr(pis<d>(i      + dist, j, k)),
                      u_2 = arr(pis<d>(i_next + dist, j, k));

          arr_t c_1 = this->c_j(pis<d>(i,      j, k)),
                c_2 = this->c_j(pis<d>(i_next, j, k)),
                d_1 = this->d_j(pis<d>(i,      j, k)),
                d_2 = this->d_j(pis<d>(i_next, j, k)),
                f_1 = this->f_j(pis<d>(i,      j, k)),
                f_2 = this->f_j(pis<d>(i_next, j, k));




          // Eqs. (5) and (6) Akinlabi et al.
          c_1 = (u_1 - u_0) / (2.) - d_1 * (u_2 - u_0) / (2.);
          c_2 = (u_2 - u_1) / (2.) - d_2 * (u_2 - u_0) / (2.);
          f_1 = (x[2] * u_0 - x[0] * u_1) / (2.) - d_1 * (x[2] * u_0 - x[0] * u_2) / (2.);
          f_2 = (x[2] * u_1 - x[0] * u_2) / (2.) - d_2 * (x[2] * u_0 - x[0] * u_2) / (2.);
//          c_1 = (u_1 - u_0) / (4 * dist) - d_1 * (u_2 - u_0) / (4 * dist);
//          c_2 = (u_2 - u_1) / (4 * dist) - d_2 * (u_2 - u_0) / (4 * dist);
//          f_1 = (x[2] * u_0 - x[0] * u_1) / (4 * dist) - d_1 * (x[2] * u_0 - x[0] * u_2) / (4 * dist);
//          f_2 = (x[2] * u_1 - x[0] * u_2) / (4 * dist) - d_2 * (x[2] * u_0 - x[0] * u_2) / (4 * dist);

          // new u, at j=1, between i=0 and i=1, result of w_j=1(u_i=1), Eq. (1) in Akinlabi et al.
          arr(pis<d>(i    , j, k)) = c_1 * 1 + d_1 * u_1 + f_1; 
          // new u, at j=2, between i=1 and i=2, result of w_j=2(u_i=1), Eq. (1) in Akinlabi et al.
          arr(pis<d>(i_next, j, k)) = c_2 * 1 + d_2 * u_1 + f_2;
          // not sure about * 2*dist here (is x_1 == 2*dist?)
          

/*
            c_j(pis<d>(i, j, k)) * x_j + 
            d_j[0] * arr(pis<d>(i + dist, j, k)) + // coarse-grained u at i=1
            f_j(pis<d>(i, j, k)); 
            */


          // new u, at j=2, between i=1 and i=2, result of w2(u1)
          // TODO

/*
          arr(pis<d>(i, j, k)) = real_t(.5) * (
            arr(pis<d>(i - dist, j, k)) +
            arr(pis<d>(i + dist, j, k))
          );

          arr(pis<d>(i_next, j, k)) = real_t(.5) * (
            arr(pis<d>(i_next - dist, j, k)) +
            arr(pis<d>(i_next + dist, j, k))
          );
          */
        }

        // fill distmem halos of refinee
        // TODO: move to bcond or sth? would be filled only by remote bcond
        void fill_refinee_distmem_halos(const int e, const int halo_size)
        {
          // TODO: we only need to xchng along distmem direction (x)
          this->xchng(e);

          // TODO: assert halo_size <= solver halo size

          switch(halo_size)
          {
            case 2:
              this->mem->psi_ref[e](
                this->mem->grid_size_ref[0].last() + 1 + this->mem->n_ref ,
                this->ijk_r2r[1],
                this->ijk_r2r[2]
              ) = 
                this->mem->psi[e][0](
                  this->ijk[0].last()+2,
                  this->ijk[1],
                  this->ijk[2]
                );
              // no break intentionally
            case 1:
              this->mem->psi_ref[e](
                this->mem->grid_size_ref[0].first() - this->mem->n_ref,
                this->ijk_r2r[1],
                this->ijk_r2r[2]
              ) = 
                this->mem->psi[e][0](
                  this->ijk[0].first()-1,
                  this->ijk[1],
                  this->ijk[2]
                );

              this->mem->psi_ref[e](
                this->mem->grid_size_ref[0].last() + 1,
                this->ijk_r2r[1],
                this->ijk_r2r[2]
              ) = 
                this->mem->psi[e][0](
                  this->ijk[0].last()+1,
                  this->ijk[1],
                  this->ijk[2]
                );
              break;
            default:
              assert(false);
              break;
          }
        }

        void refinement_ranges(const int iter, const int halo_size)
        {
          // messy, because in domain decomposition (sharedmem and distmem) some refined scalars are on the edge of the subdomain...
          if(iter==0)
          {
            mid_ijk_r2r_0 = this->rng_midpoints(this->ijk_r2r[0], this->mem->distmem.rank(), this->mem->distmem.size());
            mid_ijk_r2r_1 = this->rng_midpoints(this->ijk_r2r[1], this->rank, this->mem->size); 
            mid_ijk_r2r_2 = this->rng_midpoints(this->ijk_r2r[2]);

            ijk_r2r_1_h = this->ijk_r2r[1];
            ijk_r2r_2_h = this->ijk_r2r[2];
          }
          else
          {
            mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0);
            mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1);
            mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);

            ijk_r2r_1_h = this->rng_half_stride(ijk_r2r_1_h, this->rank, this->mem->size);
            ijk_r2r_2_h = this->rng_half_stride(ijk_r2r_2_h);
          }

          stride = ijk_r2r_1_h.stride();
          assert(ijk_r2r_1_h.stride() == ijk_r2r_2_h.stride());
          assert(stride % 2 == 0);
          hstride = stride / 2;

          ijk_r2r_0_h_with_halo = this->mem->distmem.rank() < this->mem->distmem.size() - 1 ?
            rng_t(this->ijk_r2r[0].first(), this->ijk_r2r[0].last() + halo_size * this->mem->n_ref, hstride) :
            rng_t(this->ijk_r2r[0].first(), this->ijk_r2r[0].last(),                                hstride);
        }

        // TODO: move some formulas
        // TODO: move to common?
        void generate_stretching_parameters(const int rng_seed = 44)
        {
          std::mt19937 gen(rng_seed); 
          std::uniform_real_distribution<> dis(-1, 1); // [-1,1), but whatever
          auto rand = std::bind(dis, gen);
          std::generate(this->d_j(this->ijk_ref).begin(), this->d_j(this->ijk_ref).end(), rand);
          this->d_j(this->ijk_ref) = d_of_CDF_fctr{}(this->d_j(this->ijk_ref));
          // TODO: randomly convert to negative
        }

        void interpolate_refinee(const int e = 0)
        {
          const int halo_size = 1;

          // TEMPORARY
          this->mem->psi_ref[e] = -1000;
          this->mem->barrier();

          fill_refinee_distmem_halos(e, halo_size);

          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

          generate_stretching_parameters();

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            refinement_ranges(i, halo_size);

            intrp<0>(this->mem->psi_ref[e], mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            this->mem->barrier();
            intrp<1>(this->mem->psi_ref[e], mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h_with_halo, hstride);
            intrp<2>(this->mem->psi_ref[e], mid_ijk_r2r_2, ijk_r2r_0_h_with_halo, this->rng_merge(ijk_r2r_1_h, mid_ijk_r2r_1), hstride);
          }
          this->mem->barrier();
        }

        void reconstruct_refinee(const int e = 0)
        {
          const int halo_size = 2;

          // TEMPORARY
          this->mem->psi_ref[e] = -1000;
          this->mem->barrier();

          fill_refinee_distmem_halos(e, halo_size);

          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

          generate_stretching_parameters();

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            refinement_ranges(i, halo_size);

            // if MPI domain starts with a point in the middle of a triple - we need to start calculating from a point to the left (in halo)
            const int offset = (i==0 && this->mem->grid_size_ref[0].first() > 0 && (this->mem->grid_size_ref[0].first() / stride) % 2 != 0) ? stride : 0;

            rcnstrct<0>(this->mem->psi_ref[e], this->rng_dbl_stride(mid_ijk_r2r_0, offset), ijk_r2r_1_h,           ijk_r2r_2_h,                                 hstride);
            this->mem->barrier();
            rcnstrct<1>(this->mem->psi_ref[e], this->rng_dbl_stride(mid_ijk_r2r_1)        , ijk_r2r_2_h,           ijk_r2r_0_h_with_halo,                       hstride);
            rcnstrct<2>(this->mem->psi_ref[e], this->rng_dbl_stride(mid_ijk_r2r_2)        , ijk_r2r_0_h_with_halo, this->rng_merge(ijk_r2r_1_h, mid_ijk_r2r_1), hstride);
          }
          this->mem->barrier();
        }

        public:

        // helper method to allocate n_arr refined scalar temporary arrays
        static void alloc_tmp_sclr_ref(
          typename parent_t::mem_t *mem,
          const char * __file__, const int n_arr,
          std::string name = "",
          bool srfc = false
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());

          if (!name.empty()) mem->avail_tmp[name] = std::make_pair(__file__, mem->tmp[__file__].size() - 1);

          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
              mem->grid_size_ref[0]^(mem->n_ref+1), // reconstruction based on 3 points, we need up to 2 resolved points outside of the domain
              mem->grid_size_ref[1],
              srfc ? rng_t(0, 0) : mem->grid_size_ref[2],
              arr3D_storage
            )));
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
