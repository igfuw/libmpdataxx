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

        protected:

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
//          using namespace arakawa_c;
          using real_t = typename ct_params_t::real_t;

//          if(d==0)
//            std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
//          if(d==1)
//            std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
//          if(d==2)
//            std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;

//          std::cerr << "range - dist: " << i - dist << " " << j << " " << k << std::endl;
//          std::cerr << "range + dist: " << i + dist << " " << j << " " << k << std::endl;
//
//          std::cerr << "arr range: " << arr(i, j, k) << std::endl;
//          std::cerr << "arr range - dist: " << arr(i - dist, j, k) << std::endl;
//          std::cerr << "arr range + dist: " << arr(i + dist, j, k) << std::endl;
//
//          std::cerr << "arr range pi: " << arr(pis<d>(i, j, k)) << std::endl;
//          std::cerr << "arr range - dist pi: " << arr(pis<d>(i - dist, j, k)) << std::endl;
//          std::cerr << "arr range + dist pi: " << arr(pis<d>(i + dist, j, k)) << std::endl;

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
//          using namespace arakawa_c;
          using real_t = typename ct_params_t::real_t;

          // stretching parameter, TODO: draw it from the distribution
          const real_t d_1 = -pow(2, -1./3.),
                       d_2 =  pow(2, -1./3.);

          const int x[3] = {0, 1, 2};

          // second interpolated position (j=2, between i=1 and i=2)
          const rng_t i_next = i + 2*dist;

          // helper references, follows Akinlabi et al.
          const arr_t u_0 = arr(pis<d>(i      - dist, j, k)),
                      u_1 = arr(pis<d>(i      + dist, j, k)),
                      u_2 = arr(pis<d>(i_next + dist, j, k));

          arr_t c_1 = this->c_j(pis<d>(i,      j, k)),
                c_2 = this->c_j(pis<d>(i_next, j, k)),
                f_1 = this->f_j(pis<d>(i,      j, k)),
                f_2 = this->f_j(pis<d>(i_next, j, k));

          // Eqs. (5) and (6) Akinlabi et al.
          c_1 = (u_1 - u_0) / (4 * dist) - d_1 * (u_2 - u_0) / (4 * dist);
          c_2 = (u_2 - u_1) / (4 * dist) - d_2 * (u_2 - u_0) / (4 * dist);
          f_1 = (x[2] * u_0 - x[0] * u_1) / (4 * dist) - d_1 * (x[2] * u_0 - x[0] * u_2) / (4 * dist);
          f_2 = (x[2] * u_1 - x[0] * u_2) / (4 * dist) - d_2 * (x[2] * u_0 - x[0] * u_2) / (4 * dist);

          // new u, at j=1, between i=0 and i=1, result of w_j=1(u_i=1), Eq. (1) in Akinlabi et al.
          arr(pis<d>(i    , j, k)) = c_1 * 2*dist + d_1 * u_1 + f_1; 
          // new u, at j=2, between i=1 and i=2, result of w_j=2(u_i=1), Eq. (1) in Akinlabi et al.
          arr(pis<d>(i_next, j, k)) = c_2 * 2*dist + d_2 * u_1 + f_2;
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

        void interpolate_refinee(const int e = 0)
        {
          using namespace arakawa_c; // for rng_t operator^

          rng_t mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2; // positions between already known values (to be filled during given iteration)
          rng_t ijk_r2r_0_h, ijk_r2r_1_h, ijk_r2r_2_h;       // all positions at resolution of given iteration
          int stride, hstride;

          // TEMPORARY
          this->mem->psi_ref[e] = -1000;
            this->mem->barrier();

          // fill distmem halos of refinee
          // TODO: move to bcond or sth? would be filled only by remote bcond
          this->mem->psi_ref[e](
            this->mem->grid_size_ref[0].first() - this->mem->n_ref/2,
            this->ijk_r2r[1],
            this->ijk_r2r[2]
          ) = 
          this->mem->psi[e][0](
            this->ijk[0].first()-1,
            this->ijk[1],
            this->ijk[2]
          );
          this->mem->psi_ref[e](
            this->mem->grid_size_ref[0].last() + this->mem->n_ref/2,
            this->ijk_r2r[1],
            this->ijk_r2r[2]
          ) = 
          this->mem->psi[e][0](
            this->ijk[0].last()+1,
            this->ijk[1],
            this->ijk[2]
          );

//          std::cerr << "post left halo fill psi_ref: " << this->mem->psi_ref[e];
//          std::cerr << "post left halo fill psi: " << this->mem->psi[e][0];
//          std::cerr << this->mem->psi_ref[e](
//            this->ijk_r2r[0].first()-1,
//            this->ijk_r2r[1],
//            this->ijk_r2r[2]
//          );
//          std::cerr <<           this->mem->psi[e][0](
//            this->ijk[0].first()-1,
//            this->ijk[1],
//            this->ijk[2]
//          );


          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            // messy, because in domain decomposition (sharedmem and distmem) some refined scalars are on the edge of the subdomain...
            if(i==0)
            {
              mid_ijk_r2r_0 = this->rng_midpoints(this->ijk_r2r[0], this->mem->distmem.rank(), this->mem->distmem.size());
              mid_ijk_r2r_1 = this->rng_midpoints(this->ijk_r2r[1], this->rank, this->mem->size, false); 
              mid_ijk_r2r_2 = this->rng_midpoints(this->ijk_r2r[2]);

              ijk_r2r_0_h = this->ijk_r2r[0];
              ijk_r2r_1_h = this->ijk_r2r[1];
              ijk_r2r_2_h = this->ijk_r2r[2];
            }
            else
            {
              if(i==1)
              {
                mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0, this->mem->distmem.rank(), this->mem->distmem.size());
                if(this->rank > 0)
                  mid_ijk_r2r_1 = rng_t(mid_ijk_r2r_1.first() - mid_ijk_r2r_1.stride(), mid_ijk_r2r_1.last(), mid_ijk_r2r_1.stride()); // shift back to an overlapping range along y
                mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1, this->rank, this->mem->size);
                mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);
              }
              else
              {
                mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0);
                mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1);
                mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);
              }

              ijk_r2r_0_h = this->rng_half_stride(ijk_r2r_0_h, this->mem->distmem.rank(), this->mem->distmem.size());
              ijk_r2r_1_h = this->rng_half_stride(ijk_r2r_1_h, this->rank, this->mem->size, false);
              ijk_r2r_2_h = this->rng_half_stride(ijk_r2r_2_h);
            }

            stride = ijk_r2r_0_h.stride();
            assert(ijk_r2r_0_h.stride() == ijk_r2r_1_h.stride() && ijk_r2r_1_h.stride() == ijk_r2r_2_h.stride());
            assert(stride % 2 == 0);
            hstride = stride / 2;


            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            this->mem->barrier();
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, this->rng_merge(ijk_r2r_0_h, mid_ijk_r2r_0), hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, this->rng_merge(ijk_r2r_0_h, mid_ijk_r2r_0), this->rng_merge(ijk_r2r_1_h, mid_ijk_r2r_1), hstride);


/*

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, ijk_r2r_1_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, mid_ijk_r2r_1, hstride);
            this->mem->barrier(); // necessary before interpolation along sharedmem y direction
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, mid_ijk_r2r_0, hstride);

            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);
            */

/*
            if(this->rank == 0)
            {

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, ijk_r2r_1_h, hstride);
            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, mid_ijk_r2r_1, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, mid_ijk_r2r_0, hstride);

            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, mid_ijk_r2r_0, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);

              std::cerr << std::endl;
            }
            this->mem->barrier();
            if(this->rank == 1)
            {

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, ijk_r2r_1_h, hstride);
            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, mid_ijk_r2r_1, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, mid_ijk_r2r_0, hstride);

            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, mid_ijk_r2r_0, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);

              std::cerr << std::endl;
            }
            */

            this->mem->barrier();
          }
        }

        // TODO: very similar to interpolate refinee
        void reconstruct_refinee(const int e = 0)
        {
          using namespace arakawa_c; // for rng_t operator^

          rng_t mid_ijk_r2r_0, mid_ijk_r2r_1, mid_ijk_r2r_2; // positions between already known values (to be filled during given iteration)
          rng_t ijk_r2r_0_h, ijk_r2r_1_h, ijk_r2r_2_h;       // all positions at resolution of given iteration
          int stride, hstride;

          // TEMPORARY
          this->mem->psi_ref[e] = -1000;
            this->mem->barrier();

          // fill distmem halos of refinee
          // TODO: move to bcond or sth? would be filled only by remote bcond
          this->mem->psi_ref[e](
            this->mem->grid_size_ref[0].first() - this->mem->n_ref/2,
            this->ijk_r2r[1],
            this->ijk_r2r[2]
          ) = 
          this->mem->psi[e][0](
            this->ijk[0].first()-1,
            this->ijk[1],
            this->ijk[2]
          );
          this->mem->psi_ref[e](
            this->mem->grid_size_ref[0].last() + this->mem->n_ref/2,
            this->ijk_r2r[1],
            this->ijk_r2r[2]
          ) = 
          this->mem->psi[e][0](
            this->ijk[0].last()+1,
            this->ijk[1],
            this->ijk[2]
          );

//          std::cerr << "post left halo fill psi_ref: " << this->mem->psi_ref[e];
//          std::cerr << "post left halo fill psi: " << this->mem->psi[e][0];
//          std::cerr << this->mem->psi_ref[e](
//            this->ijk_r2r[0].first()-1,
//            this->ijk_r2r[1],
//            this->ijk_r2r[2]
//          );
//          std::cerr <<           this->mem->psi[e][0](
//            this->ijk[0].first()-1,
//            this->ijk[1],
//            this->ijk[2]
//          );


          // fill refined array at position where it overlaps with the resolved array
          this->mem->refinee(e)(this->ijk_r2r) = this->mem->advectee(e)(this->ijk);

          for(int i=0; i<this->n_fra_iter; ++i)
          {
            // messy, because in domain decomposition (sharedmem and distmem) some refined scalars are on the edge of the subdomain...
            if(i==0)
            {
              mid_ijk_r2r_0 = this->rng_midpoints(this->ijk_r2r[0], this->mem->distmem.rank(), this->mem->distmem.size());
              mid_ijk_r2r_1 = this->rng_midpoints(this->ijk_r2r[1], this->rank, this->mem->size, false); 
              mid_ijk_r2r_2 = this->rng_midpoints(this->ijk_r2r[2]);

              ijk_r2r_0_h = this->ijk_r2r[0];
              ijk_r2r_1_h = this->ijk_r2r[1];
              ijk_r2r_2_h = this->ijk_r2r[2];
            }
            else
            {
              if(i==1)
              {
                mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0, this->mem->distmem.rank(), this->mem->distmem.size());
                if(this->rank > 0)
                  mid_ijk_r2r_1 = rng_t(mid_ijk_r2r_1.first() - mid_ijk_r2r_1.stride(), mid_ijk_r2r_1.last(), mid_ijk_r2r_1.stride()); // shift back to an overlapping range along y
                mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1, this->rank, this->mem->size);
                mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);
              }
              else
              {
                mid_ijk_r2r_0 = this->rng_midpoints_out(mid_ijk_r2r_0);
                mid_ijk_r2r_1 = this->rng_midpoints_out(mid_ijk_r2r_1);
                mid_ijk_r2r_2 = this->rng_midpoints_out(mid_ijk_r2r_2);
              }

              ijk_r2r_0_h = this->rng_half_stride(ijk_r2r_0_h, this->mem->distmem.rank(), this->mem->distmem.size());
              ijk_r2r_1_h = this->rng_half_stride(ijk_r2r_1_h, this->rank, this->mem->size, false);
              ijk_r2r_2_h = this->rng_half_stride(ijk_r2r_2_h);
            }

            stride = ijk_r2r_0_h.stride();
            assert(ijk_r2r_0_h.stride() == ijk_r2r_1_h.stride() && ijk_r2r_1_h.stride() == ijk_r2r_2_h.stride());
            assert(stride % 2 == 0);
            hstride = stride / 2;


            rcnstrct<0>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_0), ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            this->mem->barrier();
            rcnstrct<1>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_1), ijk_r2r_2_h, this->rng_merge(ijk_r2r_0_h, mid_ijk_r2r_0), hstride);
            rcnstrct<2>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_2), this->rng_merge(ijk_r2r_0_h, mid_ijk_r2r_0), this->rng_merge(ijk_r2r_1_h, mid_ijk_r2r_1), hstride);


//            rcnstrct<0>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_0), ijk_r2r_1_h, ijk_r2r_2_h, hstride);
//            rcnstrct<1>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_1), ijk_r2r_2_h, ijk_r2r_0_h, hstride);
//            rcnstrct<2>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_2), ijk_r2r_0_h, ijk_r2r_1_h, hstride);
//
//            rcnstrct<2>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_2), mid_ijk_r2r_0, ijk_r2r_1_h, hstride);
//            rcnstrct<2>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_2), ijk_r2r_0_h, mid_ijk_r2r_1, hstride);
//            this->mem->barrier(); // necessary before interpolation along sharedmem y direction
//            rcnstrct<1>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_1), ijk_r2r_2_h, mid_ijk_r2r_0, hstride);
//
//            rcnstrct<2>(this->mem->refinee(e), this->rng_dbl_stride(mid_ijk_r2r_2), mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);

/*
            if(this->rank == 0)
            {

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, ijk_r2r_1_h, hstride);
            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, mid_ijk_r2r_1, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, mid_ijk_r2r_0, hstride);

            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, mid_ijk_r2r_0, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);

              std::cerr << std::endl;
            }
            this->mem->barrier();
            if(this->rank == 1)
            {

            intrp<0>(this->mem->refinee(e), mid_ijk_r2r_0, ijk_r2r_1_h, ijk_r2r_2_h, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, ijk_r2r_1_h, hstride);

            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, ijk_r2r_1_h, hstride);
            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, ijk_r2r_0_h, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, ijk_r2r_0_h, mid_ijk_r2r_1, hstride);
            intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, ijk_r2r_2_h, mid_ijk_r2r_0, hstride);

            //intrp<1>(this->mem->refinee(e), mid_ijk_r2r_1, mid_ijk_r2r_2, mid_ijk_r2r_0, hstride);
            intrp<2>(this->mem->refinee(e), mid_ijk_r2r_2, mid_ijk_r2r_0, mid_ijk_r2r_1, hstride);

              std::cerr << std::endl;
            }
            */

            this->mem->barrier();
          }
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
              mem->grid_size_ref[0]^(mem->n_ref/2), // NOTE: halo size 1 along x, because in distmem reconstructed point is at the edge of a domain
              mem->grid_size_ref[1],
              srfc ? rng_t(0, 0) : mem->grid_size_ref[2],
              arr3D_storage
            )));
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
