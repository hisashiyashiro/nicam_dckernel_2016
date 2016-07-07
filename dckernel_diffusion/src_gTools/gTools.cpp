#define PEDANTIC_DISABLED
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "tools/verifier.hpp"
#include <stencil-composition/stencil-composition.hpp>
#include "problem_size.hpp"
#ifdef __CUDACC__
#include <boundary-conditions/apply_gpu.hpp>
#else
#include <boundary-conditions/apply.hpp>
#endif

using namespace gridtools;
using namespace expressions;
using gridtools::direction;
using gridtools::sign;
using gridtools::minus_;
using gridtools::zero_;
using gridtools::plus_;
  
 
 // single interval for vertical direction
 typedef gridtools::interval< level< 0, -1 >, level< 1, -1 > > x_interval;
 
 // only one axis, but different interval for different equations   
 typedef gridtools::interval< level< 0, -2 >, level< 1, 1 > > axis;


struct laplace_function {
      
  // vxt --> [I, J, T]            
  typedef accessor<0, enumtype::inout, extent<>, 4 > vxt;
  // vyt --> [I, J, T]            
  typedef accessor<1, enumtype::inout, extent<>, 4 > vyt;
  // vzt --> [I, J, T]            
  typedef accessor<2, enumtype::inout, extent<>, 4 > vzt;
  // scl --> [I, J, K]+[L] 
  typedef accessor<3, enumtype::in, extent<0, 1, 0, 1>, 5 > scl;
  // cinterp_TN  --> [I, J, A, C]+[L] 
  typedef accessor<4, enumtype::in, extent<0, 1, 0, 1>, 7 > cinterp_TN;
  // cinterp_TRA [I, J, T]+[L] 
  typedef accessor<5, enumtype::in, extent<0, 0, 0, 0>, 6 > cinterp_TRA;
    
  typedef boost::mpl::vector<vxt, vyt, vzt, scl, cinterp_TN, cinterp_TRA> arg_list;

  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<4>::Index t;
        dimension<4>::Index a;
        dimension<5>::Index c;
        
        // TI
        auto smean = eval((scl{}        + scl{i+1} + scl{i+1,j+1} ) / 3.0);
        auto u1    = eval((scl{}        + scl{i+1}))*0.5     - smean;
        auto u2    = eval((scl{i+1}     + scl{i+1,j+1}))*0.5 - smean;
        auto u3    = eval((scl{i+1,j+1} + scl{}))*0.5        - smean;
        
        eval(vxt{}) =  (- u1 * eval(cinterp_TN{})
                        + u2 * eval(cinterp_TN{i+1,a+2}) 
                        + u3 * eval(cinterp_TN{a+1} ))
			                       * eval(cinterp_TRA{});

        eval(vyt{}) =  (- u1 * eval(cinterp_TN{c+1})
                        + u2 * eval(cinterp_TN{i+1,a+2,c+1}) 
                        + u3 * eval(cinterp_TN{a+1,c+1} )) 
                             * eval(cinterp_TRA{});
             
        eval(vzt{}) =  (- u1 * eval(cinterp_TN{c+2})
                        + u2 * eval(cinterp_TN{i+1,a+2,c+2}) 
                        + u3 * eval(cinterp_TN{a+1,c+2} )) 
                              * eval(cinterp_TRA{});

        // TJ    
        smean = eval((scl{}        + scl{i+1,j+1} + scl{j+1} ) / 3.0);
        u1    = eval((scl{}        + scl{i+1,j+1}))*0.5 - smean;
        u2    = eval((scl{i+1,j+1} + scl{j+1}))*0.5     - smean;
        u3    = eval((scl{j+1}     + scl{}))*0.5        - smean;

        eval(vxt{t+1}) =  (- u1 * eval(cinterp_TN{a+1})
                           + u2 * eval(cinterp_TN{j+1}) 
                           + u3 * eval(cinterp_TN{a+2} ))   
                                * eval(cinterp_TRA{t+1});

        eval(vyt{t+1}) =  (- u1 * eval(cinterp_TN{a+1,c+1})
                           + u2 * eval(cinterp_TN{i+1,c+1}) 
                           + u3 * eval(cinterp_TN{a+2,c+1} )) 
                                * eval(cinterp_TRA{t+1});
             
        eval(vzt{t+1}) =  (- u1 * eval(cinterp_TN{a+1,c+2})
                           + u2 * eval(cinterp_TN{i+1,c+2}) 
                           + u3 * eval(cinterp_TN{a+2,c+2} )) 
                                * eval(cinterp_TRA{t+1});
        /*
        // TI
        auto smean = eval((scl{i,j,k}     + scl{i+1,j,k} + scl{i+1,j+1,k} ) / 3.0);
	    auto u1    = eval((scl{i,j,k}     + scl{i+1,j,k}))*0.5   - smean;
	    auto u2    = eval((scl{i+1,j,k}   + scl{i+1,j+1,k}))*0.5 - smean;
	    auto u3    = eval((scl{i+1,j+1,k} + scl{i,j,k}))*0.5     - smean;
	    
        eval(vxt{i,j,t}) =  (- u1 * eval(cinterp_TN{i,j,a,c})
		                         + u2 * eval(cinterp_TN{i+1,j,a+2,c}) 
		                         + u3 * eval(cinterp_TN{i,j,a+1,c} ))   
                                  * eval(cinterp_TRA{i,j,t});

        eval(vyt{i,j,t}) =  (- u1 * eval(cinterp_TN{i,j,a,c+1})
                             + u2 * eval(cinterp_TN{i+1,j,a+2,c+1}) 
                             + u3 * eval(cinterp_TN{i,j,a+1,c+1} )) 
                                  * eval(cinterp_TRA{i,j,t});
             
        eval(vzt{i,j,t}) =  (- u1 * eval(cinterp_TN{i,j,a,c+2})
                             + u2 * eval(cinterp_TN{i+1,j,a+2,c+2}) 
                             + u3 * eval(cinterp_TN{i,j,a+1,c+2} )) 
                                  * eval(cinterp_TRA{i,j,t});

        // TJ    
        smean = eval((scl{i,j,k}     + scl{i+1,j+1,k} + scl{i,j+1,k} ) / 3.0);
        u1    = eval((scl{i,j,k}     + scl{i+1,j+1,k}))*0.5 - smean;
        u2    = eval((scl{i+1,j+1,k} + scl{i,j+1,k}))*0.5   - smean;
        u3    = eval((scl{i,j+1,k}   + scl{i,j,k}))*0.5     - smean;

        eval(vxt{i,j,t+1}) =  (- u1 * eval(cinterp_TN{i,j,a+1,c})
                               + u2 * eval(cinterp_TN{i,j+1,a,c}) 
                               + u3 * eval(cinterp_TN{i,j,a+2,c} ))  
                                    * eval(cinterp_TRA{i,j,t+1});

        eval(vyt{i,j,t+1}) =  (- u1 * eval(cinterp_TN{i,j,a+1,c+1})
                               + u2 * eval(cinterp_TN{i+1,j,a,c+1}) 
                               + u3 * eval(cinterp_TN{i,j,a+2,c+1} )) 
                                    * eval(cinterp_TRA{i,j,t+1});
             
        eval(vzt{i,j,t+1}) =  (- u1 * eval(cinterp_TN{i,j,a+1,c+2})
                               + u2 * eval(cinterp_TN{i+1,j,a,c+2}) 
                               + u3 * eval(cinterp_TN{i,j,a+2,c+2} ))
                                    * eval(cinterp_TRA{i,j,t+1});
        */
    }
}; // end laplace_function

struct flux_function {
    
  // vxt  --> [I, J, T]            
  typedef accessor<0, enumtype::in, extent<-1, 0, -1, 0>, 4 > vxt;
  // vyt  --> [I, J, T]            
  typedef accessor<1, enumtype::in, extent<-1, 0, -1, 0>, 4 > vyt;
  // vzt  --> [I, J, T]              
  typedef accessor<2, enumtype::in, extent<-1, 0, -1, 0>, 4 > vzt;
  // flux --> [I, J, A] 
  typedef accessor<3, enumtype::inout, extent<>, 4 > flux;
  // kh   --> [I, J, K]+[L]
  typedef accessor<4, enumtype::in,    extent<0, 1, 0, 1>, 5 > kh;
  // cinterp_HN  --> [I, J, A, C]+[L]
  typedef accessor<5, enumtype::in, extent<>, 7 > cinterp_HN;
    
  typedef boost::mpl::vector<vxt, vyt, vzt, flux, kh, cinterp_HN> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<4>::Index a;
        dimension<4>::Index t;
        dimension<5>::Index c;
        
        //AI
        eval(flux{}) =  ( (eval(vxt{j-1,t+1}) + eval(vxt{})) * eval(cinterp_HN{}) 
                        + (eval(vyt{j-1,t+1}) + eval(vyt{})) * eval(cinterp_HN{c+1}) 
                        + (eval(vzt{j-1,t+1}) + eval(vzt{})) * eval(cinterp_HN{c+2}) ) 
                 * 0.25 * (eval(kh{}) + eval(kh{i+1}));
        
        // AIJ    
        eval(flux{a+1}) =  ((eval(vxt{}) + eval(vxt{t+1})) * eval(cinterp_HN{a+1})
                          + (eval(vyt{}) + eval(vyt{t+1})) * eval(cinterp_HN{a+1,c+1}) 
                          + (eval(vzt{}) + eval(vzt{t+1})) * eval(cinterp_HN{a+1,c+2}) )
                   * 0.25 * (eval(kh{})  + eval(kh{i+1,j+1}));

        // AIJ    
        eval(flux{ a+2 }) =   ((eval(vxt{t+1}) + eval(vxt{i-1})) * eval(cinterp_HN{ a+2 })
 		                         + (eval(vyt{t+1}) + eval(vyt{i-1})) * eval(cinterp_HN{a+2,c+1})
	  		                     + (eval(vzt{t+1}) + eval(vzt{i-1})) * eval(cinterp_HN{a+2,c+2}) )
                      * 0.25 * (eval(kh{})     + eval(kh{j+1}));

        /*
        //AI
        eval(flux{i,j,a}) =  ( (eval(vxt{i,j-1,t+1}) + eval(vxt{i,j,t})) * eval(cinterp_HN{i,j,a,c}) 
                             + (eval(vyt{i,j-1,t+1}) + eval(vyt{i,j,t})) * eval(cinterp_HN{i,j,a,c+1}) 
                             + (eval(vzt{i,j-1,t+1}) + eval(vzt{i,j,t})) * eval(cinterp_HN{i,j,a,c+2}) ) 
                      * 0.25 * (eval(kh{i,j,k})+eval(kh{i+1,j,k}));
        
        // AIJ    
        eval(flux{i,j,a+1}) =  ((eval(vxt{i,j,t}) + eval(vxt{i,j,t+1})) * eval(cinterp_HN{i,j,a+1,c})
                              + (eval(vyt{i,j,t}) + eval(vyt{i,j,t+1})) * eval(cinterp_HN{i,j,a+1,c+1}) 
                              + (eval(vzt{i,j,t}) + eval(vzt{i,j,t+1})) * eval(cinterp_HN{i,j,a+1,c+2}) )
                       * 0.25 * (eval(kh{i,j,k})  + eval(kh{i+1,j+1,k}));

        // AIJ    
        eval(flux{i,j,a+2}) = ( (eval(vxt{i,j,t+1}) + eval(vxt{i-1,j,t})  * eval(cinterp_HN{i,j,a+2,c})
                              + (eval(vyt{i,j,t+1}) + eval(vyt{i-1,j,t})) * eval(cinterp_HN{i,j,a+2,c+1})
                              + (eval(vzt{i,j,t+1}) + eval(vzt{i-1,j,t})) * eval(cinterp_HN{i,j,a+2,c+2}) )
                       * 0.25 * (eval(kh{i,j,k})    + eval(kh{i,j+1,k}));
        */
        
    }
}; // end flux_function

struct consume_function {
    
  // flux --> [I, J, A]   
  typedef accessor<0, enumtype::in, extent<-1, 0, -1, 0>, 4 > flux;
  // cinterp_PRA --> [I, J]+[L]
  typedef accessor<1, enumtype::in, extent<0, 0, 0, 0>, 4 > cinterp_PRA;
  // dscl --> [I, J, K]+[L]
  typedef accessor<2, enumtype::inout, extent<0, 0, 0, 0>, 5 > dscl;
    
  typedef boost::mpl::vector<flux, cinterp_PRA, dscl> arg_list;

  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<4>::Index a;
        
        eval(dscl{}) = ( eval(flux{})    - eval(flux{i-1})
                       + eval(flux{a+1}) - eval(flux{i-1,j-1,a+1})
			                 + eval(flux{a+2}) - eval(flux{j-1,a+2}) )
                       * eval(cinterp_PRA{});


        /*
        eval(dscl{i,j,k}) = ( eval(flux{i,j,a})   - eval(flux{i-1,j,a})
                            + eval(flux{i,j,a+1}) - eval(flux{i-1,j-1,a+1})
                            + eval(flux{i,j,a+2}) - eval(flux{i,j-1,a+2}) )
                            * eval(cinterp_PRA{i,j});          
        */
        
    }
}; // end consume_function

// Boundary condition functor
template <typename T>
struct direction_bc_input {
    T value;

    GT_FUNCTION
    direction_bc_input()
        : value(1)
    {}

    GT_FUNCTION
    direction_bc_input(T v)
        : value(v)
    {}

    // relative coordinates
  template <sign I, sign J, sign K, typename DataField0>
    GT_FUNCTION
  void operator()(direction<I,J,K>,
                    DataField0 & data_field0,
                    uint_t i, uint_t j, uint_t k) const {
        data_field0(i,j,k) = value;
    }

}; // end bounrady condition functor

int main(int argc, char const *argv[])
{
  
    using namespace enumtype;

    // [layout_map]
    typedef gridtools::layout_map<1,2,-1,0>   layout_3d_noK;
    typedef gridtools::layout_map<0,1,2>      layout_3d_K;
    typedef gridtools::layout_map<2,3,-1,0,1> layout_4d;
    typedef gridtools::layout_map<0,1>        layout_2d;
    //CUDA case
    typedef gridtools::layout_map<3,0,-1,1>   layout_3d_noK_cuda;
    typedef gridtools::layout_map<3,0,1>      layout_3d_K_cuda;
    typedef gridtools::layout_map<3,0,1,2>    layout_4d_cuda;
    typedef gridtools::layout_map<1,0>        layout_2d_cuda;

    // CUDA backend: use block strategy
#ifdef __CUDACC__
    typedef backend<Cuda,structured,Block> backend_t;
#else
    typedef backend<Host,structured,Block> backend_t;
#endif
    
    typedef backend_t::storage_info<5, layout_3d_noK> storage_info_3d_noK_flux_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_noK_flux_t>::type storage_type_3d_noK_flux;

    typedef backend_t::storage_info<0, layout_3d_noK> storage_info_3d_noK_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_noK_t>::type storage_type_3d_noK;
    
    typedef backend_t::storage_info<1, layout_3d_K> storage_info_3d_K_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_K_t>::type storage_type_3d_K;

    typedef backend_t::storage_info<6, layout_3d_K> storage_info_3d_K_scl_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_K_scl_t>::type storage_type_3d_K_scl;

    typedef backend_t::storage_info<2, layout_4d> storage_info_4d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_4d_t>::type storage_type_4d;

    typedef backend_t::storage_info<7, layout_4d> storage_info_4d_TN_t;
    typedef typename backend_t::storage_type<float_type, storage_info_4d_TN_t>::type storage_type_4d_TN;
     
    typedef backend_t::storage_info<3, layout_2d> storage_info_2d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_2d_t>::type storage_type_2d;

    storage_info_3d_noK_t      metadata_3d_noK(ADM_iall,ADM_jall,1,2);
    storage_info_3d_K_scl_t    metadata_3d_k_scl(ADM_iall+1,ADM_jall+1,ADM_kall);
    storage_info_3d_K_t        metadata_3d_k(ADM_iall,ADM_jall,ADM_kall);
    storage_info_4d_TN_t       metadata_4d_TN(ADM_iall+1,ADM_jall+1,1,3,3);
    storage_info_4d_t          metadata_4d(ADM_iall,ADM_jall,1,3,3);
    storage_info_3d_noK_flux_t metadata_3d_noK_flux(ADM_iall,ADM_jall,1,3);
    storage_info_2d_t          metadata_2d(ADM_iall,ADM_jall);

    
    // laplace functor
    storage_type_3d_noK                            vxt(metadata_3d_noK, 0.0, "vxt");
    storage_type_3d_noK                            vyt(metadata_3d_noK, 0.0, "vyt");
    storage_type_3d_noK                            vzt(metadata_3d_noK, 0.0, "vzt");
    field <storage_type_3d_K_scl,ADM_lall>::type   scl(metadata_3d_k_scl, 0.0, "scl");
    field <storage_type_4d_TN,ADM_lall>::type      cinterp_TN(metadata_4d_TN, 0.0, "cinterp_TN");
    field <storage_type_3d_noK,ADM_lall>::type     cinterp_TRA(metadata_3d_noK, 0.0, "cinterp_TRA");
    // flux functor
    float_type d[16*16*7];
    storage_type_3d_noK_flux                   flux(metadata_3d_noK_flux, (float_type*) d, "flux");
    field <storage_type_3d_K,ADM_lall>::type   kh(metadata_3d_k, 0.0, "kh");
    field <storage_type_4d,ADM_lall>::type     cinterp_HN(metadata_4d, 0.0, "cinterp_HN");
    // consume functor
    field <storage_type_2d,ADM_lall>::type     cinterp_PRA(metadata_2d, 0.0, "cinterp_PRA");
    field <storage_type_3d_K,ADM_lall>::type   dscl(metadata_3d_k, 0.0, "dscl");

    // set intial values to observe data layout
    for(int i=0; i<metadata_3d_k.template dims<0>(); ++i)
      for(int j=0; j<metadata_3d_k.template dims<1>(); ++j)
  	{
  	for(int k=0; k<metadata_3d_k.template dims<2>(); ++k)
  	  {
  	    scl.get_value<0>(i,j,k)=(i+j+k)/2.;
  	    kh.get_value<0>(i,j,k)=i+j+k;
  	  }
  	for(int a=0; a<metadata_4d.template dims<3>(); ++a)
  	  for(int c=0; c<metadata_4d.template dims<4>(); ++c)
  	    {
  	      cinterp_HN.get_value<0>(i,j,0,a,c)=(i+j+c+a)*2.;
  	      cinterp_TN.get_value<0>(i,j,0,a,c)=(i+j+c+a)/2.;
  	    }
  	for(int t=0; t<metadata_3d_noK.template dims<3>(); ++t)
  	  cinterp_TRA.get_value<0>(i,j,0,t)=t+1;
  	cinterp_PRA.get_value<0>(i,j)=i+j+1;
  	}

    scl.print();

    // laplace functor
    typedef arg<0, storage_type_3d_noK>                          p_vxt;
    typedef arg<1, storage_type_3d_noK>                          p_vyt;
    typedef arg<2, storage_type_3d_noK>                          p_vzt;
    typedef arg<3, field <storage_type_3d_K_scl,ADM_lall>::type> p_scl;
    typedef arg<4, field <storage_type_4d_TN,ADM_lall>::type>    p_cinterp_TN;
    typedef arg<5, field <storage_type_3d_noK,ADM_lall>::type>   p_cinterp_TRA;
    // flux functor
    typedef arg<6, storage_type_3d_noK_flux>                     p_flux;
    typedef arg<7, field <storage_type_3d_K,ADM_lall>::type>     p_kh;
    typedef arg<8, field <storage_type_4d,ADM_lall>::type>       p_cinterp_HN;
    // consume fucntor
    typedef arg<9, field <storage_type_2d,ADM_lall>::type >      p_cinterp_PRA;  
    typedef arg<10, field <storage_type_3d_K,ADM_lall>::type>    p_dscl;
    
    // accessor list of all place holders used
    typedef boost::mpl::vector<p_vxt, p_vyt, p_vzt, p_scl, p_cinterp_TN, p_cinterp_TRA,p_flux,p_kh,p_cinterp_HN, p_cinterp_PRA, p_dscl> accessor_list;

    gridtools::domain_type<accessor_list> domain ((p_vxt() = vxt), (p_vyt() = vyt), (p_vzt() = vzt), (p_scl() = scl ), (p_cinterp_TN() = cinterp_TN), (p_cinterp_TRA() = cinterp_TRA), (p_flux() = flux), (p_kh() = kh),(p_cinterp_HN() = cinterp_HN), (p_cinterp_PRA() = cinterp_PRA), (p_dscl() = dscl) );

    // i-2,j-2 ---- i+1, j+1
    // extra lower halo added to account for operaion dependecy
    int halo_size = 1; 
    uint_t di[5] = {halo_size*2, halo_size, halo_size*2, ADM_iall-halo_size-1, ADM_iall};
    uint_t dj[5] = {halo_size*2, halo_size, halo_size*2, ADM_jall-halo_size-1, ADM_jall};
    gridtools::grid<axis> grid(di,dj);
    
    // single splitter
    grid.value_list[0] = 0;
    grid.value_list[1] = ADM_kall-1;

    auto diffusion = make_computation<backend_t>
        (
         domain, grid,
         make_mss
         (
          execute<forward>(), // parallel means no vertical dependency, forward from 0 to K-1, otherwise it is backward
          make_esf<laplace_function>(p_vxt(), p_vyt(), p_vzt(), p_scl(), p_cinterp_TN(), p_cinterp_TRA() ),//define in same order as accessors unique identifier          
          make_esf<flux_function>(p_vxt(), p_vyt(), p_vzt(), p_flux(), p_kh(), p_cinterp_HN()),// in case of several DO functions due to dependinceis, add make_esf for each functor
          make_esf<consume_function>(p_flux(), p_cinterp_PRA(), p_dscl())// in case of several DO functions due to dependinceis, add make_esf for each functor
          )
         );

    // allocate and init storages
    diffusion->ready();

    // copy to device
    diffusion->steady();

    // TODO: ADD LOOP FOR L WHEN CYCLE IS IMPLEMENTED
    //for(int l=0: l < ADM_lall; ++l) {  
        // SWAP HERE
    diffusion->run(); 
    // run boundary condition
    gridtools::array<gridtools::halo_descriptor, 3> halos;
    halos[0] = gridtools::halo_descriptor(0,0,0,ADM_iall-1,ADM_iall);
    halos[1] = gridtools::halo_descriptor(0,0,0,ADM_jall-1,ADM_jall);
    halos[2] = gridtools::halo_descriptor(0,0,0,ADM_kall-1,ADM_kall);
#ifdef __CUDACC__
    gridtools::boundary_apply_gpu<direction_bc_input<float_type> >(halos, direction_bc_input<float_type>(0.0)).apply(dscl);
#else
    gridtools::boundary_apply<direction_bc_input<float_type> >(halos, direction_bc_input<float_type>(0.0)).apply(dscl);
#endif

    // copy data back and deallocate
    diffusion->finalize();

    dscl.print();

    // a CPU side storage to compare results
    // can pass pointer to compare to another data
    //storage_type_6d ref(metadata_vt, 0.0, "ref");
    /*
// compare ref 

    verifier verif(1e-13);
    array<array<uint_t, 2>, 3> halos{{ {halo_size,halo_size}, {halo_size,halo_size}, {halo_size,halo_size} }};
    bool result = verif.verify(grid, vt, ref, halos);

// report timing
#ifdef BENCHMARK
        std::cout << diffusion->print_meter() << std::endl;
#endif

    ASSERT_TRUE(result);
}
*/

    return 0;
}


/**
@}
*/
