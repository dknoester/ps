//
//  stripes.h
//  ealife
//
//  Created by Heather Goldsby on 11/22/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//


#ifndef _EALIFE_STRIPES_H_
#define _EALIFE_STRIPES_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/environment.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/subpopulation_founder.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/mutation.h>
#include <ea/recombination.h>

using namespace ealib;


LIBEA_MD_DECL(STRIPE_FIT, "ea.stripes.fit", int); // count the number of organisms that have the right color stripe
LIBEA_MD_DECL(ANCESTOR, "ea.stripes.ancestor", int);
LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);


//! Stripe fitness.
struct permute_stripes : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    int eval_permute_stripes(EA& ea) {
        // vert stripes
        double one_fit_not = 0;
        double one_fit_nand = 0;
        double two_fit_not = 0;
        double two_fit_nand = 0;
        // horizontal stripes
        double three_fit_not = 0;
        double three_fit_nand = 0;
        double four_fit_not = 0;
        double four_fit_nand = 0;
        // diagonal stripes
        double five_fit_not = 0;
        double five_fit_nand = 0;
        double six_fit_not = 0;
        double six_fit_nand = 0;
        int num_org = 0;
        
        for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
            for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                typename EA::location_type& l = ea.env().location(x,y);
                if (!l.occupied()) {
                    continue;
                }
                
                std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
                // Vertical stripes!
                if((y % 2) == 0) {
                    if (lt == "nand") { ++one_fit_nand; }
                    if (lt == "not") { ++two_fit_not; }
                } else {
                    if(lt == "not") { ++one_fit_not; }
                    if (lt == "nand") { ++two_fit_nand; }
                }
                
                // Horizontal stripes
                if ((x % 2) == 0) {
                    if (lt == "nand") { ++three_fit_nand; }
                    if (lt == "not") { ++four_fit_not; }
                } else {
                    if(lt == "not") { ++three_fit_not; }
                    if (lt == "nand") { ++four_fit_nand; }
                }
                
                // Diagonal stripes
                if (((x % 2) == 0) && ((y % 2) == 0 )) {
                    if(lt == "not") { ++five_fit_not; }
                    if (lt == "nand") { ++six_fit_nand; }
                } else {
                    if(lt == "nand") { ++five_fit_nand; }
                    if (lt == "not") { ++six_fit_not; }
                }
            }
        }
        
        double tmp_one_fit = (one_fit_not + 1)  * (one_fit_nand + 1);
        double tmp_two_fit = (two_fit_not + 1)  * (two_fit_nand + 1);
        double tmp_three_fit = (three_fit_not + 1)  * (three_fit_nand + 1);
        double tmp_four_fit = (four_fit_not + 1)  * (four_fit_nand + 1);
        double tmp_five_fit = (five_fit_not + 1)  * (five_fit_nand + 1);
        double tmp_six_fit = (six_fit_not + 1)  * (six_fit_nand + 1);
        double tmp_fit = std::max(tmp_one_fit, tmp_two_fit);
        tmp_fit = std::max(tmp_fit, tmp_three_fit);
        tmp_fit = std::max(tmp_fit, tmp_four_fit);
        tmp_fit = std::max(tmp_fit, tmp_five_fit);
        tmp_fit = std::max(tmp_fit, tmp_six_fit);
        
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_stripes(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};






/*! Compete to evolve stripes -- even number rows nand; odd number rows not
 */
//template <typename EA>
//struct permute_stripes : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
//    permute_stripes(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("permute_stripes.dat") {
//        _df.add_field("update")
//        .add_field("mean_fitness")
//        .add_field("max_fitness");
//    }
//    
//    virtual ~permute_stripes() {
//    }
//    
//    virtual void operator()(EA& ea) {
//        using namespace boost::accumulators;
//        accumulator_set<double, stats<tag::mean, tag::max> > fit;
//
//        
//        // calculate "fitness":
//        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
//           
//            eval_permute_stripes(i);
//            
//            // copy the stripe fit to the accumulator and also the subpop
//            double sf =get<STRIPE_FIT>(i);
//            fit(sf);
//            put<STRIPE_FIT>(sf, *i);
//            
//        }
//        
//        
//        _df.write(ea.current_update())
//        .write(mean(fit))
//        .write(max(fit))
//        .endl();
//        
//        std::size_t n=get<METAPOPULATION_SIZE>(ea);
//        typename EA::population_type offspring; // container of (pointers to) subpopulations
//        recombine_n(ea.population(), offspring,
//                    selection::tournament < access::metadata<STRIPE_FIT> > (n, ea.population(), ea),
//                    recombination::propagule_without_replacement(),
//                    n, ea);
//        
//        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
//        int s = get<POPULATION_SIZE>(ea);
//
//        // Mutate and fill each offspring group.
//        for(typename EA::population_type::iterator i=offspring.begin(); i!=offspring.end(); ++i) {
////            // for notational ease, get a ref to the subpopulation:
////            typename EA::subpopulation_type& sea = **i;
////            assert(sea.size() == 1);
////
////            // clear founders...
////            founder(sea).clear();
////            
////            // mutate it:
////            mutate(sea.population().begin(), m, sea)
////            //            mutate(**((*i)->population().begin()),m,*i);
////            typename EA::subpopulation_type::individual_type& g = *sea.begin();
////            //            typename EA::subpopulation_type::individual_type g = (**((*i)->population().begin()));
////            
////            // add first org as founder
////            (*i)->founder().insert((*i)->founder().end(), (*i)->copy_individual(g));
////            
////            // and fill up the offspring population with copies of the germ:
////            for (int k=1; k<get<NUM_PROPAGULE_CELL>(ea); ++k) {
////                typename EA::individual_type::ea_type::individual_ptr_type o = (*i)->ea().copy_individual(g);
////                (*i)->insert((*i)->end(), o);
////                
////                // move to random location
////                std::size_t pos = (*i)->ea().rng()(s);
////                (*i)->ea().env().move_ind(k, pos);
////                
////                // add org as founders
////                (*i)->ea().founder().insert((*i)->ea().founder().end(), (*i)->ea().copy_individual(*o));
////            }
//        }
//                
//        // swap populations
//        std::swap(ea.population(), offspring);
//        
//        
//    }
//    
//    datafile _df;
//};

#endif
