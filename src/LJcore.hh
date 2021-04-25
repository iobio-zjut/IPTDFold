
#ifndef INCLUDED_protocols_abinitio_LJcore_hh
#define INCLUDED_protocols_abinitio_LJcore_hh

#include<iostream>
#include<vector>
#include<string>
#include<map>

#include<fstream>
#include<sstream>


#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <protocols/abinitio/LJAngleRotation.fwd.hh>
#include <protocols/abinitio/LJAngleRotation.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>

//#include <protocols/abinitio/ClassicAbinitio.fwd.hh>
//#include <protocols/abinitio/ClassicAbinitio.hh>

using namespace std;
using namespace core;


namespace protocols {
namespace abinitio {

Real random_uniform();
Real random_range(Size const a, Size const b);


///@brief set parameters
extern Size NP_;
extern Size G_;
extern Size Refine_G_;
extern Size frag_;
extern Real KTc_;
extern Real KTd_;
extern Real KTs_;
extern bool Only_Global_;
extern bool use_best_;

extern Real F_helix;
extern Real F_beta;
extern Real F_loop;


extern bool frag_assembly_;
extern bool add_loop_sampling_;


///@brief need object from ClassicAbinitio
extern core::pose::Pose initPose;
extern pair<core::pose::Pose, Real>  best_FA_mover_pose;
extern pair<core::pose::Pose, Real>  best_FM_mover_pose;
extern pair<core::pose::Pose, Real>  best_SS_mover_pose;
extern pair<core::pose::Pose, Real>  best_LS_mover_pose;

extern core::scoring::ScoreFunctionOP score2; //score_stage3a_; //score2
extern core::scoring::ScoreFunctionOP energy_function;

extern simple_moves::FragmentMoverOP FragAssem_;
extern simple_moves::FragmentMoverOP FragAssem_9;
extern simple_moves::FragmentMoverOP FragAssem_3;

extern string sequence;
extern Size proteinLength;


extern vector<core::pose::Pose> Initial_population;
extern vector<core::pose::Pose> population_pose;
extern vector<Real> population_energy;
extern vector<Real> population_Dscore;
extern vector<Real> population_rmsd;

///@brief population_Dscore_map[i][r1-1][r2-1]; 
extern vector<vector<vector<Real> > > population_Dscore_map;
extern vector<vector<Real> > Dscore_map;

extern map<string,Real> parametersMap;

extern multimap<Real, Size> rank_Dscore;



struct SCORE
{
	Real energy_;
	Real dscore_;
	Real rmsd_;
	
	Real Rmsd(){	return rmsd_;	}
	Real Energy(){	return energy_;	}
	Real Dscore(){	return dscore_;	}
	Real total_score(){	return energy_ + dscore_;	}
	
	void Rmsd(Real value){	rmsd_ = value;	}
	void Energy(Real value){	energy_ = value;	}
	void Dscore(Real value){	dscore_ = value;	}
	
	SCORE() {energy_ = 0.0; dscore_ = 0.0; rmsd_ = 0.0;}
	SCORE(Real v1, Real v2) : energy_(v1), dscore_(v2) {rmsd_ = 0.0;}
	SCORE(Real v1, Real v2, Real v3) : energy_(v1), dscore_(v2), rmsd_(v3) {}
};
extern vector<SCORE> population_score;
SCORE cal_all_score(core::pose::Pose &pose);
SCORE cal_Dscore(core::pose::Pose &pose);
Real Energy_score(core::pose::Pose &pose);
Real Distance_score(core::pose::Pose &pose);
void add_FM_refinement();

///@brief save distance pair
struct Dist_unit 
{
	Size r1;
	Size r2;
	Real peak_dist;
	Real peak_dist_prob;
	Real mean_dist;		///Guassian
	Real mean_dist_var;
	
	Dist_unit() {r1 = 0; r2 = 0; peak_dist = 0; peak_dist_prob = 0; mean_dist = 0; mean_dist_var = 0;}
	Dist_unit(Size res1, Size res2, Real peak_distance = 0.0, Real peak_distance_prob = 0.0, Real mean_distance = 0.0, Real mean_distance_var = 0.0) : r1(res1), r2(res2), 
		  peak_dist(peak_distance), peak_dist_prob(peak_distance_prob), mean_dist(mean_distance), mean_dist_var(mean_distance_var) {}
};
extern vector<Dist_unit> Distance_map;


///@brief save distance pair
struct Residue_unit 
{
	numeric::xyzVector<Real> CA_;
	numeric::xyzVector<Real> CB_;
	numeric::xyzVector<Real> N_;	
	
	numeric::xyzVector<Real> CA() {return CA_;}
	numeric::xyzVector<Real> CB() {return CB_;}
	numeric::xyzVector<Real> N()  {return N_; }	
	
	Residue_unit(numeric::xyzVector<Real> ca, numeric::xyzVector<Real> cb, numeric::xyzVector<Real> n) : CA_(ca), CB_(cb), N_(n) {}
};


struct SS_local
{
      char ss_type;
      Size begin;
      Size end;
      Size best_index;	/// the index of pose with best local SS
      Real best_score;
      
      SS_local() {ss_type = ' '; begin = 0; end = 0; best_index = 0; best_score = 0;}
      SS_local(char type, Size be, Size en, Size index = 0, Real score = 0.0) : ss_type(type), begin(be), end(en), best_index(index), best_score(score) {}
};

void read_fasta(string const &route);
void read_distance(string const &route);
vector<SS_local> divide_SS(string SS);
void read_parameters(string const &route);
void get_parameters();

extern LJAngleRotationOP LJAR;


void run_IPTDFold();
void FM_Optimization();
std::pair<bool, vector<core::pose::Pose> > geometry_optimization(core::pose::Pose &pose, Size Begin, Size End);
void Loop_Specific_Geometry_Optimization();


double DM_score( core::pose::Pose &pose, core::pose::Pose &tempPose );
vector< std::vector< Real > > Similarity_Matrix(vector< core::pose::Pose > &population);
vector<pair<Size, vector<Size> > > Kmediods(Size K, std::vector<std::vector<Real> > &Distance_Matrix);


bool boltzmann_accept(const Real targetEnergy, const Real trialEnergy, Real recipocal_KT = 2);


}
}
#endif