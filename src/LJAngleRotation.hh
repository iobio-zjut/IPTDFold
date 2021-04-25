
#ifndef INCLUDED_protocols_abinitio_LJAngleRotation_hh
#define INCLUDED_protocols_abinitio_LJAngleRotation_hh

#include<vector>
#include<map>
#include<string>
#include <core/types.hh>
#include <math.h>
#include <numeric/xyzVector.hh>

using namespace std;
using namespace core;


namespace protocols {
namespace abinitio {

class LJAngleRotation {
public:
	LJAngleRotation(){ NP2 = 50; G2 = 50; F = 0.5; CR = 0.5; KTl = 1; KT_reciprocal = 1; Max_Disturbance = 5; Use_best = 1; N_Candidate = 10; Greedy_strategy = 1; }
		
	void get_parameters(map<string,Real> &parametersMap);
		
	Real score(vector<Real> &rotationAngles);
	
	Real difference(vector<Real> &Angles_1, vector<Real> &Angles_2);
	
	pair<bool, vector<vector<Real> > > angle_learning();
	
	Size dimension(){
		return rotationAxis.size();
	};
	
	void rotation_parameters(vector<pair<numeric::xyzVector<Real>, numeric::xyzVector<Real> > > &rotation_axis, vector<numeric::xyzVector<Real> > &fixed_points, vector<numeric::xyzVector<Real> > &rotation_points, vector<pair<pair<Size, Size>, pair<Real, Real> > > &restraint_dist_map); //, vector<pair<pair<Size, Size>, Real > > &L_contact_index_conf);
  
private:

	///@note set of rotation porint
	vector<numeric::xyzVector<Real> > rotationPoints;
	
	///@note set of fixed porint
	vector<numeric::xyzVector<Real> > fixedPoints;
	
	///@note distance matrix of objective function
	vector<pair<pair<Size, Size>, pair<Real, Real> > > restraint_distance_map;
	
	
	///@note <base_point, rotation_axis>
	vector<pair<numeric::xyzVector<Real>, numeric::xyzVector<Real> > > rotationAxis;
	
private:
	Size NP2;
	Size G2;
	Real F;
	Real CR;
	Real KTl;
	Real KT_reciprocal;
//	Real Driver_Angle;
	Real Max_Disturbance;
	Real Use_best;
	Size N_Candidate;
	Real Greedy_strategy;
//	Size Convergence_G;
	
public:
	static Real TrialScore;
	static bool isSimilar(Real &score){
		return ( fabs(score - TrialScore) < 0.00001 );
	}
};

} //abinitio
} //protocols

#endif
