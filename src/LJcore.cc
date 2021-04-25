#include <protocols/abinitio/LJcore.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>

#include<vector>
#include<algorithm>
#include<math.h>
#include<cmath>
#include<iomanip>
#include<stdlib.h>
#include<time.h>
#include<fstream>
#include<iomanip>

#include<sys/stat.h>
#include<sys/types.h>
#include<unistd.h>

#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/simple_moves/SymmetricFragmentMover.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

using namespace std;
using namespace core::scoring::dssp;

namespace protocols {
namespace abinitio {

Real random_uniform(){	return double(rand() % 10001) / 10000;}

Real random_range(Size const a, Size const b){	return (rand() % (b-a+1)) + a;}

///@brief set parameters
Size NP_;
Size G_;
Size Refine_G_;
Size frag_;
Real KTd_ = 2;

Real F_helix = 0.3;
Real F_beta = 0.1;
Real F_loop = 0.5;

bool frag_assembly_ = false;
bool add_loop_sampling_ = false;
bool use_best_ = false;


///@brief need object from ClassicAbinitio
core::pose::Pose initPose;
pair<core::pose::Pose, Real>  best_FA_mover_pose;
pair<core::pose::Pose, Real>  best_FM_mover_pose;
pair<core::pose::Pose, Real>  best_SS_mover_pose;
pair<core::pose::Pose, Real>  best_LS_mover_pose;

core::scoring::ScoreFunctionOP score2;
core::scoring::ScoreFunctionOP energy_function;

simple_moves::FragmentMoverOP FragAssem_;
simple_moves::FragmentMoverOP FragAssem_9;
simple_moves::FragmentMoverOP FragAssem_3;

string sequence;
Size proteinLength;

vector<core::pose::Pose> Initial_population;
vector<core::pose::Pose> population_pose;
vector<SCORE> population_score;
vector<Real> population_energy;
vector<Real> population_Dscore;
vector<Real> population_rmsd;

vector<vector<vector<Real> > > population_Dscore_map;
vector<vector<Real> > Dscore_map(proteinLength, vector<Real>(proteinLength, 0.0));

vector<Dist_unit> Distance_map;

map<string,Real> parametersMap;


multimap<Real, Size> rank_Dscore;

LJAngleRotationOP LJAR( new LJAngleRotation );


///==========================================================================================================================

void read_fasta(string const &route){
	ifstream sequence_file( route.c_str() );
	if ( sequence_file.is_open() ){
		string line;
		while( getline(sequence_file, line) ){
			if ( line[0] != '>' )
				sequence += line;
		}
		sequence_file.close();
		
		for (Size i = 0; i < sequence.size(); ++i){
		    if ( sequence[i] == ' ' ){
			    sequence.erase(i,1);
			    --i;
		    }
		}
		cout << "Sequence: " << sequence.size() << endl;
		cout << sequence << endl;
	}
	else{
		cout << "ERROR: core::getFasta(): can not open sequence file. " << endl;
		exit(0);
	}
}

void read_distance(string const &route){
	ifstream read_distance( route.c_str() );
	if ( !read_distance ){
		cout << "ERROR: cannot open 'distance map' in route './filtered_dmap.txt'" << endl;
		exit(0);
	}
	string distance_line;
	while (getline(read_distance, distance_line)){
		istringstream Data(distance_line);
		
		Size r1, r2;
		Real peak_dist, peak_dist_prob, mean_dist, mean_dist_var;
		Data >> r1 >> r2 >> peak_dist >> peak_dist_prob >> mean_dist >> mean_dist_var;
		
		if (r1 < r2){
			Dist_unit distance(r1, r2, peak_dist, peak_dist_prob, mean_dist, mean_dist_var);
			Distance_map.push_back(distance);
		}
	}
	read_distance.close();
	if ( Distance_map.size() == 0 ){
		cout << "ERROR!!!\n" << "ERROR : Read Distance failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}else
		cout << "# Read distace pair: " << Distance_map.size() << endl;
}

void read_parameters(string const &route){
	ifstream informParam( route.c_str() );
	if ( !informParam ){
		cout << "ERROR: cannot open 'parameters file' in route '" << route << "'" << endl;
		exit(0);
	}
	
	string line;
	while (getline(informParam,line)){
		if ( line[0] != '#' && line[0] != ' ' ){
			istringstream line_data( line );
			string parameterName;
			Real parameterValue;
			line_data >> parameterName >> parameterValue;
			parametersMap.insert(make_pair(parameterName, parameterValue));
		}
	}
	informParam.close();
}

void get_parameters(){
	if (parametersMap.find("NP") != parametersMap.end()){
		 NP_ = parametersMap["NP"];
		 cout << "# Using user defined NP: " << NP_ << endl;
	}
	else{
		 cout << "# Using default NP: " << NP_ << endl;
	}
	
	if (parametersMap.find("G") != parametersMap.end()){
		 G_ = parametersMap["G"];
		 cout << "# Using user defined G: " << G_ << endl;
	}
	else{
		cout << "# Using default G: " << G_ << endl;
	}
	
	if (parametersMap.find("FA") != parametersMap.end()){
		frag_assembly_ = parametersMap["FA"];
		cout << "# Whether using Fragment_Assembly sampling: " << frag_assembly_ << endl;
	}
	
	if (parametersMap.find("LS") != parametersMap.end()){
		add_loop_sampling_ = parametersMap["LS"];
		cout << "# Whether using Loop_Specific sampling" << add_loop_sampling_ << endl;
	}
	
	if (parametersMap.find("Refine_G") != parametersMap.end()){
		Refine_G_ = parametersMap["Refine_G"];
		cout << "# Using refinement generation: " << Refine_G_ << endl;
	}
	
	if (parametersMap.find("use_best") != parametersMap.end()){
		use_best_ = parametersMap["use_best"];
		cout << "# Whether using Loop_Specific sampling" << add_loop_sampling_ << endl;
	}
}

SCORE cal_all_score(core::pose::Pose &pose){
	SCORE Score(Energy_score( pose ), Distance_score( pose ));
	return Score;
}

SCORE cal_Dscore(core::pose::Pose &pose){
//	SCORE Score(Energy_score( pose ), Distance_score( pose ));
	SCORE Score(0.0, Distance_score( pose ));
	return Score;
}

Real Energy_score(core::pose::Pose &pose){
	return (*energy_function)(pose);
}

Real Distance_score(core::pose::Pose &pose){
	Real Dscore( 0 );
	vector<vector<Real> > score_map(proteinLength, vector<Real>(proteinLength, 0.0));
	for (Size i = 0; i < Distance_map.size(); ++i){
		Size residue1( Distance_map[i].r1 );
		Size residue2( Distance_map[i].r2 );
		Real pre_dist ( Distance_map[i].mean_dist );
		Real pre_var ( Distance_map[i].mean_dist_var );
		
		Real distance( 0 );
		if ( pose.residue(residue1).name3() == "GLY" || pose.residue(residue2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(residue1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(residue2).xyz("CA");
			distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(residue1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(residue2).xyz("CB");
			distance = CB1.distance(CB2);
		}
		
		Real score = log( (distance - pre_dist)*(distance - pre_dist)  + 1 ) / sqrt(pre_var);
		Real clash_score( 0 );
		if (distance < 3.8)
			clash_score = (distance - 3.8)*(distance - 3.8);
		Real new_score = score + clash_score;
		
		score_map[residue1-1][residue2-1] = score;
		score_map[residue2-1][residue1-1] = score;
		
		Dscore += score;
	//	Dscore += log( (distance - pre_dist)*(distance - pre_dist)/pre_var  + 1 );
	}
	Dscore_map = score_map;
	return Dscore;
}




vector<SS_local> divide_SS(string sec_struct){
	vector<SS_local> divide_region;
	Size loop_length(0);
	Size H_length(0);
	Size E_length(0);
	for (Size i = 0; i < sec_struct.size(); ++i){
		if ( sec_struct[i] == 'L' )
			++loop_length;
		else if (loop_length != 0){
			SS_local ss_('L', i+1-loop_length, i, 0, 0.0);
			divide_region.push_back(ss_);
			loop_length = 0;
		}
		if ( sec_struct[i] == 'H' )
			++H_length;
		else if (H_length != 0){
			SS_local ss_('H', i+1-H_length, i, 0, 0.0);
			divide_region.push_back(ss_);
			H_length = 0;
		}
		if ( sec_struct[i] == 'E' )
			++E_length;
		else if (E_length != 0){
			SS_local ss_('E', i+1-E_length, i, 0, 0.0);
			divide_region.push_back(ss_);
			E_length = 0;
		}
		if ( i == sec_struct.size()-1 ){
			if (loop_length != 0){
				SS_local ss_('L', i+2-loop_length, i+1, 0, 0.0);
				divide_region.push_back(ss_);
				loop_length = 0;
			}
			else if(H_length != 0){
				SS_local ss_('H', i+2-H_length, i+1, 0, 0.0);
				divide_region.push_back(ss_);
				H_length = 0;
			}
			else if(E_length != 0){
				SS_local ss_('E', i+2-E_length, i+1, 0, 0.0);
				divide_region.push_back(ss_);
				E_length = 0;
			}
		}
	}
	return divide_region;
}
	
	
bool boltzmann_accept(const Real targetEnergy, const Real trialEnergy, Real recipocal_KT){
	if ( trialEnergy <= targetEnergy )
		return true;
	else{
		Real prob = exp( -(trialEnergy - targetEnergy) / recipocal_KT );
		if ( prob >= random_uniform() )
			return true;
		else
			return false;
	}
}


void Loop_Specific_Geometry_Optimization(){
	for (Size i = 0; i < NP_; ++i){
	//	cout << "   pose: " << i << " |   ";
		///@brief get secondary structure regions
		Dssp dssp( population_pose[i] );
		string SS_pose( dssp.get_dssp_secstruct() );
		vector<SS_local> SS_elememts = divide_SS( SS_pose );
				
		///@brief Iteratively select each loop region
		for (Size ss_index = 0; ss_index < SS_elememts.size(); ++ss_index){
			Size loopIndex( 0 );
			if ( SS_elememts[ss_index].ss_type == 'L' && (SS_elememts[ss_index].end - SS_elememts[ss_index].begin) >= 2 )
				loopIndex = ss_index;
			else
				continue;
			
			///@brief split loop to several fragment if it is too long
			vector<std::pair<Size, Size> > vec_B_E;
			if ( SS_elememts[loopIndex].end - SS_elememts[loopIndex].begin + 1 > 15 ){
				for (Size ss = SS_elememts[loopIndex].begin; ss <= SS_elememts[loopIndex].end; ss+=3){
					Size fragBegin = ss;
					Size fragEnd = fragBegin + 8;
					if ( fragEnd > SS_elememts[loopIndex].end )
						fragEnd = SS_elememts[loopIndex].end;
					if (fragBegin <= 1)
						fragBegin = 2;
					vec_B_E.push_back( make_pair(fragBegin, fragEnd) );
					
					if ( fragEnd >= SS_elememts[loopIndex].end )
						break;
				}
			}else{
				Size fragBegin = SS_elememts[loopIndex].begin;
				if (fragBegin <= 1)		
					fragBegin = 2;
				Size fragEnd = SS_elememts[loopIndex].end;
				vec_B_E.push_back( make_pair(fragBegin, fragEnd) );
			}
			
			for (Size f = 0; f < vec_B_E.size(); ++f){
			  
				Size loopBegin( vec_B_E[f].first );
				Size loopEnd( vec_B_E[f].second );
								
				std::pair<bool, vector<core::pose::Pose> > trial_pose_set = geometry_optimization(population_pose[i], loopBegin, loopEnd);
				if ( trial_pose_set.first ){
					///@note comformation update
					for (Size k = 0; k < trial_pose_set.second.size(); ++k){
						core::pose::Pose trial_pose = trial_pose_set.second[k];
						SCORE trial_score = cal_all_score( trial_pose );
						if ( boltzmann_accept(population_score[i].total_score(), trial_score.total_score(), KTd_) ){
							population_pose[i] = trial_pose;
							population_score[i] = trial_score;
							
							
							if (population_score[i].total_score() < best_LS_mover_pose.second){
								best_LS_mover_pose.first = population_pose[i];
								best_LS_mover_pose.second = population_score[i].total_score();
							}
							for ( multimap<Real, Size>::iterator iter = rank_Dscore.begin(); iter != rank_Dscore.end(); ++iter){
								if (iter->second == i){
									rank_Dscore.erase( iter );
									break;
								}
							}
							rank_Dscore.insert( make_pair(population_score[i].total_score(), i) );
							
							break;
						}
					}
				}
			}
		}
	//	cout << endl;
	}
}

void FM_Optimization(){
	///@note clustering
	cout << "      Model clustering ..." << endl;
	///@note calculate population similarity
	vector< std::vector< Real > > pose_similarity_matrix = Similarity_Matrix(population_pose);

	///@note clustring
	vector<pair<Size, vector<Size> > > cluster_mediod_elements = Kmediods(5, pose_similarity_matrix);
	///@note selecte represent poses to FM
	vector<std::pair<Size, core::pose::Pose> > selected_poses; 
	cout << "      selected representive pose by cluster: ";
	for (Size i = 0; i < cluster_mediod_elements.size(); ++i){
		Size index_mediod = cluster_mediod_elements[i].first;
		Size index_best_score = 0;
		Real best_score = 1000000.0;
		for (Size j = 0; j < cluster_mediod_elements[i].second.size(); ++j){
			Size n = cluster_mediod_elements[i].second[j];
			if (n != index_mediod && population_score[n].total_score() < best_score){
				index_best_score = n;
				best_score = population_score[n].total_score();
			}
		}
		cout << index_mediod << " " << index_best_score << " ";
		
		selected_poses.push_back( make_pair(index_mediod, population_pose[index_mediod]) );
		selected_poses.push_back( make_pair(index_best_score, population_pose[index_best_score]) );
	}
	cout << endl;

	///===========================================================================================================
	///@note Refinement ......
	cout << "      Residue dihedral angle optimization for the selected representive..." << endl;
	for (Size i = 0; i < selected_poses.size(); ++i){
		core::pose::Pose target_pose = selected_poses[i].second;
		SCORE target_score = cal_all_score( target_pose );
		vector<vector<Real> > target_Dscore_map = Dscore_map;
		for (Size g = 0; g < Refine_G_; ++g){
			Size Begin, End;
			if (g < 0.5*Refine_G_){
				///@note calculate each residues score that indicate it need to be move
				multimap<Real, Size, greater<double> > Dscore_Index;
				for (Size r = 0; r < proteinLength; ++r){
					Real ave_score( 0 );
					Size count_score( 0 );
					if (r != 0 && r != proteinLength-1){
						for (Size i = 0; i <= r; ++i){
							for (Size j = r; j < proteinLength; ++j){
								if ( target_Dscore_map[i][j] != 0.0 ){
									count_score += 1;
									ave_score += target_Dscore_map[i][j];
								}
							}
						}
						ave_score /= count_score;
					}
					Dscore_Index.insert( make_pair(ave_score, r+1) );
				}
				///@note selecte Top 20 residues to adjust
				vector<std::pair<Size, Real> > top20_residues;
				Size count( 0 );
				for ( multimap<Real, Size>::iterator iter = Dscore_Index.begin(); count < 20; ++iter, ++count){
					top20_residues.push_back( make_pair(iter->second, iter->first) );
				}
				///@note calculate adjust probability and build roulette
				vector<std::pair<Size, Real> > top20_residue_move_prob;
				Real sum = 0;
				for (Size k = 0; k < top20_residues.size(); ++k)
					sum += top20_residues[k].second;
				for (Size k = 0; k < top20_residues.size(); ++k){
					  Real local_sum = 0;
					  for (Size m = 0; m <= k; ++m){
						  local_sum += top20_residues[m].second;
					  }
					  top20_residue_move_prob.push_back( make_pair(top20_residues[k].first, local_sum / sum) );
				}
				
				///@note selecte a residue to adjust according to probability
				Size perturb_index(0);
				for (Size m = 0; m < 1000; ++m){
					Real prob = random_uniform();
					for (Size k = 0; k < top20_residue_move_prob.size(); ++k){
						if (k == 0 && prob < top20_residue_move_prob[k].second){
							perturb_index = top20_residue_move_prob[k].first;
							break;
						}
						else if (top20_residue_move_prob[k-1].second <= prob && prob <= top20_residue_move_prob[k].second){
							perturb_index = top20_residue_move_prob[k].first;
							break;
						}
					}
					if (perturb_index != 0)
						break;
				}		
				if ( perturb_index == 0 ){
				//	cout << "Warning: perturb index error" << endl;
					continue;
				}
			
				///=====================================================================================================
				///@note define move region (Begin, End)
				if (g < 0.25*Refine_G_){
					if (perturb_index <= 3){
						Begin = 2; End = 6;
					}else if (perturb_index >= proteinLength - 1){
						Begin = proteinLength - 4; End = proteinLength;
					}else{
						Begin = perturb_index - 2; End = perturb_index + 2;
					}
				}else{
					if (perturb_index <= 2){
						Begin = 2; End = 4;
					}else if (perturb_index >= proteinLength){
						Begin = proteinLength - 2; End = proteinLength;
					}else{
						Begin = perturb_index - 1; End = perturb_index + 1;
					}
				}
			}else{
				if (g < 0.75*Refine_G_){
					Begin = random_range(2, proteinLength-4);
					End = Begin+4;
				}else{
					Begin = random_range(2, proteinLength-2);
					End = Begin+2;
				}
			}
			
			///@note adustment...
			std::pair<bool, vector<core::pose::Pose> > trial_pose_set = geometry_optimization(target_pose, Begin, End);
			
			if ( trial_pose_set.first ){
				///@note comformation update
				for (Size k = 0; k < trial_pose_set.second.size(); ++k){
					core::pose::Pose trial_pose = trial_pose_set.second[k];
					SCORE trial_score = cal_all_score( trial_pose );
					
					if (target_score.total_score() > trial_score.total_score()){
						target_pose = trial_pose;
						target_score = trial_score;
						target_Dscore_map = Dscore_map;
						
						break;
					}
				}	
			}
		}
		
		Size pose_index = selected_poses[i].first;
		population_pose[pose_index] = target_pose;
		target_score = cal_all_score( target_pose );
		population_score[pose_index] = target_score;
		
		
		if (population_score[pose_index].total_score() < best_FM_mover_pose.second){
			best_FM_mover_pose.first = population_pose[pose_index];
			best_FM_mover_pose.second = population_score[pose_index].total_score();
		}
		for ( multimap<Real, Size>::iterator iter = rank_Dscore.begin(); iter != rank_Dscore.end(); ++iter){
			if (iter->second == pose_index){
				rank_Dscore.erase( iter );
				break;
			}
		}
		rank_Dscore.insert( make_pair(population_score[pose_index].total_score(), pose_index) );
	}
}

std::pair<bool, vector<core::pose::Pose> > geometry_optimization(core::pose::Pose &pose, Size Begin, Size End){
	///@brief extract partial distance
	vector<std::pair<std::pair<Size, Size>, std::pair<Real, Real> > > L_distance_constrant;
	vector<numeric::xyzVector<Real> > res_left_xyz, res_right_xyz;
	vector<Size> res_left, res_right;
	Size first_right_resi_index( 0 );
	for (Size c = 0; c < Distance_map.size(); ++c){
		if (Distance_map[c].r1 <= End && Distance_map[c].r2 >= Begin){
			Size r1 = Distance_map[c].r1;
			Size r2 = Distance_map[c].r2;
			std::pair<Size, Size> pair_res_index;
			std::pair<Real, Real> dist_mean_var(Distance_map[c].mean_dist, Distance_map[c].mean_dist_var);
			
			vector<Size>::iterator find_iter = find(res_left.begin(), res_left.end(), r1);
			if ( find_iter == res_left.end() ){
				res_left.push_back( r1 );
				
				if ( pose.residue(r1).name3() == "GLY" ){
					numeric::xyzVector<Real> CA = pose.residue(r1).xyz("CA");
					res_left_xyz.push_back( CA );
				}else{
					numeric::xyzVector<Real> CB = pose.residue(r1).xyz("CB");
					res_left_xyz.push_back( CB );
				}
				
				pair_res_index.first = res_left.size() - 1;
			}else{
				pair_res_index.first = find_iter - res_left.begin();
			}
			
			vector<Size>::iterator find_iter2 = find(res_right.begin(), res_right.end(), r2);
			if ( find_iter2 == res_right.end() ){
				res_right.push_back( r2 );
				
				if ( pose.residue(r2).name3() == "GLY" ){
					numeric::xyzVector<Real> CA = pose.residue(r2).xyz("CA");
					res_right_xyz.push_back( CA );
				}else{
					numeric::xyzVector<Real> CB = pose.residue(r2).xyz("CB");
					res_right_xyz.push_back( CB );
				}
				
				pair_res_index.second = res_right.size() - 1;
			}else{
				pair_res_index.second = find_iter2 - res_right.begin();
			}
			
			L_distance_constrant.push_back( make_pair(pair_res_index, dist_mean_var) );
		}
	}
	
	///@brief extract rotation axis set
	vector<std::pair<numeric::xyzVector<Real>, numeric::xyzVector<Real> > > rotationAxis;
	bool add_a_angle = true;
	if ( add_a_angle ){
		numeric::xyzVector<Real> CA = pose.residue(Begin - 1).xyz("CA");
		numeric::xyzVector<Real> C = pose.residue(Begin - 1).atom("C").xyz();
		
		double x = (C.x() - CA.x()) / sqrt( pow((C.x() - CA.x()), 2) + pow((C.y() - CA.y()), 2) + pow((C.z() - CA.z()), 2));
		double y = (C.y() - CA.y()) / sqrt( pow((C.x() - CA.x()), 2) + pow((C.y() - CA.y()), 2) + pow((C.z() - CA.z()), 2));
		double z = (C.z() - CA.z()) / sqrt( pow((C.x() - CA.x()), 2) + pow((C.y() - CA.y()), 2) + pow((C.z() - CA.z()), 2));
		
		numeric::xyzVector<Real> axis(x, y, z);
		rotationAxis.push_back( make_pair(CA, axis) );
	}
	for (Size r = Begin; r <= End; ++r){
		numeric::xyzVector<Real> N = pose.residue(r).atom("N").xyz();
		numeric::xyzVector<Real> CA = pose.residue(r).xyz("CA");
		numeric::xyzVector<Real> C = pose.residue(r).atom("C").xyz();
		
		double x_1 = (CA.x() - N.x()) / sqrt( pow((CA.x() - N.x()), 2) + pow((CA.y() - N.y()), 2) + pow((CA.z() - N.z()), 2));
		double y_1 = (CA.y() - N.y()) / sqrt( pow((CA.x() - N.x()), 2) + pow((CA.y() - N.y()), 2) + pow((CA.z() - N.z()), 2));
		double z_1 = (CA.z() - N.z()) / sqrt( pow((CA.x() - N.x()), 2) + pow((CA.y() - N.y()), 2) + pow((CA.z() - N.z()), 2));
		
		numeric::xyzVector<Real> axis_1(x_1, y_1, z_1);
		rotationAxis.push_back( make_pair(N, axis_1) );
		
		double x_2 = (C.x() - CA.x()) / sqrt( pow((C.x() - CA.x()), 2) + pow((C.y() - CA.y()), 2) + pow((C.z() - CA.z()), 2));
		double y_2 = (C.y() - CA.y()) / sqrt( pow((C.x() - CA.x()), 2) + pow((C.y() - CA.y()), 2) + pow((C.z() - CA.z()), 2));
		double z_2 = (C.z() - CA.z()) / sqrt( pow((C.x() - CA.x()), 2) + pow((C.y() - CA.y()), 2) + pow((C.z() - CA.z()), 2));
		
		numeric::xyzVector<Real> axis_2(x_2, y_2, z_2);
		rotationAxis.push_back( make_pair(CA, axis_2) );
	}
	
	///@brief solving the rotation angle by differential evolution algorithm
	LJAR->rotation_parameters(rotationAxis, res_left_xyz, res_right_xyz, L_distance_constrant);
	std::pair<bool, vector<vector<Real> > > output = LJAR->angle_learning();
	vector<vector<Real> > Top10_Disturbance_Angles = output.second;
	
	///@brief generate candidate poses
	std::pair<bool, vector<core::pose::Pose> > trial_pose_set;
	vector<core::pose::Pose> candidate_pose;
	if ( output.first ){
		for (Size a = 0; a < Top10_Disturbance_Angles.size(); ++a){
			core::pose::Pose trial_pose( pose );
			Size angle(0);
			trial_pose.set_psi(Begin-1, trial_pose.psi(Begin-1) + Top10_Disturbance_Angles[a][angle]);
			++angle;
			for (Size b = Begin; b <= End; ++b){
				trial_pose.set_phi(b, trial_pose.phi(b) + Top10_Disturbance_Angles[a][angle]);
				++angle;
				trial_pose.set_psi(b, trial_pose.psi(b) + Top10_Disturbance_Angles[a][angle]);
				++angle;
			}
			candidate_pose.push_back( trial_pose );
		}
		trial_pose_set.first = true;
		trial_pose_set.second = candidate_pose;
	}else{
		trial_pose_set.first = false;
		trial_pose_set.second = candidate_pose;
	}
	
	return trial_pose_set;
}


void run_IPTDFold(){
	cout << "=========================================================================" << endl;
	cout << "=========================================================================" << endl;
	
	if (proteinLength > 100){
		G_ = 10 * sqrt(proteinLength);
		Refine_G_ = 50 * sqrt(proteinLength);		
	}
	energy_function = score2;

	population_pose = Initial_population;
	
	for (Size i = 0; i < NP_; ++i){
		SCORE Score = cal_all_score( population_pose[i] );
		population_score.push_back( Score );
		rank_Dscore.insert( make_pair(population_score[i].total_score(), i) );
	}
	
	best_FA_mover_pose.first = population_pose[0];
	best_FA_mover_pose.second = population_score[0].total_score();
	best_SS_mover_pose.first = population_pose[0];
	best_SS_mover_pose.second = population_score[0].total_score();
	best_LS_mover_pose.first = population_pose[0];
	best_LS_mover_pose.second = population_score[0].total_score();
	best_FM_mover_pose.first = population_pose[0];
	best_FM_mover_pose.second = population_score[0].total_score();
	
	
	for (Size g = 1; g <= G_; ++g){
		cout << "Generation: " << g << "   lowest-potential pose index: " << rank_Dscore.begin()->second << endl;
		
		cout << "   Partition sampling..." << endl;
		for (Size i = 0; i < NP_; ++i){
			///@brief ================================ Fragment insertion ================================
			for (Size t = 0; t < 100; ++t){
				core::pose::Pose trial_pose = population_pose[i];
				if (g <= 0.5 * G_)
					FragAssem_9->apply( trial_pose );
				else
					FragAssem_3->apply( trial_pose );
				
				SCORE trial_score = cal_all_score( trial_pose );
				if ( boltzmann_accept(population_score[i].total_score(), trial_score.total_score(), KTd_) ){
					population_pose[i] = trial_pose;
					population_score[i] = trial_score;
					
					if (population_score[i].total_score() < best_FA_mover_pose.second){
						best_FA_mover_pose.first = trial_pose;
						best_FA_mover_pose.second = population_score[i].total_score();
					}
					for ( multimap<Real, Size>::iterator iter = rank_Dscore.begin(); iter != rank_Dscore.end(); ++iter){
						if (iter->second == i){
							rank_Dscore.erase( iter );
							break;
						}
					}
					rank_Dscore.insert( make_pair(population_score[i].total_score(), i) );
				}
			}
			
			///@brief =================== Local dihedrall angle Mutation & Crossover ========================
			core::pose::Pose target_pose( population_pose[i] );
			SCORE target_score( population_score[i] );
			bool local_mover_success( false );
			Dssp dssp( population_pose[i] );
			string SS_pose( dssp.get_dssp_secstruct() );
			vector<SS_local> SS_elememts = divide_SS( SS_pose );
			for (Size s = 0; s < SS_elememts.size(); ++s){
				vector<pair<Size, Size> > vec_B_E;
				if (SS_elememts[s].ss_type == 'L' && (SS_elememts[s].end - SS_elememts[s].begin + 1 > 15)){
					cout << "SSmover: loop is too long, split it for ";
					for (Size ss = SS_elememts[s].begin; ss <= SS_elememts[s].end; ss+=3){
						Size fragBegin = ss;
						Size fragEnd = fragBegin + 8;
						if ( fragEnd > SS_elememts[s].end )
							fragEnd = SS_elememts[s].end;
						vec_B_E.push_back( make_pair(fragBegin, fragEnd) );
						
						if ( fragEnd >= SS_elememts[s].end )
							break;
					}
					cout << vec_B_E.size() << " fragments with length of 9." << endl;
				}else{
					vec_B_E.push_back( make_pair(SS_elememts[s].begin, SS_elememts[s].end) );
				}
				for (Size f = 0; f < vec_B_E.size(); ++f){
					Size ss_begin = vec_B_E[f].first;
					Size ss_end = vec_B_E[f].second;
					
					///@note mutation					
					vector<Size> Top15_D;
					Size count = 0;
					for ( multimap<Real, Size>::iterator iter = rank_Dscore.begin(); count < 15; ++iter, ++count){
						Top15_D.push_back( iter->second );
					}
					Size rand( random_range(1, Top15_D.size() - 1) );
					Size rand_index( Top15_D[rand] );
					while (rand_index == i){
						rand = random_range(1, Top15_D.size() - 1);
						rand_index = Top15_D[rand];
					}
					core::pose::Pose base_pose( target_pose ); 
					core::pose::Pose rand_pose( population_pose[rand_index] );
					core::pose::Pose best_pose( population_pose[rank_Dscore.begin()->second] );
					
					core::pose::Pose mutant_pose( base_pose );
					for (Size r = ss_begin; r <= ss_end; ++r){
						mutant_pose.set_phi( r, base_pose.phi(r) + F_loop * (best_pose.phi(r) - rand_pose.phi(r)) );
						mutant_pose.set_psi( r, base_pose.psi(r) + F_loop * (best_pose.psi(r) - rand_pose.psi(r)) );
					}
							
					///@note Crossover
					core::pose::Pose cross_pose( base_pose );
					Size rand_r( random_range(ss_begin, ss_end) );
					for (Size r = ss_begin; r <= ss_end; ++r){
						if ( random_uniform() > 0.5 || r == rand_r ){
							cross_pose.set_phi( r, mutant_pose.phi(r) );
							cross_pose.set_psi( r, mutant_pose.psi(r) );
						}
					}
					///@note Selection
					SCORE trial_score = cal_all_score( cross_pose );
					if ( boltzmann_accept(target_score.total_score(), trial_score.total_score()) ){
						target_pose = cross_pose;
						target_score = trial_score;
						local_mover_success = true;
					}
				}
			}
			if ( local_mover_success ){
				population_pose[i] = target_pose;
				population_score[i] = target_score;
				
			  				
				if (population_score[i].total_score() < best_SS_mover_pose.second){
					best_SS_mover_pose.first = target_pose;
					best_SS_mover_pose.second = population_score[i].total_score();
				}
				for ( multimap<Real, Size>::iterator iter = rank_Dscore.begin(); iter != rank_Dscore.end(); ++iter){
					if (iter->second == i){
						rank_Dscore.erase( iter );
						break;
					}
				}
				rank_Dscore.insert( make_pair(population_score[i].total_score(), i) );
			}
		}
		
		///================================= Topology adjustment =================================
		if ( g % 10 == 0 || g == G_ ){
			cout << "   Topology adjustment..." << endl;
			Loop_Specific_Geometry_Optimization();
		}
	
		///==================== Residue-level distance deviation optimization ====================
		if ( g % 5 == 0 || g == G_){
			cout << "   Residue-level distance deviation optimization..." << endl;
			if (g == G_){
				Refine_G_ = 2 * Refine_G_;
			}
			FM_Optimization();
		}
		cout << endl;
	}


	///@note output Five lowest-potential models
	multimap<Real, core::pose::Pose> score_poses;
	for (Size i = 0; i < NP_; ++i){
		score_poses.insert( make_pair(population_score[i].total_score(), population_pose[i]) );
	}
	score_poses.insert( make_pair(best_FA_mover_pose.second, best_FA_mover_pose.first) );
	score_poses.insert( make_pair(best_SS_mover_pose.second, best_SS_mover_pose.first) );
	score_poses.insert( make_pair(best_LS_mover_pose.second, best_LS_mover_pose.first) );
	score_poses.insert( make_pair(best_FM_mover_pose.second, best_FM_mover_pose.first) );
	
	Size n = 1;
	for ( multimap<Real, core::pose::Pose>::iterator iter = score_poses.begin(); n <= 5; ++iter, ++n){
		stringstream IO;
		IO << "./model_" << n << ".pdb";
		string pdb_route;
		IO >> pdb_route;
		iter->second.dump_pdb( pdb_route.c_str() );
	}
}


Real DM_score( core::pose::Pose &pose, core::pose::Pose &tempPose ){
	Real score_sum = 0;
	Real dm_score = 0;
	Size count = 0;
	for(Size i = 0; i < proteinLength; i++){
		for(Size j = i+1; j < proteinLength; j++){
			
			numeric::xyzVector<Real> a1;
			numeric::xyzVector<Real> a2;
			numeric::xyzVector<Real> b1;
			numeric::xyzVector<Real> b2;
			
			if ( pose.residue(i+1).name3() == "GLY" || pose.residue(j+1).name3() == "GLY" ){
				a1 = pose.residue(i+1).xyz("CA");
				a2 = pose.residue(j+1).xyz("CA");
				
				b1 = tempPose.residue(i+1).xyz("CA");
				b2 = tempPose.residue(j+1).xyz("CA");
			}else{
				a1 = pose.residue(i+1).xyz("CB");
				a2 = pose.residue(j+1).xyz("CA");
				
				b1 = tempPose.residue(i+1).xyz("CA");
				b2 = tempPose.residue(j+1).xyz("CA");
			}
			
			Real d0 = log(fabs(i-j));
			Real di = fabs( a1.distance(a2) - b1.distance(b2) );
			
			if( fabs(i-j) >= 3 ){
				score_sum += 1/(1 + pow(di/d0, 2));
				++count;
			}
		}
	}
	dm_score = score_sum/count;

	return dm_score;
}

vector< std::vector< Real > > Similarity_Matrix(vector< core::pose::Pose > &population){
	vector<vector<Real> > distance_Matrix( population.size(), vector<Real>(population.size(), 0) );
	
	for ( Size i = 0; i < population.size(); ++i ){
		for ( Size j = i+1; j < population.size(); ++j ){
			 distance_Matrix[i][j] = 1.0 - DM_score( population[i], population[j] );
			 distance_Matrix[j][i] = distance_Matrix[i][j];
		//	 if ( distance_Matrix[i][j] == 0 )
		//		cout <<  "The RMSD of " << i << " and " << j << "equal to 0. (" << i <<"," << j << ")" << endl;
		}
	}
	return distance_Matrix;
}

vector<pair<Size, vector<Size> > > Kmediods(Size K, std::vector<std::vector<Real> > &Distance_Matrix){
	///@note the number of clustring points.
	Size N( Distance_Matrix.size() );
	///@note store K clustring centers.
	vector<Size> mediods( K, 0 );
	///@brief randomly generates the cluster centers.
	for ( Size i = 0; i < mediods.size(); ++i ){
		Size rand_mediod( 0 );
		while ( 1 ){
			rand_mediod = numeric::random::rg().random_range(0, N - 1);
			Size j( 0 );
			for (; j < mediods.size(); ++j){
				if ( rand_mediod == mediods[j] )
					break;
			}
			if ( j == mediods.size() ){
				mediods[i] = rand_mediod;
				break;
			}
		}
	}
//	cout << "Initial mediods: " << mediods << endl;
	
	///@note  store the cluster to which the point belong for each clustring points.
	vector< Size > clusterAssement( N, 0 );
	///@brief cluster n points according to the cluster centers.
	for ( Size i = 0; i < N; ++i ){
		///@note store closest cluster and the distance.
		pair<Size, Real> cluster_distance( make_pair(0, 0) );
		for ( Size k = 0; k < mediods.size(); ++k ){
			if ( Distance_Matrix[i][mediods[k]] < cluster_distance.second || k == 0 ){
				cluster_distance.first  = k;
				cluster_distance.second = Distance_Matrix[i][mediods[k]];
			}
		}
		clusterAssement[i] = cluster_distance.first;
	}
//	cout << "Initial cluster: " << clusterAssement << endl << endl;
	
	///@brief avoiding have cluster with no points.
	while ( 1 ){
		Size k( 0 );
		for ( ; k < mediods.size(); ++k ){
			///@note extract all points belongs to cluster k.
			Size count_points( 0 );
			vector<Size> points_in_cluster;
			for ( Size i = 0; i < N; ++i ){
				if ( clusterAssement[i] == k ){
					points_in_cluster.push_back( i );
					++count_points;
				}
			}
			if ( count_points <= 2 ){
				Size rand_mediod( 0 );
				while ( 1 ){
					rand_mediod = numeric::random::rg().random_range(0, N - 1);
					Size j( 0 );
					for (; j < mediods.size(); ++j){
						if ( rand_mediod == mediods[j] )
							break;
					}
					if ( j == mediods.size() ){
						mediods[k] = rand_mediod;
						break;
					}
				}
				break;
			}
		}
		if ( k == mediods.size() )
			  break;
		else{
			///@brief cluster n points according to the cluster centers.
			for ( Size i = 0; i < N; ++i ){
				///@note store closest cluster and the distance.
				pair<Size, Real> cluster_distance( make_pair(0, 0) );
				for ( Size k = 0; k < mediods.size(); ++k ){
					if ( Distance_Matrix[i][mediods[k]] < cluster_distance.second || k == 0 ){
						cluster_distance.first  = k;
						cluster_distance.second = Distance_Matrix[i][mediods[k]];
					}
				}
				clusterAssement[i] = cluster_distance.first;
			}
		}
	}
//	cout << "Initial mediods: " << mediods << endl;
//	cout << "Initial cluster: " << clusterAssement << endl << endl;
	
	
	Size num( 0 );
	///@brief iteration updata
	while ( num < 100)
	{
		++num;
	//	cout << "========================= clustring cycles: " << num << "   cluster centers: ";
	//	for (Size i = 0; i < mediods.size(); ++i)
	//	    cout << mediods[i] << " ";
	//	cout << endl;
	  
		///@brief update cluster centers
		pair<Size, Size> max_cluster( make_pair(0, 0) );
		for ( Size k = 0; k < mediods.size(); ++k ){
			///@note extract all points belongs to cluster k.
			Size count_points( 0 );
			vector<Size> points_in_cluster;
			for ( Size i = 0; i < N; ++i ){
				if ( clusterAssement[i] == k ){		
					points_in_cluster.push_back( i );
					++count_points;
				}
			}
			///@note store cluster center and the total distance to all points.
			pair<Size, Real> center_distance(0, 0);
			///@note update
			for ( Size m = 0; m < points_in_cluster.size(); ++m ){
				Real total_distance( 0 );
				for ( Size n = 0; n < points_in_cluster.size(); ++n )
					total_distance += Distance_Matrix[points_in_cluster[m]][points_in_cluster[n]];
				if ( m == 0 || total_distance < center_distance.second ){
					center_distance.first = points_in_cluster[m];
					center_distance.second = total_distance;
				}
			}
			mediods[k] = center_distance.first;
		//	cout << mediods[k] << "(" << count_points << ")   ";
			
			if ( k == 0 || count_points > max_cluster.second ){
				max_cluster.first = k;
				max_cluster.second = count_points;
			}
				
		}
	//	cout << endl;
		
		Size count_of_state_transition( 0 );
		///@brief update cluster (reclassification).
		for ( Size i = 0; i < N; ++i ){
			///@note store closest cluster and the distance.
			pair<Size, Real> cluster_distance( make_pair(0, 0) );
			for ( Size k = 0; k < mediods.size(); ++k ){
				if ( k == 0 || Distance_Matrix[i][mediods[k]] < cluster_distance.second ){
					cluster_distance.first = k;
					cluster_distance.second = Distance_Matrix[i][mediods[k]];
				}
			}
			if ( cluster_distance.first != clusterAssement[i] ){
				++count_of_state_transition;
				clusterAssement[i] = cluster_distance.first;
			}
		}
		
		///@brief avoiding have cluster with no points.
		bool redistribure_medoid( false );
		while ( 1 ){
			Size k( 0 );
			for ( ; k < mediods.size(); ++k ){
				///@note extract all points belongs to cluster k.
				Size count_points( 0 );
				vector<Size> points_in_cluster;
				for ( Size i = 0; i < N; ++i ){
					if ( clusterAssement[i] == k ){
						points_in_cluster.push_back( i );
						++count_points;
					}
				}
				if ( count_points <= 2 ){
					redistribure_medoid = true;
				//	cout << k << "-th cluster " << mediods[k] << " not have element, reselect a new cluster center!!!" << endl;
					Size rand_mediod( 0 );
					while ( 1 ){
						rand_mediod = numeric::random::rg().random_range(0, N - 1);
						Size j( 0 );
						for (; j < mediods.size(); ++j){
							if ( rand_mediod == mediods[j] )
								break;
						}
						if ( j == mediods.size() ){
							mediods[k] = rand_mediod;
							break;
						}
					}
					break;
				}
			}
			if ( k == mediods.size() )
				  break;
			else{
				///@brief cluster n points according to the cluster centers.
				for ( Size i = 0; i < N; ++i ){
					///@note store closest cluster and the distance.
					pair<Size, Real> cluster_distance( make_pair(0, 0) );
					for ( Size k = 0; k < mediods.size(); ++k ){
						if ( Distance_Matrix[i][mediods[k]] < cluster_distance.second || k == 0 ){
							cluster_distance.first  = k;
							cluster_distance.second = Distance_Matrix[i][mediods[k]];
						}
					}
					clusterAssement[i] = cluster_distance.first;
				}
			}
		}
	//	cout << " mediods: " << mediods << endl;
	//	cout << " cluster: " << clusterAssement << endl << endl;
		
		
		
		///@note if cluster is not change, break.
		if ( (count_of_state_transition == 0 && !redistribure_medoid) || num == 150 ){
			Size center_of_max_cluster = mediods[max_cluster.first];
			
		//	cout << "cluster centers : ";
		//	for ( Size k = 0; k < mediods.size(); ++k ){
		//		cout << mediods[k] << " ";
		//	}
		//	cout << " center of max cluster : " << center_of_max_cluster << endl;
			break;
		}
	}
	
	vector<pair<Size, vector<Size> > > vec_mediod_elements;
	for (Size i = 0; i < mediods.size(); ++i){
		Size mediod = mediods[i];
		vector<Size> elements;
		for(Size j = 0; j < clusterAssement.size(); ++j){
			if (clusterAssement[j] == i)
				elements.push_back( j );
		}
		vec_mediod_elements.push_back( make_pair(mediod, elements) );
	}
	
	return vec_mediod_elements;
}

}
}