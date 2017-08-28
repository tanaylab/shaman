#include "ContactShuffler.h"
#include "Random.h"
#include <ctime>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
int shaman_hic_matrix_shuffler_cpp(Rcpp::DataFrame& raw_contacts,
		std::string shuf_contacts,
		int shuffle_factor,
		int proposal_from_contacts,
		double proposal_iterations_d,
		int dist_resolution,
		int decay_smooth,
		int decay_regularization,
		double proposal_correction_factor_d,
		int max_dist,
		int min_dist,
		int grid_switch_bin_dist,
		int grid_x_min_bin,
		int grid_x_max_bin,
		int grid_x_increase,
		int grid_x_increase_iter,
		int input_symmetric_mat,
		int output_symmetric_mat)
{
	Rcpp::Rcerr << "running new shuffler" << endl;

	Random::reset(-1);
	clock_t begin = clock();
	//Rcpp::Rcout << "raw contacts file = " << raw_contacts << endl;


	int dist_log_scale = 2;
	int output_header = 1;
	float proposal_iterations = (float) proposal_iterations_d;
	int proposal_from_constant = 0;
	float proposal_correction_factor = (float)proposal_correction_factor_d;

	Rcpp::Rcerr << "proposal iterations = " << proposal_iterations << "\t"
			<< " proposal_correction_factor = " << proposal_correction_factor << endl;

	Rcerr << "shuffle factor = " << shuffle_factor << endl;
	ContactShuffler shuffler(dist_log_scale,dist_resolution,
			grid_x_min_bin, grid_switch_bin_dist,
			proposal_correction_factor,  decay_smooth,
			decay_regularization, min_dist, max_dist);

	vector< vector <int > > r_contacts = Rcpp::as<vector<vector<int> > >(raw_contacts);
	int contacts = shuffler.load_contacts(r_contacts, input_symmetric_mat);
	Rcerr << "finished loading from contacts" << endl;

	shuffler.init_exp_decay_from_obs();
	cerr << "finished init exp decay from observed" << endl;

	if (proposal_from_constant) {
		shuffler.init_proposal_const();
		cerr << "finished init proposal to constant" << endl;
	} else {
	 if (proposal_from_contacts) {
		shuffler.init_proposal_from_contacts(floor(proposal_iterations * contacts));
		cerr << "finished init proposal from contacts" << endl;
	 } else {
		shuffler.init_proposal_from_area();
		cerr << "finished init proposal from area" << endl;
	 }
	}
	int iter = 0;
	while(iter < shuffle_factor) {
		shuffler.shuffle_contacts(grid_x_increase_iter, 0.0001, 0.01, 0);
		iter += grid_x_increase_iter;
		if (grid_x_min_bin != grid_x_max_bin) {
			grid_x_min_bin += grid_x_increase;
			shuffler.reset_grid(grid_x_min_bin);
		}
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cerr << "run time = " << elapsed_secs/60/60 << " hours" << endl;
	}
	shuffler.save_contacts(shuf_contacts.c_str(), output_symmetric_mat, output_header);
	return(0);
}

