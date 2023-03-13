#include "utils.h"

cl_args::cl_args(int argc, char* argv[]) : switches(), options(), trailing() {
	bool last_value_name = false;
	bool l = false;
	char name;
	for (int i = 0; i < argc; i++) {
		if (l) {
			trailing.push_back(string(argv[i]));	
		} else if (argv[i][0] == ':') {
			l = true;
		} else if (argv[i][0] == '-') {
			if (last_value_name) {
				if (switches.count(name) > 0) throw std::invalid_argument("Duplicate argument given!");
				switches.insert(name);
			}
			last_value_name = true;
			name = argv[i][1];
		} else {
			if (!last_value_name) throw std::invalid_argument("Stray argument value!");
			last_value_name = false;
			if (options.count(name) > 0) throw std::invalid_argument("Duplicate argument given!");
			options[name] = string(argv[i]);
		}
	}
	if (last_value_name) { 
		if (switches.count(name) > 0) throw std::invalid_argument("Duplicate argument given!");
		switches.insert(name);
	}
}

uint32_t ac_thres(vector<double>& lfacs, uint32_t an, double cond_af, double prob_thres) {
	double prob = 0.0;
	if (an >= lfacs.size()) throw std::invalid_argument("AN out of range");
	uint32_t ac = 0;
	while (ac <= an) {
		prob += exp(lfacs[an]-lfacs[ac]-lfacs[an-ac] + ac * log(cond_af) + (an-ac)*log(1-cond_af));
		if (prob > (1-prob_thres)) {
			if (ac > 0) ac-=1;
			break;
		}
		ac++;
	}
	return ac; 
}

vector<uint32_t> ac_thres_vec(vector<double>& lfacs, double cond_af, double prob_thres) {
	auto upto = lfacs.size();
	vector<uint32_t> r(upto);
	for (uint32_t i = 1; i < upto; i++) {
		r[i] = ac_thres(lfacs, i, cond_af, prob_thres);
	}
	return r; 
}

rsvr_pop::rsvr_pop(int64_t in_id) : id(in_id), pop_ac_ans() {}
void rsvr_pop::judge(vector<uint32_t>& pmaf1k, vector<uint32_t>& pmaf10k) {
	bool onek = true;
	bool tenk = true;
	uint32_t tot = 0;
	uint32_t hon = 0;
	for (const auto &keyval : pop_ac_ans) {
		if (keyval.first.compare("raw") == 0) { continue; } 
		auto ac = keyval.second.first;
		auto an = keyval.second.second;
		if (ac > pmaf1k[an] && ac > 1) { onek=false; }
		if (ac > pmaf10k[an] && ac > 1) { tenk=false; }
		tot += ac;
	}
	if (tot == 0) hon = 3;
	else if (tenk) hon = 2;
	else if (onek) hon = 1;
	else hon = 0;

	cout << id << '\t' << hon << endl;
}

