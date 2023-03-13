#ifndef PMAF_H
#define PMAF_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <getopt.h>

using namespace std;

struct cl_args {
	set<char> switches;
	map<char, string> options;
	vector<string> trailing;
	cl_args(int argc, char* argv[]);
};

uint32_t ac_thres(vector<double>& lfacs, uint32_t an, double cond_af, double prob_thres);
vector<uint32_t> ac_thres_vec(vector<double>& lfacs, double cond_af, double prob_thres);

struct rsvr_pop {
	int64_t id;
	map<string, pair<uint32_t, uint32_t> > pop_ac_ans;
	rsvr_pop(int64_t in_id);
	void judge(vector<uint32_t>& pmaf1k, vector<uint32_t>& pmaf10k);
};

#endif
