#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <string>
#include <map> 
#include <set> 
#include <memory> 
#include <array> 
#include <vector> 
#include <algorithm>
#include <stdexcept>
#include "utils.h"
#include "enc.h"
#include "idsets.h"
#include "idtypes.h"
#include "seqfx.h"

#define VERSION "1.0"

using namespace std;

struct gtmap {
	int64_t curv;
	map<uint32_t, uint32_t> m;
	void dump() {
		for (const auto &keyval : m) {
			cout << curv << ',' << keyval.first << ',' << keyval.second << endl;
		}
		m.clear();
	}
	void push(int64_t v, uint32_t s, uint32_t gt) {
		if (v != curv) {
			dump();
			curv = v;
		}
		if (m.count(s) == 0) m[s] = gt;
		else m[s] = max(m[s], gt);
	}
	gtmap() : curv(0), m({}) {}
};

bool ambiguous_rsvr_id(int64_t v) {
	int64_t refl = (v >> 24) % 64;
	int64_t altl = (v >> 18) % 64;
	return (refl >= 63) || (altl >= 10);
}

int enc(cl_args& args, istream& inp) {
	bool sorted = false;
	bool normalised = true;
	bool prepend = false;
	bool allow_dups = true;
	bool print_once = false;
	bool record = false;
	string record_file_name;
	char delim = '\t';
	int j = 0;
	int check = 0;
	int tolerance = 2;

	array<char, 4> fields({'c','p','r','a'});
	array<int, 4> indices({1,2,3,4});

	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Assign RSVR IDs to input CHROM, POS, REF, ALT tuples.\n";
				cout << "-d\tDelimiter.\n";
				cout << "-f\tOrder in which CHROM (c), POS (p), REF (r) and ALT(a) appear in lines of input (defaults to 'cpra').\n";
				cout << "-i\tComma-separated, one-based indices of fields in each line corresponding to properties given with -f option (defaults to 1,2,3,4).\n";
				cout << "-h\tPrint this message (again).\n";
				cout << "-k\tIf '-s' is set then terminate the program if a duplicated RSVR ID is encountered\n";
				cout << "-n\tDon't normalise.\n";
				cout << "-p\tPrepend the RSVR ID to each line.\n";
				cout << "-r\tFile in which to record rows with ambiguous RSVR IDs.\n";
				cout << "-s\tSort output by RSVR ID.\n";
				cout << "-x\tIf '-s' is set then exclude duplicated RSVR IDs (i.e. keeping the first one)\n";
				return 0;
			case 'n':
				normalised = false;
				break;
			case 's':
				sorted = true;
				break;
			case 'x':
				print_once = true;
				break;
			case 'k':
				allow_dups = false;
				break;
			case 'p':
				prepend = true;
				break;
		}
	}

	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'd':
				delim = keyval.second[0];
				break;
			case 't':
				indstr = istringstream(keyval.second);
				indstr >> tolerance;
				break;
			case 'f':
				indstr = istringstream(keyval.second);
				while (indstr >> fields[j++] && j < 4);
				for (int i = 0; i < 4; i++) 
					check += 1 << ((int)fields[i] - 97);
				if (check != 163845)
					throw std::invalid_argument("When specifying fields must have c (CHROM), p (POS), r (REF) and a (ALT)!");
				break;
			case 'r':
				indstr = istringstream(keyval.second);
				record = true;
				record_file_name = indstr.str();
				break;
			case 'i': 
				indstr = istringstream(keyval.second);
				for (std::string num; std::getline(indstr, num, ',') && j < 4; j++) {
					indices[j] = stoi(num);
				}
				if (j != 4)
					throw std::invalid_argument("Must specify four indices for c (CHROM), p (POS), r (REF) and a (ALT)!");
				for (int i = 0; i < 3; i++) 
					if (indices[i] >= indices[i+1])
						throw std::invalid_argument("Given indices must be strictly increasing!");
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	ofstream record_file(record_file_name);

	if (!sorted && (!allow_dups || print_once)) throw std::invalid_argument("Only specify '-x' or '-k' when sorting!");

	sorting_printer prms(fields, indices, normalised, delim, prepend);
	read_sorted<sp_entry, sorting_printer> sp(prms, allow_dups, print_once, 30, tolerance);
	if (sorted) { sp.print_sorted(inp); }
	else { sp.print(inp); }
	return 0;
}

int gtnorm(cl_args& args, istream& inp) {
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Normalise input genotype triples formatted as <RSVR ID><tab><sample ID><tab><GT> by selecting max GT per RSVR ID per sample.\n";
				cout << "-h\tPrint this message (again).\n";
				return 0;
			}
	}

	int64_t var;
	uint32_t samp;
	uint32_t gt;
	gtmap x;
	while (inp >> var >> samp >> gt) {
		x.push(var, samp, gt);
	}
	x.dump();
	return 0;
}

int amb(cl_args& args, istream& inp) {
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Filter input lines to retain those prepended with an ambiguous RSVR ID (i.e. having an ALT allele length greater than 9 or a REF allele length greater than 62.\n";
				cout << "-h\tPrint this message (again).\n";
				return 0;
			}
	}

	for (std::string line; std::getline(inp, line);) {
		int64_t var;
		istringstream ss(line);
		ss >> var;
		if (ambiguous_rsvr_id(var)) { cout << line << endl; }
	}
	return 0;
}

void bank(typename id_tabr<rsvr<int64_t>, rsvr_record>::Toutput& o, vector<int64_t>& sample_ids, vector<size_t>& heights, bool sort_sample_ids=false) {
	if (sort_sample_ids) sort(sample_ids.begin(), sample_ids.end());
	auto N = sample_ids.size();
	if (N > 1) for (size_t i = 0; i < N-1; i++) {
		if (sample_ids[i] >= sample_ids[i+1]) throw std::invalid_argument("Input must contain no duplicates/be strictly increasing!");
	}
	if (o.first.size() == 0) {
		o.first = sample_ids;
		o.second = heights;
	} else {
		id_tab<rsvr<int64_t> > t_sample(sample_ids, heights);
		id_tab<rsvr<int64_t> > t_track(o.first, o.second);
		o = iio<id_tabr<rsvr<int64_t>, rsvr_record> >(t_track, t_sample);
	}
	sample_ids.clear();
	heights.clear();
}

void bank1(typename id_tabr<rsvr<int64_t>, rsvr_record>::Toutput& o, vector<int64_t>& sample_ids, bool sort_sample_ids=false) {
	vector<size_t> x(sample_ids.size(), 1);
	bank(o, sample_ids, x, sort_sample_ids);
}

int depth(cl_args& args, istream& inp) {
	typedef id_tabr<rsvr<int64_t>, rsvr_record> T;
	typedef vector<int64_t> vll;

	size_t cut_at = 0;
	bool do_cut = false;
	bool sort_breakpoints = false;
	char delim = '\t';
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Evaluate depth of piled sets of ordered non-overlapping intervals supplied as <start><delimiter><end>\\n, with <end> exclusive and interval sets separated by empty lines.\n";
				cout << "Output: <region start><delimiter><depth> (final line should indicate a depth of 0).\n";
				cout << "Example use case: finding depth of PASSing regions from a set of gVCF files.\n";
				cout << "-c\tCut interval piles at given height and yield output as intervals achieving at least given height <start><delimiter><end>.\n";
				cout << "-d\tDelimiter.\n";
				cout << "-h\tPrint this message (again).\n";
				cout << "-s\tSort inputs.\n";
				return 0;
			case 's':
				sort_breakpoints = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'c':
				do_cut = true;
				indstr = istringstream(keyval.second);
				indstr >> cut_at;
				break;
			case 'd':
				delim = keyval.second[0];
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	typename T::Toutput s(
		vector<int64_t>(0),
		vector<size_t>(0)
	);
	typename T::Toutput e(
		vector<int64_t>(0),
		vector<size_t>(0)
	);
	vll sids(0);
	vll eids(0);

	for (std::string line; std::getline(inp, line);) {
		bool line_read = false;
		int64_t start = 0;
		std::istringstream l(line);
		for (std::string field; std::getline(l, field, delim); ) {
			int64_t i = stoll(field);
			if (line_read)
				if (i > start)
					eids.push_back(i);
				else
					throw std::invalid_argument("Interval end points must be strictly greater than starts!");
			else {
				start = i;
				sids.push_back(i);
			}
			line_read = true;
		}
		if (!line_read) {
			bank1(s, sids, sort_breakpoints);
			bank1(e, eids, sort_breakpoints);
		}
	}
	bank1(s, sids, sort_breakpoints);
	bank1(e, eids, sort_breakpoints);

	id_tab<rsvr<int64_t> > stab(s.first, s.second);
	id_tab<rsvr<int64_t> > etab(e.first, e.second);

	typename T::Toutput o = iio<cum_id_tabr<rsvr<int64_t> , rsvr_record, false> >(stab, etab);

	if (!do_cut) {
		size_t n = o.first.size();
		for (int j = 0; j < n; j++) {
			cout << o.first[j] << delim << o.second[j] << endl;
		}
	} else {
		pair<vll, vll> ints = cut_cum_tab<vll>(o.first, o.second, cut_at);
		size_t n = ints.first.size();
		for (int j = 0; j < n; j++) {
			cout << ints.first[j] << delim << ints.second[j] << endl;
		}
	}
	return 0;
}

int ann(cl_args& args, istream& inp) {
	istringstream indstr;
	bool one = false;
	bool file_given = false;
	bool m_given = false;
	string mng = "";
	string file_name;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Prepend lines prefixed RSVR IDs with values attached to RSVR IDs stored in a tab delimited file.\n";
				cout << "Output: Appends <tab><value> to lines with matching RSVR IDs.\n";
				cout << "-1\tPrint exactly one annotated line per input RSVR ID, thus lines not matched in annotated RSVR IDs file.\n";
				cout << "-f\tFile of annotated RSVR IDs.\n";
				cout << "-h\tPrint this message (again).\n";
				cout << "-m\tString to print when RSVR ID is missing from annotation.\n";
				return 0;
			case '1':
				one = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
		}
	}
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'f':
				indstr = istringstream(keyval.second);
				file_name = indstr.str();
				file_given = true;
				break;
			case 'm':
				indstr = istringstream(keyval.second);
				m_given = true;
				mng = indstr.str();
				break;
			default:
				throw std::invalid_argument("Invalid option!");
		}
	}

	if (m_given && !one) 
		throw std::invalid_argument("Cannot specify 'm' without '1' option set!");
	if (!file_given)
		throw std::invalid_argument("Must specify file!");

	ids_file<id_prep_cin> f1(0, inp);
	ifstream infile(file_name);
	ids_file<id_prep_cin> f2(0, infile);

	unique_ptr<T2> o;
	if (one) o = make_unique<anno_one>(f1, f2, mng);
	else o = make_unique<anno_many>(f1, f2);
	o->process();
	return 0;
}

int mix(cl_args& args, istream& inp) {
	istringstream indstr;
	bool as_intervals = false;
	bool anno = false;
	bool exclude = false;
	bool file_given = false;
	string file_name;
	int move_bits = 0;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Intersect/exclude input RSVR IDs with intervals/RSVR IDs from a file.\n";
				cout << "-a\tInstead of printing/not-printing, annotate 1 or 0 respectively for included/excluded.\n";
				cout << "-b\tBits to shift IDs supplied to the right.\n";
				cout << "-f\tFile of RSVR IDs/intervals.\n";
				cout << "-h\tPrint this message (again).\n";
				cout << "-i\tTreat file as file of intervals (exclusive end points) rather than points (RSVRs).\n";
				cout << "-x\tExclude points/intervals from file.\n";
				return 0;
			case 'i':
				as_intervals = true;
				break;
			case 'x':
				exclude = true;
				break;
			case 'a':
				anno = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
		}
	}
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'f':
				indstr = istringstream(keyval.second);
				file_name = indstr.str();
				file_given = true;
				break;
			case 'b':
				indstr = istringstream(keyval.second);
				indstr >> move_bits;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
		}
	}

	if (!file_given)
		throw std::invalid_argument("Must specify file!");

	ids_file<id_prep_cin> f1(move_bits, inp);
	ifstream infile(file_name);
	ids_file<id_file> f2(0, infile);

	unique_ptr<T> o;
	if (as_intervals) o = make_unique<f_ints>(f1, f2, !exclude, anno);
	else if (anno) { if (exclude) throw invalid_argument("Can't specify -a and -x"); o = make_unique<anno_ids>(f1, f2, anno_type::anno); }
	else if (exclude) o = make_unique<anno_ids>(f1, f2, anno_type::exclude);
	else o = make_unique<anno_ids>(f1, f2, anno_type::intersect);
	o->process();
	return 0;
}

int merge(cl_args& args, istream& inp) {
	istringstream indstr;
	bool file_given = false;
	string file_name;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Merge two sets of optionally annotated, sorted RSVR IDs into a single sorted list.\n";
				cout << "-f\tFile of RSVRs IDs/intervals.\n";
				cout << "-h\tPrint this message (again).\n";
				return 0;
			default:
				throw std::invalid_argument("Invalid option!");
		}
	}
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'f':
				indstr = istringstream(keyval.second);
				file_name = indstr.str();
				file_given = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
		}
	}

	if (!file_given)
		throw std::invalid_argument("Must specify file!");

	ids_file<id_prep_cin> f1(0, inp);
	ifstream infile(file_name);
	ids_file<id_prep_cin> f2(0, infile);

	merge_lines m(f1, f2);
	m.process();
	return 0;
}


int get_chrom_num(string& chrom) {
	int chrom_num = -1;
	if (chrom.compare("MT")==0) chrom_num = 25;
	else if (chrom.compare("Y")==0) chrom_num = 24;
	else if (chrom.compare("X")==0) chrom_num = 23;
	else chrom_num = std::stoi(chrom);
	if (chrom_num > 25 || chrom_num < 1)
		throw std::invalid_argument("Invalid chromosome: must belong to {1-22, X, Y, MT}!");
	return chrom_num;
}

int dec(cl_args& args, istream& inp) {
	istringstream indstr;
	string ref;
	bool unique = false;
	bool ref_spec = false;
	bool chr = false;
	bool vcf_style = true;
	bool prep_rsvr = false;
	bool ref_len = false;
	string mt = "MT";
	char delim = '\t';
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Decode input RSVR IDs to CHROM, POS, REF, ALT tuples against a given reference fasta file\n";
				cout << "-c\tAdd 'chr's to chromosome identifiers.\n";
				cout << "-d\tDelimiter to separate fields of output.\n";
				cout << "-f\tLocation of reference fasta file.\n";
				cout << "-h\tPrint this message again.\n";
				cout << "-l\tDecode to reference allele length instead of reference allele (reference genome not required).\n";
				cout << "-m\tString to use for mitochondrial genome (defaults to M).\n";
				cout << "-n\tPrint normalised alleles rather than VCF style (default) with enforced non-empty REF/ALT strings.\n";
				cout << "-p\tPrepend the RSVR ID.\n";
				cout << "-u\tFail when an ambiguous RSVR ID is encountered.\n";
				return 0;
			case 'n':
				vcf_style = false;
				break;
			case 'p':
				prep_rsvr = true;
				break;
			case 'l':
				ref_len = true;
				break;
			case 'c':
				chr = true;
				break;
			case 'u':
				unique = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'd':
				delim = keyval.second[0];
				break;
			case 'm':
				indstr = istringstream(keyval.second);
				mt = indstr.str();
				break;
			case 'f':
				indstr = istringstream(keyval.second);
				ref = indstr.str();
				ref_spec = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	unique_ptr<id_dec> d;
	std::ifstream refstr;
	
	if (ref_spec && ref_len) throw invalid_argument("Cannot set -l and -f options!");
	if (ref_spec) {
		refstr.open(ref);
		d = make_unique<dec_ref>(
			delim,
			chr,
			mt,
			vcf_style,
			refstr
		);
	} else {
		if (ref_len) 
			d = make_unique<dec_noref>(
				delim,
				chr,
				mt
			);
		else throw std::invalid_argument("Reference genome not specified with -f!");
	}

	int64_t v;
	string line;
	while (getline(inp, line)) {
		std::istringstream l(line);
		l >> v;

		if (unique && (ambiguous_rsvr_id(v)))
			throw invalid_argument("Ambiguous RSVR ID (" + to_string(v) + ") detected!");

		if (prep_rsvr)
			cout << line << delim;

		d->operator()(v);
	}

	return 0;
}

int seqfx(cl_args& args, istream& inp) {
	istringstream indstr;
	string gtf;
	string ref;
	AGGREGATE agg = AGGREGATE::GENE;
	bool unique = false;
	bool vep = false;
	bool one_per_var = false;
	bool append = false;
	char delim = '\t';
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Annotate input RSVR IDs with transcript consequences stored in a tab delimited table\n";
				cout << "Usage: rsvr seqfx -f <reference fasta file> -g <>\n";
				cout << "-1\tExpect one annotation per variant.\n";
				cout << "-a\tAppend aggregation to lines of input (collapsing with comma-separated lists of pipe-delimited fields).\n";
				cout << "-d\tDelimiter to use for output (defaults to tab).\n";
				cout << "-f\tLocation of reference fasta file.\n";
				cout << "-g\tLocation of gtf file.\n";
				cout << "-h\tPrint this message again.\n";
				cout << "-t\tAggregate consequences by transcript (rather than gene by default).\n";
				cout << "-u\tFail when non-unique RSVR ID is encountered.\n";
				cout << "-v\tPrint consequences used in VEP.\n";
				return 0;
			case 'u':
				unique = true;
				break;
			case 'v':
				vep = true;
				break;
			case 't':
				agg = AGGREGATE::TRANSCRIPT;
				break;
			case 'a':
				append = true;
				break;
			case '1':
				one_per_var = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'd':
				delim = keyval.second[0];
				break;
			case 'f':
				indstr = istringstream(keyval.second);
				ref = indstr.str();
				break;
			case 'g':
				indstr = istringstream(keyval.second);
				gtf = indstr.str();
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	std::ifstream gtf_file(gtf);
	feature_stack<transcript> features(gtf_file);

	std::ifstream refstr(ref);
	ref_fasta reff(refstr);

	tater2 t1(features, 63, reff, agg);

	int64_t v;
	int64_t lastv = 0;

	string line;

	unique_ptr<seqfx_print> o;
	if (append)
		o = make_unique<seqfx_line_print>(vep, '|',',',delim);
	else 
		o = make_unique<seqfx_db_print>(vep, delim);

	while (getline(inp, line)) {
		std::istringstream l(line);
		l >> v;

		if (v < lastv) {
			throw invalid_argument("Input variants not sorted by RSVR ID!");
		}

		if (unique && ambiguous_rsvr_id(v))
			throw invalid_argument("Ambiguous variant RSVR ID (" + to_string(v) + ") detected!");

		VV av(v);
		auto ann = t1.anno_gene(av);
		if (one_per_var && (ann.size() > 1))
			throw invalid_argument("Multiple annotations found for variant with RSVR ID " + to_string(v));
		o->print(line, v, ann);

		lastv = v;
	}

	return 0;
}

int norm(cl_args& args, istream& inp) {
	istringstream indstr;
	string ref;
	bool must_norm = false;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Normalise RSVR IDs with respect to a reference genome\n";
				cout << "Usage: rsvr norm -f <reference fasta file>\n";
				cout << "-f\tLocation of reference fasta file.\n";
				cout << "-h\tPrint this message again.\n";
				cout << "-x\tFail if unnormalised RSVR ID encountered.\n";
				return 0;
			case 'x':
				must_norm = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'f':
				indstr = istringstream(keyval.second);
				ref = indstr.str();
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	std::ifstream refstr(ref);
	ref_fasta reff(refstr);

	int64_t v;
	int64_t lastv = 0;
	for (std::string line; std::getline(inp, line);) {
		istringstream l(line);

		l >> v;

		if (v < lastv) {
			throw invalid_argument("Input variants not sorted by RSVR ID!");
		}
		
		int64_t nv = make_norm_rsvr(v, reff, true);
		if (must_norm && (v != nv))
			throw invalid_argument("Unnormalised RSVR ID encountered!");

		string rest = "";
		getline(l, rest);

		cout << nv << rest << endl;
	}
	return 0;
}

int gref(cl_args& args, istream& inp) {
	istringstream indstr;
	string ref;
	bool use_ids = false;
	uint32_t store_bases = 1000;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Get sequences from reference genome within input intervals (specified with exclusive end points)\n";
				cout << "Usage: rsvr gref -f <reference fasta>\n";
				cout << "-f\tLocation of reference fasta file.\n";
				cout << "-h\tPrint this message again.\n";
				cout << "-i\tIntervals are prefixed with IDs to be printed too.\n";
				cout << "-n\tNumber of bases to store, defaults to " << store_bases << ".\n";
				return 0;
			case 'i':
				use_ids = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'f':
				indstr = istringstream(keyval.second);
				ref = indstr.str();
				break;
			case 'n':
				indstr = istringstream(keyval.second);
				indstr >> store_bases;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	std::ifstream refstr(ref);
	ref_fasta reff(refstr, store_bases);

	int64_t from;
	int64_t to;
	int64_t id;

	if (use_ids) {
		while (inp >> id >> from >> to) {
			int64_t chrnum = from >> 28;
			cout << id << '\t';
			reff.print_ref(chrnum, (from % (1 << 28))-1, to-from);

			cout << endl;
		}
	} else {
		while (inp >> from >> to) {
			int64_t chrnum = from >> 28;

			reff.print_ref(chrnum, (from % (1 << 28))-1, to-from);

			cout << endl;
		}
	}
	return 0;
}

int sort_at_cp(cl_args& args, istream& inp) {
	int tolerance = 2;
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'h':
				cout << "Sort input lines with leading RSVR IDs by RSVR ID, assuming already sorted to chromosome and position.\n";
				cout << "-h\tPrint this message again.\n";
				return 0;
			case 't':
				indstr = istringstream(keyval.second);
				indstr >> tolerance;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}
	anno_rsvr prms;
	read_sorted<sp_entry, anno_rsvr> rs1(prms, true, false, 30, tolerance);
	rs1.print_sorted(inp);
	return 0;
}

int filtgt(cl_args& args, istream& inp) {
	typedef tva<rsvr<gt_rec>, ids_file<id_file>, vector<gt_rec> > S;
	int had = 1;
	istringstream indstr;

	int tolerance = 2;
	string excludes;
	string intervals;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Filter input genotype information against sets of RSVR IDs or intervals stored in files.\n";
				cout << "Input genotype information should be provided as lines formatted as: <chromosome><tab><position><tab><ref allele><tab><alt allele><tab><genotype><tab><pass>\n";
				cout << "Usage: rsvr filtgt [-i <y/n intervals>] [-e <y/n exclude>] : file1 file2 ...\n";
				cout << "-e\tString of 'y's and 'n's determining whether to exclude points/intervals from respective files.\n";
				cout << "-i\tString of 'y's and 'n's determining whether to treat input files respectively as intervals/points.\n";
				cout << "-h\tPrint this message again.\n";
				return 0;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'e':
				indstr = istringstream(keyval.second);
				excludes = indstr.str();
				had += 2;
				break;
			case 't':
				indstr = istringstream(keyval.second);
				indstr >> tolerance;
				break;
			case 'i':
				indstr = istringstream(keyval.second);
				intervals = indstr.str();
				had += 2;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	if ((args.trailing.size() != excludes.length()) || (excludes.length() != intervals.length()))
		throw invalid_argument("The lengths of the strings specified with -i and -e must match the number of given files!");

	gt_vec gtv;
	read_sorted<gt_rec, gt_vec> gts(gtv, true, false, 30, tolerance);

	gts.print_sorted(inp);

	vector<gt_rec> v(gts.ac.gts);
	for (int i = 0; i < args.trailing.size(); i++) {
		bool ints = intervals[i] == 'y';
		bool excl = excludes[i] == 'y';
		rsvr<gt_rec> u(v, ints ? 30 : 0);
		string filename(args.trailing[i]);
		unique_ptr<S> o;
		ifstream infile(filename);
		ids_file<id_file> f2(0, infile);
		if (ints)
			o = make_unique<ss_vec_file_ints<gt_rec> >(u, f2, excl);
		else 
			o = make_unique<ss_vec_file_pts<gt_rec> >(u, f2, excl);
		o->process();
		v = o->yield();
	}

	for (int i = 0; i < v.size(); i++) {
		cout << v[i].rsvr << '\t' << (int)v[i].gt << '\t' << (int)v[i].pass << endl;
	}

	return 0;
}

int tabulate(cl_args& args, istream& inp) {
	typedef id_tabr<rsvr<int64_t>, rsvr_record> T;
	typedef vector<int64_t> vll;

	size_t cut_at = 0;
	bool do_cut = false;
	bool sort_breakpoints = false;
	char delim = '\t';
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Determine count of each input RSVR ID from a set of ordered sets of distinct RSVR IDs formatted as <RSVR ID><delimiter><frequency>\\n and separated by empty lines.\n";
				cout << "-c\tSpecify frequency and output list of RSVR IDs at least given frequency.\n";
				cout << "-d\tDelimiter.\n";
				cout << "-h\tPrint this message (again).\n";
				cout << "-s\tSort inputs within each set of ordered RSVR IDs.\n";
				return 0;
			case 's':
				sort_breakpoints = true;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'c':
				indstr = istringstream(keyval.second);
				indstr >> cut_at;
				do_cut = true;
				break;
			case 'd':
				delim = keyval.second[0];
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	typename T::Toutput s(
		vector<int64_t>(0),
		vector<size_t>(0)
	);
	vll sids(0);
	vector<size_t> freqs(0);

	for (std::string line; std::getline(inp, line);) {
		bool line_read = false;
		size_t freq_i = 1;
		std::istringstream l(line);
		for (std::string field; std::getline(l, field, delim); ) {
			int64_t i = stoll(field);
			if (line_read)
				freq_i = (size_t)i;
			else
				sids.push_back(i);
			line_read = true;
		}
		if (!line_read) {
			bank(s, sids, freqs, sort_breakpoints);
		} else {
			freqs.push_back(freq_i);
		}
	}

	for (int j = 0; j < s.first.size(); j++) {
		if (do_cut) {
			if (s.second[j] >= cut_at) 
				cout << s.first[j] << endl;
		} else {
			cout << s.first[j] << delim << s.second[j] << endl;
		}
	}
	return 0;
}

int pmaf(cl_args& args, istream& inp) {
	uint32_t maxAN = 300000;
	double prob_thres = 0.05;
	for (auto it = args.switches.begin(); it != args.switches.end(); ++it) {
		switch((*it)) {
			case 'h':
				cout << "Annotate RSVR IDs with PMAF scores.\n";
				cout << "-h\tPrint this message (again).\n";
				cout << "-N\tMax AC (default " << maxAN << ").\n";
				return 0;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}
	for (const auto &keyval : args.options) {
		istringstream indstr;
		switch(keyval.first) {
			case 'N':
				indstr = istringstream(keyval.second);
				indstr >> maxAN;
				break;
			default:
				throw std::invalid_argument("Invalid option!");
				break;
		}
	}

	vector<double> lfacs(maxAN);
	lfacs[0] = 0.0;
	for (uint32_t i = 1; i < maxAN; i++) {
		lfacs[i] += lfacs[i-1] + log((double)i);
	}
	auto pmaf1k = ac_thres_vec(lfacs, 0.001, prob_thres);
	auto pmaf10k = ac_thres_vec(lfacs, 0.0001, prob_thres);

	string field;
	rsvr_pop rp(0);
	bool first = true;
	while (std::getline(inp, field)) {
		std::istringstream l(field);
		string s;
		std::getline(l, s, ' ');
		int64_t id = stoll(s);
		if (first) {
			rp = rsvr_pop(id);
			first = false;
		} else if (id != rp.id) {
			rp.judge(pmaf1k, pmaf10k); 
			rp = rsvr_pop(id); 
		}
		pair<uint32_t, uint32_t> acn;
		std::getline(l, s, ' ');
		acn.first = stoul(s);
		std::getline(l, s, ' ');
		acn.second = stoul(s);
		std::getline(l, s, ' ');
		string pop = s;
		if (rp.pop_ac_ans.count(pop) == 0) {
			rp.pop_ac_ans[pop] = acn;
		} else {
			auto newac = rp.pop_ac_ans[pop].first+acn.first;
			auto newan = max(newac, max(rp.pop_ac_ans[pop].second, acn.second));
			rp.pop_ac_ans[pop] = make_pair(newac, newan);
		}
	}
	rp.judge(pmaf1k, pmaf10k);
	return 0;
}

void usage() {
	cerr << "rsvr [-v|--version] <command> [-S <input string> | -F <input file>] <command specific options>\n"; 
	cerr << "\tNote: For all commands STDIN is treated as input if neither the -S or -F option is specified.\n\n";
	cerr << "\tThe -h option may be passed to any command for command specific help and parameters\n";
	cerr << "\n";
	cerr << "Command options:\n";
	cerr << "\tamb: Filter input lines to retain those starting with an ambiguous RSVR ID.\n";
	cerr << "\tann: Annotate input RSVR IDs with annotation assigned to RSVR IDs in a file.\n";
	cerr << "\tdec: Decode input RSVR IDs to CHROM, POS, REF, ALT.\n";
	cerr << "\tdepth: Depth of overlaps of a set of sets of intervals.\n";
	cerr << "\tenc: Encode input variants specified as CHROM, POS, REF ALT as RSVR IDs.\n";
	cerr << "\tfiltgt: Filter input genotype information against sets of RSVR IDs or intervals from files\n";
	cerr << "\tgtnorm: Normalise genotype triples (RSVR ID, sample ID, GT) by selecting max GT per RSVR ID per sample.\n";
	cerr << "\tpmaf: Compute PMAF scores for input RSVR IDs.\n";
	cerr << "\tmix: Set operations on input RSVR IDs against those stored in a file.\n";
	cerr << "\tmerge: Merge sets of lines prefixed by RSVR IDs in two files sorted by RSVR ID.\n";
	cerr << "\nnorm: Normalise input RSVR IDs with respect to reference genome.\n";
	cerr << "\tseqfx: Annotate input RSVR IDs with transcript consequences.\n";
	cerr << "\tsort: Sort input lines with leading RSVR IDs by RSVR ID (assuming already sorted to chromosome and position).\n";
	cerr << "\ttabulate: Tabulate input RSVR IDs.\n";
	cerr << "\n\tcite: Print citation information.\n";
}

void version() {
	cout << "rsvr version " << VERSION << endl;
}

int main(int argc, char *argv[]) {
	try {
		int (*program)(cl_args&, istream&);
		if (argc < 2) {
			usage();
			return 1;
		}
		if (strcmp(argv[1], "enc") == 0)
			program = enc;
		else if (strcmp(argv[1], "depth") == 0)
			program = depth;
		else if (strcmp(argv[1], "ann") == 0)
			program = ann;
		else if (strcmp(argv[1], "merge") == 0)
			program = merge;
		else if (strcmp(argv[1], "mix") == 0)
			program = mix;
		else if (strcmp(argv[1], "tabulate") == 0)
			program = tabulate;
		else if (strcmp(argv[1], "filtgt") == 0)
			program = filtgt;
		else if (strcmp(argv[1], "dec") == 0)
			program = dec;
		else if (strcmp(argv[1], "seqfx") == 0)
			program = seqfx;
		else if (strcmp(argv[1], "norm") == 0)
			program = norm;
		else if (strcmp(argv[1], "gref") == 0)
			program = gref;
		else if (strcmp(argv[1], "sort") == 0)
			program = sort_at_cp;
		else if (strcmp(argv[1], "gtnorm") == 0)
			program = gtnorm;
		else if (strcmp(argv[1], "pmaf") == 0)
			program = pmaf;
		else if (strcmp(argv[1], "amb") == 0)
			program = amb;
		else if ((strcmp(argv[1], "--version") == 0) || (strcmp(argv[1], "-v") == 0)) {
			version();
			return 0;
		}
		else if (strcmp(argv[1], "cite") == 0) {
			cout << "Greene et al., Nature Medicine (2023), https://doi.org/10.1038/s41591-023-02211-z\n";
			return 0;
		}
		else {
			usage();
			return 1;
		}

		cl_args args(argc-2, argv+2);
		bool use_file = false;
		bool use_string = false;
		if (args.options.count('F') > 0) use_file = true;
		if (args.options.count('S') > 0) use_string = true;
		bool input_set = use_file || use_string;
		if (!(!input_set & !use_file & !use_string) && !(input_set && (use_string != use_file))) {
			throw std::invalid_argument("Specify at most one of 'F' (file) and 'S' as input.");
		}
		ifstream input_file_stream; 
		istringstream input_string_stream;
		istream& inp  = use_file ? input_file_stream : (input_set ? input_string_stream : cin);
		
		if (use_file) input_file_stream.open(args.options['F']);
		if (use_string) input_string_stream = istringstream(args.options['S']);

		if (input_set) {
			if (use_file) {
				args.options.erase(args.options.find('F'));
			} else if (use_string) {
				args.options.erase(args.options.find('S'));
			}
		}

		return program(args, inp);
	}

	catch ( const std::invalid_argument& ex ) {
		cerr << ex.what() << endl;
		return -1;
	}
	return 0;
}
