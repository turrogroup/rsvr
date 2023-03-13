#ifndef ENC_H
#define ENC_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <list> 
#include <array> 
#include <vector> 
#include <algorithm>
#include <stdexcept>
#include <memory> 
#include <limits> 

using namespace std;

uint8_t gt_code(string&, bool);

const char bases[4] = {'A', 'C', 'G', 'T'}; 

template <class T>
int64_t allele_num(T& str, int max_bases, int max_front, int from, int to) {
	int chop = max<int>(0, to-from-max_bases);
	int bases_to_read = min<int>(to-from, max_bases);
	int64_t result = 0;
	int64_t four = 1;
	for (int base_num = 0; base_num < bases_to_read; base_num++) {
		char base = str[(base_num < max_front) ? (base_num+from): (base_num+from+chop)];
		int64_t bn(0);
		switch (base) {
			case 'A':
				bn = 0;
				break;
			case 'C':
				bn = 1;
				break;
			case 'G':
				bn = 2;
				break;
			case 'T':
				bn = 3;
				break;
			case '*':
				if (result > 0) throw std::invalid_argument("Invalid ALT allele!");
				return 0;
			default:
				throw std::invalid_argument("Invalid base character in ALT field!");
		}
		result += bn * four;
		four *= 4;
	}
	return result;
}

struct variant_record {
	string chrom;
	int64_t pos;
	string ref;
	string alt;
};

int64_t get_rsvr(string&, int64_t, string&, string&, bool);
int64_t strsvr(string&, char, array<char, 4>&, array<int, 4>&, bool);
void complain_dups(int64_t);

struct sp_entry {
	int64_t rsvr;
	string str;
	sp_entry(int64_t in_rsvr, string& s);
};

template <class Tentry, class Taction>
struct read_sorted {
	struct comp {
		inline bool operator() (const Tentry& a, const Tentry& b) {
			return a.rsvr < b.rsvr;
		}
	};
	Taction& ac;
	list<Tentry> strings;
	int64_t last_rsvr_printed;
	bool first;
	int64_t chrpos;
	bool allow_dups;
	bool print_once;
	uint8_t sort_at_shift;
	int64_t toler;

	read_sorted(Taction& in_ac, bool in_allow_dups, bool in_print_once, uint8_t in_sort_at_shift = 30, int64_t in_toler = 2) : 
		ac(in_ac),
		strings(),
		last_rsvr_printed(0),
		first(true),
		chrpos(0),
		allow_dups(in_allow_dups),
		print_once(in_print_once),
		sort_at_shift(in_sort_at_shift),
		toler(in_toler)
	{
		if (!allow_dups && print_once) throw std::invalid_argument("Cannot specify '-x' and '-k'!");
	}
	void sort_and_print(int64_t upto) {
		if (!strings.empty()) {
			while (!strings.empty()) {
				auto i = strings.front();
				if ((i.rsvr == last_rsvr_printed) && !first) {
					if (!allow_dups) complain_dups(i.rsvr);
					else if (print_once) { strings.pop_front(); continue; }
				}
				if (i.rsvr < last_rsvr_printed) {
					std::ostringstream ss;
					ss << "The input is not sorted to chromosome/position: ";
					ss << (i.rsvr >> 58) << ":" << ((i.rsvr >> 30) % (1 << 28)) << " comes after "; 
					ss << (last_rsvr_printed >> 58) << ":" << ((last_rsvr_printed >> 30) % (1 << 28)); 
					std::string s = ss.str();
					throw std::invalid_argument(s);
				}
				if ((i.rsvr >> sort_at_shift) > upto) break;
				else {
					first = false;
					last_rsvr_printed = i.rsvr;
					ac.action(i);
					strings.pop_front();
				}
			}
		}
	}
	void sort_and_print() {
		sort_and_print(26LL << 58);
	}
	void add_entry(Tentry& e) {
		int64_t newchrpos = e.rsvr >> sort_at_shift;
		if (chrpos != newchrpos) {
			int64_t put = newchrpos-toler;
			sort_and_print(put);
			chrpos = newchrpos;
		}

		auto iter = strings.rbegin();
		auto end = strings.rend(); 
		while ((iter != end) && ((*iter).rsvr > e.rsvr)) {
			iter++;
		}
		strings.insert(iter.base(), e);
	}
	void print_sorted(istream& inp) {
		for (std::string line; std::getline(inp, line);) {
			Tentry e = ac.create(line);
			add_entry(e);
		}
		sort_and_print();
	}
	void print(istream& inp) {
		for (std::string line; std::getline(inp, line);) {
			Tentry e = ac.create(line);
			ac.action(e);
		}
	}
};

struct print_rec {
	void action(sp_entry&);
};

struct anno_rsvr : print_rec {
	anno_rsvr();
	sp_entry create(string& s);
};

struct sorting_printer : print_rec {
	array<char, 4> fields;
	array<int, 4> indices;
	bool normalised;
	char delim;
	bool prepend;
	sorting_printer(array<char, 4>, array<int, 4>, bool, char, bool);
	sorting_printer();
	sp_entry create(string& s);
};

struct gt_rec {
	int64_t rsvr;
	uint8_t gt;
	uint8_t pass;
	operator int64_t() const;
	gt_rec();
	gt_rec(int64_t, uint8_t, uint8_t);
	gt_rec(string&);
};	

struct gt_vec {
	gt_vec();
	vector<gt_rec> gts;
	long int cur;
	gt_rec create(string& s);
	void action(gt_rec&);
};

struct refstore {
	ifstream& sin;
	virtual void get(char&) = 0;
	virtual void ignore() = 0;
	refstore(ifstream& in_sin);
	virtual ~refstore() = default;
};

struct refstore_text : refstore {
	using refstore::refstore;
	void get(char&) override;
	void ignore() override;
};

unique_ptr<refstore> get_refstore(ifstream& sin, bool text = true);

struct ref_fasta {
	unique_ptr<refstore> ref;
	size_t chrom_num;
	size_t max_pos;
	size_t cursor;
	bool first;
	uint32_t store_bases;
	vector<char> letters;
	void next();
	ref_fasta(ifstream&, uint32_t);
	ref_fasta(ifstream&);
	void gotocoord(size_t, size_t);
	void print_ref(size_t, size_t, size_t);
	char read_base(size_t chr, size_t pos);
	char read_base(int64_t pos33);
	uint64_t refstring(size_t chr, size_t pos, size_t len);
	uint64_t refstring(int64_t pos33, size_t len);
};

struct id_dec {
	char delim;
	bool chr = false;
	string mt_str;
	id_dec(char, bool, string&);
	virtual tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> parts(int64_t);
	virtual void print_ref_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool>) = 0;
	virtual void print_alt_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool>);
	virtual void operator()(int64_t);
	virtual ~id_dec() = default;
};

struct dec_ref : id_dec {
	bool vcf_style;
	int64_t last;
	ref_fasta ref;
	dec_ref(char, bool, string&, bool, ifstream& in_ref);
	tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> parts(int64_t) override;
	void print_ref_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool>) override;
	void print_alt_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool>) override;
	void operator()(int64_t) override;
};

struct dec_noref : id_dec {
	using id_dec::id_dec;
	void print_ref_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool>) override;
};

#endif
