#ifndef SEQFX_H
#define SEQFX_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <memory> 
#include <array> 
#include <vector> 
#include <map>
#include <variant> 
#include <algorithm>
#include <stdexcept>
#include "enc.h"

using namespace std;

struct interval {
	int64_t start33;
	int64_t exclusive_end33;
	int64_t width();
	interval(int64_t,int64_t);  
};

const char AA[64] = {'K','Q','E','X','T','P','A','S','R','R','G','X','I','L','V','L','N','H','D','Y','T','P','A','S','S','R','G','C','I','L','V','F','K','Q','E','X','T','P','A','S','R','R','G','W','M','L','V','L','N','H','D','Y','T','P','A','S','S','R','G','C','I','L','V','F'};
const uint8_t AA_type[64] = { 5,6,6,7,2,3,1,2,5,5,1,7,1,1,1,1,6,5,6,4,2,3,1,2,2,5,1,2,1,1,1,4,5,6,6,7,2,3,1,2,5,5,1,4,2,1,1,1,6,5,6,4,2,3,1,2,2,5,1,2,1,1,1,4 };

const char nucs[4] = {'A','C','G','T'};

enum class TX_CSQ_TYPE {
transcript_variant=0,
intron_variant=1,
exon_variant=2,
five_prime_UTR_variant=3,
three_prime_UTR_variant=4,
coding_sequence_variant=5,
synonymous_variant=6,
splice_region_variant=7,
extended_intronic_splice_region_variant=8,
exonic_splice_region_variant=9,
stop_retained_variant=10,
missense_variant=11,
conservative_missense_variant=12,
non_conservative_missense_variant=13,
inframe_deletion=14,
conservative_inframe_deletion=15,
disruptive_inframe_deletion=16,
inframe_insertion=17,
conservative_inframe_insertion=18,
disruptive_inframe_insertion=19,
start_lost=20,
stop_lost=21,
frameshift_variant=22,
stop_gained=23,
splice_donor_variant=24,
splice_acceptor_variant=25
};

enum class FEATURE_TYPE {
CDS=0,
exon=1,
five_prime_utr=2,
three_prime_utr=3,
transcript=4,
intron=5
};

enum class AGGREGATE { GENE=0, TRANSCRIPT=1 };

const TX_CSQ_TYPE feature_overlaps[6] = { 
	TX_CSQ_TYPE::coding_sequence_variant,
	TX_CSQ_TYPE::exon_variant,
	TX_CSQ_TYPE::five_prime_UTR_variant,
	TX_CSQ_TYPE::three_prime_UTR_variant,
	TX_CSQ_TYPE::transcript_variant,
	TX_CSQ_TYPE::intron_variant
};

const uint32_t N_GTF_FIELDS = 17;

int64_t make_norm_rsvr(int64_t, ref_fasta&, bool);
int64_t as_vep(int64_t);

struct VV {
	int64_t rsvr;
	VV(int64_t);
	bool ol(int64_t, int64_t);
	uint8_t chr();
	uint8_t refl();
	uint8_t altl();
	int64_t alts();
	int64_t palts();
	int64_t pos();
	int64_t pos33();
	vector<uint8_t> alt_bases();
};

struct AV : VV {
	int64_t orsvr;
	AV(int64_t, ref_fasta&, bool);
	
};

template <class T>
inline bool overlaps(T a_first, T a_second, T b_first, T b_second) {
	return (a_second > b_first) && (b_second > a_first);
}

template <class T>
inline bool touches(T a_first, T a_second, T b_first, T b_second) {
	return (a_second >= b_first) && (b_second >= a_first);
}

template <class T>
inline bool contains(T a_first, T a_second, T b_first, T b_second) {
	
	return (a_first <= b_first) && (a_second >= b_second);
}

struct transcript_record {
	int64_t csq;
	uint32_t txid;
	variant<monostate, int32_t> cdna;
	variant<monostate, int32_t> cds;
	int32_t cdna_len;
	variant<monostate, int32_t> cds_len;
	transcript_record();
	transcript_record(int64_t, uint32_t, variant<monostate, int32_t>, variant<monostate, int32_t>, int32_t, variant<monostate, int32_t>);  
};

struct seqfx_record : transcript_record {
	uint32_t record_id;
	int64_t allcsqs;
	int64_t allcsqscan;
	seqfx_record(transcript_record&, uint32_t, int64_t, int64_t);
};

struct feature : interval {
	uint32_t id;
	TX_CSQ_TYPE basic_overlap_type;
	bool strand;
	uint32_t gene_id;
	feature(uint32_t, FEATURE_TYPE, int64_t, int64_t, bool, uint32_t); 
	feature(array<string, N_GTF_FIELDS>&);
	virtual int64_t process(VV&, ref_fasta&);
	virtual void printID(VV&);
};

struct transcript : feature {
	uint32_t transcript_id;
	int32_t cdna_len;
	bool canonical;
	transcript(uint32_t, FEATURE_TYPE, int64_t, int64_t, bool, uint32_t, uint32_t, int32_t, bool); 
	transcript(array<string, N_GTF_FIELDS>&);
	
	
	virtual transcript_record rec(VV&, ref_fasta&);
};

struct exon : transcript {
	int64_t cdna_finish;
	int64_t cdnapos(VV&);
	exon(array<string, N_GTF_FIELDS>&);
	void printID(VV&) override;
	int64_t process(VV&, ref_fasta&) override;
	
	transcript_record rec(VV&, ref_fasta&) override;
};

struct intron : transcript {
	array<interval, 2> ex_spr;
	array<interval, 2> flank8;
	array<interval, 2> flank2;
	intron(array<string, N_GTF_FIELDS>&);
	interval proper_intron();
	int64_t process(VV&, ref_fasta&) override;
};

struct cds : exon {
	uint8_t phase;
	int64_t cds_finish;
	int32_t cds_len;
	bool first_cds;
	bool last_cds;
	char first_codon[3];
	char last_codon[3];
	cds(array<string, N_GTF_FIELDS>&);
	interval containing_orf();
	uint8_t cds_base(int64_t, ref_fasta&);
	variant<monostate, pair<int8_t, char> > read_codon(int64_t codon, ref_fasta& ref);
	variant<monostate, pair<int8_t, char> > read_codon_with_variant(int64_t codon, interval& vloc, vector<uint8_t>& alt_bases, ref_fasta& ref);
	int64_t length();
	int64_t orfs();
	int64_t cdspos(VV&);
	void printID(VV&) override;
	bool affects_cds(VV&);
	int64_t process(VV&, ref_fasta&) override;
	
	transcript_record rec(VV&, ref_fasta&) override;
};

shared_ptr<transcript> string2txfeature(string& str, char delim);

template <class Tbasefeature>
struct feature_stack {
	ifstream& gtf;
	shared_ptr<Tbasefeature> cur_feature;
	vector<shared_ptr<Tbasefeature> > features;

	feature_stack(ifstream& in_gtf) :
		gtf(in_gtf),
		cur_feature(nullptr),
		features(0)
	{
		pull_feature();
	}

	bool pull_feature() {
		std::string line; 
		bool result = (bool)std::getline(gtf, line);
		if (result) { 
			cur_feature = string2txfeature(line, '\t'); 
		}
		return result;
	}

	void push_features(int64_t loc33, int64_t lpad) {
		while (loc33 >= (cur_feature->start33-lpad)) {
			features.push_back(cur_feature);
			if (!pull_feature())
				break;
		}

		auto end = std::remove_if(
			features.begin(), 
			features.end(), 
			[loc33,lpad](const shared_ptr<Tbasefeature>& f) { 
				return !(((f->start33-lpad) <= loc33) && (f->exclusive_end33 >= loc33)); 
			}
		);

		features.erase(end, features.end());
	}
};

struct tater2 {
	feature_stack<transcript>& fs;
	int64_t left_padding;
	AGGREGATE agg;
	ref_fasta& ref;
	tater2(feature_stack<transcript>&, int64_t, ref_fasta&, AGGREGATE);
	tater2(feature_stack<transcript>&, ref_fasta&);
	void anno(VV&);
	vector<seqfx_record> anno_gene(VV&);
};

struct seqfx_print {
	bool vep;
	char delim;
	void recstr(seqfx_record& rec);
	seqfx_print(bool, char);
	virtual void print(string&, int64_t rsvr, vector<seqfx_record>& sfx) = 0;
	virtual ~seqfx_print() = default;
};

struct seqfx_line_print : seqfx_print {
	char sep_sm;
	char sep_bg;
	seqfx_line_print(bool, char, char, char);
	void print(string&, int64_t rsvr, vector<seqfx_record>& sfx) override;
};

struct seqfx_db_print : seqfx_print {
	using seqfx_print::seqfx_print;
	void print(string&, int64_t rsvr, vector<seqfx_record>& sfx) override;
};

#endif

