#include "seqfx.h"

int64_t make_norm_rsvr(int64_t orsvr, ref_fasta& ref, bool right_norm) {
	VV v(orsvr);
	int64_t orefl = min<int64_t>(v.refl(), 63);
	int64_t oaltl = min<int64_t>(v.altl(), 63);
	int s = 0;
	int e = 0;
	int64_t chr = v.chr();
	auto pos = v.pos();
	vector<char> refb(orefl);
	for (int64_t i = 0; i < orefl; i++) {
		refb[i]=ref.read_base(chr, pos+i-1);
	}
	int64_t alts = v.alts();
	vector<char> altb(oaltl);
	for (uint8_t i = 0; i < oaltl; i++) {
		altb[i] = nucs[alts % 4];
		alts /= 4LL;
	};

	while ((s < min<int64_t>(orefl, oaltl)) && (refb[s] == altb[s])) s++;
	if (right_norm) { while ((e < (min<int64_t>(orefl, oaltl)-s)) && (refb[orefl-1-e] == altb[oaltl-1-e])) e++; }
	
	int64_t aln = allele_num<vector<char> >(altb, 9, 5, s, oaltl-e);
	if (chr > 25 || chr < 1)
		throw std::invalid_argument("invalid chromosome: must belong to {1-22, X, Y, MT}!");

	return (chr << (28+6+6+18)) + ((pos+s) << (6+6+18)) + ((orefl-s-e) << (6+18)) + ((oaltl-s-e) << 18) + aln;
}

int64_t as_vep(int64_t csqs) {
	int64_t ii = (1ULL << (uint8_t)TX_CSQ_TYPE::conservative_inframe_insertion) + (1ULL << (uint8_t)TX_CSQ_TYPE::disruptive_inframe_insertion);
	int64_t id = (1ULL << (uint8_t)TX_CSQ_TYPE::conservative_inframe_deletion) + (1ULL << (uint8_t)TX_CSQ_TYPE::disruptive_inframe_deletion);
	int64_t miss = (1ULL << (uint8_t)TX_CSQ_TYPE::conservative_missense_variant) + (1ULL << (uint8_t)TX_CSQ_TYPE::non_conservative_missense_variant);
	int64_t spr = (1ULL << (uint8_t)TX_CSQ_TYPE::extended_intronic_splice_region_variant) + (1ULL << (uint8_t)TX_CSQ_TYPE::exonic_splice_region_variant);
	int64_t excl = ii | id | miss | spr;
	return (csqs | (((csqs & ii) > 0) ? (1ULL << (uint8_t)TX_CSQ_TYPE::inframe_insertion) : 0) | (((csqs & id) > 0) ? (1ULL << (uint8_t)TX_CSQ_TYPE::inframe_deletion) : 0) | (((csqs & miss) > 0) ? (1ULL << (uint8_t)TX_CSQ_TYPE::missense_variant) : 0) | (((csqs & spr) > 0) ? (1ULL << (uint8_t)TX_CSQ_TYPE::splice_region_variant) : 0)) & ~excl;
}

interval::interval(int64_t s, int64_t e) : start33(s), exclusive_end33(e) {}

int64_t interval::width() {
	return exclusive_end33-start33;
}

VV::VV(int64_t in_rsvr) : rsvr(in_rsvr) {}

bool VV::ol(int64_t s33, int64_t e33) {
	return true;
}

uint8_t VV::chr() {
	return (rsvr >> 58);
}

uint8_t VV::refl() {
	return (rsvr >> 24) % 64;
}

uint8_t VV::altl() {
	return (rsvr >> 18) % 64;
}

int64_t VV::alts() {
	return rsvr % (1 << 18);
}

int64_t VV::palts() {
	return rsvr % (1 << 18);
}

int64_t VV::pos33() { return rsvr >> 30; }
int64_t VV::pos() { return pos33() % (1 << 28); }

vector<uint8_t> VV::alt_bases() {
	auto altl1 = altl();
	vector<uint8_t> bases(altl1);
	auto alts1 = alts();
	for (uint8_t i = 0; i < altl1; i++) {
		bases[i] = alts1 % 4;
		alts1 /= 4LL;
	}
	return bases;
}

AV::AV(int64_t in_rsvr, ref_fasta& ref, bool rnorm) : VV(make_norm_rsvr(in_rsvr, ref, ((uint8_t)((in_rsvr >> 18) % 64) >= 9) ? false : rnorm)), orsvr(in_rsvr) {}

shared_ptr<transcript> string2txfeature(string& str, char delim) {
	std::istringstream l(str);
	array<string, N_GTF_FIELDS> rec ;
	for (int field_index = 0; field_index < N_GTF_FIELDS; field_index++) {
		getline(l, rec[field_index], '\t');
	}

	switch ((FEATURE_TYPE)stoul(rec[0])) {
		case FEATURE_TYPE::exon:
		case FEATURE_TYPE::five_prime_utr:
		case FEATURE_TYPE::three_prime_utr:
			return make_shared<exon>(rec);
		case FEATURE_TYPE::transcript:
			return make_shared<transcript>(rec);
		case FEATURE_TYPE::intron:
			return make_shared<intron>(rec);
		case FEATURE_TYPE::CDS:
			return make_shared<cds>(rec);
		default:
			throw invalid_argument("unrecognised feature type");
	}
}

transcript_record::transcript_record(int64_t in_csq, uint32_t in_txid, variant<monostate, int32_t> in_cdna, variant<monostate, int32_t> in_cds, int32_t in_cdna_len, variant<monostate, int32_t> in_cds_len) : 
	csq(in_csq),
	txid(in_txid),
	cdna(in_cdna),
	cds(in_cds),
	cdna_len(in_cdna_len),
	cds_len(in_cds_len)
{}

seqfx_record::seqfx_record(transcript_record& tr, uint32_t in_record_id, int64_t in_allcsqs, int64_t in_allcsqscan) : transcript_record(tr), record_id(in_record_id), allcsqs(in_allcsqs), allcsqscan(in_allcsqscan) {}

transcript_record::transcript_record() : csq(0), txid(0), cdna(monostate()), cds(monostate()), cdna_len(0), cds_len(monostate()) {}

transcript_record transcript::rec(VV& v, ref_fasta& ref) {
	return transcript_record(process(v, ref), transcript_id, monostate(), monostate(), cdna_len, monostate());
}

feature::feature(uint32_t in_id, FEATURE_TYPE f, int64_t s, int64_t e, bool str, uint32_t id) :
	interval(s, e),
	id(in_id),
	basic_overlap_type((TX_CSQ_TYPE)feature_overlaps[(uint8_t)f]),
	strand(str),
	gene_id(id)
{}

feature::feature(array<string, N_GTF_FIELDS>& rec) :
	id(stoul(rec[7])),
	interval(stoll(rec[1]), stoll(rec[2])),
	basic_overlap_type((TX_CSQ_TYPE)feature_overlaps[stol(rec[0])]),
	strand((bool)stoul(rec[3])),
	gene_id(stoul(rec[5]))
{}

int64_t feature::process(VV& v, ref_fasta& ref) {
	interval vloc(v.pos33(), v.pos33()+v.refl());
	return overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, start33, exclusive_end33) ? (1ULL << (uint8_t)basic_overlap_type) : 0;
}

void feature::printID(VV& v) {
	cout << "\t" << gene_id << "\t\t";
}

transcript::transcript(uint32_t id, FEATURE_TYPE f, int64_t s, int64_t e, bool str, uint32_t gene_id, uint32_t tx_id, int32_t in_cdna_len, bool in_canonical) :
	feature(id, f, s, e, str, gene_id),
	transcript_id(tx_id),
	cdna_len(in_cdna_len),
	canonical(in_canonical)
{}

transcript::transcript(array<string, N_GTF_FIELDS>& rec) :
	feature(rec),
	transcript_id(stoul(rec[6])),
	cdna_len(stoul(rec[13])),
	canonical(stoul(rec[16])==1)
{}

exon::exon(array<string, N_GTF_FIELDS>& rec) :
	transcript(rec),
	cdna_finish(stoll(rec[11]))
{}

int64_t exon::cdnapos(VV& v) {
	return strand ? (cdna_finish-width()+v.pos33()-start33+1) : (1+cdna_finish-width()+exclusive_end33-(v.pos33()+v.refl()));
}

void exon::printID(VV& v) {
	cout << transcript_id << "\t" << gene_id << "\t" << cdnapos(v) << "\t";
}

int64_t exon::process(VV& v, ref_fasta& ref) {
	interval vloc(v.pos33(), v.pos33()+v.refl());
	return touches<int64_t>(vloc.start33, vloc.exclusive_end33, start33, exclusive_end33) ? (1ULL << (uint8_t)basic_overlap_type) : 0;
}

transcript_record exon::rec(VV& v, ref_fasta& ref) {
	return transcript_record(process(v, ref), transcript_id, (int32_t)cdnapos(v), monostate(), cdna_len, monostate());
}

intron::intron(array<string, N_GTF_FIELDS>& rec) :
	transcript(rec),
	ex_spr({ interval(start33, start33+3), interval(exclusive_end33-3, exclusive_end33) }),
	flank8({ interval(start33+3, min(start33+11, exclusive_end33)), interval(max(exclusive_end33-11, start33), exclusive_end33-3) }),
	flank2({ interval(start33+3, min(start33+5, exclusive_end33)), interval(max(exclusive_end33-5, start33), exclusive_end33-3) })
{
	if ((exclusive_end33 - start33) <= 6)
		throw invalid_argument("intron must have greater than 0 length!");
}

interval intron::proper_intron() {
	return interval(start33+3, exclusive_end33-3);
}

int64_t intron::process(VV& v, ref_fasta& ref) {
	interval vloc(v.pos33(), v.pos33()+v.refl());
	
	bool splice_reg = overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, flank8[0].start33, flank8[0].exclusive_end33) || overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, flank8[1].start33, flank8[1].exclusive_end33);
	bool exonic_sp = overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, ex_spr[0].start33, ex_spr[0].exclusive_end33) || overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, ex_spr[1].start33, ex_spr[1].exclusive_end33);
	bool splice_start = overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, flank2[0].start33, flank2[0].exclusive_end33);
	bool splice_end = overlaps<int64_t>(vloc.start33, vloc.exclusive_end33, flank2[1].start33, flank2[1].exclusive_end33);

	if (splice_start) {
		return 1ULL << (uint8_t)(strand ? TX_CSQ_TYPE::splice_donor_variant : TX_CSQ_TYPE::splice_acceptor_variant);
	} else if (splice_end) {
		return 1ULL << (uint8_t)(!strand ? TX_CSQ_TYPE::splice_donor_variant : TX_CSQ_TYPE::splice_acceptor_variant);
	} else if (exonic_sp) {
		return 1ULL << (uint8_t)TX_CSQ_TYPE::exonic_splice_region_variant;
	} else if (splice_reg) {
		return 1ULL << (uint8_t)TX_CSQ_TYPE::splice_region_variant;
	} else {
		return overlaps(vloc.start33, vloc.exclusive_end33, start33+3, exclusive_end33-3) ? (1ULL << (uint8_t)TX_CSQ_TYPE::intron_variant) : 0;
	}
}

cds::cds(array<string, N_GTF_FIELDS>& rec) :
	exon(rec),
	phase(stol(rec[4])),
	cds_finish(stoll(rec[12])),
	cds_len(stol(rec[14])),
	first_cds((bool)(stol(rec[10]))),
	last_cds((bool)(stol(rec[15]))),
	first_codon{ rec[8][0], rec[8][1], rec[8][2] },
	last_codon{ rec[9][0], rec[9][1], rec[9][2] }
{}

interval cds::containing_orf() {
	if (strand) {
		int64_t s = start33 - ((phase == 0) ? 0 : (3-phase));
		int64_t e = s + ((exclusive_end33-s)/3)*3;
		if (e < exclusive_end33) e+=3;
		return interval(s, e);
	} else {
		int64_t e = exclusive_end33 + ((phase == 0) ? 0 : (3-phase));
		int64_t s = e - ((e-start33)/3)*3;
		if (s > start33) s-=3;
		return interval(s, e);
	}
}

uint8_t cds::cds_base(int64_t base33, ref_fasta& ref) {
	uint8_t result(0);
	interval corf = containing_orf();

	if ((base33 < corf.start33) || (base33 >= corf.exclusive_end33))
		throw invalid_argument("trying to read base outside the containing ORF!");

	if (base33 < start33) {
		int64_t ind = base33-corf.start33;
		if (abs(ind) >= 3) throw invalid_argument("trying to read outside ORF!");
		result = first_codon[ind];
	} else if (base33 >= exclusive_end33) {
		int64_t ind = 3-(corf.exclusive_end33-base33);
		if (abs(ind) >= 3) throw invalid_argument("trying to read outside ORF!");
		result = last_codon[ind];
	} else {
		result = ref.read_base(base33-1);
	}
	switch (result) {
		case 'A':
			result =  0;
			break;
		case 'C':
			result =  1;
			break;
		case 'G':
			result =  2;
			break;
		case 'T':
			result =  3;
			break;
		default:
			result = 4;
			break;
	}
	return result;
}

variant<monostate, pair<int8_t, char> > cds::read_codon(int64_t codon, ref_fasta& ref) {
	interval corf = containing_orf();
	if (codon >= corf.width()/3) {
		return monostate();
	} else {
		uint8_t pow = 1;
		int8_t cod_val = 0;
		for (int64_t i = 0; i < 3; i++) {
			int64_t pos = strand ? (corf.start33 + codon * 3 + i) : (corf.exclusive_end33 - 1 - codon * 3 - i);
			uint8_t base = cds_base(pos, ref);

			if (base > 3) return make_pair(-1, 'U');
			cod_val += pow * (strand ? base : (3-base));
			pow *= 4;
		}
		return make_pair(cod_val, (first_cds && (codon == 0) && AA[cod_val] == 'M') ? 'O' : AA[cod_val]);
	}
}

variant<monostate, pair<int8_t, char> > cds::read_codon_with_variant(int64_t codon, interval& vloc, vector<uint8_t>& alt_bases, ref_fasta& ref) {
	
	interval corf = containing_orf();
	uint8_t pow = 1;
	int8_t cod_val = 0;
	int64_t nalt = alt_bases.size();
	int64_t m = ((corf.width()+nalt-vloc.width())/3);
	if (codon >= m) {
		return monostate();
	} else {
		for (int64_t i = 0; i < 3; i++) {
			if (strand) {
				interval use_alt(vloc.start33, vloc.start33 + nalt);
				int64_t pos = corf.start33 + codon * 3 + i;
				uint8_t base(0);
				if (pos < use_alt.start33) {
					base = cds_base(pos, ref);
				}
				else if ((pos >= use_alt.start33) && (pos < use_alt.exclusive_end33)) {
					base = alt_bases[pos-use_alt.start33];
				} else if (pos >= use_alt.exclusive_end33) {
					base = cds_base(pos+vloc.exclusive_end33-vloc.start33-use_alt.exclusive_end33+use_alt.start33, ref);
				}

				if (base > 3) return make_pair(-1, 'U');
				cod_val += pow * base;
			} else {
				interval use_alt(vloc.exclusive_end33-nalt, vloc.exclusive_end33);

				int64_t pos = corf.exclusive_end33 - 1 - codon * 3 - i;
				uint8_t base(0);
				if (pos >= use_alt.exclusive_end33) {
					base = cds_base(pos, ref);
				} else if ((pos >= use_alt.start33) && (pos < use_alt.exclusive_end33)) {
					base = alt_bases[pos-use_alt.start33];
				} else if (pos < use_alt.start33) {
					base = cds_base(pos-vloc.exclusive_end33+vloc.start33+use_alt.exclusive_end33-use_alt.start33, ref);
				}
				if (base > 3) return make_pair(-1, 'U');
				cod_val += pow * (3-base);
			}
			pow *= 4;
		}

		return make_pair(cod_val, (first_cds && (codon == 0) && AA[cod_val] == 'M') ? 'O' : AA[cod_val]);
	}
}

int64_t cds::length() {
	return ((exclusive_end33-phase-start33)/3);
}

int64_t cds::cdspos(VV& v) {
	return strand ? (cds_finish-width()+v.pos33()-start33+1) : (1+cds_finish-width()+exclusive_end33-(v.pos33()+v.refl()));
}

void cds::printID(VV& v) {
	cout << transcript_id << "\t" << gene_id << "\t" << cdnapos(v) << "\t" << cdspos(v);
}

bool cds::affects_cds(VV& v) {
	interval vloc(v.pos33(), v.pos33() + v.refl());
	bool admissible = false;
	if (touches<int64_t>(start33, exclusive_end33, v.pos33(), v.pos33()+v.refl())) {
		if (first_cds) {
			if ((strand && (vloc.exclusive_end33 > start33)) || (!strand && (vloc.start33 < exclusive_end33))) {
				admissible = true;
			}
		} else if (last_cds) {
			if ((strand && (vloc.start33 < exclusive_end33)) || (!strand && (vloc.exclusive_end33 > start33))) {
				admissible = true;
			}
		} else {
			admissible = true;
		}
	}	
	return admissible;
}

int64_t cds::process(VV& v, ref_fasta& ref) {
	
	int64_t result(0);
	interval corf = containing_orf();
	interval vloc(v.pos33(), v.pos33() + v.refl());

	bool admissible = affects_cds(v);

	if (admissible) {
		if (!contains<int64_t>(start33, exclusive_end33, v.pos33(), v.pos33()+v.refl())) {
			result |= 1ULL << (uint8_t)TX_CSQ_TYPE::coding_sequence_variant;
		} else {
			int8_t d = v.refl() - v.altl();
			if (d != 0) {
				bool disrup = abs(v.pos33() - corf.start33) % 3 != 0;
				result |= (1ULL << (uint8_t)(((d % 3) != 0) ? TX_CSQ_TYPE::frameshift_variant : ((d < 0) ? (disrup ? TX_CSQ_TYPE::disruptive_inframe_insertion : TX_CSQ_TYPE::conservative_inframe_insertion) : (disrup ? TX_CSQ_TYPE::disruptive_inframe_deletion : TX_CSQ_TYPE::conservative_inframe_deletion))));
			}
		
			int8_t sign = strand ? 1 : -1;

			auto bases = v.alt_bases();
			int64_t first_codon_to_read = max<int64_t>(0, strand ? ((v.pos33()-corf.start33)/3) : ((corf.exclusive_end33-v.pos33()-v.refl())/3));
			int max_codons = (max(v.refl(), v.altl())+4)/3;

			if (((v.rsvr >> 18) % 64) <= 9) {
				int codons_read = 0;

				array<int, 2> started = {false, false};
				array<int, 2> stopped = {false, false};
				array<int, 2> read = {true, true};

				bool exists_unknown = false;

				int64_t csq_if_length_preserved(0);

				while ((codons_read < max_codons) && (read[0] || read[1])) {
					
					int codon = first_codon_to_read + codons_read;

					array<variant<monostate, pair<int8_t, char> >, 2> cods = {
						read[0] ? read_codon(codon, ref) : monostate(),
						read[1] ? read_codon_with_variant(codon, vloc, bases, ref) : monostate()
					};

					for (int g = 0; g < 2; g++) {
						read[g] = cods[g].index()!=0;

						if (read[g]) {
							char b = get<1>(cods[g]).second;
							started[g] |= b == 'O';
							stopped[g] |= b == 'X';
						}
					}

					if ((d == 0) && (read[0] != read[1]))
						throw invalid_argument("misaligned missense!");

					if (read[0] && read[1]) {
						auto c0 = get<1>(cods[0]);
						auto c1 = get<1>(cods[1]);
						auto c0b = c0.second;
						auto c1b = c1.second;
						if ((c0b == 'U') || (c1b == 'U')) {
							exists_unknown |= true;
						} else if (c0.first != c1.first) {
							
							if (c0b != c1b) {
								csq_if_length_preserved |= (1ULL << (uint8_t)((AA_type[c0.first] != AA_type[c1.first]) ? TX_CSQ_TYPE::non_conservative_missense_variant : TX_CSQ_TYPE::conservative_missense_variant));
							} else if (c0b == 'X' && c1b == 'X') {
								csq_if_length_preserved	|= (1ULL << (uint8_t)(TX_CSQ_TYPE::stop_retained_variant));
							} else {
								csq_if_length_preserved	|= (1ULL << (uint8_t)(TX_CSQ_TYPE::synonymous_variant));
							}
						}
					}

					codons_read++;
				}

				if (!stopped[0] && stopped[1]) 
					result |= (1ULL << (uint8_t)(TX_CSQ_TYPE::stop_gained));
				else if (!stopped[1] && stopped[0]) 
					result |= (1ULL << (uint8_t)(TX_CSQ_TYPE::stop_lost));
				else if (!started[1] && started[0]) 
					result |= (1ULL << (uint8_t)(TX_CSQ_TYPE::start_lost));
				else {
					
					if (d==0) {
						result |= csq_if_length_preserved;
					}
					if (exists_unknown)
						result |= (1ULL << (uint8_t)(TX_CSQ_TYPE::coding_sequence_variant));

				}
			} else if (d == 0) {
				
				result |= (1ULL << (uint8_t)(TX_CSQ_TYPE::coding_sequence_variant));
			}
		}
	}

	return result;
}

transcript_record cds::rec(VV& v, ref_fasta& ref) {
	return transcript_record(process(v, ref), transcript_id, (int32_t)cdnapos(v), (int32_t)cdspos(v), cdna_len, cds_len);
}


tater2::tater2(feature_stack<transcript>& in_fs, int64_t in_left_padding, ref_fasta& in_ref, AGGREGATE in_agg) :
	fs(in_fs),
	left_padding(in_left_padding),
	ref(in_ref),
	agg(in_agg)
{}

tater2::tater2(feature_stack<transcript>& in_fs, ref_fasta& in_ref) : tater2(in_fs, 63, in_ref, AGGREGATE::GENE) {}

void tater2::anno(VV& av) {
	fs.push_features(av.pos33(), left_padding);
	for (uint32_t i = 0; i < fs.features.size(); i++) {
		shared_ptr<feature> f = fs.features[i];
		cout << av.rsvr << '\t' << f->id << '\t'; 
		f->printID(av);
		cout << '\t' << f->process(av, ref) << endl;
	}
}

vector<seqfx_record> tater2::anno_gene(VV& av) {
	fs.push_features(av.pos33(), left_padding);
	map<uint32_t, vector<pair<bool, transcript_record> > > txmap;
	for (uint32_t i = 0; i < fs.features.size(); i++) {
		shared_ptr<transcript> f = fs.features[i];
		uint32_t grp = (agg == AGGREGATE::GENE) ? f->gene_id : f->transcript_id;
		transcript_record tr = f->rec(av, ref);
		if (tr.csq > 0) {
			if (txmap.count(grp) == 0)
				txmap[grp] = vector<pair<bool, transcript_record> >(0);
			txmap[grp].push_back(make_pair(f->canonical, tr));
		}
	}

	vector<seqfx_record> result;
	for (const auto &keyval : txmap) {
		int64_t worst_rec_can = 0;
		int64_t worst_rec_ind = 0;
		int64_t worst_feat_csq = 0;
		
		int32_t worst_len = 0;
		int64_t allcsqs = 0;
		int64_t allcsqscan = 0;
		for (int i = 0; i < keyval.second.size(); i++) {
			auto fcsq = keyval.second[i].second;
			auto newcsq = fcsq.csq;
			auto newlen = fcsq.cds_len.index() != 0 ? get<int32_t>(fcsq.cds_len) : fcsq.cdna_len;
			auto newcan = keyval.second[i].first;
			if (newcan > worst_rec_can) {
				worst_rec_ind = i;
				worst_rec_can = newcan;
				worst_feat_csq = newcsq;
				worst_len = newlen;
			} else if (newcan == worst_rec_can) {
				if (newcsq > worst_feat_csq) {
					worst_rec_ind = i;
					worst_feat_csq = newcsq;
					worst_len = newlen;
				} else if (newcsq == worst_feat_csq) {
					if (newlen > worst_len) {
						worst_rec_ind = i;
						worst_len = newlen;
					}
				}
			}
			allcsqs |= newcsq;
			if (newcan) allcsqscan |= newcsq;
		}
		transcript_record tr = keyval.second[worst_rec_ind].second;
		result.push_back(seqfx_record(tr, keyval.first, allcsqs, allcsqscan));		
	}
	return result;
}

seqfx_print::seqfx_print(bool in_vep, char in_delim) : vep(in_vep), delim(in_delim) {}

void seqfx_print::recstr(seqfx_record& rec) {
	cout << rec.record_id << delim;
	cout << rec.txid << delim;
	cout << ((rec.cdna.index()==0) ? "" : to_string(get<int32_t>(rec.cdna))) << delim;
	cout << ((rec.cds.index()==0) ? "" : to_string(get<int32_t>(rec.cds))) << delim;
	cout << (vep ? as_vep(rec.csq) : rec.csq) << delim;
	cout << (vep ? as_vep(rec.allcsqs) : rec.allcsqs) << delim;
	cout << (vep ? as_vep(rec.allcsqscan) : rec.allcsqscan);
}

seqfx_line_print::seqfx_line_print(bool in_vep, char in_delim, char in_sep_sm, char in_sep_bg) : seqfx_print(in_vep, in_delim), sep_sm(in_sep_sm), sep_bg(in_sep_bg) {}

void seqfx_line_print::print(string& l, int64_t rsvr, vector<seqfx_record>& sfx) {
	cout << l << sep_bg;
	for (size_t i = 0; i < sfx.size(); i++) {
		if (i > 0) {
			cout << sep_sm;
		}
		recstr(sfx[i]);
	}
	cout << endl;
};

void seqfx_db_print::print(string& l, int64_t rsvr, vector<seqfx_record>& sfx) {
	for (auto &i : sfx) {
		cout << rsvr << delim;
		recstr(i);
		cout << endl;
	}
};

