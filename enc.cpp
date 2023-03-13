#include "enc.h"

uint8_t gt_code(string& gtstr, bool break_if_unrecognised = true) {
	if ((gtstr.length() == 3) && ((gtstr[1] == '/' || gtstr[1] == '|'))) {
		return (gtstr[0] == '1') + (gtstr[2] == '1');
	} else if (gtstr.length() == 1 && (gtstr[0] == '1')) {
		return 1;
	} else if (break_if_unrecognised) {
		throw invalid_argument("Unrecognised genotype! " + gtstr);
	} else {
		return 0;
	}
}

int64_t get_rsvr(string& chrom, int64_t pos, string& ref, string& alt, bool normalise) {
	int64_t orefl = ref.length();
	int64_t oaltl = alt.length();
	int s = 0;
	int e = 0;
	if (normalise) {
		while ((s < min<int64_t>(orefl, oaltl)) && (ref[s] == alt[s])) s++;
		while ((e < (min<int64_t>(orefl, oaltl)-s)) && (ref[orefl-1-e] == alt[oaltl-1-e])) e++;
	}
	
	int64_t aln = allele_num<string>(alt, 9, 5, s, oaltl-e);
	int64_t chrom_num = 0;
	if (chrom.compare("MT")==0) chrom_num = 25;
	else if (chrom.compare("M")==0) chrom_num = 25;
	else if (chrom.compare("Y")==0) chrom_num = 24;
	else if (chrom.compare("X")==0) chrom_num = 23;
	else chrom_num = std::stoi(chrom);
	if (chrom_num > 25 || chrom_num < 1)
		throw std::invalid_argument("Invalid chromosome: must belong to {1-22, X, Y, MT/M}!");

	return (chrom_num << (28+6+6+18)) + ((pos+s) << (6+6+18)) + (min<int64_t>(orefl-s-e, 63) << (6+18)) + (min<int64_t>(oaltl-s-e, 63) << 18) + aln;
}

int64_t strsvr(string& str, char delim, array<char, 4>& prop, array<int, 4>& prop_ind, bool normalise) {
	std::istringstream l(str);
	int field_index = 1;
	int cur_prop = 0;
	variant_record v;
	string field;

	while (std::getline(l, field, delim) || field.empty()) {
		if (field_index == prop_ind[cur_prop]) {
			switch (prop[cur_prop]) {
				case 'c':
					v.chrom = field;
					break;
				case 'p':
					v.pos = stoi(field);
					break;
				case 'r':
					v.ref = field;
					break;
				case 'a':
					v.alt = field;
					break;
			}
			cur_prop++;
		}
		field_index++;
		if (!l && field.empty())
			break;
	}
	if (cur_prop != 4) {
		throw std::invalid_argument("Could not find all of CHROM (c), POS (p), REF (r) and ALT (a)!");
	}

	return get_rsvr(v.chrom, v.pos, v.ref, v.alt, normalise);
}

void complain_dups(int64_t rsvr) {
	std::ostringstream ss;
	ss << "Duplicate RSVR ID detected! " << rsvr << endl;
	std::string s = ss.str();
	throw std::invalid_argument(s);
}

sp_entry::sp_entry(int64_t in_rsvr, string& s) : rsvr(in_rsvr), str(s) {}

void print_rec::action(sp_entry& e) {
	std::cout << e.str << endl;
}

anno_rsvr::anno_rsvr() {}

sp_entry anno_rsvr::create(string& s) {
	std::istringstream l(s);
	int64_t id;
	l >> id;
	return sp_entry(id, s);
}

sorting_printer::sorting_printer(
	array<char, 4> in_fields, 
	array<int, 4> in_indices, 
	bool in_normalised, 
	char in_delim,
	bool in_prepend
) :
	fields(in_fields),
	indices(in_indices),
	normalised(in_normalised),
	delim(in_delim),
	prepend(in_prepend)
{}

sorting_printer::sorting_printer() :
	fields({'c','p','r','a'}),
	indices({1,2,3,4}),
	normalised(true),
	delim('\t')
{}

sp_entry sorting_printer::create(string& s) {
	int64_t RSVR = strsvr(s, delim, fields, indices, normalised);
	string to_print = std::to_string(RSVR) + (prepend ? delim + s : "");
	return sp_entry(RSVR, to_print);
}

gt_rec::gt_rec() {}

gt_rec::gt_rec(int64_t in_rsvr, uint8_t in_gt, uint8_t in_pass) : rsvr(in_rsvr), gt(in_gt), pass(in_pass) {}

gt_rec::gt_rec(string& s) {
	std::istringstream l(s);
	variant_record v;
	std::getline(l, v.chrom, '\t');
	string postr; std::getline(l, postr, '\t'); v.pos = stoi(postr);
	std::getline(l, v.ref, '\t');
	std::getline(l, v.alt, '\t');
	string gtstr; std::getline(l, gtstr, '\t');
	uint8_t gtc = gt_code(gtstr, false);
	string passstr; std::getline(l, passstr, '\t');
	pass = (uint8_t)(stoi(passstr));
	rsvr = get_rsvr(v.chrom, v.pos, v.ref, v.alt, true);
	if (gtc == 0) throw std::invalid_argument("Attempting to add empty genotype!");
	gt = gtc;
}

gt_rec::operator int64_t() const { return rsvr; }

gt_vec::gt_vec() : gts(0), cur(-1) {}

void gt_vec::action(gt_rec& x) {
	if (x.rsvr < cur) {
		throw std::invalid_argument("Input not sorted!");
	} else if (x.rsvr == cur) {
		auto i = gts.size()-1;
		gts[i].gt = max(gts[i].gt, x.gt);	
		gts[i].pass = max(gts[i].pass, x.pass);	
	} else {
		gts.push_back(gt_rec(x.rsvr, x.gt, x.pass));
		cur = x.rsvr;
	}
}

gt_rec gt_vec::create(string& s) {
	return gt_rec(s);
}

refstore::refstore(ifstream& in_sin) : sin(in_sin) {}

void refstore_text::get(char& a) {
	sin.get(a);
	if (sin.eof()) throw std::invalid_argument("Read beyond end of reference!");
}

void refstore_text::ignore() {
	sin.ignore(numeric_limits<streamsize>::max(), '\n');
}

unique_ptr<refstore> get_refstore(ifstream& sin, bool text) {
	return make_unique<refstore_text>(sin);
}

ref_fasta::ref_fasta(ifstream& in_ref, uint32_t in_store_bases) : 
	ref(get_refstore(in_ref)), 
	chrom_num(1), 
	max_pos(0), 
	cursor(0),
	first(true),
	store_bases(in_store_bases),
	letters(store_bases)
{}

ref_fasta::ref_fasta(ifstream& in_ref) : ref_fasta(in_ref, 1000) {}

void ref_fasta::next() {
	ref->get(letters[cursor]);
	if (letters[cursor] == '>') {
		ref->ignore();
		if (!first) chrom_num++;
		max_pos = 0;
		first = false;
	} else if (first) {
		throw invalid_argument("Unexpected '>'!");
	} else if (!isspace(letters[cursor])) {
		cursor = (cursor + 1) % store_bases;
		max_pos++;
	}
}

void ref_fasta::gotocoord(size_t chr, size_t pos) {
	auto start = max_pos > (store_bases - 1) ? (max_pos + 1 - store_bases) : 0;
	while (chr != chrom_num) next();
	while (pos > max_pos) next();
}

void ref_fasta::print_ref(size_t chr, size_t pos, size_t len) {
	gotocoord(chr, pos + len);
	for (size_t i = 0; i < len; i++) {
		size_t ind = (store_bases + cursor-(max_pos-pos)+i) % store_bases;
		cout << letters[ind];
	}
}

char ref_fasta::read_base(size_t chr, size_t pos) {
	gotocoord(chr, pos + 1);
	return letters[(store_bases + cursor-(max_pos-pos)) % store_bases];
}

char ref_fasta::read_base(int64_t pos33) {
	return read_base(pos33 >> 28, pos33 % (1 << 28));
}

uint64_t ref_fasta::refstring(size_t chr, size_t pos, size_t len) {
	gotocoord(chr, pos + len);

	uint64_t result = 0ULL;
	uint64_t pow = 1ULL;
	for (size_t i = 0; i < len; i++) {
		size_t ind = (store_bases + cursor-(max_pos-pos)+i) % store_bases;
		uint64_t bn = 0ULL;
		switch (letters[ind]) {
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
		}
		result += pow * bn;
		pow *= 4ULL;
	}
	return result;
}

uint64_t ref_fasta::refstring(int64_t pos33, size_t len) {
	return refstring(pos33 >> 28, pos33 % (1 << 28), len);
}

id_dec::id_dec(char in_delim, bool in_chr, string& in_mt_str) : delim(in_delim), chr(in_chr), mt_str(in_mt_str) {} 

void id_dec::operator()(int64_t v) {
	auto p = parts(v);
	auto chrnum = get<0>(p);
	if (chr) cout << "chr";
	switch (chrnum) {
		case 25:
			cout << mt_str;
			break;
		case 24:
			cout << "Y";
			break;
		case 23:
			cout << "X";
			break;
		default:
			if ((chrnum < 1) || (chrnum > 25))
				throw invalid_argument("Unrecognised chromosome!");
			else {
				cout << chrnum;
				break;
			}
	}

	cout << delim << get<1>(p) << delim;

	print_ref_data(p);
	cout << delim;

	print_alt_data(p);
	cout << endl;
}

tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> id_dec::parts(int64_t v) {
	int64_t chrnum = v >> 58;
	int64_t pos = (v >> 30) % 268435456;
	int64_t refl = (v >> 24) % 64;
	int64_t altl = (v >> 18) % 64;
	int64_t alt = v % (1 << 18);
	return make_tuple(chrnum, pos, refl, altl, alt, false);
}

void id_dec::print_alt_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> p) {
	int done = 0;
	auto altl = get<3>(p);
	auto alt = get<4>(p);
	while (altl > 0) {
		if ((altl <= 4) || (done < 5)) {
			cout << bases[alt % 4];
			alt /= 4LL;
		} else {
			cout << '?';
		} 
		altl--;
		done++;
	};
}

dec_ref::dec_ref(char in_delim, bool in_chr, string& in_mt_str, bool in_vcf_style, ifstream& in_ref) : id_dec(in_delim, in_chr, in_mt_str), vcf_style(in_vcf_style), ref(in_ref), last(0) {}

tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> dec_ref::parts(int64_t v) {
	auto p = id_dec::parts(v);
	get<5>(p) = vcf_style && ((get<2>(p) == 0) || (get<3>(p) == 0));
	if (get<5>(p)) {
		get<1>(p) -= 1;
		get<2>(p) += 1;
	}
	return p;
}

void dec_ref::print_ref_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> p) {
	ref.print_ref(get<0>(p), get<1>(p)-1, get<2>(p));
}

void dec_ref::print_alt_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> p) {
	if (get<5>(p)) cout << bases[ref.refstring(get<0>(p), get<1>(p)-1, 1)];
	id_dec::print_alt_data(p);
}

void dec_ref::operator()(int64_t v) {
	if (v < last) {
		throw invalid_argument("Input variants not sorted by RSVR!");
	}
	id_dec::operator()(v);
	last = v;
}

void dec_noref::print_ref_data(tuple<int64_t, int64_t, int64_t, int64_t, int64_t, bool> p) {
	cout << get<2>(p);
}
