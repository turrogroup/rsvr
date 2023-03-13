#include "idtypes.h"

rsvr_record::rsvr_record(size_t rows) : n(rows), record(n), cur(0) {}
void rsvr_record::set_record_value(int64_t v) {
	record[cur++] = v;
}

id_file::id_file(istream& in_f) : f(in_f) { }
bool id_file::load(int64_t& val) { return !(f >> val); };

id_prep_cin::id_prep_cin(istream& in_f) : f(in_f), line("") { }
id_prep_cin::id_prep_cin() : id_prep_cin(cin) { }

bool id_prep_cin::load(int64_t& val) { 
	bool fin = !(getline(f, line)); 
	if (!fin) {
		istringstream l(line);
		l >> val;
	}
	return fin;
};

size_t f_ids::yield() { return 0; }
void f_ids::stamp() { cout << T::x.line << endl; }

void exclude_ids::fx(int64_t v) { stamp(); }

void intersect_ids::fxy(int64_t v) { stamp(); }

void merge_lines::fx(int64_t v) { cout << x.line << endl; }
void merge_lines::fxy(int64_t v) { cout << x.line << endl; cout << y.line << endl; }
void merge_lines::fy(int64_t v) { cout << y.line << endl; }
size_t merge_lines::yield() { return 0; }

anno_ids::anno_ids(
	ids_file<id_prep_cin>& x, 
	ids_file<id_file>& y,
	anno_type in_type
) :
	f_ids(x, y),
	seen_y(false),
	last_y(0LL),
	type(in_type)
{}

void anno_ids::anno(char c) { 
	switch (type) {
		case anno_type::anno:
			cout << x.line << '\t' << c << endl; 
			break;
		case anno_type::intersect:
			if (c == '1') stamp();
			break;
		case anno_type::exclude:
			if (c == '0') stamp();
			break;
	}
}
void anno_ids::fx(int64_t v) { if (seen_y && (v == last_y)) anno('1'); else anno('0'); }
void anno_ids::fxy(int64_t v) { seen_y = true; last_y = v; anno('1'); }

size_t anno_with::yield() { return 0; }

anno_one::anno_one(
	ids_file<id_prep_cin>& x, 
	ids_file<id_prep_cin>& y,
	string& in_mng
) :
	anno_with(x, y),
	mng(in_mng)
{}

void anno_one::fx(int64_t v) { cout << T2::x.line << '\t' << mng << endl; }
void anno_one::fxy(int64_t v) { 
	int64_t dummy;
	istringstream l(T2::y.line);
	l >> dummy;
	string s;
	l >> s;
	cout << T2::x.line << '\t' << s << endl; 
}

void anno_many::print_anno(int64_t v) { 
	if (v == last_x) cout << T2::y.line << endl; 
}
void anno_many::fx(int64_t v) { last_x = v; }
void anno_many::fy(int64_t v) { print_anno(v); }
void anno_many::fxy(int64_t v) { last_x = v; print_anno(v); }

f_ints::f_ints(
	ids_file<id_prep_cin>& x, 
	ids_file<id_file>& y,
	bool in_intersect,
	bool in_anno
) :
	y_ints<T>::y_ints(x, y),
	intersect(in_intersect),
	anno(in_anno)
{}

void f_ints::anyx() {
	if (anno) { cout << T::x.line << '\t' << (in_y == intersect ? 1 : 0) << endl; }
	else if (in_y == intersect) { cout << T::x.line << endl; }
}

size_t f_ints::yield() { return 0; }

