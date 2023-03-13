#ifndef IDTYPES_H
#define IDTYPES_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector> 
#include <stdexcept>
#include "idsets.h"

using namespace std;

template <class T>
struct rsvr {
	typedef vector<T> Tini;
	Tini& v;
	size_t n;
	size_t cur;
	int move_bits;
	rsvr(vector<T>& in_v, int in_move_bits) : v(in_v), n(v.size()), cur(0), move_bits(in_move_bits) {}
	rsvr(vector<T>& v) : rsvr(v, 0) {}
	int64_t operator[](size_t i) const {
		return (int64_t)v[i];
	}
	size_t size() { return n; }
	bool finished() { return cur >= n; }
	void push() {
		cur++;
	}
	void reset() {
		cur = 0;
	}
	void terminate() {}
	T get_el() {
		return v[cur];
	}
	int64_t get() {
		return (int64_t)v[cur] >> move_bits;
	}
};

struct rsvr_record {
	typedef vector<int64_t> record_type;
	size_t n;
	record_type record;
	size_t cur;
	rsvr_record(size_t);
	void set_record_value(int64_t);
};

struct id_file {
	istream& f;
	id_file(istream&);
	bool load(int64_t&);
};

struct id_prep_cin {
	string line;
	istream& f;
	id_prep_cin(istream&);
	id_prep_cin();
	bool load(int64_t&);
};

template <class F>
struct ids_file : F {
	bool fin;
	int64_t val;
	int move_bits;
	bool done;
	template <class ... ctypes>
	ids_file(
		int in_move_bits,
		ctypes&& ... cargs
	) : 
		F(cargs...),
		fin(false), 
		val(0),
		move_bits(in_move_bits),
		done(false)
	{ 
		fin = F::load(val); 
	}
	bool finished() { return fin; }
	void push() {
		if (fin) throw invalid_argument("already finished!");
		fin = F::load(val);
	}
	void reset() {
		if (done)
			throw invalid_argument("cannot re-read the file!");
	}
	void terminate() { done = true; }
	int64_t get() {
		return val >> move_bits;
	}
};

typedef tva<ids_file<id_prep_cin>, ids_file<id_file>, size_t> T;
typedef tva<ids_file<id_prep_cin>, ids_file<id_prep_cin>, size_t> T2;

struct f_ids : T {
	using T::tva;
	void stamp();
	size_t yield() override;
};

struct exclude_ids : f_ids {
	using f_ids::f_ids;
	void fx(int64_t) override;
	
};

struct intersect_ids : f_ids {
	using f_ids::f_ids;
	void fxy(int64_t) override;
	
};

struct anno_with : T2 {
	using T2::tva;
	size_t yield() override;
	
};

struct anno_one : anno_with {
	string& mng;
	anno_one(
		ids_file<id_prep_cin>&, 
		ids_file<id_prep_cin>&,
		string&
	);
	void fx(int64_t) override;
	void fxy(int64_t) override;
};

struct anno_many : anno_with {
	using anno_with::anno_with;
	int64_t last_x;
	void print_anno(int64_t);
	void fx(int64_t) override;
	void fy(int64_t) override;
	void fxy(int64_t) override;
};

enum class anno_type { anno, intersect, exclude };

struct anno_ids : f_ids {
	bool seen_y;
	anno_type type;
	int64_t last_y;
	void anno(char);
	void fx(int64_t) override;
	void fxy(int64_t) override;
	
	anno_ids(
		ids_file<id_prep_cin>&, 
		ids_file<id_file>&,
		anno_type
	);
};

struct merge_lines : T2 {
	void fx(int64_t) override;
	void fxy(int64_t) override;
	void fy(int64_t) override;
	size_t yield() override;
	using T2::tva;
};

template <class Ttva>
struct y_ints : Ttva {
	bool in_y;
	y_ints(
		typename Ttva::Tinput1& x, 
		typename Ttva::Tinput2& y
	) :
		Ttva(x, y),
		in_y(false)
	{}
	virtual void anyx() = 0;
	void anyy() { in_y = !in_y; }
	void fx(int64_t v) override { anyx(); }
	void fy(int64_t v) override { anyy(); }
	void fxy(int64_t v) override { anyy(); anyx(); }
	virtual ~y_ints() = default;
};

struct f_ints : y_ints<T> {
	bool intersect;
	bool anno;
	f_ints(
		ids_file<id_prep_cin>& x, 
		ids_file<id_file>& y,
		bool in_intersect,
		bool in_anno
	);
	void anyx() override;
	size_t yield() override;
};

template <class T>
struct ss_vec_file_pts : tva<rsvr<T>, ids_file<id_file>, vector<T> > {
	bool exclude;
	vector<T> result;
	ss_vec_file_pts(
		rsvr<T>& x, 
		ids_file<id_file>& y,
		bool in_exclude
	) : 
		tva<rsvr<T>, ids_file<id_file>, vector<T> >(x, y),
		exclude(in_exclude) 
	{}
	void fx(int64_t v) override {
		if (exclude) result.push_back(tva<rsvr<T>, ids_file<id_file>, vector<T> >::x.get_el());
	}
	void fxy(int64_t v) override {
		if (!exclude) result.push_back(tva<rsvr<T>, ids_file<id_file>, vector<T> >::x.get_el());
	}
	vector<T> yield() override { return result; }
};

template <class T>
struct ss_vec_file_ints : y_ints<tva<rsvr<T>, ids_file<id_file>, vector<T> > > {
	bool exclude;
	vector<T> result;
	ss_vec_file_ints(
		rsvr<T>& x, 
		ids_file<id_file>& y,
		bool in_exclude
	) :
		y_ints<tva<rsvr<T>, ids_file<id_file>, vector<T> > >(x, y),
		exclude(in_exclude)
	{}
	void anyx() override {
		if (exclude == !y_ints<tva<rsvr<T>, ids_file<id_file>, vector<T> > >::in_y) result.push_back(y_ints<tva<rsvr<T>, ids_file<id_file>, vector<T> > >::x.get_el());
	}
	vector<T> yield() override { return result; }
};

#endif
