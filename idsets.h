#ifndef IDSETS_H
#define IDSETS_H

#include <variant>
#include <stdexcept>
#include <memory>
#include <vector>
#include <array>
#include <limits>

using namespace std;

template <class Ti1, class Ti2, class To>
struct tva {
	typedef Ti1 Tinput1;
	typedef Ti2 Tinput2;
	typedef To Toutput;
	Tinput1& x;
	Tinput2& y;
	int64_t last_x;
	int64_t last_y;
	tva(Tinput1& in_x, Tinput2& in_y) : 
		x(in_x), 
		y(in_y),
		last_x(numeric_limits<int64_t>::min()),
		last_y(numeric_limits<int64_t>::min())
	{}
	virtual void fx(int64_t) { };
	virtual void fy(int64_t) { };
	virtual void fxy(int64_t) { };
	virtual Toutput yield() = 0;
	void reset() {
		x.reset();
		y.reset();
	}
	virtual bool finished() {
		return x.finished() && y.finished();
	}
	int64_t get() {
		return min<int64_t>(x.get(), y.get());
	}
	void push() {
		bool xf = x.finished();
		bool yf = y.finished();
		if (xf && yf)
			throw invalid_argument("can't push when x and y are both finished!");
		int64_t xv = xf ? last_x : x.get();
		int64_t yv = yf ? last_y : y.get();
		if ((last_x > xv) || (last_y > yv)) {
			throw invalid_argument("input not sorted!");
		}
		last_x = xv;
		last_y = yv;
		switch (xf ? 1 : (yf ? -1 : ((0LL < (xv-yv)) - ((xv-yv) < 0LL)))) {
			case 0:
				fxy(xv);
				x.push();
				y.push();
				break;
			case -1:
				fx(xv);
				x.push();
				break;
			case 1:
				fy(yv);
				y.push();
				break;
		}
	}
	void terminate() {
		x.terminate();
		y.terminate();
	}
	void process() { 
		reset();
		while (!finished()) {
			push();
		}
		terminate();
	}
	virtual ~tva() = default;
};

template <class Ti1, class Ti2, class To>
struct tva_int : tva<Ti1, Ti2, To> {
	bool on_x;
	bool on_y;
	bool is_ok;
	int64_t start;
	tva_int(typename tva<Ti1, Ti2, To>::Tinput1& in_x, typename tva<Ti1, Ti2, To>::Tinput2& in_y) : 
		tva<Ti1, Ti2, To>(in_x, in_y),
		on_x(false),
		on_y(false),
		is_ok(false),
		start(0)
	{ }
	virtual bool ok() = 0; 
	virtual void do_reg(int64_t, int64_t) = 0;
	void make_reg(int64_t v) {
		bool now_ok = ok();
		if (is_ok != now_ok) {
			if (now_ok) {
				start = v;
			}
			else if (v > start) do_reg(start, v);
		}
		is_ok = ok();	
	}
	void fx(int64_t v) override {
		on_x = !on_x;
		make_reg(v);
	};
	void fy(int64_t v) override {
		on_y = !on_y;
		make_reg(v);
	};
	void fxy(int64_t v) override {
		on_x = !on_x;
		on_y = !on_y;
		make_reg(v);
	};
	virtual To yield() = 0;
};

template <class T>
typename T::Toutput iio(typename T::Tinput1& x, typename T::Tinput2& y) {
	T t(x, y);
	t.process();
	return t.yield();
}

template <class T, class Targ>
typename T::Toutput iio2(typename T::Tinput1& x, typename T::Tinput2& y, Targ X) {
	T t(x, y, X);
	t.process();
	return t.yield();
}

enum set_op_type { SETDIFF = 0, INTERSECT = 1, UNION = 2 };

template <class Tin1, class Tin2, class Tout>
struct set_op_int : tva_int<Tin1, Tin2, Tout> {
	set_op_type so;
	set_op_int(Tin1& x, Tin2& y, int soi) : 
		tva_int<Tin1, Tin2, Tout>(x, y), 
		so((set_op_type)soi)
	{ }
	bool ok() override {
		switch (so) {
			case SETDIFF:
				return tva_int<Tin1, Tin2, Tout>::on_x && !tva_int<Tin1, Tin2, Tout>::on_y;
			case INTERSECT:
				return tva_int<Tin1, Tin2, Tout>::on_x && tva_int<Tin1, Tin2, Tout>::on_y;
			case UNION:
				return tva_int<Tin1, Tin2, Tout>::on_x || tva_int<Tin1, Tin2, Tout>::on_y;
			default:
				throw std::invalid_argument("unrecognised operation type!");
		}
	};
}; 

template <class Tid>
struct counted_int : set_op_int<Tid, Tid, size_t> {
	size_t n;
	counted_int(typename tva_int<Tid, Tid, size_t>::Tinput1& x, typename tva_int<Tid, Tid, size_t>::Tinput2& y, int soi) : 
		set_op_int<Tid, Tid, size_t>(x, y, soi), 
		n(0)
	{}
	void do_reg(int64_t s, int64_t e) override { n++; };
	typename set_op_int<Tid, Tid, size_t>::Toutput yield() override { return n; }
};

template <class Tid, class Trec>
struct get_ints : set_op_int<Tid, Tid, typename Trec::record_type>, Trec {
	get_ints(Tid& x, Tid& y, int soi) : 
		set_op_int<Tid, Tid, typename Trec::record_type>(x, y, soi),
		Trec(2*iio2<counted_int<Tid>, int>(x, y, soi))
	{}
	void do_reg(int64_t s, int64_t e) override { 
		Trec::set_record_value(s);
		Trec::set_record_value(e);
       	};
	typename Trec::record_type yield() override { return Trec::record; }
};

template <class Tfuncs>
typename Tfuncs::Tbase::Toutput tiio(
	set_op_type op_type, 
	typename Tfuncs::Tbase::Tinput1& x, 
	typename Tfuncs::Tbase::Tinput2& y
) {
	unique_ptr<typename Tfuncs::Tbase> t;
	switch (op_type) {
		case SETDIFF:
			t = make_unique<typename Tfuncs::Tsetdiff>(x, y);
			break;
		case INTERSECT:
			t = make_unique<typename Tfuncs::Tintersect>(x, y);
			break;
		case UNION:
			t = make_unique<typename Tfuncs::Tunion>(x, y);
			break;
	}
	t->process();
	return t->yield();
}

template <class Tid>
struct counted : tva<Tid, Tid, size_t> {
	size_t n;
	counted(
		typename tva<Tid, Tid, size_t>::Tinput1& x, 
		typename tva<Tid, Tid, size_t>::Tinput2& y
	) : 
		tva<Tid, Tid, size_t>(x, y), 
		n(0) 
	{ }
	size_t yield() override { return n; }
}; 

template <class Tid>
struct cexcl : counted<Tid> {
	using counted<Tid>::counted;
	void fx(int64_t v) override { counted<Tid>::n++; }
};

template <class Tid>
struct cboth : counted<Tid> {
	using counted<Tid>::counted;
	void fxy(int64_t v) override { counted<Tid>::n++; }
};

template <class Tid>
struct ceither : counted<Tid> {
	using counted<Tid>::counted;
	void fx(int64_t v) override { counted<Tid>::n++; }
	void fy(int64_t v) override { counted<Tid>::n++; }
	void fxy(int64_t v) override { counted<Tid>::n++; }
};

template <class Tid>
struct lexcl : tva<Tid, Tid, vector<bool> > {
	vector<bool> x_keep;
	lexcl(
		Tid& x, 
		Tid& y
	) :
		tva<Tid, Tid, vector<bool> >::tva(x, y),
		x_keep(x.size(), false)
	{}
	void fx(int64_t v) override {
		x_keep[tva<Tid, Tid, vector<bool> >::x.cur] = true;
	}
	vector<bool> yield() override { 
		return x_keep; 
	}
};

template <class Tpts, class Tints>
struct subset_points_by_int : tva<Tpts, Tints, vector<bool> > {
	vector<bool> x_keep;
	bool intersect;
	bool in_y;
	subset_points_by_int(
		Tpts& x, 
		Tints& y,
		bool in_intersect
	) :
		tva<Tpts, Tints, vector<bool> >(x, y),
		x_keep(tva<Tpts, Tints, vector<bool> >::x.n, false),
		intersect(in_intersect),
		in_y(false)
	{}
	void anyx() {
		x_keep[tva<Tpts, Tints, vector<bool> >::x.cur] = in_y;
	}
	void anyy() {
		in_y = !in_y;
	}
	void fx(int64_t v) override {
		anyx();
	}
	void fy(int64_t v) override {
		anyy();
	}
	void fxy(int64_t v) override {
		anyy();
		anyx();
	}
	vector<bool> yield() override { 
		return x_keep; 
	}
};

template <class Tid, class Trec>
struct returning_set : tva<Tid, Tid, typename Trec::record_type>, Trec {
	returning_set(Tid& x, Tid& y, size_t n) :
		tva<Tid, Tid, typename Trec::record_type>::tva(x, y),
		Trec(n)
	{}
	typename Trec::record_type yield() override { return Trec::record; }
};

template <class Tid, class Trec>
struct excl : returning_set<Tid, Trec> {
	excl(Tid& x, Tid& y) : returning_set<Tid, Trec>(x, y, iio<cexcl<Tid> >(x, y)) {}
	void fx(int64_t v) override { returning_set<Tid, Trec>::set_record_value(v); }
};

template <class Tid, class Trec>
struct both : returning_set<Tid, Trec> {
	both(Tid& x, Tid& y) : returning_set<Tid, Trec>::returning_set(x, y, iio<cboth<Tid> >(x, y)) {}
	void fxy(int64_t v) override { returning_set<Tid, Trec>::set_record_value(v); }
};

template <class Tid, class Trec>
struct either : returning_set<Tid, Trec> {
	either(Tid& x, Tid& y) : returning_set<Tid, Trec>::returning_set(x, y, iio<ceither<Tid> >(x, y)) {}
	void fx(int64_t v) override { returning_set<Tid, Trec>::set_record_value(v); }
	void fy(int64_t v) override { returning_set<Tid, Trec>::set_record_value(v); }
	void fxy(int64_t v) override { returning_set<Tid, Trec>::set_record_value(v); }
};

template <class Tid>
struct count_funcs {
	typedef counted<Tid> Tbase;
	typedef cexcl<Tid> Tsetdiff;
	typedef cboth<Tid> Tintersect;
	typedef ceither<Tid> Tunion;
};

template <class Tid, class Trec>
struct set_funcs {
	typedef returning_set<Tid, Trec> Tbase;
	typedef excl<Tid, Trec> Tsetdiff;
	typedef both<Tid, Trec> Tintersect;
	typedef either<Tid, Trec> Tunion;
};


template <class Tid>
struct id_tab : Tid {
	vector<size_t>& count;
	id_tab(typename Tid::Tini& in_v, vector<size_t>& in_c) :
		Tid(in_v),
		count(in_c)
	{}
};

template <class Tid, class Trec>
struct id_tab_base : tva<id_tab<Tid>, id_tab<Tid>, pair<typename Trec::record_type, vector<size_t> > >, Trec {
	vector<size_t> tab;
	id_tab_base(
		typename tva<id_tab<Tid>, id_tab<Tid>, pair<typename Trec::record_type, vector<size_t> > >::Tinput1& x, 
		typename tva<id_tab<Tid>, id_tab<Tid>, pair<typename Trec::record_type, vector<size_t> > >::Tinput2& y
	) :
		tva<id_tab<Tid>, id_tab<Tid>, pair<typename Trec::record_type, vector<size_t> > >::tva(x, y),
		Trec(iio<ceither<Tid> >(x, y)),
		tab(Trec::n, 0)
	{ }
	void set_val(int64_t v, size_t tab_val) {
		tab[Trec::cur] = tab_val;
		Trec::set_record_value(v);
	}
	typename tva<id_tab<Tid>, id_tab<Tid>, pair<typename Trec::record_type, vector<size_t> > >::Toutput yield() override {
		typename tva<id_tab<Tid>, id_tab<Tid>, pair<typename Trec::record_type, vector<size_t> > >::Toutput r(Trec::record, tab);
		return r; 
	}
};

template <class Tid, class Trec>
struct id_tabr : id_tab_base<Tid, Trec> {
	using id_tab_base<Tid, Trec>::id_tab_base;
	void fx(int64_t v) override {
		id_tab_base<Tid, Trec>::set_val(v, id_tab_base<Tid, Trec>::x.count[id_tab_base<Tid, Trec>::x.cur]);
	}
	void fy(int64_t v) override {
		id_tab_base<Tid, Trec>::set_val(v, id_tab_base<Tid, Trec>::y.count[id_tab_base<Tid, Trec>::y.cur]);
	}
	void fxy(int64_t v) override {
		id_tab_base<Tid, Trec>::set_val(v, id_tab_base<Tid, Trec>::x.count[id_tab_base<Tid, Trec>::x.cur] + id_tab_base<Tid, Trec>::y.count[id_tab_base<Tid, Trec>::y.cur]);
	}
};

template <class Tid, class Trec, bool Tadd>
struct cum_id_tabr : id_tab_base<Tid, Trec> {
	int64_t curval;
	bool add;
	cum_id_tabr(
		typename id_tab_base<Tid, Trec>::Tinput1& x, 
		typename id_tab_base<Tid, Trec>::Tinput2& y
	) :
		id_tab_base<Tid, Trec>::id_tab_base(x, y),
		curval(0),
		add(Tadd)
	{}
	void fx(int64_t v) override {
		curval += id_tab_base<Tid, Trec>::x.count[id_tab_base<Tid, Trec>::x.cur];
		id_tab_base<Tid, Trec>::set_val(v, curval);
	}
	void fy(int64_t v) override {
		int64_t d = id_tab_base<Tid, Trec>::y.count[id_tab_base<Tid, Trec>::y.cur];
		if (add)
			curval += d;
		else {
			if (d > curval)
				throw invalid_argument("can't have depth < 0!");
			curval -= d;
		}
		id_tab_base<Tid, Trec>::set_val(v, curval);
	}
	void fxy(int64_t v) override {
		curval += id_tab_base<Tid, Trec>::x.count[id_tab_base<Tid, Trec>::x.cur];
		int64_t d = id_tab_base<Tid, Trec>::y.count[id_tab_base<Tid, Trec>::y.cur];
		if (add)
			curval += d;
		else {
			if (d > curval)
				throw invalid_argument("can't have depth < 0!");
			curval -= d;
		}
		id_tab_base<Tid, Trec>::set_val(v, curval);
	}
};

template <class T>
auto cut_cum_tab(T& values, vector<size_t>& height, size_t threshold) {
	int peaks = 0;
	int n = values.size();
	bool in_peak = false;
	for (int i = 0; i < n; i++) {
		if (height[i] >= threshold) {
			if (!in_peak) {
				in_peak = true;
				peaks++;
			}
		} else {
			in_peak = false;
		}
	}
	if (in_peak) throw invalid_argument("must end below threshold!");
	pair<T, T> result = pair<T, T>(T(peaks), T(peaks));
	int cur = 0;
	for (int i = 0; i < n; i++) {
		if (height[i] >= threshold) {
			if (!in_peak) {
				in_peak = true;
				result.first[cur] = values[i];
			}
		} else {
			if (in_peak) result.second[cur++] = values[i];
			in_peak = false;
		}
	}
	return result;
}

#endif
