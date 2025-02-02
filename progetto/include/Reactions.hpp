#include<vector>
#include<array>
#include<cstdint>
#include<algorithm>
#include<exception>
#include<functional>
#include<random>
#include<stdexcept>
#include<limits>
#include<cmath>
#include<memory>
#include<string>
#ifndef REACTIONS_REACTION_AUTOMATA
#define REACTIONS_REACTION_AUTOMATA

constexpr char literals[] = "zabcdefghijklmnopqrtsuvwxyZABCDEFGHIJKLMNOPQRSTUVWXY";

template<typename N, std::size_t SIZE, class OP>
std::array<N, SIZE> multi_operation(const std::array<N, SIZE> & m1, const std::array<N, SIZE> & m2){
  std::array<N, SIZE> res;
  for (std::size_t i = 0; i < SIZE; i++){
    res[i] = OP{}(m1[i], m2[i]);
  }
  return res;
  }
template<typename N, std::size_t SIZE, class OP>
void multi_operation_inplace(std::array<N, SIZE> & m1, const std::array<N, SIZE> & m2){
  for (std::size_t i = 0; i < SIZE; i++){
    m1[i] = OP{}(m1[i], m2[i]);
  }
  return;
}


template<typename N, std::size_t SIZE>
std::array<N, SIZE> multi_sum (const std::array<N, SIZE> & m1, const std::array<N, SIZE> & m2){}
template<std::size_t SIZE>
std::array<bool, SIZE> multi_sum (const std::array<bool, SIZE> & m1, const std::array<bool, SIZE> & m2){return multi_operation<bool,SIZE,std::logical_or<bool>>(m1, m2);}

template<typename N, std::size_t SIZE>
void multi_sum_inplace (std::array<N, SIZE> & m1, const std::array<N, SIZE> & m2){multi_operation_inplace<N,SIZE,std::plus<N>>(m1, m2);}
template<std::size_t SIZE>
void multi_sum_inplace (std::array<bool, SIZE> & m1, const std::array<bool, SIZE> & m2){multi_operation_inplace<bool,SIZE,std::logical_or<bool>>(m1, m2);}


template<typename N, std::size_t SIZE>
class Multiset : public std::array<N, SIZE>{
    friend std::array<N, SIZE> multi_sum<N, SIZE>(const std::array<N, SIZE> &, const std::array<N, SIZE> &);
    friend void multi_sum_inplace<N, SIZE>(std::array<N, SIZE> & m1, const std::array<N, SIZE> & m2);
public:
    Multiset() : std::array<N, SIZE>() {for (std::size_t i = 0; i < SIZE; i++) (*this)[i] = 0; };
    Multiset(const std::array<N, SIZE> & in_values) : std::array<N, SIZE>{in_values} {}

    Multiset<N, SIZE> operator+ (const Multiset<N, SIZE> & other) const {return multi_operation<N,SIZE,std::plus<N>>(*this, other);}
    void operator+= (const Multiset<N, SIZE> & other) {multi_sum_inplace<N,SIZE>(*this, other);}
    Multiset<N, SIZE> operator* (const N & coeff) const {Multiset<N, SIZE> result;
                                                          std::transform(this->cbegin(), this->cend(), result.begin(), std::bind(std::multiplies<>{}, std::placeholders::_1, coeff));
                                                          return result;};
    Multiset<N, SIZE> operator- (const Multiset<N, SIZE> & other) const;
    Multiset<bool, SIZE> to_set() const{Multiset<bool, SIZE> res = *this; return res;};
    bool empty_intersection(const Multiset<bool, SIZE> & m) const;
    bool operator< (const Multiset<N, SIZE> & other) const;
    bool operator<= (const Multiset<N, SIZE> & other) const;
    bool operator== (const Multiset<N, SIZE> & other) const;
    bool is_empty() const{ return std::all_of(this->cbegin(), this->cend(), [](N x){return x==0;});}
    void print() const;
    std::string to_string(bool as_chars = false) const;
};
template<std::size_t SIZE>
using Set = Multiset<bool, SIZE>;



template<typename N, std::size_t SIZE>
class  CReaction{
public:
    CReaction() = default;
    CReaction(const Multiset<N, SIZE> & reactants, const Multiset<N, SIZE> & products)
        : r(reactants), p(products) {}

    CReaction<N, SIZE> operator+(const CReaction<N, SIZE> & other) const {return CReaction(r + other.r, p + other.p);}
    CReaction<N, SIZE> operator=(const CReaction<N, SIZE> & other){return *this;}
    const Multiset<N, SIZE> & get_reactants() const {return r;}
    const Multiset<N, SIZE> & get_products() const {return p;}
    void set_reactants(const Multiset<N, SIZE> & reactants) {r = reactants;}
    void set_products(const Multiset<N, SIZE> & products) {p = products;}
    virtual bool enabled(const Multiset<N, SIZE> & state) const {return r <= state;}
    virtual N enabled_n(const Multiset<N, SIZE> & state) const;
    virtual void print() const{this->r.print();this->p.print();};
    virtual std::string to_string(bool as_chars = false) const {return this->r.to_string(as_chars) + '\n' + this->p.to_string(as_chars);};

protected:
    Multiset<N, SIZE> r, p;
};


template<typename N, std::size_t SIZE>
class  Reaction : public CReaction<N, SIZE>{
public:
  	Reaction() = default;
    Reaction(const Multiset<N, SIZE> & reactants, const Set<SIZE> & inhibitors, const Multiset<N, SIZE> & products)
        : CReaction<N, SIZE>::CReaction(reactants, products), i(inhibitors) {}
    Reaction<N, SIZE> operator+(const Reaction<N, SIZE> & other) const {return Reaction<N, SIZE>(this->r + other.r, this->i + other.i, this->p + other.p);}
    const Set<SIZE> & get_inhibitors() const {return this->i;}
    void set_inhibitors(const Set<SIZE> & i_) {this->i = i_;}
    bool enabled(const Multiset<N, SIZE> & state) const override {return (this->r <= state) && state.empty_intersection(i);}
    N enabled_n(const Multiset<N, SIZE> & state) const override;
    void print() const override {this->r.print();this->i.print();this->p.print();};
    std::string to_string(bool as_chars = false) const override {return this->r.to_string(as_chars) + '\n' + this->i.to_string(as_chars) + '\n' + this->p.to_string(as_chars);};
private:
    Set<SIZE> i;
};

#endif