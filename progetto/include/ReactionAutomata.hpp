#include"Reactions.hpp"
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
#ifndef REACTION_AUTOMATA
#define REACTION_AUTOMATA

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
class FuncAutomaton{
public:
  	FuncAutomaton() = default;
    virtual ~FuncAutomaton() = default;
    FuncAutomaton(const Multiset<N, SIZE> & _initial_state, const std::vector<R<N, SIZE>> & _reactions)
        : state(_initial_state), initial_state(_initial_state), reactions(_reactions) {}

    const Multiset<N, SIZE> & get_state() const {return state;}

    virtual const N & result_indicator() const {return this->state[0];}

    const bool & is_halted() const {return this->halted;}

    virtual bool deterministic() const = 0;
	virtual void react();
    virtual void verbose_react();
    virtual void reset(){this->state = this->initial_state; this->halted = false; this->time = 0;}
    void read_input(const Multiset<N, SIZE> & input){this->state += input;};
    void read_single_input(const N & input){this->state[input] += 1;};

//protected:
    virtual std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> verbose_get_products__residue() const = 0;
    virtual std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> get_products__residue() const = 0;
    std::vector<R<N, SIZE>> reactions;
    Multiset<N, SIZE> initial_state;
    Multiset<N, SIZE> state;
    uint32_t time = 0;
    bool halted = false;
};

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
class AccAutomaton : public FuncAutomaton<N, R, SIZE, Pure>{
public:
  	AccAutomaton() = default;
    ~AccAutomaton() = default;
    AccAutomaton(const Multiset<N, SIZE> & _initial_state, const std::vector<R<N, SIZE>> & _reactions, const std::vector<Multiset<N, SIZE>> & _accepting_states)
        : FuncAutomaton<N,R,SIZE,Pure>::FuncAutomaton(_initial_state, _reactions), accepting_states(_accepting_states) {}
    bool result() const{return this->accepted;};
    void reset() override {this->accepted = false; FuncAutomaton<N,R,SIZE,Pure>::reset();};
	void react() override;
    void verbose_react() override;
//protected:
    bool accepted;
    std::vector<Multiset<N,SIZE>> accepting_states;
};

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
class Automaton_sq : public AUT<N,R,SIZE,Pure>{
public:
    using AUT<N, R, SIZE, Pure>::AUT;
    bool deterministic() const override;
    //protected:
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> get_products__residue() const override;
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> verbose_get_products__residue() const override;
};

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
using FuncAutomaton_sq = Automaton_sq<FuncAutomaton,N,R,SIZE,Pure>;
template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
using AccAutomaton_sq = Automaton_sq<AccAutomaton,N,R,SIZE,Pure>;


template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
class Automaton_mp : public AUT<N,R,SIZE,Pure>{
public:
    using AUT<N, R, SIZE, Pure>::AUT;
    bool deterministic() const override;
    //protected:
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> get_products__residue() const override;
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> verbose_get_products__residue() const override;
};

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
using FuncAutomaton_mp = Automaton_mp<FuncAutomaton,N,R,SIZE,Pure>;
template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
using AccAutomaton_mp = Automaton_mp<AccAutomaton,N,R,SIZE,Pure>;

#include "ReactionAutomata.tpl.hpp"

#endif