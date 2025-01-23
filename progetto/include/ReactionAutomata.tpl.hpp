#include"ReactionAutomata.hpp"
#include<vector>
#include<array>
#include<cstdint>
#include<algorithm>
#include<exception>
#include<functional>
#include<random>
#include<stdexcept>
#include<iostream>
#ifndef REACTION_AUTOMATA_TEMPL
#define REACTION_AUTOMATA_TEMPL


#include"ReactionAutomata.hpp"
#include<vector>
#include<array>
#include<cstdint>
#include<algorithm>
#include<exception>
#include<functional>
#include<random>
#include<stdexcept>



template<typename N, std::size_t SIZE>
bool Multiset<N, SIZE>::operator<= (const Multiset<N, SIZE> & other) const{
    for(std::size_t i = 0; i < SIZE; i++){
      if(this->operator[](i) > other[i]){
        return false;
      }
    }
    return true;
};

template<typename N, std::size_t SIZE>
bool Multiset<N, SIZE>::operator< (const Multiset<N, SIZE> & other) const{
    for(std::size_t i = 0; i < SIZE; i++){
        if(this[i] > other[i]){
            return false;
        }
    }
    return true;
}


template<typename N, std::size_t SIZE>
bool Multiset<N, SIZE>::operator== (const Multiset<N, SIZE> & other) const{
    for(std::size_t i = 0; i < SIZE; i++){
        if(this->operator[](i) != other[i]){
            return false;
        }
    }
    return true;
}

template<typename N, std::size_t SIZE>
Multiset<N, SIZE> Multiset<N, SIZE>::operator- (const Multiset<N, SIZE> & other) const{
    if (!(other <= *this)){
        throw std::runtime_error("Negative result in subtraction is not allowed.");
    }
    return Multiset<N, SIZE>(multi_operation<N, SIZE, std::minus<N>>(*this, other));
}

template<typename N, std::size_t SIZE>
bool Multiset<N, SIZE>::empty_intersection(const Multiset<bool, SIZE> & inhibitor) const{
    for(std::size_t i = 0; i < SIZE; i++){
      if((this->operator[](i)) && inhibitor[i]){
        return false;
      }
    }
    return true;
}

template<typename N, std::size_t SIZE>
void Multiset<N, SIZE>::print() const{
  for(std::size_t i = 0; i < SIZE; i++){
    std::cout << unsigned(this->operator[](i)) << " ";
  }
  std::cout << std::endl;
}

template<typename N, std::size_t SIZE>
N CReaction<N, SIZE>::enabled_n(const Multiset<N, SIZE> & state) const{
    if(!(this->r <= state) || (this->r == Multiset<N, SIZE>())){
        return 0;
    }
    N coeff = std::numeric_limits<N>::max();
    for(std::size_t i = 0; i < SIZE; i++){
      if(this->r[i] != 0){
        coeff = std::min<N>(coeff, state[i]/this->r[i]);
      }
    }
    return coeff;
}


template<typename N, std::size_t SIZE>
N Reaction<N, SIZE>::enabled_n(const Multiset<N, SIZE> & state) const{
    if(!(this->r <= state) || !(state.empty_intersection(this->i)) || (this->r == Multiset<N, SIZE>())){
        return 0;
    }
    N coeff = std::numeric_limits<N>::max();
    for(std::size_t i = 0; i < SIZE; i++){
      if(this->r[i] != 0){
        coeff = std::min<N>(coeff, state[i]/this->r[i]);
      }
    }
    return coeff;
}



template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
void FuncAutomaton<N,R,SIZE,Pure>::react(){
  	(this->time)++;
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> prod_res = this->get_products__residue();
    Multiset<N, SIZE> result;
    if constexpr (Pure){
		result = prod_res.first;
        this->halted = result == this->state;
    }else{
        result = prod_res.first + prod_res.second;
        this->halted = result == this->state;
    }
    this->state = result;
}

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
void FuncAutomaton<N,R,SIZE,Pure>::verbose_react(){
  	(this->time)++;
    std::cout << "state:" << std::endl;
    this->state.print();
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> prod_res = this->verbose_get_products__residue();

    Multiset<N, SIZE> result;
    if constexpr (Pure){
		result = prod_res.first;
        this->halted = result == this->state;
    }else{
        result = prod_res.first + prod_res.second;
        this->halted = result == this->state;
    }
    std::cout << "new state:" << std::endl << "--->";
    result.print();
    this->state = result;
}

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
void AccAutomaton<N,R,SIZE,Pure>::react(){
  	(this->time)++;
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> prod_res = this->get_products__residue();
    Multiset<N, SIZE> result;
    if constexpr (Pure){
		result = prod_res.first;
    }else{
        result = prod_res.first + prod_res.second;
    }
    this->halted = (result == this->state) || (std::find(this->accepting_states.cbegin(), this->accepting_states.cend(), result) != this->accepting_states.cend());
    this->state = result;
    this->accepted = std::find(this->accepting_states.cbegin(), this->accepting_states.cend(), this->state) != this->accepting_states.cend();
}

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
void AccAutomaton<N,R,SIZE,Pure>::verbose_react(){
  	(this->time)++;
    std::cout << "state:" << std::endl;
    this->state.print();
    std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> prod_res = this->verbose_get_products__residue();
    Multiset<N, SIZE> result;
    if constexpr (Pure){
		result = prod_res.first;
    }else{
        result = prod_res.first + prod_res.second;
    }
    std::cout << "new state:" << std::endl << "--->";
    result.print();
    this->halted = (result == this->state) || (std::find(this->accepting_states.cbegin(), this->accepting_states.cend(), result) != this->accepting_states.cend());
    this->state = result;
    this->accepted = std::find(this->accepting_states.cbegin(), this->accepting_states.cend(), this->state) != this->accepting_states.cend();
}

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> Automaton_sq<AUT,N,R,SIZE,Pure>::get_products__residue() const {
    for(std::size_t i = 0; i < this->reactions.size(); ++i){
      if(this->reactions[i].enabled(this->state)){
        return std::make_pair(this->reactions[i].get_products(), this->state - this->reactions[i].get_reactants());
      }
    }
    return std::make_pair(Multiset<N,SIZE>{}, this->state);
}

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
bool Automaton_sq<AUT,N,R,SIZE,Pure>::deterministic() const {
    std::vector<Multiset<N,SIZE>> result;
    std::transform(this->reactions.cbegin(), this->reactions.cend(), std::back_inserter(result), [this](const R<N,SIZE> & x){if(x.enabled(this->state)){return x.p;} return Multiset<N,SIZE>{};});
    bool active = false;
    Multiset<N,SIZE> next_state;
    for (Multiset<N,SIZE> x : result) {
        if(x == Multiset<N,SIZE>{}){
          break;
        }
        if(!active){
            active = true;
            next_state = x;
        }
        if(next_state == x){
          break;
        }
        return false;
    }
    return true;
}



template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> Automaton_mp<AUT,N,R,SIZE,Pure>::get_products__residue() const {
    Multiset<N,SIZE> products, residue = this->state;
    for(std::size_t i = 0; i < this->reactions.size(); ++i){
    	N coeff = this->reactions[i].enabled_n(residue);
		products += this->reactions[i].get_products()*coeff;
        residue = residue - this->reactions[i].get_reactants()*coeff;
    }
    return std::make_pair(products, residue);
}

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
std::pair<Multiset<N,SIZE>,Multiset<N,SIZE>> Automaton_mp<AUT,N,R,SIZE,Pure>::verbose_get_products__residue() const {
    Multiset<N,SIZE> products, residue = this->state;
    for(std::size_t i = 0; i < this->reactions.size(); ++i){
    	N coeff = this->reactions[i].enabled_n(residue);
        if(coeff > 0){
            std::cout  << coeff << " * "<< "r" << static_cast<unsigned>(i) << ", ";
        }
		products += this->reactions[i].get_products()*coeff;
        residue = residue - this->reactions[i].get_reactants()*coeff;
    }
    std::cout << std::endl;
    return std::make_pair(products, residue);
}

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
bool Automaton_mp<AUT,N,R,SIZE,Pure>::deterministic() const {
    return true;
    }


#endif