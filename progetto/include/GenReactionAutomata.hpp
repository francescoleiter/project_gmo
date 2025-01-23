#include"Reactions.hpp"
#include"ReactionAutomata.hpp"
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
#ifndef GENETIC_REACTION_AUTOMATA
#define GENETIC_REACTION_AUTOMATA


template<typename N, template<typename, std::size_t> class R, std::size_t SIZE>
struct GeneSize;

template<typename N, std::size_t SIZE>
struct GeneSize<N,CReaction,SIZE>{
    static constexpr std::size_t size = SIZE*2;
};
template<typename N, std::size_t SIZE>
struct GeneSize<N,Reaction,SIZE>{
    static constexpr std::size_t size = SIZE*3;
};

template<typename N, std::size_t SIZE>
CReaction<N,SIZE> reaction_from_gene(const std::array<N, SIZE*2> & gene);

template<typename N, std::size_t SIZE>
std::array<N, GeneSize<N,CReaction,SIZE>::size>  gene_from_reaction(const CReaction<N,SIZE> &reaction);

template<typename N, std::size_t SIZE>
Reaction<N,SIZE> reaction_from_gene(const std::array<N, SIZE*3> & genome);

template<typename N, std::size_t SIZE>
std::array<N, GeneSize<N,Reaction,SIZE>::size> gene_from_reaction(const Reaction<N,SIZE> & reaction);



template<typename N, template<typename, std::size_t> class R, std::size_t SIZE>
class FAGenome : public std::vector<N>{
public:
    FAGenome() = default;
    FAGenome(const std::vector<N> & genome) : std::vector<N>(genome){
      if(((genome.size() - SIZE)%GeneSize<N,R,SIZE>::size)){
          throw std::invalid_argument("genome size incompatible");
      }
    }
    const N & operator[](std::size_t i) const{return std::vector<N>::operator[](i);}
    N & operator[](std::size_t i) {return std::vector<N>::operator[](i);}

    std::size_t reactions_size() const{return (this->size() - SIZE)/GeneSize<N,R,SIZE>::size;};
    void mutate_reaction(std::size_t idx);
    std::array<N,GeneSize<N,R,SIZE>::size> * access_reaction_gene(std::size_t i) {return &(this->operator[](i*GeneSize<N,R,SIZE>::size+SIZE));};
    void print() const;
};

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
FAGenome<N,R,SIZE> genome_from_automaton(const FuncAutomaton<N, R, SIZE, Pure> & aut);

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class FAUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
FAUT<N,R,SIZE,Pure> funcautomaton_from_genome(const FAGenome<N,R,SIZE> & genome);


template<typename N, template<typename, std::size_t> class R, std::size_t SIZE>
class AAGenome : public std::vector<N>{
public:
    AAGenome() = default;
    AAGenome(const std::vector<N> & genome, const std::size_t & _n_reactions) : std::vector<N>(genome), n_reactions(_n_reactions){
      if((genome.size() - SIZE - GeneSize<N,R,SIZE>::size*(n_reactions))%SIZE){
          throw std::invalid_argument("genome size incompatible");
      }
    }
    const N & operator[](std::size_t i) const{return std::vector<N>::operator[](i);}
    N & operator[](std::size_t i) {return std::vector<N>::operator[](i);};

    std::array<N,GeneSize<N,R,SIZE>::size> * access_reaction_gene(std::size_t i) {return &(this->operator[](i*GeneSize<N,R,SIZE>::size+SIZE));};
    const std::size_t & reactions_size() const{return n_reactions;};

    void print() const;
protected:
  std::size_t n_reactions;
};

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
AAGenome<N,R,SIZE> genome_from_automaton(const AccAutomaton<N, R, SIZE, Pure> & aut);

template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AAUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
AAUT<N,R,SIZE,Pure> accautomaton_from_genome(const AAGenome<N,R,SIZE> & genome);

template<typename N,  std::size_t SIZE>
class TestFunction{
public:
  TestFunction() = default;
  TestFunction(const std::function<N(std::vector<N>)> & function) : func(function) {}
  N operator()(std::vector<N> input) const {return func(input);}
  std::function<N(std::vector<N>)> func;
};

template<typename N, std::size_t SIZE, bool CHEMICAL>
void mutate_sequence(N * gene_ptr, std::mt19937 random_gen, double m_prob, double rp_prob ,double i_prob){
    std::geometric_distribution ms_distribution(m_prob);
    std::discrete_distribution i_distribution{1.0-i_prob,i_prob}, rp_distribution{1.0-rp_prob,rp_prob};
    if constexpr(CHEMICAL){
        for(std::size_t i=0;i<SIZE;++i){
            gene_ptr[i] = (1+ms_distribution(random_gen))*rp_distribution(random_gen);
            gene_ptr[i+SIZE] = (1+ms_distribution(random_gen))*rp_distribution(random_gen);
        }
    }else{
        for(std::size_t i=0;i<SIZE;++i){
            gene_ptr[i] = (1 + ms_distribution(random_gen))*rp_distribution(random_gen);
            gene_ptr[i+SIZE] = i_distribution(random_gen) && (gene_ptr[i]==0);
            gene_ptr[i+SIZE*2] = (1 + ms_distribution(random_gen))*rp_distribution(random_gen);
        }
    }
    return;
}
template<typename N, std::size_t SIZE>
void mutate_short_sequence(N * gene_ptr, std::mt19937 random_gen, double m_prob, double rp_prob){
    std::geometric_distribution ms_distribution(m_prob);
    std::discrete_distribution rp_distribution{1.0-rp_prob,rp_prob};
    for(std::size_t i=0;i<SIZE;++i){
        gene_ptr[i] = (1+ms_distribution(random_gen))*rp_distribution(random_gen);
    }
    return;
}


template<typename N, std::size_t SIZE>
void mutate_reaction(CReaction<N,SIZE> & reaction, std::mt19937 random_gen, double m_prob, double rp_prob ,double i_prob){
    std::geometric_distribution ms_distribution(m_prob);
    std::discrete_distribution rp_distribution{1.0-rp_prob,rp_prob};
    std::array<N,SIZE*2> reaction_gene;
    for(std::size_t i = 0; i < SIZE; ++i){
        reaction_gene[i] = ms_distribution(random_gen)*rp_distribution(random_gen);
        reaction_gene[i+SIZE] = ms_distribution(random_gen)*rp_distribution(random_gen);
    }
    reaction = reaction_from_gene(reaction_gene);
    return;
}


template<template<typename, template<typename, std::size_t> class, std::size_t> class GEN, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Acc>
void mutation(GEN<N,R,SIZE> & genome, const std::array<double,3> & probs, double m_prob, double rp_prob, double i_prob){
    std::random_device rd{};
    std::array<std::uint32_t, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(rd));
    std::seed_seq seq(seed_data.begin(), seed_data.end());
    std::mt19937 random_gen(seq);
    std::discrete_distribution<std::size_t> mutation_type(probs.begin(), probs.end());
    std::uniform_int_distribution<std::size_t> reaction_idx(0, genome.reactions_size()-1);
    std::uniform_int_distribution<std::size_t> _idx(1,(genome.size() - SIZE - genome.reactions_size()*GeneSize<N,R,SIZE>::size)/SIZE);
    switch(mutation_type(random_gen)){
        case 0:
            mutate_short_sequence<N,SIZE>(&genome[0], random_gen, m_prob, rp_prob);
            break;
        case 1:
            mutate_sequence<N,SIZE,!(GeneSize<N,R,SIZE>::size-2*SIZE)>(&genome[SIZE + reaction_idx(random_gen)*GeneSize<N,R,SIZE>::size], random_gen, m_prob, rp_prob,i_prob);
            break;
        case 2:
            if constexpr(Acc){
                std::size_t idx = _idx(random_gen);
                mutate_short_sequence<N,SIZE>(&genome[genome.size()-idx*SIZE], random_gen, m_prob, rp_prob);
            }
    return;
    }
}
template<template<typename, template<typename, std::size_t> class, std::size_t> class GEN, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Acc>
void mutation_vt(std::vector<N> & genome_vt, const std::array<double,3> & probs, double m_prob, double rp_prob, double i_prob){
    GEN<N,R,SIZE> * genome;
    genome = static_cast<GEN<N,R,SIZE> *>(&genome_vt);
    mutation<GEN,N,R,SIZE,Acc>(*genome, probs, m_prob, i_prob);
}
template<template<typename, template<typename, std::size_t> class, std::size_t> class GEN, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Acc>
void mutation_raw(N * ptr, std::size_t len, std::size_t n_reactions, const std::array<double,3>& probs, double m_prob, double rp_prob, double i_prob){
    std::random_device rd{};
    std::array<std::uint32_t, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(rd));
    std::seed_seq seq(seed_data.begin(), seed_data.end());
    std::mt19937 random_gen(seq);
    std::discrete_distribution<std::size_t> mutation_type(probs.begin(), probs.end());
    std::uniform_int_distribution<std::size_t> reaction_idx(0, n_reactions-1);
    std::uniform_int_distribution<std::size_t> _idx(1, len/SIZE - n_reactions*GeneSize<N,R,SIZE>::size/SIZE - 1);
    switch(mutation_type(random_gen)){
        case 0:
            mutate_short_sequence<N,SIZE>(ptr, random_gen, m_prob, rp_prob);
            break;
        case 1:
            mutate_sequence<N,SIZE,!(GeneSize<N,R,SIZE>::size-2*SIZE)>(&ptr[SIZE + reaction_idx(random_gen)*GeneSize<N,R,SIZE>::size], random_gen, m_prob, rp_prob, i_prob);
            break;
        case 2:
            if constexpr(Acc){
                std::size_t idx = _idx(random_gen);
                mutate_short_sequence<N,SIZE>(&ptr[len-idx*SIZE], random_gen, m_prob, rp_prob);
            }
    return;
    }
}

template<template<typename, template<typename, std::size_t> class, std::size_t> class GEN, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Acc>
void mutation_raw_vctr(N * ptr, std::size_t len, std::size_t population, std::size_t n_reactions, double mutation_prob, const std::array<double,3>& probs, double m_prob, double rp_prob,double i_prob){
    std::random_device rd{};
    std::array<std::uint32_t, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(rd));
    std::seed_seq seq(seed_data.begin(), seed_data.end());
    std::mt19937 random_gen(seq);
    std::discrete_distribution<std::size_t> mutation_type(probs.begin(), probs.end()), is_mutation({1.0-mutation_prob, mutation_prob});
    std::uniform_int_distribution<std::size_t> reaction_idx(0, n_reactions-1);
    std::uniform_int_distribution<std::size_t> _idx(1, len/SIZE - n_reactions*GeneSize<N,R,SIZE>::size/SIZE);
    int idx_r;
    for (int i = 0; i < population; i++){
        if(is_mutation(random_gen)){
            switch(mutation_type(random_gen)){
            case 0:
                mutate_short_sequence<N,SIZE>(&ptr[i*len], random_gen, m_prob, rp_prob);
            case 1:
                idx_r = reaction_idx(random_gen)*GeneSize<N,R,SIZE>::size;
                mutate_sequence<N,SIZE,!(GeneSize<N,R,SIZE>::size-2*SIZE)>(&ptr[i*len + SIZE + idx_r ], random_gen, m_prob, rp_prob, i_prob);
            case 2:
                if constexpr(Acc){
                    std::size_t idx = _idx(random_gen);
                    mutate_short_sequence<N,SIZE>(&ptr[i*len + len -idx*SIZE], random_gen ,m_prob, rp_prob);
            }
        }
    }
    }
}


#include "GenReactionAutomata.tpl.hpp"

#endif