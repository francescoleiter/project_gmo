#include"Reactions.hpp"
#include"ReactionAutomata.hpp"
#include"GenReactionAutomata.hpp"
#include<vector>
#include<array>
#include<cstdint>
#include<algorithm>
#include<exception>
#include<functional>
#include<random>
#include<stdexcept>
#include<iostream>
#ifndef GENETIC_REACTION_AUTOMATA_TEMPL
#define GENETIC_REACTION_AUTOMATA_TEMPL


template<typename N, std::size_t SIZE>
CReaction<N,SIZE> reaction_from_gene(const std::array<N, SIZE*2> & gene){
    std::array<N, SIZE> reactants, products;
    std::copy(gene.cbegin(), gene.cend() - SIZE, reactants.begin());
    std::copy(gene.cbegin() + SIZE, gene.cend(), products.begin());
    return CReaction<N,SIZE>(reactants, products);
}

template<typename N, std::size_t SIZE>
std::array<N,GeneSize<N,CReaction,SIZE>::size> gene_from_reaction(const CReaction<N,SIZE> &reaction){
    std::array<N, SIZE*2> gene;
    auto iter = gene.begin();
    std::copy(reaction.get_reactants().cbegin(), reaction.get_reactants().cend(), iter);
    std::copy(reaction.get_products().cbegin(), reaction.get_products().cend(), iter);
    return gene;
}

template<typename N, std::size_t SIZE>
Reaction<N,SIZE> reaction_from_gene(const std::array<N, SIZE*3> & gene){
  std::array<N, SIZE> reactants, products;
  std::array<bool, SIZE> inhibitors;
  std::copy(gene.cbegin(), gene.cend() - SIZE*2, reactants.begin());
  std::copy(gene.cbegin() + SIZE, gene.cend() - SIZE, inhibitors.begin());
  std::copy(gene.cbegin() + SIZE*2, gene.cend(), products.begin());
  return Reaction<N,SIZE>(reactants, inhibitors, products);
}

template<typename N, std::size_t SIZE>
std::array<N,GeneSize<N,Reaction,SIZE>::size> gene_from_reaction(const Reaction<N,SIZE> & reaction){
    std::array<N, SIZE*3> gene;
    auto iter = gene.begin();
    iter = std::copy(reaction.get_reactants().cbegin(), reaction.get_reactants().cend(), iter);
    iter = std::copy(reaction.get_inhibitors().cbegin(), reaction.get_inhibitors().cend(), iter);
    std::copy(reaction.get_products().cbegin(), reaction.get_products().cend(), iter);
    return gene;
}

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
FAGenome<N,R,SIZE> genome_from_automaton(const FuncAutomaton<N, R, SIZE, Pure> & aut){
    std::vector<N> genome;
    auto iter = std::back_inserter(genome);
    genome.reserve(SIZE + GeneSize<N,R,SIZE>::size*aut.reactions.size());
    std::copy(aut.initial_state.cbegin(), aut.initial_state.cend(), iter);
    for(auto & reaction : aut.reactions){
        std::array<N,GeneSize<N,R,SIZE>::size> gene = gene_from_reaction(reaction);
        std::copy(gene.cbegin(), gene.cend(), iter);
    }
    return FAGenome<N,R,SIZE>(genome);
}
template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class FAUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
FAUT<N,R,SIZE,Pure> funcautomaton_from_genome(const FAGenome<N,R,SIZE> & genome){
    Multiset<N, SIZE> initial_state;
    std::copy(genome.cbegin(), genome.cbegin()+SIZE,initial_state.begin());
    std::vector<R<N, SIZE>> reactions;
    for(std::size_t i = SIZE; i < genome.size(); i+=GeneSize<N,R,SIZE>::size){
        std::array<N, GeneSize<N,R,SIZE>::size> gene;
        std::copy(genome.cbegin() + i, genome.cbegin() + i + GeneSize<N,R,SIZE>::size, gene.begin());
        reactions.push_back(reaction_from_gene<N, SIZE>(gene));
    }
    return FAUT<N,R,SIZE,Pure>(initial_state, reactions);
}
template<template<typename , template<typename, std::size_t> class , std::size_t, bool> class AAUT, typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
AAUT<N,R,SIZE,Pure> accautomaton_from_genome(const AAGenome<N,R,SIZE> & genome){
    Multiset<N, SIZE> initial_state;
    std::copy(genome.cbegin(), genome.cbegin()+SIZE,initial_state.begin());
    std::vector<R<N, SIZE>> reactions;
    for(std::size_t i = SIZE; i < genome.reactions_size()*GeneSize<N,R,SIZE>::size + SIZE; i+=GeneSize<N,R,SIZE>::size){
        std::array<N, GeneSize<N,R,SIZE>::size> gene;
        std::copy(genome.cbegin() + i, genome.cbegin() + i + GeneSize<N,R,SIZE>::size, gene.begin());
        reactions.push_back(reaction_from_gene<N, SIZE>(gene));
    }
    std::vector<Multiset<N, SIZE>> final_states;
    for(std::size_t i = genome.reactions_size()*GeneSize<N,R,SIZE>::size + SIZE; i<genome.size(); i+=SIZE){
        Multiset<N, SIZE> state;
        std::copy(genome.cbegin() + i, genome.cbegin() + i + SIZE,state.begin());
        final_states.push_back(state);
    }
    return AAUT<N,R,SIZE,Pure>(initial_state, reactions, final_states);
}
template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
AAGenome<N,R,SIZE> genome_from_automaton(const AccAutomaton<N, R, SIZE, Pure> & aut){
    std::vector<N> genome;
    auto iter = std::back_inserter(genome);
    std::copy(aut.initial_state.cbegin(), aut.initial_state.cend(), iter);
    for(auto & r : aut.reactions){
        std::array<N,GeneSize<N,R,SIZE>::size> gene = gene_from_reaction(r);
        std::copy(gene.cbegin(), gene.cend(), iter);
    }
    for(auto & state : aut.accepting_states){
        std::copy(state.cbegin(), state.cend(), iter);
    }
    return AAGenome<N,R,SIZE>(genome, aut.reactions.size());
}

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE>
void FAGenome<N,R,SIZE>::print() const{
    std::cout << "init_state:" << std::endl;
    for(int i = 0; i < SIZE; ++i){
        std::cout << static_cast<unsigned>((*this)[i]) << " ";
    }
    std::cout << std::endl;
    std::cout << "reactions:" << std::endl;
    for(int i = SIZE; i < this->size(); i+=GeneSize<N,R,SIZE>::size){
        std::cout << "r" << (i-SIZE)/GeneSize<N,R,SIZE>::size << ":"<< std::endl;
        for(int j = 0; j < GeneSize<N,R,SIZE>::size; j+=SIZE){
            for(int k = 0; k < SIZE; ++k){
                std::cout << static_cast<unsigned>((*this)[i+j+k]) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}
template<typename N, template<typename, std::size_t> class R, std::size_t SIZE>
void AAGenome<N,R,SIZE>::print() const{
    std::cout << "init_state:" << std::endl;
    for(int i = 0; i < SIZE; ++i){
        std::cout << static_cast<unsigned>((*this)[i]) << " ";
    }
    std::cout << std::endl;
    std::cout << "reactions:" << std::endl;
    for(int i = SIZE; i < SIZE + GeneSize<N,R,SIZE>::size*this->n_reactions; i+=GeneSize<N,R,SIZE>::size){
        std::cout << "r" << (i-SIZE)/GeneSize<N,R,SIZE>::size << ":"<< std::endl;
        for(int j = 0; j < GeneSize<N,R,SIZE>::size; j+=SIZE){
            for(int k = 0; k < SIZE; ++k){
                std::cout << static_cast<unsigned>((*this)[i+j+k]) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "final_states:" << std::endl;
    for(int i = SIZE + GeneSize<N,R,SIZE>::size*this->n_reactions; i < this->size(); i+=SIZE){
        for(int j = 0; j < SIZE; ++j){
            std::cout << static_cast<unsigned>((*this)[i+j]) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}



#ifndef BYTES
#define BYTES 2
#endif

#if BYTES == 1
#define DTYPE uint8_t
#elif BYTES == 2
#define DTYPE uint16_t
#elif BYTES == 4
#define DTYPE uint32_t
#elif BYTES == 8
#define DTYPE uint64_t
#else
#define DTYPE uint16_t
#endif

#ifndef ALPHABET
#define ALPHABET 25
#endif



using FCRA_sq = Automaton_sq<FuncAutomaton,DTYPE, CReaction, ALPHABET, false>;
using FRA_sq = Automaton_sq<AccAutomaton,DTYPE, Reaction, ALPHABET, false>;
using FCRA_mp = Automaton_mp<FuncAutomaton,DTYPE, CReaction, ALPHABET, false>;
using FPCRA_mp = Automaton_mp<FuncAutomaton,DTYPE, CReaction, ALPHABET, true>;
using FRA_mp = Automaton_mp<FuncAutomaton,DTYPE, Reaction, ALPHABET, false>;
using FPRA_mp = Automaton_mp<FuncAutomaton,DTYPE, Reaction, ALPHABET, true>;

using FCRAGen = FAGenome<DTYPE, CReaction, ALPHABET>;
using FRAGen = FAGenome<DTYPE, Reaction, ALPHABET>;

using ARA_sq = AccAutomaton_sq<DTYPE, Reaction, ALPHABET, false>;
using ARA_mp = Automaton_mp<AccAutomaton, DTYPE, Reaction, ALPHABET, false>;
using APRA_mp = Automaton_mp<AccAutomaton, DTYPE, Reaction, ALPHABET, true>;

using ACRAGen = AAGenome<DTYPE, CReaction, ALPHABET>;
using ARAGen = AAGenome<DTYPE, Reaction, ALPHABET>;


using multiset = Multiset<DTYPE, ALPHABET>;
using set = Set<ALPHABET>;
using reaction = Reaction<DTYPE, ALPHABET>;
using creaction = CReaction<DTYPE, ALPHABET>;

FRAGen get_genome(const FuncAutomaton<DTYPE, Reaction, ALPHABET, false> & aut){
  return genome_from_automaton<DTYPE, Reaction, ALPHABET, false>(aut);}
FRA_mp get_FRA_mp(const FRAGen & v){
  return funcautomaton_from_genome<FuncAutomaton_mp,DTYPE, Reaction, ALPHABET, false>(v);
}
ARAGen get_genome(const AccAutomaton<DTYPE, Reaction, ALPHABET, false> & aut){
  return genome_from_automaton<DTYPE, Reaction, ALPHABET, false>(static_cast<const AccAutomaton<DTYPE,Reaction,ALPHABET,false> & >(aut));}
ARA_mp get_ARA_mp(const ARAGen & v){
  return accautomaton_from_genome<AccAutomaton_mp,DTYPE, Reaction, ALPHABET, false>(v);
}
FRAGen get_genome(const FPRA_mp & aut){
    return genome_from_automaton<DTYPE, Reaction, ALPHABET, true>(aut);}
FPRA_mp get_FPRA_mp(const FRAGen & v){
    return funcautomaton_from_genome<FuncAutomaton_mp,DTYPE, Reaction, ALPHABET, true>(v);
}
ARAGen get_genome(const APRA_mp & aut){
    return genome_from_automaton<DTYPE, Reaction, ALPHABET, true>(aut);}
APRA_mp get_APRA_mp(const ARAGen & v){
    return accautomaton_from_genome<AccAutomaton_mp,DTYPE, Reaction, ALPHABET, true>(v);
}

void fgenome_mutation_ngh(DTYPE * ptr, std::size_t len, std::size_t population, std::size_t n_reactions, double mutation_prob, const std::array<double,3>& probs, double m_prob, double rp_prob, double i_prob){
  mutation_raw_vctr<DTYPE,Reaction,ALPHABET,false, false>(ptr, len, population, n_reactions, mutation_prob, probs, m_prob, rp_prob, i_prob);
}
void fgenome_mutation_ow(DTYPE * ptr, std::size_t len, std::size_t population, std::size_t n_reactions, double mutation_prob, const std::array<double,3>& probs, double m_prob, double rp_prob, double i_prob){
  mutation_raw_vctr<DTYPE,Reaction,ALPHABET,false, true>(ptr, len, population, n_reactions, mutation_prob, probs, m_prob, rp_prob, i_prob);
}
void agenome_mutation_ngh(DTYPE * ptr, std::size_t len, std::size_t population, std::size_t n_reactions, double mutation_prob, const std::array<double,3>& probs, double m_prob, double rp_prob, double i_prob){
  mutation_raw_vctr<DTYPE,Reaction,ALPHABET,true, false>(ptr, len, population, n_reactions, mutation_prob, probs, m_prob, rp_prob, i_prob);
}
void agenome_mutation_ow(DTYPE * ptr, std::size_t len, std::size_t population, std::size_t n_reactions, double mutation_prob, const std::array<double,3>& probs, double m_prob, double rp_prob, double i_prob){
  mutation_raw_vctr<DTYPE,Reaction,ALPHABET,true, true>(ptr, len, population, n_reactions, mutation_prob, probs, m_prob, rp_prob, i_prob);
}
/*using RA_sq = Automaton_sq<DTYPE, Reaction>;
using CRA_mp = Automaton_sq<DTYPE, CReaction>;
using RA_mp = Automaton_sq<DTYPE, Reaction>;
using CRA_mr = Automaton_mr<DTYPE, CReaction>;
using RA_mr = Automaton_mr<DTYPE, Reaction>;
*/
#endif