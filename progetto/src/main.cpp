#include"Reactions.hpp"
#include"ReactionAutomata.hpp"
#include"GenReactionAutomata.hpp"
#include<iostream>
template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
double fitness_on_input(FuncAutomaton<N, R, SIZE, Pure> & automaton, const std::function<N(std::vector<N>)> & func, std::vector<N> input,const unsigned int & max_iter){
    if(input.size() > SIZE - 1){
        throw std::invalid_argument("Input size is too big");
    }
    Multiset<N, SIZE> in_ms;
    for (std::size_t i = 0; i < input.size(); ++i){
        in_ms[i+1] = input[i];
    }
    automaton.read_input(in_ms);
    for(auto & i : automaton.state){
      std::cout << i << " ";
    }
    std::cout << std::endl;
    unsigned int count = 0;
    while((count < max_iter) && !(automaton.is_halted)()){
        ++count;
        automaton.react();
        for(auto & i : automaton.state){
            std::cout << i << " ";
        }
    }
    if(!automaton.is_halted()){
        return std::numeric_limits<double>::min();
    }
    return -std::abs(double(func(input)) - double(automaton.result_indicator()));
}


template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
double fitness_on_input_range(FuncAutomaton<N, R, SIZE, Pure> & automaton, const std::function<N(std::vector<N>)> & func, const std::vector<std::vector<N>> & inputs, unsigned int max_iter){
    double fitness = 0.0;
    for(auto & input : inputs){
        fitness += fitness_on_input<N,R,SIZE,Pure>(automaton, func, input, max_iter);
        automaton.reset();
    }
    return fitness;
}

double eval_fitness(const FRAGen & genome,const std::function<DTYPE(std::vector<DTYPE>)> & func, const std::vector<std::vector<DTYPE>> & inputs, unsigned int max_iter){
    auto aut = get_FRA_mp(genome);
    return fitness_on_input_range<DTYPE,Reaction,ALPHABET,false>(aut, func, inputs, max_iter);
}

int main(){
    multiset zero{{0}}, a{{1}}, b{{0,1}}, c{{0,0,1}}, d{{0,0,0,1}};

    set empty{{0}}, set_a{{1}}, set_b{{0,1}}, set_c{{0,0,1}}, set_d{{0,0,0,1}};

    reaction r1(b*2, empty , c);
    reaction r2(c+d, empty , a*3);
    reaction r3(a, set_d , d*2);

    multiset state = zero;

    uint8_t my_int = 42;
    std::cout << "uint8_t: "<< unsigned(my_int) << std::endl;

    std::cout << std::endl;
    FRA_mp aut(state, std::vector<reaction>{r1, r2, r3});
    ARA_mp a_aut(state, std::vector<reaction>{r1, r2, r3}, std::vector<multiset>{zero});
    r2.print();
    r3.print();
    return 0;
}