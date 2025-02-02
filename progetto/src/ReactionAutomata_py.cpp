#include"Reactions.hpp"
#include"ReactionAutomata.hpp"
#include"GenReactionAutomata.hpp"
#include<vector>
#include<pybind11/pybind11.h>
#include<pybind11/stl.h>
#include<pybind11/functional.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include<string>
#include<iostream>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<DTYPE>);

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
double fitness_on_input(FuncAutomaton<N, R, SIZE, Pure> & automaton, const std::function<double(std::vector<N>)> & func, std::vector<N> input,const unsigned int & max_iter){
  if(input.size() > SIZE - 1){
    throw std::invalid_argument("Input size is too big");
  }
  Multiset<N, SIZE> in_ms;
  for (std::size_t i = 0; i < input.size(); ++i){
    in_ms[i+1] = input[i];
  }
  automaton.read_input(in_ms);
  unsigned int count = 0;
  while((count < max_iter) && !(automaton.is_halted)()){
    ++count;
    automaton.react();
  }
  if(!automaton.is_halted()){
    return -std::numeric_limits<double>::max();
  }
  return -std::abs(func(input) - automaton.result_indicator());
}

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
double fitness_on_sequence(AccAutomaton<N, R, SIZE, Pure> & automaton, const std::vector<N> & seq, const bool & accepted,const unsigned int & max_iter){
  unsigned int count = 0;
  for(const N & input : seq){
    automaton.read_single_input(input);
    automaton.react();
    double sum = 0.0;
    for(std::size_t i = 0; i<ALPHABET; ++i){
        sum += automaton.state[i];
    }
    if(sum > std::numeric_limits<N>::max()/3){
        return -std::numeric_limits<double>::max();
    }
  }
  while((count < max_iter) && !(automaton.is_halted())){
    ++count;
    automaton.react();
    double sum = 0.0;
    for(std::size_t i = 0; i<ALPHABET; ++i){
        sum += automaton.state[i];
    }
  }
  if(!automaton.is_halted()){
    return -std::numeric_limits<double>::max();
  }
  return static_cast<double>(accepted == automaton.result());
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

template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
double fitness_on_sequence_range(AccAutomaton<N, R, SIZE, Pure> & automaton,const std::vector<std::vector<N>> & inputs, const std::vector<bool> & accepted, const unsigned int & max_iter){
  double fitness = 0.0;
  for(std::size_t i = 0; i < inputs.size(); ++i){
    fitness += fitness_on_sequence<N,R,SIZE,Pure>(automaton, inputs[i], accepted[i], max_iter);
    automaton.reset();
  }
  return fitness;
}

double eval_fitness_input(const FRAGen & genome,const std::function<DTYPE(std::vector<DTYPE>)> & func, const std::vector<std::vector<DTYPE>> & inputs, unsigned int max_iter){
    auto aut = get_FRA_mp(genome);
    return fitness_on_input_range<DTYPE,Reaction,ALPHABET,false>(aut, func, inputs, max_iter);
}
double eval_fitness_input_pure(const FRAGen & genome,const std::function<DTYPE(std::vector<DTYPE>)> & func, const std::vector<std::vector<DTYPE>> & inputs, unsigned int max_iter){
    auto aut = get_FPRA_mp(genome);
    return fitness_on_input_range<DTYPE,Reaction,ALPHABET,true>(aut, func, inputs, max_iter);
}


double eval_fitness_sequence(const ARAGen & genome,const std::vector<std::vector<DTYPE>> & inputs, const std::vector<bool> & accepted, const unsigned int & max_iter){
    auto aut = get_ARA_mp(genome);
    return fitness_on_sequence_range<DTYPE,Reaction,ALPHABET,false>(aut, inputs, accepted, max_iter);
}
double eval_fitness_sequence_pure(const ARAGen & genome,const std::vector<std::vector<DTYPE>> & inputs, const std::vector<bool> & accepted, const unsigned int & max_iter){
    auto aut = get_APRA_mp(genome);
    return fitness_on_sequence_range<DTYPE,Reaction,ALPHABET,true>(aut, inputs, accepted, max_iter);
}



template<typename N, template<typename, std::size_t> class R, std::size_t SIZE, bool Pure>
FRAGen get_FRAGen(const FuncAutomaton<N, R, SIZE, Pure> & aut){
    return FRAGen();
}


multiset a{{1}}, b{{0,1}}, c{{0,0,1}}, d{{0,0,0,1}};

set empty{{0}}, set_a{{1}}, set_b{{0,1}}, set_c{{0,0,1}}, set_d{{0,0,1}};

reaction r1(b*2, empty , c);
reaction r2(c+d, empty , a*3);
reaction r3(a, set_d , d*2);

multiset state = a+d+c;

FRA_mp aut(state, std::vector<reaction>{r1, r2, r3});



PYBIND11_MODULE(pyRA, m) {
  using namespace pybind11::literals;


  py::bind_vector<std::vector<DTYPE>>(m, "CVect")
      .def("to_list", [](const std::vector<DTYPE> & vct) -> py::list{
        py::list res;
        for(auto & v : vct){
          res.append(v);
        }
        return res;
      });
  m.attr("ALPHABET") = ALPHABET;
  m.attr("RGENE_SIZE") = GeneSize<DTYPE,Reaction,ALPHABET>::size;
  py::class_<FRAGen>(m, "FGenome")
      .def(py::init<const std::vector<DTYPE> &>(), "genome"_a)
      .def("print", static_cast<void (FRAGen::*)(void) const>(&FRAGen::print))
      .def("__getitem__",static_cast<DTYPE & (FRAGen::*)(std::size_t)>(&FRAGen::operator[]), "index"_a);

  py::class_<ARAGen>(m, "AGenome")
      .def(py::init<const std::vector<DTYPE> &, const std::size_t &>(), "genome"_a, "n_reactions"_a)
      .def("print", static_cast<void (ARAGen::*)(void) const>(&ARAGen::print))
      .def("__getitem__",static_cast<DTYPE & (ARAGen::*)(std::size_t)>(&ARAGen::operator[]), "index"_a);

  py::class_<multiset>(m, "Multiset")
      .def(py::init<const std::array<DTYPE, ALPHABET> &>());
  py::class_<reaction>(m, "Reaction")
      .def(py::init<const multiset &,const set &,const multiset &>());

  py::class_<FRA_mp>(m, "FAutomaton")
      .def(py::init<const std::array<DTYPE,ALPHABET> &, const std::vector<reaction> &>())
      .def("read_input", static_cast<void (FRA_mp::*)(const multiset &)>(&FRA_mp::read_input))
      .def("react", static_cast<void (FRA_mp::*)(void)>(&FRA_mp::react))
      .def("verbose_react", static_cast<void (FRA_mp::*)(void)>(&FRA_mp::verbose_react))
      .def("is_halted", static_cast<const bool & (FRA_mp::*)(void) const>(&FRA_mp::is_halted))
      .def("get_state", static_cast<const multiset & (FRA_mp::*)(void) const>(&FRA_mp::get_state))
      .def("result_indicator", static_cast<const DTYPE & (FRA_mp::*)(void) const>(&FRA_mp::result_indicator));

  py::class_<ARA_mp>(m, "AAutomaton")
      .def(py::init<const std::array<DTYPE,ALPHABET> &, const std::vector<reaction> &, const std::vector<multiset> &>())
      .def("read_input", static_cast<void (ARA_mp::*)(const multiset &)>(&ARA_mp::read_input))
      .def("read_single_input", static_cast<void (ARA_mp::*)(const DTYPE &)>(&ARA_mp::read_single_input))
      .def("react", static_cast<void (ARA_mp::*)(void)>(&ARA_mp::react))
      .def("verbose_react", static_cast<void (ARA_mp::*)(void)>(&ARA_mp::verbose_react))
      .def("is_halted", static_cast<const bool & (ARA_mp::*)(void) const>(&ARA_mp::is_halted))
      .def("get_state", static_cast<const multiset & (ARA_mp::*)(void) const>(&ARA_mp::get_state))
      .def("accepted" ,static_cast<bool (ARA_mp::*)(void) const>(&ARA_mp::result));
  py::class_<FPRA_mp>(m, "FPAutomaton")
      .def(py::init<const std::array<DTYPE,ALPHABET> &, const std::vector<reaction> &>())
      .def("read_input", static_cast<void (FPRA_mp::*)(const multiset &)>(&FPRA_mp::read_input))
      .def("react", static_cast<void (FPRA_mp::*)(void)>(&FPRA_mp::react))
      .def("verbose_react", static_cast<void (FPRA_mp::*)(void)>(&FPRA_mp::verbose_react))
      .def("is_halted", static_cast<const bool & (FPRA_mp::*)(void) const>(&FPRA_mp::is_halted))
      .def("get_state", static_cast<const multiset & (FPRA_mp::*)(void) const>(&FPRA_mp::get_state))
      .def("result_indicator", static_cast<const DTYPE & (FPRA_mp::*)(void) const>(&FPRA_mp::result_indicator));

  py::class_<APRA_mp>(m, "APAutomaton")
      .def(py::init<const std::array<DTYPE,ALPHABET> &, const std::vector<reaction> &, const std::vector<multiset> &>())
      .def("read_input", static_cast<void (APRA_mp::*)(const multiset &)>(&APRA_mp::read_input))
      .def("read_single_input", static_cast<void (APRA_mp::*)(const DTYPE &)>(&APRA_mp::read_single_input))
      .def("react", static_cast<void (APRA_mp::*)(void)>(&APRA_mp::react))
      .def("verbose_react", static_cast<void (APRA_mp::*)(void)>(&APRA_mp::verbose_react))
      .def("is_halted", static_cast<const bool & (APRA_mp::*)(void) const>(&APRA_mp::is_halted))
      .def("get_state", static_cast<const multiset & (APRA_mp::*)(void) const>(&APRA_mp::get_state))
      .def("accepted" ,static_cast<bool (APRA_mp::*)(void) const>(&APRA_mp::result));;



  m.def("read_fgenome", &get_FRA_mp);
  m.def("read_agenome", &get_ARA_mp);
  m.def("read_fpgenome", &get_FPRA_mp);
  m.def("read_apgenome", &get_APRA_mp);

  m.def("eval_fitness_func", &eval_fitness_input);
  m.def("eval_fitness_seq", &eval_fitness_sequence);
  m.def("eval_fitness_func_pure", &eval_fitness_input_pure);
  m.def("eval_fitness_seq_pure", &eval_fitness_sequence_pure);

  m.def("fgenome_mutation_ow", &fgenome_mutation_ow);
  m.def("agenome_mutation_ow", &agenome_mutation_ow);
  m.def("fgenome_mutation_ngh", &fgenome_mutation_ngh);
  m.def("agenome_mutation_ngh", &agenome_mutation_ngh);

  m.def("fgenome_2Dmutation_ow", [](py::array_t<DTYPE, py::array::c_style> & arr, const std::array<double,3>& probs, double mutation_prob, double m_prob, double rp_prob, double i_prob ){
        auto r = arr.mutable_unchecked<2>();

        auto buffer = arr.request();
        if (buffer.ndim == 2){
            fgenome_mutation_ow(r.mutable_data(0,0), r.shape(1), r.shape(0), (r.shape(1)/(ALPHABET) - 1)/3, mutation_prob, probs, m_prob, rp_prob, i_prob);
        }
  });
  m.def("fgenome_2Dmutation_ngh", [](py::array_t<DTYPE, py::array::c_style> & arr, const std::array<double,3>& probs, double mutation_prob, double m_prob, double rp_prob, double i_prob ){
        auto r = arr.mutable_unchecked<2>();

        auto buffer = arr.request();
        if (buffer.ndim == 2){
            fgenome_mutation_ngh(r.mutable_data(0,0), r.shape(1), r.shape(0), (r.shape(1)/(ALPHABET) - 1)/3, mutation_prob, probs, m_prob, rp_prob, i_prob);
        }
  });
  m.def("agenome_2Dmutation_ow", [](py::array_t<DTYPE, py::array::c_style> & arr, std::size_t n_reactions, const std::array<double,3>& probs, double mutation_prob, double m_prob, double rp_prob, double i_prob ){
        auto r = arr.mutable_unchecked<2>();
        auto buffer = arr.request();
        if (buffer.ndim == 2){
            agenome_mutation_ow(r.mutable_data(0,0), r.shape(1), r.shape(0), n_reactions, mutation_prob, probs, m_prob, rp_prob, i_prob);
        }
  });
  m.def("agenome_2Dmutation_ngh", [](py::array_t<DTYPE, py::array::c_style> & arr, std::size_t n_reactions, const std::array<double,3>& probs, double mutation_prob, double m_prob, double rp_prob, double i_prob ){
        auto r = arr.mutable_unchecked<2>();
        auto buffer = arr.request();
        if (buffer.ndim == 2){
            agenome_mutation_ngh(r.mutable_data(0,0), r.shape(1), r.shape(0), n_reactions, mutation_prob, probs, m_prob, rp_prob, i_prob);
        }
  });
  m.def("to_vect", [](py::array_t<DTYPE, py::array::c_style> arr) {
     return std::vector<DTYPE>(arr.mutable_data(),arr.mutable_data() + arr.size());});

  py::register_exception<std::runtime_error>(m, "RuntimeError");
}