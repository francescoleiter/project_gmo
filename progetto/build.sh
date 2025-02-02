#!/bin/bash

CYAN='\033[0;36m'
RUN_CLR='\033[0;36m'
ERR_CLR='\033[0;31m'
BLD_CLR='\033[0;35m'
ENDCOLOR='\033[0m'

echo -e "${BLD_CLR}╔═════╗              ╦        ╦                            ${ENDCOLOR}"
echo -e "${BLD_CLR}║     ║           o  ║        ║  o                         ${ENDCOLOR}"
echo -e "${BLD_CLR}╠═════╝╗  ╦    ╦  ╦  ║   ╔════╣  ╦  ╗════╗  ╔════╣         ${ENDCOLOR}"
echo -e "${BLD_CLR}║      ║  ║    ║  ║  ║   ║    ║  ║  ║    ║  ║    ║         ${ENDCOLOR}"
echo -e "${BLD_CLR}╚══════╝  ╚════╝  ╩  ╚═╝ ╚════╩  ╩  ╩    ╩  ╚════╣  O  O  O${ENDCOLOR}"
echo -e "${BLD_CLR}                                                 ║         ${ENDCOLOR}"
echo -e "${BLD_CLR}╚════════════════════════════════════════════════╝         ${ENDCOLOR}"

if [ ! -d "./build/" ]; then
    mkdir build
    cmake . -B ./build/
    make --directory=./build/ all
else 
    cmake . -B ./build/
    make --directory=./build/ all
fi

if [ $? -ne 0 ]; then
  echo -e "${ERR_CLR}╔═════╗                                                     ${ENDCOLOR}"
  echo -e "${ERR_CLR}║                                                  ----     ${ENDCOLOR}"
  echo -e "${ERR_CLR}╠════          ╗═══╗       ╗═══╗            ___  | X  X |   ${ENDCOLOR}"
  echo -e "${ERR_CLR}║        ╗═══╗ ║    ╔════╗ ║               / __|  \----/ /  ${ENDCOLOR}"
  echo -e "${ERR_CLR}╩═════╝  ║     ╩    ║    ║ ╩      O  O  O  |/  | \ __ __/   ${ENDCOLOR}"
  echo -e "${ERR_CLR}         ╩          ╩════╝                     |    - -     ${ENDCOLOR}"
  echo -e "${ERR_CLR}                                               |  _|   |_   ${ENDCOLOR}"
else
  echo -e "${RUN_CLR}╔═════╗                                                    ${ENDCOLOR}"
  echo -e "${RUN_CLR}║     ║                          o                         ${ENDCOLOR}"
  echo -e "${RUN_CLR}╠═══╦═╝  ╦    ╦  ╗════╗  ╗════╗  ╦  ╗════╗  ╔════╣         ${ENDCOLOR}"
  echo -e "${RUN_CLR}║   ║    ║    ║  ║    ║  ║    ║  ║  ║    ║  ║    ║         ${ENDCOLOR}"
  echo -e "${RUN_CLR}╩   ╚═╝  ╚════╝  ╩    ╩  ╩    ╩  ╩  ╩    ╩  ╚════╣  O  O  O${ENDCOLOR}"
  echo -e "${RUN_CLR}                                                 ║         ${ENDCOLOR}"
  echo -e "${RUN_CLR}╚════════════════════════════════════════════════╝         ${ENDCOLOR}"
  ./my_program/bin/main
fi