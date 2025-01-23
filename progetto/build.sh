#!/bin/bash

echo "╔═════╗              ╦        ╦                            "
echo "║     ║           o  ║        ║  o                         "
echo "╠═════╝╗  ╦    ╦  ╦  ║   ╔════╣  ╦  ╗════╗  ╔════╣         "
echo "║      ║  ║    ║  ║  ║   ║    ║  ║  ║    ║  ║    ║         "
echo "╚══════╝  ╚════╝  ╩  ╚═╝ ╚════╩  ╩  ╩    ╩  ╚════╣  O  O  O"
echo "                                                 ║         "
echo "╚════════════════════════════════════════════════╝         "

if [ ! -d "./build/" ]; then
    mkdir build
    cmake . -B ./build/
    make --directory=./build/ all
else 
    cmake . -B ./build/
    make --directory=./build/ all
fi

if [ $? -ne 0 ]; then
  echo "╔═════╗                                                    "
  echo "║                                                          "
  echo "╠════          ╗═══╗       ╗═══╗                           "
  echo "║        ╗═══╗ ║    ╔════╗ ║                               "
  echo "╩═════╝  ║     ╩    ║    ║ ╩      O  O  O                "
  echo "         ╩          ╩════╝                                 "
  echo "                                                           "
else
  echo "╔═════╗                                                    "
  echo "║     ║                          o                         "
  echo "╠═══╦═╝  ╦    ╦  ╗════╗  ╗════╗  ╦  ╗════╗  ╔════╣         "
  echo "║   ║    ║    ║  ║    ║  ║    ║  ║  ║    ║  ║    ║         "
  echo "╩   ╚═╝  ╚════╝  ╩    ╩  ╩    ╩  ╩  ╩    ╩  ╚════╣  O  O  O"
  echo "                                                 ║         "
  echo "╚════════════════════════════════════════════════╝         "
  ./my_program/bin/main
fi