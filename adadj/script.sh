#!/bin/bash -x

tapenade -reverse                                  \
         -inputlanguage c                          \
         -outputlanguage c                         \
         -head "pressure_cost(q J)/(q J)"          \
         -I ../include                             \
         cost.c

# tapenade -tangent                  \
#          -inputlanguage c          \
#          -outputlanguage c         \
#          -head "timestep(q dt)/(q dt)"  \
#          timestep.c


tail *.msg
