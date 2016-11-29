#!/bin/bash -x

# tapenade -reverse                                    \
#          -inputlanguage  c                           \
#          -outputlanguage c                           \
#          -I ../include                               \
#          -head           "pressure_cost(q J)/(q J)"  \
#          cost.c

tapenade -reverse                      \
	 -inputlanguage  c             \
	 -outputlanguage c             \
	 -I              ../include    \
	 -head           "timestep(q xy1 xy2 xy3 dt)/(q xy1 xy2 xy3 dt)"      \
	 timestep.c

# tapenade -reverse                      \
# 	 -inputlanguage  c             \
# 	 -outputlanguage c             \
# 	 -I              ../include    \
# 	 -head           "periodic_bc(q)/(q)"      \
# 	 bc_routines.c

# tapenade -reverse                      \
# 	 -inputlanguage  c             \
# 	 -outputlanguage c             \
# 	 -I              ../include    \
# 	 -head           "wall_bc(q xy1 xy2)/(q xy1 xy2)"      \
# 	 bc_routines.c

tapenade -reverse                      \
	 -inputlanguage  c             \
	 -outputlanguage c             \
	 -I              ../include    \
	 -head           "jroeflux(q_l q_r rhs_m1 rhs xy1 xy2)/(q_l q_r rhs_m1 rhs xy1 xy2)"      \
	 roeflux.c

tapenade -reverse                      \
	 -inputlanguage  c             \
	 -outputlanguage c             \
	 -I              ../include    \
	 -head           "kroeflux(q_l q_r rhs_m1 rhs xy1 xy2)/(q_l q_r rhs_m1 rhs xy1 xy2)"      \
	 roeflux.c


# tapenade -tangent                  \
#          -inputlanguage c          \
#          -outputlanguage c         \
#          -head "timestep(q dt)/(q dt)"  \
#          timestep.c


tail *.msg
