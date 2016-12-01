TPN = ${TAPENADE_HOME}/bin/tapenade

XDIFF = timestep(q xy1 xy2 xy3 dt)/(q xy1 xy2 xy3 dt) \
	wall_bc(q xy1 xy2)/(q xy1 xy2) \
	jroeflux(q_l q_r rhs_m1 rhs xy1 xy2)/(q_l q_r rhs_m1 rhs xy1 xy2) \
	kroeflux(q_l q_r rhs_m1 rhs xy1 xy2)/(q_l q_r rhs_m1 rhs xy1 xy2)

QDIFF = timestep(q dt)/(q dt) \
	periodic_bc(q)/(q)    \
	wall_bc(q)/(q) \
	pressure_cost(q J)/(q J) \
	jroeflux(q_l q_r rhs_m1 rhs)/(q_l q_r rhs_m1 rhs) \
	kroeflux(q_l q_r rhs_m1 rhs)/(q_l q_r rhs_m1 rhs)

FILES = timestep.c bc_routines.c cost.c roeflux.c

# ALL : timestep_b.c  periodic_bc_b.c wall_bc_b.c jkroeflux_b.c kroeflux_b.c pressure_cost_b.c\
#       timestep_bx.c wall_bc_bx.c jkroeflux_bx.c kroeflux_bx.c
ALL :  pressure_cost_b.c timestep_b.c  timestep_bx.c periodic_bc_b.c wall_bc_b.c wall_bc_bx.c \
	jkroeflux_b.c jkroeflux_bx.c lift_cost_b.c lift_cost_bx.c

pressure_cost_b.c : cost.c
	${TPN} -reverse                                   \
	       -inputlanguage  c                          \
	       -outputlanguage c                          \
	       -I              ../include                 \
	       -head           "pressure_cost(q J)/(q J)" \
	       -adjfuncname    "_b"                       \
	       -o pressure_cost                           \
	       cost.c

lift_cost_b.c : cost.c
	${TPN} -reverse                                   \
	       -inputlanguage  c                          \
	       -outputlanguage c                          \
	       -I              ../include                 \
	       -head           "lift_cost(q J)/(q J)"     \
	       -adjfuncname    "_b"                       \
	       -o lift_cost                               \
	       cost.c

timestep_b.c : timestep.c
	${TPN} -reverse                                   \
	       -inputlanguage  c                          \
	       -outputlanguage c                          \
	       -I              ../include                 \
	       -head           "ad_timestep(q dt)/(q dt)" \
	       -adjfuncname    "_b"                       \
	       -o timestep                                \
	       timestep.c

periodic_bc_b.c : bc_routines.c
	${TPN} -reverse                                 \
	       -inputlanguage  c                        \
	       -outputlanguage c                        \
	       -I              ../include               \
	       -head           "periodic_bc(q)/(q)"     \
	       -adjfuncname    "_b"                     \
	       -o periodic_bc                           \
	       bc_routines.c

wall_bc_b.c : bc_routines.c
	${TPN} -reverse                      \
	       -inputlanguage  c             \
	       -outputlanguage c             \
	       -I              ../include    \
	       -head           "wall_bc(q)/(q)" \
	       -adjfuncname    "_b"                     \
	       -o wall_bc \
	       bc_routines.c

jkroeflux_b.c : roeflux.c
	${TPN} -reverse                      \
	       -inputlanguage  c             \
	       -outputlanguage c             \
	       -I              ../include    \
	       -head           "jkroeflux(q_l q_r rhs_m1 rhs)/(q_l q_r rhs_m1 rhs)" \
	       -adjfuncname    "_b"                     \
	       -o jkroeflux  \
	       roeflux.c

# =----------------------------------------------------------
timestep_bx.c : timestep.c
	${TPN} -reverse                      \
	       -inputlanguage  c             \
	       -outputlanguage c             \
	       -I              ../include    \
	       -head           "ad_timestep(q xy1 xy2 xy3 dt)/(q xy1 xy2 xy3 dt)" \
	       -adjfuncname    "_bx"                     \
	       -o timestep \
	       timestep.c

wall_bc_bx.c : bc_routines.c
	${TPN} -reverse                      \
	       -inputlanguage  c             \
	       -outputlanguage c             \
	       -I              ../include    \
	       -head           "wall_bc(q xy1 xy2)/(q xy1 xy2)" \
	       -adjfuncname    "_bx"                     \
	       -o wall_bc \
	       bc_routines.c

jkroeflux_bx.c : roeflux.c
	${TPN} -reverse                      \
	       -inputlanguage  c             \
	       -outputlanguage c             \
	       -I              ../include    \
	       -head           "jkroeflux(q_l q_r rhs_m1 rhs xy1 xy2)/(q_l q_r rhs_m1 rhs xy1 xy2)" \
	       -adjfuncname    "_bx"                     \
	       -o jkroeflux  \
	       roeflux.c

lift_cost_bx.c : cost.c
	${TPN} -reverse                                               \
	       -inputlanguage  c                                      \
	       -outputlanguage c                                      \
	       -I              ../include                             \
	       -head           "lift_cost(q xy1 xy2 J)/(q xy1 xy2 J)" \
	       -adjfuncname    "_bx"                                  \
	       -o lift_cost                                           \
	       cost.c

# qdiff_routines_b.c : ${FILES}
# 	${TPN} -reverse                      \
# 	       -inputlanguage  c             \
# 	       -outputlanguage c             \
# 	       -I              ../include    \
# 	       -head           "${QDIFF}"    \
# 	       -adjfuncname "_b"             \
# 	       -o qdiff_routines             \
# 	       ${FILES}

# xdiff_routines_bx.c : ${FILES}
# 	${TPN} -reverse                      \
# 	       -inputlanguage  c             \
# 	       -outputlanguage c             \
# 	       -I              ../include    \
# 	       -head           "${XDIFF}"    \
# 	       -adjfuncname "_bx"            \
# 	       -o xdiff_routines             \
# 	       ${FILES}

clean :
	rm *~ *_b.c *_b.h *.msg *_bx.c 
