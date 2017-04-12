TPN = ${TAPENADE_HOME}/bin/tapenade

ALL :  tmpflux_b.c dadi_b.c

dadi_b.c : dadi.c
	${TPN} -reverse                              \
	       -inputlanguage  c                     \
	       -outputlanguage c                     \
	       -I              ../include            \
	       -head           "dadi(q rhs)/(q rhs)" \
	       -adjfuncname    "_b"                  \
	       -o dadi                               \
	       dadi.c && cat *.msg


tmpflux_b.c : tmpflux.c
	${TPN} -reverse                      \
	       -inputlanguage  c             \
	       -outputlanguage c             \
	       -I              ../include    \
	       -head           "tmpflux(q_l q_r rhs_m1 rhs)/(q_l q_r rhs_m1 rhs)" \
	       -adjfuncname    "_b"                     \
	       -o tmpflux  \
	       tmpflux.c && cat *.msg

#tmpflux_b.c : tmpflux.c
#	${TPN} -reverse                      \
#	       -inputlanguage  c             \
#	       -outputlanguage c             \
#	       -I              ../include    \
#	       -head           "tmpflux(q, rrhs)/(q rrhs)" \
#	       -adjfuncname    "_b"                     \
#	       -o tmpflux  \
#	       tmpflux.c && cat *.msg

# tmpflux_b.c : tmpflux.c
# 	${TPN} -reverse                      \
# 	       -inputlanguage  c             \
# 	       -outputlanguage c             \
# 	       -I              ../include    \
# 	       -head           "tmpflux(xy1 xy2 rhs_m1 rhs)/(xy1 xy2 rhs_m1 rhs)" \
# 	       -adjfuncname    "_b"                     \
# 	       -o tmpflux  \
# 	       tmpflux.c && cat *.msg


clean :
	rm *~ *_b.c *_b.h *.msg *_bx.c 
