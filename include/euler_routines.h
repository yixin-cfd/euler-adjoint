
#ifdef __cplusplus
extern "C" {
#endif

  void timestep(double (*q)[4], double (*Sj)[2], double (*Sk)[2],  
		double *V, Dim *dim, double cfl, double *dt);

  
#ifdef __cplusplus
}
#endif
