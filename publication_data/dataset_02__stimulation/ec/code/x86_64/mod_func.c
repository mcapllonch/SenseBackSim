#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _gaines_sensory_flut_reg(void);
extern void _gaines_sensory_mysa_reg(void);
extern void _gaines_sensory_node_reg(void);
extern void _gaines_sensory_stin_reg(void);
extern void _MRG_AXNODE_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," gaines_sensory_flut.mod");
    fprintf(stderr," gaines_sensory_mysa.mod");
    fprintf(stderr," gaines_sensory_node.mod");
    fprintf(stderr," gaines_sensory_stin.mod");
    fprintf(stderr," MRG_AXNODE.mod");
    fprintf(stderr, "\n");
  }
  _gaines_sensory_flut_reg();
  _gaines_sensory_mysa_reg();
  _gaines_sensory_node_reg();
  _gaines_sensory_stin_reg();
  _MRG_AXNODE_reg();
}
