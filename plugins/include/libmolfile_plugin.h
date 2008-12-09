#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H
#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DECLARE_PLUGIN(PLUGIN) \
  extern int molfile_ ## PLUGIN ## _init(void); \
  extern int molfile_ ## PLUGIN ## _register(void *, vmdplugin_register_cb); \
  extern int molfile_ ## PLUGIN ## _fini(void);

DECLARE_PLUGIN(dcdplugin)
DECLARE_PLUGIN(jsplugin)
DECLARE_PLUGIN(pdbplugin)
DECLARE_PLUGIN(psfplugin)

#undef DECLARE_PLUGIN

#ifdef __cplusplus
}
#endif
#endif
