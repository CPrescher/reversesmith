#ifndef RSMITH_GAP_QUIP_SHIM_H
#define RSMITH_GAP_QUIP_SHIM_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct RsmithGapQuipHandle RsmithGapQuipHandle;

RsmithGapQuipHandle *rsmith_gap_quip_create(
    const char *model_path,
    const char *init_args);

void rsmith_gap_quip_destroy(RsmithGapQuipHandle *handle);

int rsmith_gap_quip_set_structure(
    RsmithGapQuipHandle *handle,
    size_t n_atoms,
    const char *const *species,
    const double *positions,
    const double *box_lengths);

int rsmith_gap_quip_set_local_cluster(
    RsmithGapQuipHandle *handle,
    size_t n_atoms,
    size_t n_central,
    const size_t *atom_ids,
    const char *const *species,
    const double *positions,
    const double *box_lengths);

int rsmith_gap_quip_move_atom(
    RsmithGapQuipHandle *handle,
    size_t atom_idx,
    const double *position);

int rsmith_gap_quip_per_atom_energy(
    RsmithGapQuipHandle *handle,
    size_t n_indices,
    const size_t *indices,
    double *out_energy);

int rsmith_gap_quip_per_atom_energies(
    RsmithGapQuipHandle *handle,
    size_t n_indices,
    const size_t *indices,
    double *out_energies);

#ifdef __cplusplus
}
#endif

#endif
