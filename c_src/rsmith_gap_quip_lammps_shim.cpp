#include "rsmith_gap_quip_shim.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <vector>

extern "C" {
void quip_lammps_wrapper(int *, int *, int *, int *, int *, int *, int *,
                         int *, int *, double *, int *, int *, double *,
                         double *, double *, double *, double *, double *);
void quip_lammps_potential_initialise(int *, int *, double *, char *, int *,
                                      char *, int *);
}

struct RsmithGapQuipHandle {
  int n_atoms = 0;
  int n_quip_potential = 0;
  double cutoff = 0.0;
  std::vector<int> quip_potential;
  std::vector<int> atomic_numbers;
  std::vector<int> tags;
  std::vector<int> ilist;
  std::vector<double> positions;
  double box_lengths[3] = {0.0, 0.0, 0.0};
};

static int atomic_number(const char *symbol) {
  if (std::strcmp(symbol, "H") == 0)
    return 1;
  if (std::strcmp(symbol, "He") == 0)
    return 2;
  if (std::strcmp(symbol, "Li") == 0)
    return 3;
  if (std::strcmp(symbol, "Be") == 0)
    return 4;
  if (std::strcmp(symbol, "B") == 0)
    return 5;
  if (std::strcmp(symbol, "C") == 0)
    return 6;
  if (std::strcmp(symbol, "N") == 0)
    return 7;
  if (std::strcmp(symbol, "O") == 0)
    return 8;
  if (std::strcmp(symbol, "F") == 0)
    return 9;
  if (std::strcmp(symbol, "Ne") == 0)
    return 10;
  if (std::strcmp(symbol, "Na") == 0)
    return 11;
  if (std::strcmp(symbol, "Mg") == 0)
    return 12;
  if (std::strcmp(symbol, "Al") == 0)
    return 13;
  if (std::strcmp(symbol, "Si") == 0)
    return 14;
  if (std::strcmp(symbol, "P") == 0)
    return 15;
  if (std::strcmp(symbol, "S") == 0)
    return 16;
  if (std::strcmp(symbol, "Cl") == 0)
    return 17;
  if (std::strcmp(symbol, "Ar") == 0)
    return 18;
  if (std::strcmp(symbol, "K") == 0)
    return 19;
  if (std::strcmp(symbol, "Ca") == 0)
    return 20;
  return 0;
}

static double pbc_delta(double dx, double box) {
  return dx - box * std::round(dx / box);
}

static double distance2(const RsmithGapQuipHandle &h, int i, int j) {
  double r2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    const double dx = pbc_delta(h.positions[3 * j + d] - h.positions[3 * i + d],
                               h.box_lengths[d]);
    r2 += dx * dx;
  }
  return r2;
}

static void build_full_neighbor_list(const RsmithGapQuipHandle &h,
                                     std::vector<int> &num_neigh,
                                     std::vector<int> &neigh) {
  const double cutoff2 = h.cutoff * h.cutoff;
  num_neigh.assign(h.n_atoms, 0);
  neigh.clear();

  for (int i = 0; i < h.n_atoms; ++i) {
    for (int j = 0; j < h.n_atoms; ++j) {
      if (i == j)
        continue;
      if (distance2(h, i, j) < cutoff2) {
        ++num_neigh[i];
        neigh.push_back(j + 1);
      }
    }
  }
}

extern "C" RsmithGapQuipHandle *
rsmith_gap_quip_create(const char *model_path, const char *init_args) {
  try {
    if (model_path == nullptr)
      return nullptr;

    auto h = std::make_unique<RsmithGapQuipHandle>();
    std::string model(model_path);
    std::string init = init_args == nullptr ? "Potential" : init_args;
    int model_len = static_cast<int>(model.size());
    int init_len = static_cast<int>(init.size());

    int dummy = 0;
    h->n_quip_potential = 0;
    quip_lammps_potential_initialise(&dummy, &h->n_quip_potential, &h->cutoff,
                                     model.data(), &model_len, init.data(),
                                     &init_len);
    if (h->n_quip_potential <= 0)
      return nullptr;

    h->quip_potential.resize(static_cast<std::size_t>(h->n_quip_potential));
    quip_lammps_potential_initialise(h->quip_potential.data(),
                                     &h->n_quip_potential, &h->cutoff,
                                     model.data(), &model_len, init.data(),
                                     &init_len);
    return h.release();
  } catch (...) {
    return nullptr;
  }
}

extern "C" void rsmith_gap_quip_destroy(RsmithGapQuipHandle *handle) {
  delete handle;
}

extern "C" int rsmith_gap_quip_set_structure(
    RsmithGapQuipHandle *handle, std::size_t n_atoms, const char *const *species,
    const double *positions, const double *box_lengths) {
  if (handle == nullptr || species == nullptr || positions == nullptr ||
      box_lengths == nullptr)
    return 1;

  handle->n_atoms = static_cast<int>(n_atoms);
  handle->atomic_numbers.resize(n_atoms);
  handle->tags.resize(n_atoms);
  handle->ilist.resize(n_atoms);
  handle->positions.assign(positions, positions + 3 * n_atoms);
  std::copy(box_lengths, box_lengths + 3, handle->box_lengths);

  for (std::size_t i = 0; i < n_atoms; ++i) {
    const int z = atomic_number(species[i]);
    if (z == 0)
      return 2;
    handle->atomic_numbers[i] = z;
    handle->tags[i] = static_cast<int>(i + 1);
    handle->ilist[i] = static_cast<int>(i);
  }
  return 0;
}

extern "C" int rsmith_gap_quip_move_atom(RsmithGapQuipHandle *handle,
                                         std::size_t atom_idx,
                                         const double *position) {
  if (handle == nullptr || position == nullptr ||
      atom_idx >= static_cast<std::size_t>(handle->n_atoms))
    return 1;
  std::copy(position, position + 3, &handle->positions[3 * atom_idx]);
  return 0;
}

extern "C" int rsmith_gap_quip_per_atom_energy(RsmithGapQuipHandle *handle,
                                               std::size_t n_indices,
                                               const std::size_t *indices,
                                               double *out_energy) {
  if (handle == nullptr || out_energy == nullptr ||
      (n_indices > 0 && indices == nullptr))
    return 1;

  std::vector<int> num_neigh;
  std::vector<int> neigh;
  build_full_neighbor_list(*handle, num_neigh, neigh);

  int nlocal = handle->n_atoms;
  int nghost = 0;
  int inum = handle->n_atoms;
  int sum_num_neigh = static_cast<int>(neigh.size());
  double lattice[9] = {handle->box_lengths[0], 0.0, 0.0,
                       0.0, handle->box_lengths[1], 0.0,
                       0.0, 0.0, handle->box_lengths[2]};
  double energy = 0.0;
  std::vector<double> local_e(static_cast<std::size_t>(handle->n_atoms), 0.0);
  std::vector<double> virial(9, 0.0);
  std::vector<double> local_virial(static_cast<std::size_t>(9 * handle->n_atoms),
                                   0.0);
  std::vector<double> force(static_cast<std::size_t>(3 * handle->n_atoms), 0.0);

  quip_lammps_wrapper(&nlocal, &nghost, handle->atomic_numbers.data(),
                      handle->tags.data(), &inum, &sum_num_neigh,
                      handle->ilist.data(), num_neigh.data(), neigh.data(),
                      lattice, handle->quip_potential.data(),
                      &handle->n_quip_potential, handle->positions.data(),
                      &energy, local_e.data(), virial.data(), local_virial.data(),
                      force.data());

  double sum = 0.0;
  for (std::size_t i = 0; i < n_indices; ++i) {
    if (indices[i] >= static_cast<std::size_t>(handle->n_atoms))
      return 2;
    sum += local_e[indices[i]];
  }
  *out_energy = sum;
  return 0;
}
