#include "rsmith_gap_quip_shim.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <vector>

extern "C" {
void rsmith_quip_local_energy_wrapper(int *, int *, int *, int *, int *, int *,
                                      int *, int *, int *, double *, int *, int *,
                                      double *, double *);
void quip_lammps_potential_initialise(int *, int *, double *, char *, int *,
                                      char *, int *);
}

struct RsmithGapQuipHandle {
  int n_atoms = 0;
  int n_central = 0;
  bool periodic = true;
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

static std::mutex potential_mutex;
static std::vector<int> shared_quip_potential;
static double shared_cutoff = 0.0;
static std::string shared_model;
static std::string shared_init;

static double pbc_delta(double dx, double box) {
  return dx - box * std::round(dx / box);
}

static double distance2(const RsmithGapQuipHandle &h, int i, int j) {
  double r2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    double dx = h.positions[3 * j + d] - h.positions[3 * i + d];
    if (h.periodic)
      dx = pbc_delta(dx, h.box_lengths[d]);
    r2 += dx * dx;
  }
  return r2;
}

static void build_neighbor_list(const RsmithGapQuipHandle &h,
                                std::vector<int> &num_neigh,
                                std::vector<int> &neigh) {
  const double cutoff2 = h.cutoff * h.cutoff;
  num_neigh.assign(h.n_central, 0);
  neigh.clear();

  for (int i = 0; i < h.n_central; ++i) {
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

static RsmithGapQuipHandle periodic_expansion(const RsmithGapQuipHandle &source) {
  RsmithGapQuipHandle expanded = source;
  expanded.periodic = false;
  for (int atom = 0; atom < source.n_atoms; ++atom) {
    int minimum_shift[3];
    int maximum_shift[3];
    for (int d = 0; d < 3; ++d) {
      const double position = source.positions[3 * atom + d];
      minimum_shift[d] = static_cast<int>(
          std::ceil((-source.cutoff - position) / source.box_lengths[d]));
      maximum_shift[d] = static_cast<int>(std::floor(
          (source.box_lengths[d] + source.cutoff - position) /
          source.box_lengths[d]));
    }
    for (int sx = minimum_shift[0]; sx <= maximum_shift[0]; ++sx) {
      for (int sy = minimum_shift[1]; sy <= maximum_shift[1]; ++sy) {
        for (int sz = minimum_shift[2]; sz <= maximum_shift[2]; ++sz) {
          if (sx == 0 && sy == 0 && sz == 0)
            continue;
          expanded.atomic_numbers.push_back(source.atomic_numbers[atom]);
          expanded.tags.push_back(source.tags[atom]);
          expanded.positions.push_back(source.positions[3 * atom] +
                                       sx * source.box_lengths[0]);
          expanded.positions.push_back(source.positions[3 * atom + 1] +
                                       sy * source.box_lengths[1]);
          expanded.positions.push_back(source.positions[3 * atom + 2] +
                                       sz * source.box_lengths[2]);
        }
      }
    }
  }
  expanded.n_atoms = static_cast<int>(expanded.atomic_numbers.size());
  return expanded;
}

static int set_atoms(RsmithGapQuipHandle *handle, std::size_t n_atoms,
                     std::size_t n_central, const char *const *species,
                     const double *positions, const double *box_lengths,
                     bool periodic, const std::size_t *atom_ids) {
  if (handle == nullptr || species == nullptr || positions == nullptr ||
      box_lengths == nullptr || n_central == 0 || n_central > n_atoms)
    return 1;

  handle->n_atoms = static_cast<int>(n_atoms);
  handle->n_central = static_cast<int>(n_central);
  handle->periodic = periodic;
  handle->atomic_numbers.resize(n_atoms);
  handle->tags.resize(n_atoms);
  handle->ilist.resize(n_central);
  handle->positions.assign(positions, positions + 3 * n_atoms);
  std::copy(box_lengths, box_lengths + 3, handle->box_lengths);

  for (std::size_t i = 0; i < n_atoms; ++i) {
    const int z = atomic_number(species[i]);
    if (z == 0)
      return 2;
    handle->atomic_numbers[i] = z;
    handle->tags[i] =
        static_cast<int>((atom_ids == nullptr ? i : atom_ids[i]) + 1);
    if (i < n_central)
      handle->ilist[i] = static_cast<int>(i);
  }
  return 0;
}

extern "C" RsmithGapQuipHandle *
rsmith_gap_quip_create(const char *model_path, const char *init_args) {
  try {
    if (model_path == nullptr)
      return nullptr;

    auto h = std::make_unique<RsmithGapQuipHandle>();
    std::string model(model_path);
    std::string init = init_args == nullptr ? "Potential" : init_args;
    std::lock_guard<std::mutex> lock(potential_mutex);
    if (!shared_quip_potential.empty()) {
      if (model != shared_model || init != shared_init)
        return nullptr;
      h->quip_potential = shared_quip_potential;
      h->n_quip_potential = static_cast<int>(shared_quip_potential.size());
      h->cutoff = shared_cutoff;
      return h.release();
    }
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
    shared_quip_potential = h->quip_potential;
    shared_cutoff = h->cutoff;
    shared_model = model;
    shared_init = init;
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
  return set_atoms(handle, n_atoms, n_atoms, species, positions, box_lengths,
                   true, nullptr);
}

extern "C" int rsmith_gap_quip_set_local_cluster(
    RsmithGapQuipHandle *handle, std::size_t n_atoms, std::size_t n_central,
    const std::size_t *atom_ids, const char *const *species,
    const double *positions, const double *box_lengths) {
  return set_atoms(handle, n_atoms, n_central, species, positions, box_lengths,
                   false, atom_ids);
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

  std::vector<double> energies(n_indices, 0.0);
  const int rc = rsmith_gap_quip_per_atom_energies(handle, n_indices, indices,
                                                   energies.data());
  if (rc != 0)
    return rc;
  *out_energy = std::accumulate(energies.begin(), energies.end(), 0.0);
  return 0;
}

extern "C" int rsmith_gap_quip_per_atom_energies(
    RsmithGapQuipHandle *handle, std::size_t n_indices,
    const std::size_t *indices, double *out_energies) {
  if (handle == nullptr || out_energies == nullptr ||
      (n_indices > 0 && indices == nullptr))
    return 1;

  std::vector<int> num_neigh;
  std::vector<int> neigh;
  RsmithGapQuipHandle expanded;
  RsmithGapQuipHandle *evaluation = handle;
  if (handle->periodic) {
    expanded = periodic_expansion(*handle);
    evaluation = &expanded;
  }
  build_neighbor_list(*evaluation, num_neigh, neigh);

  int nlocal = evaluation->n_central;
  int nghost = evaluation->n_atoms - evaluation->n_central;
  int inum = evaluation->n_central;
  int sum_num_neigh = static_cast<int>(neigh.size());
  double lattice[9] = {evaluation->box_lengths[0], 0.0, 0.0,
                       0.0, evaluation->box_lengths[1], 0.0,
                       0.0, 0.0, evaluation->box_lengths[2]};
  std::vector<double> local_e(static_cast<std::size_t>(evaluation->n_atoms), 0.0);

  rsmith_quip_local_energy_wrapper(
      &nlocal, &nghost, evaluation->atomic_numbers.data(),
      evaluation->tags.data(), &inum, &sum_num_neigh, evaluation->ilist.data(),
      num_neigh.data(), neigh.data(), lattice,
      evaluation->quip_potential.data(), &evaluation->n_quip_potential,
      evaluation->positions.data(), local_e.data());

  for (std::size_t i = 0; i < n_indices; ++i) {
    if (indices[i] >= static_cast<std::size_t>(handle->n_atoms))
      return 2;
    out_energies[i] = local_e[indices[i]];
  }
  return 0;
}
