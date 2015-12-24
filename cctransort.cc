/*
 *@BEGIN LICENSE
 *
 * cctransort by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

using namespace boost;

namespace psi{ namespace cctransort {

extern "C"
int read_options(std::string name, Options& options)
{
  if (name == "CCTRANSORT"|| options.read_globals()) {
    /*- The amount of information printed to the output file -*/
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF") ;
  }

    return true;
}

extern "C"
PsiReturnType cctransort(Options& options)
{
  int print = options.get_int("PRINT");

  shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");

  std::vector<shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);

  shared_ptr<IntegralTransform> ints;
  if(options.get_str("REFERENCE") == "RHF")
    ints = shared_ptr<IntegralTransform> (new IntegralTransform(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly));
  else if(options.get_str("REFERENCE") == "ROHF") // not sure if I want semicanonical here yet
    ints = shared_ptr<IntegralTransform> (new IntegralTransform(ref, spaces, IntegralTransform::SemiCanonical, IntegralTransform::DPDOnly));
  else if(options.get_str("REFERENCE") == "UHF")
    ints = shared_ptr<IntegralTransform> (new IntegralTransform(ref, spaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly));
  else
    throw PSIEXCEPTION("Invalid choice of reference wave function.");

  ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  double efzc = ints->get_frozen_core_energy();
  outfile->Printf("Frozen core energy = %20.12f\n", efzc);

  shared_ptr<PSIO> psio(_default_psio_lib_);
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  dpd_set_default(ints->get_dpd_id());
  dpdbuf4 K;
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
  for(int h=0; h < ref->nirrep(); h++) {
    global_dpd_->buf4_mat_irrep_row_init(&K, h);
    for(int pq=0; pq < K.params->rowtot[h]; pq++) {
      global_dpd_->buf4_mat_irrep_row_rd(&K, h);

    }
    global_dpd_->buf4_mat_irrep_row_close(&K, h);
  }
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_LIBTRANS_DPD, 1);

  return Success;
}

}} // End namespaces

