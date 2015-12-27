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

#define ID(x) ints->DPD_ID(x)

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
  spaces.push_back(MOSpace::occ);
  spaces.push_back(MOSpace::vir);

  shared_ptr<IntegralTransform> ints;
  if(options.get_str("REFERENCE") == "RHF")
    ints = shared_ptr<IntegralTransform> (new IntegralTransform(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly));
  else if(options.get_str("REFERENCE") == "ROHF") // not sure if I want semicanonical here yet
    ints = shared_ptr<IntegralTransform> (new IntegralTransform(ref, spaces, IntegralTransform::SemiCanonical, IntegralTransform::DPDOnly));
  else if(options.get_str("REFERENCE") == "UHF")
    ints = shared_ptr<IntegralTransform> (new IntegralTransform(ref, spaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly));
  else
    throw PSIEXCEPTION("Invalid choice of reference wave function.");


  shared_ptr<PSIO> psio(_default_psio_lib_);
  dpdbuf4 K;
  dpd_set_default(ints->get_dpd_id());
  ints->set_keep_dpd_so_ints(true);

  outfile->Printf("[O,O] = %d\n", ID("[O,O]"));
  outfile->Printf("[V,V] = %d\n", ID("[V,V]"));
  outfile->Printf("[O,V] = %d\n", ID("[O,V]"));
  outfile->Printf("[V,O] = %d\n", ID("[V,O]"));
  outfile->Printf("[O>=O]+ = %d\n", ID("[O>=O]+"));
  outfile->Printf("[V>=V]+ = %d\n", ID("[V>=V]+"));
  outfile->Printf("[O>O]- = %d\n", ID("[O>O]-"));
  outfile->Printf("[V>V]- = %d\n", ID("[V>V]-"));

  ints->print_dpd_lookup();

  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::MakeAndKeep);
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::MakeAndNuke);

  double efzc = ints->get_frozen_core_energy();
  outfile->Printf("Frozen core energy = %20.12f\n", efzc);

  // Transformations are complete; need to re-set the 

  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
  global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, ID("[O,O]"), ID("[O,O]"), "A <ij|kl>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_AINTS, 1);

  psio->open(PSIF_CC_BINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"), ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, ID("[V,V]"), ID("[V,V]"), "B <ab|cd>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_BINTS, 1);

  psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, ID("[O,V]"), ID("[O,V]"), "C <ia|jb>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_CINTS, 1);

  psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_DINTS, 1);

  psio->open(PSIF_CC_EINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"), ID("[O>=O]+"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, ID("[V,O]"), ID("[O,O]"), "E <ai|jk>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_EINTS, 1);

  psio->open(PSIF_CC_FINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, ID("[O,V]"), ID("[V,V]"), "D <ia|bc>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_FINTS, 1);

  psio->close(PSIF_LIBTRANS_DPD, 1);

  // Additional sorts to mimic ccsort behavior
  for(int i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio->open(i, PSIO_OPEN_OLD);
  dpdbuf4 C, D;

  // C_SORT() 
  // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba>
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia|jb>");
  global_dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <ia||jb>");
  global_dpd_->buf4_close(&C);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, ID("[O,V]"), ID("[O,V]"), "D <ij|ab> (ib,ja)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "D <ij|ab> (ib,ja)");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia||jb>");
  global_dpd_->buf4_axpy(&D, &C, -1);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&C);

  /* <ia|jb> (bi,ja) */
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia|jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, sprq, ID("[V,O]"), ID("[O,V]"), "C <ia|jb> (bi,ja)");
  global_dpd_->buf4_close(&C);

  /* <ia||jb> (bi,ja) */
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia||jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, sprq, ID("[V,O]"), ID("[O,V]"), "C <ia||jb> (bi,ja)");
  global_dpd_->buf4_close(&C);

  /* <ia|jb> (ia,bj) */
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia|jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, ID("[O,V]"), ID("[V,O]"), "C <ia|jb> (ia,bj)");
  global_dpd_->buf4_close(&C);

  /* <ia||jb> (ia,bj) (Wmbej.c) */
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia||jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, ID("[O,V]"), ID("[V,O]"), "C <ia||jb> (ia,bj)");
  global_dpd_->buf4_close(&C);

  /* <ai|bj> (cchbar/Wabei_RHF.c) */
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "C <ia|jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, qpsr, ID("[V,O]"), ID("[V,O]"), "C <ai|bj>");
  global_dpd_->buf4_close(&C);

  // Close up all DPD files
  for(int i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio->close(i,1);
  for(int i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio->close(i,0); /* delete CC_TMP files */
  for(int i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio->close(i,1);

  return Success;
}

}} // End namespaces

