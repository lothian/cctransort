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
#include <libciomr/libciomr.h>

#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

// #define ID(x) ints->DPD_ID(x)

using namespace boost;

namespace psi{ namespace cctransort {

double trans(boost::shared_ptr<Wavefunction> ref, Options &options, std::vector<boost::shared_ptr<MOSpace> > spaces);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);

void c_sort(int reference);
void d_sort(int reference);
void e_sort(int reference);
void f_sort(int reference);
void d_spinad();
void e_spinad();

extern "C"
int read_options(std::string name, Options& options)
{
  if (name == "CCTRANSORT"|| options.read_globals()) {
    /*- The amount of information printed to the output file -*/
    options.add_int("PRINT_LEVEL", 1);
    options.add_str("REFERENCE", "RHF") ;
    options.add_str("WFN", "CCSD") ;
  }

    return true;
}

extern "C"
PsiReturnType cctransort(Options& options)
{
  int print = options.get_int("PRINT_LEVEL");

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);

  // Grab all the basic MO data
  boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");

  bool semicanonical = false;
  int reference = 0;
  if(options.get_str("REFERENCE") =="RHF") reference = 0;
  else if(options.get_str("REFERENCE") =="ROHF" &&
          (options.get_str("WFN")=="MP2" || options.get_str("WFN")=="CCSD_T" || options.get_str("WFN")=="CCSD_AT" ||
           options.get_str("WFN")=="CC3" || options.get_str("WFN")=="EOM_CC3" ||
           options.get_str("WFN")=="CC2" || options.get_str("WFN")=="EOM_CC2")) {
    reference = 2;
    semicanonical = true;
  }
  else if(options.get_str("REFERENCE") == "ROHF") reference = 1;
  else if(options.get_str("REFERENCE") == "UHF") reference = 2;
  else {
    outfile->Printf("Invalid value of input keyword REFERENCE: %s\n", options.get_str("REFERENCE").c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  int nirreps = ref->nirrep();
  int nmo = ref->nmo();
  int nso = ref->nso();
  int nao = ref->basisset()->nao();
  char **labels = ref->molecule()->irrep_labels();
  double enuc = ref->molecule()->nuclear_repulsion_energy();
  double escf;
  if(ref->reference_wavefunction())
      escf = ref->reference_wavefunction()->reference_energy();
  else
      escf = ref->reference_energy();
  int *sopi = ref->nsopi();
  int *orbspi = ref->nmopi();
  int *openpi = ref->soccpi();
  int *clsdpi = init_int_array(nirreps);
  for(int h = 0; h < nirreps; ++h)
      clsdpi[h] = ref->doccpi()[h];

  int *frdocc = ref->frzcpi();
  int *fruocc = ref->frzvpi();
  int nfzc = 0; int nfzv = 0;
  for(int i=0; i < nirreps; i++) { nfzc += frdocc[i]; nfzv += fruocc[i]; }

  for(int i=0; i < nirreps; i++) clsdpi[i] -= frdocc[i];

  int *uoccpi = init_int_array(nirreps);
  for(int i=0; i < nirreps; i++)
    uoccpi[i] = orbspi[i] - clsdpi[i] - openpi[i] - fruocc[i] - frdocc[i];
  int nclsd = 0; int nopen = 0; int nuocc = 0;
  for(int i=0; i < nirreps; i++) {
    nclsd += clsdpi[i]; nopen += openpi[i]; nuocc += uoccpi[i];
  }
  int nactive = nclsd + nopen + nuocc;

  psio->open(PSIF_CC_INFO, PSIO_OPEN_OLD);

  psio->write_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(reference), sizeof(int));
  psio->write_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep", (char *) frdocc, sizeof(int)*nirreps);
  psio->write_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep", (char *) fruocc, sizeof(int)*nirreps);
  psio->write_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive), sizeof(int));

  vector<int> aoccpi, boccpi, avirtpi, bvirtpi;
  vector<int> occpi, virtpi;
  if(reference == 2) { // UHF or semicanonical
    for(int h=0; h < nirreps; h++) {
      aoccpi.push_back(clsdpi[h] + openpi[h]);
      boccpi.push_back(clsdpi[h]);
      avirtpi.push_back(uoccpi[h]);
      bvirtpi.push_back(uoccpi[h] + openpi[h]);
    }
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Occ Orbs Per Irrep", (char *) aoccpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Occ Orbs Per Irrep", (char *) boccpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Virt Orbs Per Irrep", (char *) avirtpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Virt Orbs Per Irrep", (char *) bvirtpi.data(), sizeof(int)*nirreps);
  }
  else { // RHF or ROHF
    for(int h=0; h < nirreps; h++) {
      occpi.push_back(clsdpi[h] + openpi[h]);
      virtpi.push_back(uoccpi[h] + openpi[h]);
    }
    psio->write_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep", (char *) occpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep", (char *) virtpi.data(), sizeof(int)*nirreps);
  }

  psio->close(PSIF_CC_INFO, 1);

  outfile->Printf("\n\tWfn Parameters:\n");
  outfile->Printf("\t--------------------\n");
  outfile->Printf("\tPrint Level          = %d\n",print);
  outfile->Printf("\tNumber of irreps     = %d\n",nirreps);
  outfile->Printf("\tNumber of MOs        = %d\n",nmo);
  outfile->Printf("\tNumber of active MOs = %d\n",nactive);
  outfile->Printf("\tIRREP\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
  outfile->Printf(
            "\t-----\t-----\t------\t------\t------\t------\t------\n");
  for(int i=0; i < nirreps; i++) {
    outfile->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
                    labels[i],orbspi[i],frdocc[i],clsdpi[i],openpi[i],uoccpi[i],fruocc[i]);
  }
  outfile->Printf("\tNuclear Rep. energy    =  %20.14f\n", enuc);
  outfile->Printf(  "\tSCF energy             =  %20.14f\n", escf);

  std::vector<boost::shared_ptr<MOSpace> > transspaces;
  transspaces.push_back(MOSpace::occ);
  transspaces.push_back(MOSpace::vir);

  double efzc = trans(ref, options, transspaces);

  outfile->Printf(  "\tFrozen core energy     =  %20.14f\n", efzc);

  if(nfzc && (fabs(efzc) < 1e-7)) {
    outfile->Printf( "\tCCSORT Error: Orbitals are frozen in input,\n");
    outfile->Printf( "\tbut frozen core energy is small!\n");
    outfile->Printf( "\tCalculation will be aborted...\n");
    exit(PSI_RETURN_FAILURE);
  }
  else if(!nfzc && fabs(efzc)) {
    outfile->Printf( "\tCCSORT Warning: No orbitals are frozen,\n");
    outfile->Printf( "\tbut the frozen-core energy in wfn is non-zero.\n");
    outfile->Printf( "\tCalculation will continue with zero efzc...\n");
    efzc = 0.0;
  }

  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  dpdbuf4 K;
  std::vector<DPDMOSpace> spaces;
  int *cachefiles = init_int_array(PSIO_MAXUNIT);
  int **cachelist;
  if(reference == 2) { // UHF/semicanonical
    cachelist = cacheprep_uhf(options.get_int("CACHELEVEL"), cachefiles);

    DPDMOSpace aocc('O', "IJKLMN", aoccpi);
    DPDMOSpace avir('V', "ABCDEF", avirtpi);
    DPDMOSpace bocc('o', "ijklmn", boccpi);
    DPDMOSpace bvir('v', "abcdef", bvirtpi);
    spaces.push_back(aocc);
    spaces.push_back(avir);
    spaces.push_back(bocc);
    spaces.push_back(bvir);

    if(dpd_list[0]) throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[0] = new DPD(0, nirreps, Process::environment.get_memory(), 0, cachefiles, cachelist, NULL, 4, spaces);
    dpd_default = 0;
    global_dpd_ = dpd_list[0];

    psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "KL", "I>=J+", "K>=L+", 0, "MO Ints (OO|OO)");
    global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "IJ", "KL", "A <IJ|KL>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "IJ", "KL", 0, "A <IJ|KL>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "kl", "i>=j+", "k>=l+", 0, "MO Ints (oo|oo)");
    global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "ij", "kl", "A <ij|kl>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "kl", "I>=J+", "k>=l+", 0, "MO Ints (OO|oo)");
    global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "Ij", "Kl", "A <Ij|Kl>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "Ij", "Kl", 0, "A <Ij|Kl>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_AINTS, 1);

    psio->open(PSIF_CC_BINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "AB", "CD", "A>=B+", "C>=D+", 0, "MO Ints (VV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "AB", "CD", "B <AB|CD>");
    global_dpd_->buf4_close(&K);
    if(print > 10) {
      global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "AB", "CD", 0, "B <AB|CD>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ab", "cd", "a>=b+", "c>=d+", 0, "MO Ints (vv|vv)");
    global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "ab", "cd", "B <ab|cd>");
    global_dpd_->buf4_close(&K);
    if(print > 10) {
      global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "ab", "cd", 0, "B <ab|cd>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "AB", "cd", "A>=B+", "c>=d+", 0, "MO Ints (VV|vv)");
    global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "Ab", "Cd", "B <Ab|Cd>");
    global_dpd_->buf4_close(&K);
    if(print > 10) {
      global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "Ab", "Cd", 0, "B <Ab|Cd>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_BINTS, 1);

    psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "AB", "I>=J+", "A>=B+", 0, "MO Ints (OO|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "IA", "JB", "C <IA|JB>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "IA", "JB", 0, "C <IA|JB>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ab", "i>=j+", "a>=b+", 0, "MO Ints (oo|vv)");
    global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "ia", "jb", "C <ia|jb>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "ia", "jb", 0, "C <ia|jb>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "ab", "I>=J+", "a>=b+", 0, "MO Ints (OO|vv)");
    global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "Ia", "Jb", "C <Ia|Jb>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "Ia", "Jb", 0, "C <Ia|Jb>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "AB", "ij", "A>=B+", "i>=j+", 0, "MO Ints (VV|oo)");
    global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "Ai", "Bj", "C <Ai|Bj>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "Ai", "Bj", 0, "C <Ai|Bj>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_CINTS, 1);

    psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "JB", 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "IJ", "AB", "D <IJ|AB>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "IJ", "AB", 0, "D <IJ|AB>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "jb", "ia", "jb", 0, "MO Ints (ov|ov)");
    global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "ij", "ab", "D <ij|ab>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "ij", "ab", 0, "D <ij|ab>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "jb", 0, "MO Ints (OV|ov)");
    global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "Ij", "Ab", "D <Ij|Ab>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "Ij", "Ab", 0, "D <Ij|Ab>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_DINTS, 1);

    psio->open(PSIF_CC_EINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "KA", "I>=J+", "KA", 0, "MO Ints (OO|OV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "AI", "JK", "E <AI|JK>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "AI", "JK", 0, "E <AI|JK>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ka", "i>=j+", "ka", 0, "MO Ints (oo|ov)");
    global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "ai", "jk", "E <ai|jk>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "ai", "jk", 0, "E <ai|jk>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "KA", "ij", "KA", "i>=j+", 0, "MO Ints (OV|oo)");
    global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, qspr, "Ai", "Jk", "E <Ai|Jk>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "Ai", "Jk", 0, "E <Ai|Jk>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "ka", "I>=J+", "ka", 0, "MO Ints (OO|ov)");
    global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, prqs, "Ij", "Ka", "E <Ij|Ka>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "Ij", "Ka", 0, "E <Ij|Ka>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_EINTS, 1);

    psio->open(PSIF_CC_FINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "BC", "IA", "B>=C+", 0, "MO Ints (OV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "IA", "BC", "F <IA|BC>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "IA", "BC", 0, "F <IA|BC>");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "AI", "BC", "F <AI|BC>");
    global_dpd_->buf4_close(&K);
    if(print > 8) {
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "IA", "BC", 0, "F <IA|BC>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "bc", "ia", "b>=c+", 0, "MO Ints (ov|vv)");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "ia", "bc", "F <ia|bc>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "ai", "bc", "F <ai|bc>");
    global_dpd_->buf4_close(&K);
    if(print > 8) {
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "bc", "IA", "b>=c+", 0, "MO Ints (OV|vv)");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "Ia", "Bc", "F <Ia|Bc>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ia", "Bc", 0, "F <Ia|Bc>");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "aI", "bC", "F <aI|bC>");
    global_dpd_->buf4_close(&K);
    if(print > 8) {
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ia", "Bc", 0, "F <Ia|Bc>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "aI", "bC", 0, "F <aI|bC>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "BC", "ia", "B>=C+", "ia", 0, "MO Ints (VV|ov)");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "Ai", "Bc", "F <Ai|Bc>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ai", "Bc", 0, "F <Ai|Bc>");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, rspq, "Ab", "Ci", "F <Ab|Ci>");
    global_dpd_->buf4_close(&K);
    if(print > 8) {
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ai", "Bc", 0, "F <Ai|Bc>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ab", "Ci", 0, "F <Cb|Ci>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_FINTS, 1);
  }
  else { // RHF/ROHF
    cachelist = cacheprep_rhf(options.get_int("CACHELEVEL"), cachefiles);

    DPDMOSpace occ('o', "ijklmn", occpi);
    DPDMOSpace vir('v', "abcdef", virtpi);
    spaces.push_back(occ);
    spaces.push_back(vir);

    if(dpd_list[0]) throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[0] = new DPD(0, nirreps, Process::environment.get_memory(), 0, cachefiles, cachelist, NULL, 2, spaces);
    dpd_default = 0;
    global_dpd_ = dpd_list[0];

    psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "kl", "i>=j+", "k>=l+", 0, "MO Ints (OO|OO)");
    global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "ij", "kl", "A <ij|kl>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_AINTS, 1);

    psio->open(PSIF_CC_BINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ab", "cd", "a>=b+", "c>=d+", 0, "MO Ints (VV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "ab", "cd", "B <ab|cd>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "ab", "cd", 0, "B <ab|cd>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_BINTS, 1);

    psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ab", "i>=j+", "a>=b+", 0, "MO Ints (OO|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "ia", "jb", "C <ia|jb>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "ia", "jb", 0, "C <ia|jb>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_CINTS, 1);

    psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "jb", "ia", "jb", 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "ij", "ab", "D <ij|ab>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "ij", "ab", 0, "D <ij|ab>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_DINTS, 1);

    psio->open(PSIF_CC_EINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ka", "i>=j+", "ka", 0, "MO Ints (OO|OV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "ai", "jk", "E <ai|jk>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "ai", "jk", 0, "E <ai|jk>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_EINTS, 1);

    psio->open(PSIF_CC_FINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "bc", "ia", "b>=c+", 0, "MO Ints (OV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "ia", "bc", "F <ia|bc>");
    global_dpd_->buf4_close(&K);
    if(print > 6) {
      global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
      global_dpd_->buf4_print(&K, "outfile", 1);
      global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_FINTS, 1);

  }
  psio->close(PSIF_LIBTRANS_DPD, 0);

  for(int i =PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio->open(i,1);

  c_sort(reference);
  d_sort(reference);
  e_sort(reference);
  f_sort(reference);
  if(reference == 0) {
    d_spinad();
    e_spinad();
  }

  for(int i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio->close(i,1);
  for(int i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio->close(i,0); /* delete CC_TMP files */
  for(int i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio->close(i,1);

  return Success;
}

}} // End namespaces

