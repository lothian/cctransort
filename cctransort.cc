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

void cc_memcheck(int reference);

void sort_tei_rhf(boost::shared_ptr<PSIO> psio, int print);
void sort_tei_uhf(boost::shared_ptr<PSIO> psio, int print);

void c_sort(int reference);
void d_sort(int reference);
void e_sort(int reference);
void f_sort(int reference);
void d_spinad();
void e_spinad();

void fock_rhf(boost::shared_ptr<Wavefunction> ref, vector<int> &occpi, vector<int> &openpi,
              vector<int> &virpi, vector<int> &frdocc, int print);
void fock_uhf(boost::shared_ptr<Wavefunction> ref, vector<int> &aoccpi, vector<int> &boccpi,
              vector<int> &avirpi, vector<int> &bvirpi, vector<int> &frdocc, int print);

double scf_check(int reference);

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

  vector<int> orbspi, openpi, clsdpi, frdocc, fruocc;
  for(int h = 0; h < nirreps; ++h) {
    orbspi.push_back(ref->nmopi()[h]);
    openpi.push_back(ref->soccpi()[h]);
    clsdpi.push_back(ref->doccpi()[h]);
    frdocc.push_back(ref->frzcpi()[h]); 
    fruocc.push_back(ref->frzvpi()[h]);
  }

  int nfzc = 0; int nfzv = 0;
  for(int h=0; h < nirreps; h++) { nfzc += frdocc[h]; nfzv += fruocc[h]; }
  for(int h=0; h < nirreps; h++) clsdpi[h] -= frdocc[h];

  vector<int> uoccpi;
  for(int h=0; h < nirreps; h++) uoccpi.push_back(orbspi[h] - clsdpi[h] - openpi[h] - fruocc[h] - frdocc[h]);
  int nclsd = 0; int nopen = 0; int nuocc = 0;
  for(int h=0; h < nirreps; h++) {
    nclsd += clsdpi[h]; nopen += openpi[h]; nuocc += uoccpi[h];
  }
  int nactive = nclsd + nopen + nuocc;

  psio->open(PSIF_CC_INFO, PSIO_OPEN_OLD);

  psio->write_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(reference), sizeof(int));
  psio->write_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep", (char *) frdocc.data(), sizeof(int)*nirreps);
  psio->write_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep", (char *) fruocc.data(), sizeof(int)*nirreps);
  psio->write_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive), sizeof(int));

  vector<int> aoccpi, boccpi, avirpi, bvirpi;
  vector<int> occpi, virpi;
  if(reference == 2) { // UHF or semicanonical
    for(int h=0; h < nirreps; h++) {
      aoccpi.push_back(clsdpi[h] + openpi[h]);
      boccpi.push_back(clsdpi[h]);
      avirpi.push_back(uoccpi[h]);
      bvirpi.push_back(uoccpi[h] + openpi[h]);
    }
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Occ Orbs Per Irrep", (char *) aoccpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Occ Orbs Per Irrep", (char *) boccpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Virt Orbs Per Irrep", (char *) avirpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Virt Orbs Per Irrep", (char *) bvirpi.data(), sizeof(int)*nirreps);
  }
  else { // RHF or ROHF
    for(int h=0; h < nirreps; h++) {
      occpi.push_back(clsdpi[h] + openpi[h]);
      virpi.push_back(uoccpi[h] + openpi[h]);
    }
    psio->write_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep", (char *) occpi.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep", (char *) virpi.data(), sizeof(int)*nirreps);
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

  // Transformation

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

  // Set up DPD object
  std::vector<DPDMOSpace> spaces;
  int *cachefiles = init_int_array(PSIO_MAXUNIT);
  int **cachelist;
  if(reference == 2) { // UHF/semicanonical
    cachelist = cacheprep_uhf(options.get_int("CACHELEVEL"), cachefiles);

    DPDMOSpace aocc('O', "IJKLMN", aoccpi);
    DPDMOSpace avir('V', "ABCDEF", avirpi);
    DPDMOSpace bocc('o', "ijklmn", boccpi);
    DPDMOSpace bvir('v', "abcdef", bvirpi);
    spaces.push_back(aocc);
    spaces.push_back(avir);
    spaces.push_back(bocc);
    spaces.push_back(bvir);

    if(dpd_list[0]) throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[0] = new DPD(0, nirreps, Process::environment.get_memory(), 0, cachefiles, cachelist, NULL, 4, spaces);
    dpd_default = 0;
    global_dpd_ = dpd_list[0];
  }
  else { // RHF/ROHF
    cachelist = cacheprep_rhf(options.get_int("CACHELEVEL"), cachefiles);

    DPDMOSpace occ('o', "ijklmn", occpi);
    DPDMOSpace vir('v', "abcdef", virpi);
    spaces.push_back(occ);
    spaces.push_back(vir);

    if(dpd_list[0]) throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[0] = new DPD(0, nirreps, Process::environment.get_memory(), 0, cachefiles, cachelist, NULL, 2, spaces);
    dpd_default = 0;
    global_dpd_ = dpd_list[0];
  }

  cc_memcheck(reference);

  // Sort integrals into main categories
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  if(reference == 2) sort_tei_uhf(psio, print);
  else sort_tei_rhf(psio, print);
  psio->close(PSIF_LIBTRANS_DPD, 0);

  for(int i =PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio->open(i,1);

  // Generate additional orderings of basic integrals
  c_sort(reference);
  d_sort(reference);
  e_sort(reference);
  f_sort(reference);
  if(reference == 0) {
    d_spinad();
    e_spinad();
  }

  // Organize Fock matrices
  if(reference == 2) fock_uhf(ref, aoccpi, boccpi, avirpi, bvirpi, frdocc, print);
  else fock_rhf(ref, occpi, openpi, virpi, frdocc, print);

  double eref = scf_check(reference);
  outfile->Printf("Reference energy = %20.14f\n", eref + enuc + efzc);

  for(int i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio->close(i,1);
  for(int i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio->close(i,0); /* delete CC_TMP files */
  for(int i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio->close(i,1);

  return Success;
}

}} // End namespaces

