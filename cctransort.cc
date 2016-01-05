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
void b_spinad();
void a_spinad();
void d_spinad();
void e_spinad();

void fock_rhf(boost::shared_ptr<Wavefunction> ref, vector<int> &occpi, vector<int> &openpi,
              vector<int> &virpi, vector<int> &frdocc, int print);
void fock_uhf(boost::shared_ptr<Wavefunction> ref, vector<int> &aoccpi, vector<int> &boccpi,
              vector<int> &avirpi, vector<int> &bvirpi, vector<int> &frdocc, int print);

double scf_check(int reference, vector<int> &openpi);

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

  vector<int> cc_aocc, cc_bocc, cc_avir, cc_bvir;
  vector<int> qt_aocc, qt_bocc, qt_avir, qt_bvir;
  vector<int> aocc_sym, bocc_sym, avir_sym, bvir_sym;
  vector<int> aocc_off, bocc_off, avir_off, bvir_off;
  vector<int> cc_occ, cc_vir;
  vector<int> qt_occ, qt_vir;
  vector<int> occ_sym, vir_sym;
  if(reference == 2) {
    cc_aocc.assign(nactive, -1); cc_bocc.assign(nactive, -1);
    cc_avir.assign(nactive, -1); cc_bvir.assign(nactive, -1);
    qt_aocc.assign(nactive, -1); qt_bocc.assign(nactive, -1);
    qt_avir.assign(nactive, -1); qt_bvir.assign(nactive, -1);
    aocc_sym.assign(nactive, -1); bocc_sym.assign(nactive, -1);
    avir_sym.assign(nactive, -1); bvir_sym.assign(nactive, -1);
    for(int h=0, count=0, offset=0; h < nirreps; h++) {
      if(h) offset += clsdpi[h-1] + openpi[h-1];
      for(int i=0; i < clsdpi[h] + openpi[h]; i++, count++) {
        cc_aocc[offset+i] = count;
        qt_aocc[count] = nfzc + offset + i;
        aocc_sym[count] = h;
      }
    }
    for(int h=0, count=0, offset=0; h < nirreps; h++) {
      if(h) offset += clsdpi[h-1];
      for(int i=0; i < clsdpi[h]; i++, count++) {
        cc_bocc[offset+i] = count;
        qt_bocc[count] = nfzc + offset + i;
        bocc_sym[count] = h;
      }
    }
    for(int h=0, count=0, offset=nclsd+nopen; h < nirreps; h++) {
      if(h) offset += uoccpi[h-1];
      for(int i=0; i < uoccpi[h]; i++, count++) {
        cc_avir[offset+i] = count;
        qt_avir[count] = nfzc + offset + i;
        avir_sym[count] = h;
      }
    }
    for(int h=0, count=0, offset=nclsd; h < nirreps; h++) {
      if(h) offset += uoccpi[h-1] + openpi[h-1];
      for(int i=0; i < uoccpi[h] + openpi[h]; i++, count++) {
        cc_bvir[offset+i] = count;
        qt_bvir[count] = nfzc + offset + i;
        bvir_sym[count] = h;
      }
    }
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry", (char *) aocc_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry", (char *) bocc_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry", (char *) avir_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry", (char *) bvir_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Alpha Active Occ Order", (char *) cc_aocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Beta Active Occ Order", (char *) cc_bocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Alpha Active Virt Order", (char *) cc_avir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Beta Active Virt Order", (char *) cc_bvir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Alpha Active Occ Order", (char *) qt_aocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Beta Active Occ Order", (char *) qt_bocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Alpha Active Virt Order", (char *) qt_avir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Beta Active Virt Order", (char *) qt_bvir.data(), sizeof(int)*nactive);
  }
  else {
    cc_occ.assign(nactive, -1); cc_vir.assign(nactive, -1);
    qt_occ.assign(nactive, -1); qt_vir.assign(nactive, -1);
    occ_sym.assign(nactive, -1); vir_sym.assign(nactive, -1);
    for(int h=0, count=0, cl_offset=0, op_offset=nclsd; h < nirreps; h++) {
      if(h) cl_offset += clsdpi[h-1];
      for(int i=0; i < clsdpi[h]; i++, count++) {
        cc_occ[cl_offset+i] = count;
        qt_occ[count] = nfzc + cl_offset+i;
        occ_sym[count] = h;
      }
      if(h) op_offset += openpi[h-1];
      for(int i=0; i < openpi[h]; i++, count++) {
        cc_occ[op_offset+i] = count;
        qt_occ[count] = nfzc + op_offset+i;
        occ_sym[count] = h;
      }
    }

    for(int h=0, vr_offset=nclsd+nopen, op_offset=nclsd; h < nirreps; h++) {
      if(h) vr_offset += uoccpi[h-1];
      for(int i=0; i < uoccpi[h]; i++, count++) {
        cc_vir[vr_offset+i] = count;
        qt_vir[count] = nfzc + vr_offset+i;
        vir_sym[count] = h;
      }
      if(h) op_offset += openpi[h-1];
      for(int i=0; i < openpi[h]; i++, count++) {
        cc_vir[op_offset+i] = count;
        qt_vir[count] = nfzc + op_offset+i;
        vir_sym[count] = h;
      }
    }
    psio->write_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry", (char *) occ_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry", (char *) vir_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Active Occ Order", (char *) cc_occ.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Active Virt Order", (char *) cc_vir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Active Occ Order", (char *) qt_occ.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Active Virt Order", (char *) qt_vir.data(), sizeof(int)*nactive);
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
  psio->close(PSIF_LIBTRANS_DPD, 0); // delete file

  for(int i =PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio->open(i,1);

  // Generate additional orderings of basic integrals
  c_sort(reference);
  d_sort(reference);
  e_sort(reference);
  f_sort(reference);
  if(reference == 0) {
    b_spinad();
    a_spinad();
    d_spinad();
    e_spinad();
  }

  // Organize Fock matrices
  if(reference == 2) fock_uhf(ref, aoccpi, boccpi, avirpi, bvirpi, frdocc, print);
  else fock_rhf(ref, occpi, openpi, virpi, frdocc, print);

  double eref = scf_check(reference, openpi) + enuc + efzc;
  outfile->Printf("\tReference energy = %20.14f\n", eref);
  psio->write_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(eref), sizeof(double));

  for(int i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio->close(i,1);
  for(int i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio->close(i,0); /* delete CC_TMP files */
  for(int i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio->close(i,1);

  return Success;
}

}} // End namespaces

